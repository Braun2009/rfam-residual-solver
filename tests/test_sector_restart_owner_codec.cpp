#include "framework/SectorRestartOwnerCodec.H"
#include "mesh/SectorMeshMetrics.H"
#include "thermo/MixingRules/BinaryInteractionDB.H"
#include "thermo/MixingRules/SpeciesCriticalDB.H"
#include "thermo/PhaseEquilibrium/PREquilibriumThermoModel.H"
#include "thermo/SpeciesRegistry.H"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace Foam;

namespace
{

int checks = 0;

void check(bool condition, const char* message)
{
    ++checks;
    if (!condition) throw std::runtime_error(message);
}

void putU32(std::string& bytes, std::size_t offset, std::uint32_t value)
{
    check(offset + 4U <= bytes.size(),"u32 mutation offset is in range");
    for (unsigned shift = 0; shift < 32; shift += 8)
        bytes[offset++] = static_cast<char>(value >> shift);
}

void putDouble(std::string& bytes, std::size_t offset, double value)
{
    std::uint64_t bits = 0;
    std::memcpy(&bits,&value,sizeof(bits));
    check(offset + 8U <= bytes.size(),"double mutation offset is in range");
    for (unsigned shift = 0; shift < 64; shift += 8)
        bytes[offset++] = static_cast<char>(bits >> shift);
}

std::uint32_t readU32(const std::string& bytes, std::size_t& offset)
{
    check(offset + 4U <= bytes.size(),"u32 read offset is in range");
    std::uint32_t value = 0;
    for (unsigned shift = 0; shift < 32; shift += 8)
        value |= static_cast<std::uint32_t>(
            static_cast<unsigned char>(bytes[offset++])) << shift;
    return value;
}

std::size_t firstFluidCellOffset(const std::string& bytes)
{
    // magic/version/dimensions/lengths/step counters
    std::size_t offset = 8U + 4U + 4U*4U + 8U*3U + 8U*2U;
    const std::uint32_t thermoLength = readU32(bytes,offset);
    offset += thermoLength;
    const std::uint32_t speciesCount = readU32(bytes,offset);
    for (std::uint32_t species = 0; species < speciesCount; ++species)
        offset += readU32(bytes,offset);
    offset += 8U; // cell count
    check(offset < bytes.size(),"first fluid cell offset is in range");
    return offset;
}

std::uint32_t fluidSpeciesCount(const std::string& bytes)
{
    std::size_t offset = 8U + 4U + 4U*3U;
    return readU32(bytes,offset);
}

std::size_t firstFluidTemperatureOffset(const std::string& bytes)
{
    return firstFluidCellOffset(bytes)
        + (static_cast<std::size_t>(fluidSpeciesCount(bytes)) + 7U)
            * sizeof(double);
}

std::size_t firstFluidPhaseFingerprintOffset(const std::string& bytes)
{
    std::size_t offset = firstFluidCellOffset(bytes)
        + (static_cast<std::size_t>(fluidSpeciesCount(bytes)) + 8U)
            * sizeof(double);
    const std::uint32_t length = readU32(bytes,offset);
    check(length == 64U,"fluid phase fingerprint has canonical length");
    check(offset + length <= bytes.size(),
        "fluid phase fingerprint is in range");
    return offset;
}

BinaryInteractionDB completeTestInteractions(const SpeciesRegistry& registry)
{
    BinaryInteractionDB result = BinaryInteractionDB::createDefault();
    const std::vector<std::string> names = registry.names();
    for (std::size_t i = 0; i < names.size(); ++i)
        for (std::size_t j = i + 1; j < names.size(); ++j)
            if (!result.has(names[i],names[j])) result.add(names[i],names[j],0.0);
    return result;
}

phaseEquilibrium::RhoEFlashOptions twoPhaseOptions()
{
    phaseEquilibrium::RhoEFlashOptions options;
    options.minimumPressure = 1.0e6;
    options.maximumPressure = 4.0e6;
    options.minimumTemperature = 250.0;
    options.maximumTemperature = 1500.0;
    options.initialPressure = 2.0e6;
    options.initialTemperature = 300.0;
    options.pressureScanIntervals = 64;
    options.maximumIterations = 120;
    options.finiteDifferenceStep = 1.0e-4;
    options.logStepLimit = 0.75;
    options.densityAbsoluteTolerance = 1.0e-7;
    options.densityRelativeTolerance = 2.0e-8;
    options.energyAbsoluteTolerance = 1.0e-4;
    options.energyRelativeTolerance = 2.0e-10;
    return options;
}

void initialize(
    Domain& domain,
    const SpeciesRegistry& registry,
    double pressure,
    double temperature)
{
    std::vector<double> composition(registry.size(),0.0);
    composition[registry.index("O2")] = 1.0;
    for (int k = 0; k < domain.Nz(); ++k)
        for (int j = 0; j < domain.Ny(); ++j)
            for (int i = 0; i < domain.Nx(); ++i)
            {
                PrimitiveVars& state = domain.stagedW(i,j,k);
                state.p = pressure + 100.0*i;
                state.T = temperature + 2.0*j + k;
                state.U = {0.25*i,-0.5*j,0.75*k};
                state.Y = composition;
            }
    check(domain.primitiveToConserved(),"source Domain initializes");
}

bool sameConserved(const ConservedVars& lhs, const ConservedVars& rhs)
{
    return lhs.rho == rhs.rho && lhs.rhoU == rhs.rhoU
        && lhs.rhoE == rhs.rhoE && lhs.rhoY == rhs.rhoY
        && lhs.rhoZ == rhs.rhoZ && lhs.rhoC == rhs.rhoC;
}

bool close(double lhs, double rhs, double relative = 1.0e-10,
           double absolute = 1.0e-10)
{
    if (std::isnan(lhs) || std::isnan(rhs))
        return std::isnan(lhs) && std::isnan(rhs);
    return std::abs(lhs-rhs)
        <= absolute + relative*std::max(std::abs(lhs),std::abs(rhs));
}

bool closeVector(const std::vector<double>& lhs,
                 const std::vector<double>& rhs)
{
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t index = 0; index < lhs.size(); ++index)
        if (!close(lhs[index],rhs[index])) return false;
    return true;
}

void compareRecovered(const Domain& expected, const Domain& actual)
{
    check(expected.N() == actual.N(),"recovered Domain cell count matches");
    for (int cell = 0; cell < expected.N(); ++cell)
    {
        check(sameConserved(expected.U(cell),actual.U(cell)),
            "recovered committed rho-e-Y matches exactly");
        check(sameConserved(expected.Unew(cell),actual.Unew(cell)),
            "recovered candidate rho-e-Y matches exactly");
        check(close(expected.W(cell).p,actual.W(cell).p)
                && close(expected.W(cell).T,actual.W(cell).T)
                && closeVector(expected.W(cell).U,actual.W(cell).U)
                && closeVector(expected.W(cell).Y,actual.W(cell).Y),
            "recovered primitive cache matches within restart tolerance");
        check(expected.Q(cell).region == actual.Q(cell).region
                && close(expected.Q(cell).p,actual.Q(cell).p)
                && close(expected.Q(cell).T,actual.Q(cell).T)
                && close(expected.Q(cell).rho,actual.Q(cell).rho)
                && close(expected.Q(cell).e,actual.Q(cell).e)
                && close(expected.Q(cell).vaporMassFraction,
                    actual.Q(cell).vaporMassFraction),
            "recovered equilibrium state matches within restart tolerance");
    }
}

} // namespace

int main()
{
    try
    {
        SpeciesRegistry registry;
        constexpr int nx = 8;
        constexpr int ny = 6;
        constexpr int nz = 6;
        SectorMeshMetrics mesh(vortexProfilePath(),nx,ny,nz,120.0);

        Domain source(
            nx,ny,nz,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
        source.setSpecies(registry.names());
        initialize(source,registry,4.0e6,700.0);
        check(source.beginStepTransaction(),"accepted ledger step begins");
        source.acceptStepTransaction();
        check(source.beginStepTransaction(),"rejected ledger step begins");
        source.rollbackStepTransaction();
        check(source.acceptedSteps() == 1U && source.rejectedSteps() == 1U,
            "source step ledger is nontrivial");

        SectorTurbulenceOperator sourceTurbulence(mesh,registry);
        std::vector<double> sourceK = sourceTurbulence.kField();
        std::vector<double> sourceOmega = sourceTurbulence.omegaField();
        for (double& value : sourceK) value *= 1.5;
        for (double& value : sourceOmega) value *= 1.25;
        sourceTurbulence.restoreFields(sourceK,sourceOmega);

        const SectorRestartEncodeResult fluid =
            SectorRestartOwnerCodec::encodeFluid(source);
        check(fluid.success() && !fluid.bytes.empty(),
            "fluid owner encodes at an accepted boundary");
        const SectorRestartEncodeResult turbulence =
            SectorRestartOwnerCodec::encodeTurbulence(
                source,sourceTurbulence);
        check(turbulence.success() && !turbulence.bytes.empty(),
            "turbulence owner encodes at an accepted boundary");

        Domain targetTemplate(
            nx,ny,nz,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
        targetTemplate.setSpecies(registry.names());
        SectorTurbulenceOperator targetTurbulence(mesh,registry);
        check(!targetTemplate.stateReady(),
            "detached decode target begins unpublished");

        SectorRestartFluidDecodeResult decodedFluid =
            SectorRestartOwnerCodec::decodeFluid(
                fluid.bytes,targetTemplate);
        check(decodedFluid.success()
                && decodedFluid.stagedDomain->stateReady()
                && !targetTemplate.stateReady(),
            "fluid decode constructs a fresh Domain without live publication");
        check(decodedFluid.stagedDomain->acceptedSteps() == 1U
                && decodedFluid.stagedDomain->rejectedSteps() == 1U,
            "fluid decode restores the Domain step ledger");
        compareRecovered(source,*decodedFluid.stagedDomain);
        const SectorRestartEncodeResult fluidReplay =
            SectorRestartOwnerCodec::encodeFluid(*decodedFluid.stagedDomain);
        check(fluidReplay.success()
                && fluidReplay.bytes == fluid.bytes,
            "detached fluid Domain re-encodes byte identically");

        const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
        const BinaryInteractionDB testKij = completeTestInteractions(registry);
        const auto hemModel =
            std::make_shared<const thermo::PREquilibriumThermoModel>(
                critical,testKij,registry.names(),twoPhaseOptions());
        std::vector<double> hemComposition(
            static_cast<std::size_t>(registry.size()),0.0);
        constexpr double productTrace = 1.0e-4;
        const double reactantScale = 1.0 - 2.0*productTrace;
        phaseEquilibrium::PRTPFlash feedFlash(
            critical,testKij,{"O2","C3H8"});
        const phaseEquilibrium::TPFlashResult feed =
            feedFlash.flashMoleFractions(2.0e6,300.0,{0.2,0.8});
        check(feed.valid()
                && feed.region
                    == phaseEquilibrium::EquilibriumRegion::TWO_PHASE,
            "HEM fixture has a two-phase O2-C3H8 feed");
        hemComposition[registry.index("O2")] =
            reactantScale*feed.feedMassFractions[0];
        hemComposition[registry.index("C3H8")] =
            reactantScale*feed.feedMassFractions[1];
        hemComposition[registry.index("CO")] = productTrace;
        hemComposition[registry.index("H2O")] = productTrace;
        const thermo::EquilibriumResult hemInitial =
            hemModel->initializePTY(2.0e6,300.0,hemComposition);
        check(hemInitial.success() && hemInitial.state.has_value()
                && hemInitial.state->region
                    == thermo::EquilibriumRegion::TWO_PHASE,
            "HEM model constructs a complete two-phase phase assertion");

        Domain hemSource(1,0.01,registry.size(),hemModel);
        hemSource.setSpecies(registry.names());
        PrimitiveVars hemPrimitive(registry.size());
        hemPrimitive.p = hemInitial.state->p;
        hemPrimitive.T = hemInitial.state->T;
        hemPrimitive.U = {1.0,-0.5,0.25};
        hemPrimitive.Y = hemComposition;
        check(hemSource.replaceCompleteField(
                std::vector<PrimitiveVars>{hemPrimitive},
                std::vector<thermo::EquilibriumState>{*hemInitial.state}),
            "HEM source publishes an authoritative rho-e-Y/Q pair");
        const SectorRestartEncodeResult hemFluid =
            SectorRestartOwnerCodec::encodeFluid(hemSource);
        check(hemFluid.success(),
            "HEM source encodes its phase assertion fingerprint");

        Domain hemTarget(1,0.01,registry.size(),hemModel);
        hemTarget.setSpecies(registry.names());
        const SectorRestartFluidDecodeResult decodedHem =
            SectorRestartOwnerCodec::decodeFluid(hemFluid.bytes,hemTarget);
        if (!decodedHem.success())
            std::cerr << "HEM decode status="
                      << sectorRestartOwnerCodecStatusName(decodedHem.status)
                      << " detail=" << decodedHem.detail << '\n';
        check(decodedHem.success()
                && decodedHem.stagedDomain->Q(0).region
                    == thermo::EquilibriumRegion::TWO_PHASE,
            "HEM rho-e-Y recovery reproduces the encoded two-phase assertion");
        compareRecovered(hemSource,*decodedHem.stagedDomain);
        const SectorRestartEncodeResult hemReplay =
            SectorRestartOwnerCodec::encodeFluid(*decodedHem.stagedDomain);
        check(hemReplay.success() && hemReplay.bytes == hemFluid.bytes,
            "HEM phase assertion survives a byte-identical detached replay");

        const std::size_t hemFingerprintOffset =
            firstFluidPhaseFingerprintOffset(hemFluid.bytes);
        const std::size_t legacyFingerprintOffset =
            firstFluidPhaseFingerprintOffset(fluid.bytes);
        std::string semanticMismatch = hemFluid.bytes;
        semanticMismatch.replace(
            hemFingerprintOffset,64U,fluid.bytes,legacyFingerprintOffset,64U);
        check(semanticMismatch.substr(hemFingerprintOffset,64U)
                != hemFluid.bytes.substr(hemFingerprintOffset,64U),
            "semantic mismatch uses a different canonical phase assertion");
        const SectorRestartEncodeResult hemBeforeMismatch =
            SectorRestartOwnerCodec::encodeFluid(hemSource);
        const SectorRestartFluidDecodeResult rejectedHem =
            SectorRestartOwnerCodec::decodeFluid(semanticMismatch,hemSource);
        check(rejectedHem.status
                    == SectorRestartOwnerCodecStatus::PHASE_ASSERTION_MISMATCH
                && !rejectedHem.stagedDomain,
            "valid but wrong HEM phase assertion rejects detached publication");
        const SectorRestartEncodeResult hemAfterMismatch =
            SectorRestartOwnerCodec::encodeFluid(hemSource);
        check(hemBeforeMismatch.success() && hemAfterMismatch.success()
                && hemAfterMismatch.bytes == hemBeforeMismatch.bytes,
            "phase assertion rejection leaves the live Domain byte-exact");

        std::string malformedFingerprint = hemFluid.bytes;
        malformedFingerprint[hemFingerprintOffset] = 'g';
        check(SectorRestartOwnerCodec::decodeFluid(
                malformedFingerprint,hemTarget).status
                == SectorRestartOwnerCodecStatus::INVALID_VALUE,
            "noncanonical phase fingerprint fails before recovery publication");
        malformedFingerprint = hemFluid.bytes;
        putU32(malformedFingerprint,hemFingerprintOffset-4U,63U);
        check(SectorRestartOwnerCodec::decodeFluid(
                malformedFingerprint,hemTarget).status
                == SectorRestartOwnerCodecStatus::INVALID_VALUE,
            "short phase fingerprint fails before recovery publication");

        SectorRestartTurbulenceDecodeResult decodedTurbulence =
            SectorRestartOwnerCodec::decodeTurbulence(
                turbulence.bytes,targetTemplate,targetTurbulence);
        check(decodedTurbulence.success()
                && targetTurbulence.kField() != decodedTurbulence.stagedState->k,
            "turbulence decode leaves the live target unpublished");
        SectorTurbulenceOperator stagedTurbulence(mesh,registry);
        stagedTurbulence.restoreFields(
            std::move(decodedTurbulence.stagedState->k),
            std::move(decodedTurbulence.stagedState->omega));
        check(stagedTurbulence.kField() == sourceK
                && stagedTurbulence.omegaField() == sourceOmega,
            "detached turbulence fields restore exactly");
        const SectorRestartEncodeResult turbulenceReplay =
            SectorRestartOwnerCodec::encodeTurbulence(
                source,stagedTurbulence);
        check(turbulenceReplay.success()
                && turbulenceReplay.bytes == turbulence.bytes,
            "detached turbulence state re-encodes byte identically");

        Domain divergent = source;
        divergent.Unew(0).rhoE += 1.0;
        check(SectorRestartOwnerCodec::encodeFluid(divergent).status
                == SectorRestartOwnerCodecStatus::INVALID_SOURCE,
            "fluid encode rejects a noncanonical candidate field");

        std::string mutation = fluid.bytes;
        mutation[0] ^= 1;
        check(SectorRestartOwnerCodec::decodeFluid(mutation,targetTemplate).status
                == SectorRestartOwnerCodecStatus::INVALID_MAGIC,
            "fluid magic mutation fails closed");
        mutation = fluid.bytes;
        mutation[6] = '1';
        putU32(mutation,8U,1U);
        check(SectorRestartOwnerCodec::decodeFluid(mutation,targetTemplate).status
                == SectorRestartOwnerCodecStatus::INVALID_MAGIC,
            "legacy fluid v1 payload cannot decode as phase-bound v2");
        mutation = fluid.bytes;
        putU32(mutation,8U,3U);
        check(SectorRestartOwnerCodec::decodeFluid(mutation,targetTemplate).status
                == SectorRestartOwnerCodecStatus::UNSUPPORTED_VERSION,
            "fluid version mutation fails closed");
        mutation = fluid.bytes;
        putU32(mutation,12U,static_cast<std::uint32_t>(nx+1));
        check(SectorRestartOwnerCodec::decodeFluid(mutation,targetTemplate).status
                == SectorRestartOwnerCodecStatus::STRUCTURE_MISMATCH,
            "fluid dimension mutation fails closed");
        mutation = fluid.bytes;
        mutation.push_back('x');
        check(SectorRestartOwnerCodec::decodeFluid(mutation,targetTemplate).status
                == SectorRestartOwnerCodecStatus::TRAILING_DATA,
            "fluid trailing data fails closed");
        mutation.assign(fluid.bytes.begin(),fluid.bytes.end()-1);
        check(SectorRestartOwnerCodec::decodeFluid(mutation,targetTemplate).status
                == SectorRestartOwnerCodecStatus::TRUNCATED,
            "fluid truncation fails closed");
        mutation = fluid.bytes;
        putDouble(mutation,firstFluidCellOffset(mutation),-1.0);
        check(SectorRestartOwnerCodec::decodeFluid(mutation,targetTemplate).status
                == SectorRestartOwnerCodecStatus::STATE_RECOVERY_FAILED,
            "invalid fluid conserved state fails detached recovery");
        mutation = fluid.bytes;
        putDouble(mutation,firstFluidTemperatureOffset(mutation),
            std::numeric_limits<double>::quiet_NaN());
        check(SectorRestartOwnerCodec::decodeFluid(mutation,targetTemplate).status
                == SectorRestartOwnerCodecStatus::INVALID_VALUE,
            "nonfinite fluid temperature seed fails closed");

        Domain wrongShape(
            nx+1,ny,nz,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
        wrongShape.setSpecies(registry.names());
        check(SectorRestartOwnerCodec::decodeFluid(fluid.bytes,wrongShape).status
                == SectorRestartOwnerCodecStatus::STRUCTURE_MISMATCH,
            "fluid target shape mismatch fails closed");

        mutation = turbulence.bytes;
        mutation[0] ^= 1;
        check(SectorRestartOwnerCodec::decodeTurbulence(
                mutation,targetTemplate,targetTurbulence).status
                == SectorRestartOwnerCodecStatus::INVALID_MAGIC,
            "turbulence magic mutation fails closed");
        mutation = turbulence.bytes;
        putU32(mutation,8U,2U);
        check(SectorRestartOwnerCodec::decodeTurbulence(
                mutation,targetTemplate,targetTurbulence).status
                == SectorRestartOwnerCodecStatus::UNSUPPORTED_VERSION,
            "turbulence version mutation fails closed");
        mutation = turbulence.bytes;
        mutation.push_back('x');
        check(SectorRestartOwnerCodec::decodeTurbulence(
                mutation,targetTemplate,targetTurbulence).status
                == SectorRestartOwnerCodecStatus::TRAILING_DATA,
            "turbulence trailing data fails closed");
        mutation.assign(turbulence.bytes.begin(),turbulence.bytes.end()-1);
        check(SectorRestartOwnerCodec::decodeTurbulence(
                mutation,targetTemplate,targetTurbulence).status
                == SectorRestartOwnerCodecStatus::TRUNCATED,
            "turbulence truncation fails closed");
        mutation = turbulence.bytes;
        putDouble(mutation,mutation.size()-8U,
            std::numeric_limits<double>::quiet_NaN());
        check(SectorRestartOwnerCodec::decodeTurbulence(
                mutation,targetTemplate,targetTurbulence).status
                == SectorRestartOwnerCodecStatus::TURBULENCE_REJECTED,
            "nonfinite turbulence field fails target validation");

        targetTurbulence.beginStepTransaction(targetTemplate);
        check(SectorRestartOwnerCodec::decodeTurbulence(
                turbulence.bytes,targetTemplate,targetTurbulence).status
                == SectorRestartOwnerCodecStatus::INVALID_TARGET,
            "active turbulence transaction blocks detached decode");
        targetTurbulence.rollbackStepTransaction();

        std::cout << "SECTOR_RESTART_OWNER_CODEC_PASS checks=" << checks
                  << " scope=detached-no-live-publication\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "SECTOR_RESTART_OWNER_CODEC_FAIL " << error.what()
                  << '\n';
        return 1;
    }
}
