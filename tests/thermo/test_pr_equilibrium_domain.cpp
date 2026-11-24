#include "PhaseEquilibrium/ConservedPhaseMixture.H"
#include "PhaseEquilibrium/PREquilibriumThermoModel.H"
#include "PhaseEquilibrium/PRPhaseCaloric.H"
#include "PhaseEquilibrium/PRTPFlash.H"
#include "PhaseEquilibrium/TabulatedIdealCaloric.H"
#include "NASA7Thermo.H"
#include "SectorCheckpoint.H"
#include "SolverBase.H"
#include "SpeciesRegistry.H"
#include "ThermoDataProvenance.H"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#ifdef RFAM_USE_OPENMP
#include <omp.h>
#endif

using namespace Foam;
using namespace Foam::phaseEquilibrium;

namespace
{

int failures = 0;

void check(bool condition, const std::string& label)
{
    std::cout << "  [" << (condition ? "PASS" : "FAIL") << "] "
              << label << '\n';
    if (!condition) ++failures;
}

bool close(double lhs, double rhs, double relativeTolerance = 1.0e-8)
{
    if (std::isnan(lhs) && std::isnan(rhs)) return true;
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) return lhs == rhs;
    return std::abs(lhs - rhs) <= relativeTolerance
        *std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

bool sameVector(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs,
    double tolerance = 1.0e-8)
{
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t i = 0; i < lhs.size(); ++i)
        if (!close(lhs[i], rhs[i], tolerance)) return false;
    return true;
}

bool samePhase(
    const thermo::EquilibriumPhaseProperties& lhs,
    const thermo::EquilibriumPhaseProperties& rhs)
{
    return sameVector(lhs.massFractions, rhs.massFractions)
        && sameVector(lhs.moleFractions, rhs.moleFractions)
        && close(lhs.molecularWeight, rhs.molecularWeight)
        && close(lhs.rho, rhs.rho)
        && close(lhs.h, rhs.h)
        && close(lhs.e, rhs.e)
        && close(lhs.s, rhs.s)
        && close(lhs.compressibility, rhs.compressibility)
        && lhs.entropyAvailable == rhs.entropyAvailable;
}

bool sameOptionalPhase(
    const std::optional<thermo::EquilibriumPhaseProperties>& lhs,
    const std::optional<thermo::EquilibriumPhaseProperties>& rhs)
{
    return lhs.has_value() == rhs.has_value()
        && (!lhs.has_value() || samePhase(*lhs, *rhs));
}

bool sameOptionalThermodynamics(
    const std::optional<thermo::EquilibriumThermodynamicPayload>& lhs,
    const std::optional<thermo::EquilibriumThermodynamicPayload>& rhs)
{
    if (lhs.has_value() != rhs.has_value()) return false;
    if (!lhs.has_value()) return true;
    return lhs->version == rhs->version
        && lhs->closure == rhs->closure
        && close(lhs->Cp, rhs->Cp)
        && close(lhs->Cv, rhs->Cv)
        && close(lhs->dpdrho_T, rhs->dpdrho_T)
        && close(lhs->dpdT_rho, rhs->dpdT_rho)
        && close(lhs->dedrho_T, rhs->dedrho_T)
        && close(lhs->dedT_rho, rhs->dedT_rho);
}

bool sameEquilibrium(
    const thermo::EquilibriumState& lhs,
    const thermo::EquilibriumState& rhs)
{
    return lhs.region == rhs.region
        && sameVector(lhs.overallMassFractions, rhs.overallMassFractions)
        && close(lhs.p, rhs.p)
        && close(lhs.T, rhs.T)
        && close(lhs.rho, rhs.rho)
        && close(lhs.h, rhs.h)
        && close(lhs.e, rhs.e)
        && close(lhs.s, rhs.s)
        && close(lhs.vaporMoleFraction, rhs.vaporMoleFraction)
        && close(lhs.vaporMassFraction, rhs.vaporMassFraction)
        && close(lhs.vaporVolumeFraction, rhs.vaporVolumeFraction)
        && lhs.acoustic.closure == rhs.acoustic.closure
        && lhs.acoustic.evidence == rhs.acoustic.evidence
        && close(lhs.acoustic.pressureDensityDerivative,
            rhs.acoustic.pressureDensityDerivative)
        && lhs.entropyAvailable == rhs.entropyAvailable
        && sameOptionalThermodynamics(
            lhs.thermodynamics, rhs.thermodynamics)
        && sameOptionalPhase(lhs.liquid, rhs.liquid)
        && sameOptionalPhase(lhs.vapor, rhs.vapor)
        && sameOptionalPhase(lhs.singlePhase, rhs.singlePhase);
}

bool sameConserved(const ConservedVars& lhs, const ConservedVars& rhs)
{
    return close(lhs.rho, rhs.rho)
        && sameVector(lhs.rhoU, rhs.rhoU)
        && close(lhs.rhoE, rhs.rhoE)
        && sameVector(lhs.rhoY, rhs.rhoY)
        && close(lhs.rhoZ, rhs.rhoZ)
        && close(lhs.rhoC, rhs.rhoC);
}

bool samePrimitive(const PrimitiveVars& lhs, const PrimitiveVars& rhs)
{
    return close(lhs.p, rhs.p)
        && close(lhs.T, rhs.T)
        && sameVector(lhs.U, rhs.U)
        && sameVector(lhs.Y, rhs.Y)
        && close(lhs.Z, rhs.Z)
        && close(lhs.C, rhs.C);
}

struct DomainSnapshot
{
    std::vector<ConservedVars> current;
    std::vector<ConservedVars> candidate;
    std::vector<PrimitiveVars> primitive;
    std::vector<thermo::EquilibriumState> equilibrium;
};

DomainSnapshot capture(const Domain& domain)
{
    DomainSnapshot result;
    for (int cell = 0; cell < domain.N(); ++cell)
    {
        result.current.push_back(domain.U(cell));
        result.candidate.push_back(domain.Unew(cell));
        result.primitive.push_back(domain.W(cell));
        result.equilibrium.push_back(domain.Q(cell));
    }
    return result;
}

bool sameDomain(const Domain& domain, const DomainSnapshot& expected)
{
    if (domain.N() != static_cast<int>(expected.current.size())) return false;
    for (int cell = 0; cell < domain.N(); ++cell)
    {
        if (!sameConserved(domain.U(cell), expected.current[cell])
            || !sameConserved(domain.Unew(cell), expected.candidate[cell])
            || !samePrimitive(domain.W(cell), expected.primitive[cell])
            || !sameEquilibrium(domain.Q(cell), expected.equilibrium[cell]))
            return false;
    }
    return true;
}

ConservedPhaseProperties conservedPhase(
    const PRPhaseCaloricProperties& phase)
{
    ConservedPhaseProperties result;
    result.molecularWeight = phase.molecularWeight;
    result.rho = phase.rho;
    result.h = phase.h;
    result.e = phase.e;
    result.s = phase.s;
    result.massFractions = phase.massFractions;
    return result;
}

struct ReferenceState
{
    bool valid = false;
    double pressure = 2.0e6;
    double temperature = 300.0;
    double rho = std::numeric_limits<double>::quiet_NaN();
    double e = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> massFractions;
};

ReferenceState makeReference(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& kij)
{
    ReferenceState reference;
    PRTPFlash tp(critical, kij, {"O2", "C3H8"});
    const TPFlashResult split = tp.flashMoleFractions(
        reference.pressure, reference.temperature, {0.2, 0.8});
    if (!split.valid() || split.region != EquilibriumRegion::TWO_PHASE
        || !split.phaseSplitReady())
        return reference;

    PRPhaseCaloric caloric(critical, kij, {"O2", "C3H8"});
    const PRPhaseCaloricResult liquid = caloric.evaluate(
        reference.pressure, reference.temperature,
        split.liquidMoleFractions, PhaseRoot::LIQUID);
    const PRPhaseCaloricResult vapor = caloric.evaluate(
        reference.pressure, reference.temperature,
        split.vaporMoleFractions, PhaseRoot::VAPOR);
    if (!liquid.success() || !vapor.success()) return reference;

    reference.massFractions = split.feedMassFractions;
    const ConservedMixtureResult mixture = assembleConservedPhaseMixture(
        reference.pressure, split.vaporMoleFraction,
        conservedPhase(*liquid.state), conservedPhase(*vapor.state),
        reference.massFractions);
    if (!mixture.success()) return reference;

    reference.rho = mixture.state->rho;
    reference.e = mixture.state->e;
    reference.valid = true;
    return reference;
}

RhoEFlashOptions flashOptions(const ReferenceState& reference)
{
    RhoEFlashOptions options;
    options.minimumPressure = 1.0e6;
    options.maximumPressure = 3.0e6;
    options.minimumTemperature = 280.0;
    options.maximumTemperature = 320.0;
    options.initialPressure = reference.pressure;
    options.initialTemperature = reference.temperature;
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

RhoEFlashOptions stateOptions(
    double minimumPressure,
    double maximumPressure,
    double minimumTemperature,
    double maximumTemperature,
    double initialPressure,
    double initialTemperature)
{
    RhoEFlashOptions options;
    options.minimumPressure = minimumPressure;
    options.maximumPressure = maximumPressure;
    options.minimumTemperature = minimumTemperature;
    options.maximumTemperature = maximumTemperature;
    options.initialPressure = initialPressure;
    options.initialTemperature = initialTemperature;
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

ConservedVars makeConserved(const ReferenceState& reference)
{
    const std::vector<double> velocity = {12.0, -3.0, 1.5};
    double kineticEnergy = 0.0;
    for (double component : velocity)
        kineticEnergy += 0.5*component*component;

    ConservedVars result(2);
    result.rho = reference.rho;
    for (std::size_t direction = 0; direction < velocity.size(); ++direction)
        result.rhoU[direction] = reference.rho*velocity[direction];
    result.rhoE = reference.rho*(reference.e + kineticEnergy);
    for (std::size_t species = 0;
         species < reference.massFractions.size(); ++species)
        result.rhoY[species] =
            reference.rho*reference.massFractions[species];
    result.rhoZ = 0.2*reference.rho;
    result.rhoC = 0.1*reference.rho;
    return result;
}

class StageTwoRejectOperator final : public Operator
{
    mutable int calls_ = 0;

public:
    void apply(Domain& domain, double) const override
    {
        ++calls_;
        if (calls_ == 2) domain.Unew(0).rho = -1.0;
    }

    int calls() const { return calls_; }
    std::string name() const override { return "PR stage-two rejection"; }
};

class ModelIdOverride final : public thermo::EquilibriumThermoModel
{
    std::shared_ptr<const thermo::EquilibriumThermoModel> delegate_;
    std::string modelId_;

public:
    ModelIdOverride(
        std::shared_ptr<const thermo::EquilibriumThermoModel> delegate,
        std::string modelId)
    :
        delegate_(std::move(delegate)),
        modelId_(std::move(modelId))
    {}

    thermo::EquilibriumResult recoverRhoE(
        double density,
        double internalEnergy,
        const std::vector<double>& massFractions,
        const thermo::EquilibriumSeed* seed = nullptr) const override
    {
        return delegate_->recoverRhoE(
            density, internalEnergy, massFractions, seed);
    }

    thermo::EquilibriumResult initializePTY(
        double pressure,
        double temperature,
        const std::vector<double>& massFractions) const override
    {
        return delegate_->initializePTY(
            pressure, temperature, massFractions);
    }

    const std::string& modelId() const override { return modelId_; }
    const std::vector<std::string>& speciesNames() const override
    {
        return delegate_->speciesNames();
    }
    EnergyConvention energyConvention() const override
    {
        return delegate_->energyConvention();
    }
    bool stateValidForModel(
        const thermo::EquilibriumState& state) const override
    {
        return delegate_->stateValidForModel(state);
    }
    bool supportsPrimitivePTYInitialization() const override
    {
        return delegate_->supportsPrimitivePTYInitialization();
    }
};

} // namespace

int main()
{
    std::cout << "=== PR equilibrium facade and Domain integration ===\n";
    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB kij = BinaryInteractionDB::createDefault();
    const std::vector<std::string> species = {"O2", "C3H8"};
    const ReferenceState reference = makeReference(critical, kij);
    check(reference.valid, "two-phase O2/C3H8 reference is constructible");
    if (!reference.valid) return 1;

    const auto model = std::make_shared<const thermo::PREquilibriumThermoModel>(
        critical, kij, species, flashOptions(reference));
    check(
        thermoProvenance::sha256Hex("abc")
            == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
        "thermodynamic provenance uses the standard SHA-256 digest");
    check(
        thermoProvenance::sha256Hex("")
            == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
        "thermodynamic provenance hashes an empty manifest correctly");
    check(
        thermoProvenance::sha256Hex(
            "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq")
            == "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1",
        "thermodynamic provenance hashes a multi-block manifest correctly");

    const NASA7Thermo& nasa = nasa7Thermo();
    std::vector<NASA7Thermo::Poly> changedNasaPolynomials =
        nasa.polynomials();
    changedNasaPolynomials.front().lo[0] += 1.0e-12;
    const std::string changedNasaSha = NASA7Thermo::provenanceSha256(
        nasa.speciesNames(),
        changedNasaPolynomials,
        nasa.speciesRegistry().entries());
    check(changedNasaSha != nasa.provenanceSha256(),
        "one NASA7 coefficient changes the authoritative data fingerprint");

    const phaseEquilibrium::TabulatedIdealCaloric& tabulated =
        phaseEquilibrium::tabulatedIdealCaloric();
    std::vector<phaseEquilibrium::TabulatedIdealCaloricTableData>
        changedTabulatedData = tabulated.data();
    changedTabulatedData.front().nodes.front().sensibleEnthalpy += 1.0e-9;
    const std::string changedTabulatedSha =
        phaseEquilibrium::TabulatedIdealCaloric::provenanceSha256(
            changedTabulatedData);
    check(changedTabulatedSha != tabulated.provenanceSha256(),
        "one low-temperature caloric point changes its data fingerprint");

    SpeciesRegistry registry;
    std::vector<SpeciesRegistry::Entry> changedRegistry = registry.entries();
    changedRegistry.front().MW += 1.0e-12;
    const std::string changedRegistrySha =
        SpeciesRegistry::provenanceSha256(changedRegistry);
    check(changedRegistrySha != registry.provenanceSha256(),
        "one SpeciesRegistry molecular weight changes its data fingerprint");

    thermo::PREquilibriumThermoDataFingerprints changedCaloricFingerprints =
        thermo::authoritativePREquilibriumThermoDataFingerprints();
    changedCaloricFingerprints.nasa7Sha256 = changedNasaSha;
    const std::string changedCaloricModelId =
        thermo::makePREquilibriumThermoModelId(
            critical,
            kij,
            species,
            flashOptions(reference),
            PhaseSoundSpeedOptions(),
            changedCaloricFingerprints);
    check(model->modelId() != changedCaloricModelId,
        "caloric data content participates in PR model provenance");
    const auto changedCaloricIdentityModel =
        std::make_shared<const ModelIdOverride>(
            model, changedCaloricModelId);
    const auto equivalentModel =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical, kij, species, flashOptions(reference));
    RhoEFlashOptions changedOptions = flashOptions(reference);
    ++changedOptions.maximumIterations;
    const auto changedOptionsModel =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical, kij, species, changedOptions);
    SpeciesCriticalDB changedCritical = critical;
    const CriticalProperties oxygen = changedCritical.get("O2");
    changedCritical.add(
        "O2", oxygen.Tc + 1.0e-6, oxygen.Pc, oxygen.omega, oxygen.MW);
    const auto changedDataModel =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            changedCritical, kij, species, flashOptions(reference));
    BinaryInteractionDB changedKij = kij;
    changedKij.add("O2", "C3H8", kij.get("O2", "C3H8") + 1.0e-9);
    const auto changedKijModel =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical, changedKij, species, flashOptions(reference));
    BinaryInteractionDB changedDiagonalKij = kij;
    changedDiagonalKij.add(
        "O2", "O2",
        (kij.has("O2", "O2") ? kij.get("O2", "O2") : 0.0) + 1.0e-9);
    const auto changedDiagonalKijModel =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical, changedDiagonalKij, species, flashOptions(reference));
    BinaryInteractionMetadata metadataA;
    metadataA.evidence = BinaryInteractionEvidenceClass::VALIDATED;
    metadataA.zeroClass = BinaryInteractionZeroClass::NOT_ZERO;
    metadataA.sourceId = "rfam-test-source";
    metadataA.datasetId = "dataset-a";
    metadataA.applicability.minimumTemperature = 80.0;
    metadataA.applicability.maximumTemperature = 3000.0;
    metadataA.applicability.minimumPressure = 1.0e5;
    metadataA.applicability.maximumPressure = 3.0e7;
    metadataA.absoluteUncertainty = 1.0e-4;
    BinaryInteractionDB metadataKijA = kij;
    metadataKijA.add(
        "O2", "C3H8", kij.get("O2", "C3H8"), metadataA);
    BinaryInteractionMetadata metadataB = metadataA;
    metadataB.datasetId = "dataset-b";
    BinaryInteractionDB metadataKijB = kij;
    metadataKijB.add(
        "O2", "C3H8", kij.get("O2", "C3H8"), metadataB);
    const auto metadataModelA =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical, metadataKijA, species, flashOptions(reference));
    const auto metadataModelB =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical, metadataKijB, species, flashOptions(reference));
    PhaseSoundSpeedOptions changedAcousticOptions;
    changedAcousticOptions.relativeDensityPerturbation *= 2.0;
    const auto changedAcousticModel =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical, kij, species, flashOptions(reference),
            changedAcousticOptions);
    check(model->modelId() == equivalentModel->modelId(),
        "equivalent PR configurations have identical provenance");
    check(model->modelId() != changedOptionsModel->modelId(),
        "rho-e options participate in PR provenance");
    check(model->modelId() != changedDataModel->modelId(),
        "critical-property data participate in PR provenance");
    check(model->modelId() != changedKijModel->modelId(),
        "binary-interaction data participate in PR provenance");
    check(model->modelId() != changedDiagonalKijModel->modelId(),
        "diagonal binary-interaction data participate in PR provenance");
    check(metadataModelA->modelId() != metadataModelB->modelId(),
        "binary-interaction dataset metadata participates in PR provenance");
    check(model->modelId() != changedAcousticModel->modelId(),
        "acoustic options participate in PR provenance");
    check(model->supportsPrimitivePTYInitialization(),
        "complete PR PTY initialization is exposed through the facade");

    const thermo::EquilibriumResult initializedTwoPhase =
        model->initializePTY(
            reference.pressure,
            reference.temperature,
            reference.massFractions);
    check(initializedTwoPhase.success()
            && initializedTwoPhase.state->region
                == thermo::EquilibriumRegion::TWO_PHASE
            && close(initializedTwoPhase.state->rho, reference.rho, 2.0e-8)
            && close(initializedTwoPhase.state->e, reference.e, 2.0e-8),
        "PTY initializes the determined O2/C3H8 two-phase split");

    const auto oxygenModel =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical,
            kij,
            std::vector<std::string>{"O2"},
            stateOptions(2.0e6, 6.0e6, 100.0, 220.0, 4.0e6, 130.0));
    const thermo::EquilibriumResult liquidOxygen =
        oxygenModel->initializePTY(4.0e6, 130.0, {1.0});
    check(liquidOxygen.success()
            && liquidOxygen.state->region
                == thermo::EquilibriumRegion::SINGLE_LIQUID
            && liquidOxygen.state->singlePhase.has_value()
            && close(liquidOxygen.state->p, 4.0e6)
            && close(liquidOxygen.state->T, 130.0),
        "PTY initializes the 4 MPa, 130 K liquid-oxygen state");
    if (liquidOxygen.success())
    {
        thermo::EquilibriumSeed oxygenSeed;
        oxygenSeed.pressure = liquidOxygen.state->p;
        oxygenSeed.temperature = liquidOxygen.state->T;
        const thermo::EquilibriumResult oxygenRoundTrip =
            oxygenModel->recoverRhoE(
                liquidOxygen.state->rho,
                liquidOxygen.state->e,
                {1.0},
                &oxygenSeed);
        check(oxygenRoundTrip.success()
                && close(oxygenRoundTrip.state->p, 4.0e6, 2.0e-8)
                && close(oxygenRoundTrip.state->T, 130.0, 2.0e-8),
            "liquid-oxygen PTY state round trips through rho-e");
    }

    const auto propaneModel =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical,
            kij,
            std::vector<std::string>{"C3H8"},
            stateOptions(2.0e6, 6.0e6, 350.0, 600.0, 4.0e6, 450.0));
    const thermo::EquilibriumResult hotPropane =
        propaneModel->initializePTY(4.0e6, 450.0, {1.0});
    check(hotPropane.success()
            && hotPropane.state->region
                != thermo::EquilibriumRegion::TWO_PHASE
            && hotPropane.state->singlePhase.has_value()
            && close(hotPropane.state->p, 4.0e6)
            && close(hotPropane.state->T, 450.0),
        "PTY initializes the 4 MPa, 450 K propane state");

    const auto injectionModel =
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical,
            kij,
            species,
            stateOptions(1.0e6, 6.0e6, 100.0, 3000.0, 4.0e6, 700.0));
    const thermo::EquilibriumResult hotTorch =
        injectionModel->initializePTY(4.0e6, 2400.0, {0.784, 0.216});
    check(hotTorch.success()
            && hotTorch.state->region
                != thermo::EquilibriumRegion::TWO_PHASE
            && hotTorch.state->singlePhase.has_value(),
        "PTY initializes the hot O2/C3H8 torch mixture");

    const auto checkConfiguredRangeRejection =
        [&](double pressure,double temperature,const std::string& label)
    {
        const thermo::EquilibriumResult rejected =
            injectionModel->initializePTY(
                pressure,temperature,{0.784,0.216});
        check(rejected.status
                    == thermo::EquilibriumStatus::UNSUPPORTED_PROPERTY_RANGE
                && !rejected.state.has_value()
                && rejected.diagnostics.backendStatus
                    == "OUT_OF_CONFIGURED_RANGE",
            label);
    };
    checkConfiguredRangeRejection(
        0.5e6,700.0,"PTY distinguishes pressure below configured range");
    checkConfiguredRangeRejection(
        6.5e6,700.0,"PTY distinguishes pressure above configured range");
    checkConfiguredRangeRejection(
        4.0e6,90.0,"PTY distinguishes temperature below configured range");
    checkConfiguredRangeRejection(
        4.0e6,3100.0,"PTY distinguishes temperature above configured range");

    PRTPFlash oxygenTP(critical, kij, {"O2"});
    const PureSaturationResult oxygenSaturation =
        oxygenTP.solvePureSaturationTemperature(4.0e6, 120.0, 154.5);
    check(oxygenSaturation.valid(),
        "pure O2 saturation reference is available");
    if (oxygenSaturation.valid())
    {
        const thermo::EquilibriumResult underdetermined =
            oxygenModel->initializePTY(
                4.0e6, oxygenSaturation.temperature, {1.0});
        check(underdetermined.status
                    == thermo::EquilibriumStatus::
                        UNDERDETERMINED_PHASE_FRACTION
                && !underdetermined.state.has_value(),
            "pure saturation PTY rejects the undetermined phase fraction");

        Domain atomicInitialization(2, 0.02, 1, oxygenModel);
        for (int cell = 0; cell < atomicInitialization.N(); ++cell)
        {
            PrimitiveVars& primitive = atomicInitialization.stagedW(cell);
            primitive.p = 4.0e6;
            primitive.T = 130.0;
            primitive.U = {0.0, 0.0, 0.0};
            primitive.Y = {1.0};
        }
        check(atomicInitialization.primitiveToConserved(),
            "valid liquid-oxygen primitive field initializes atomically");
        const DomainSnapshot beforeFailedInitialization =
            capture(atomicInitialization);
        atomicInitialization.stagedW(1).T = oxygenSaturation.temperature;
        check(!atomicInitialization.primitiveToConserved(),
            "one underdetermined PTY cell rejects the primitive field");
        check(sameDomain(
                atomicInitialization, beforeFailedInitialization),
            "failed PTY initialization preserves committed U/Unew/W/Q");
        check(atomicInitialization.failure().cell == 1
                && atomicInitialization.failure().equilibriumStatus
                    == thermo::EquilibriumStatus::
                        UNDERDETERMINED_PHASE_FRACTION,
            "failed PTY initialization reports the first bad cell");
    }

    thermo::EquilibriumSeed seed;
    seed.pressure = reference.pressure;
    seed.temperature = reference.temperature;
    const thermo::EquilibriumResult facade = model->recoverRhoE(
        reference.rho, reference.e, reference.massFractions, &seed);
    if (!facade.success())
    {
        std::cerr << "    facade status="
                  << thermo::equilibriumStatusName(facade.status)
                  << " backend=" << facade.diagnostics.backendStatus
                  << " densityResidual=" << facade.diagnostics.densityResidual
                  << " energyResidual=" << facade.diagnostics.energyResidual
                  << '\n';
    }
    check(facade.success(), "PREquilibriumThermoModel recovers rho-e-Y");
    check(facade.success()
            && facade.state->region == thermo::EquilibriumRegion::TWO_PHASE
            && facade.state->valid(),
        "facade publishes a closed generic two-phase Q");
    check(facade.success()
            && facade.state->acoustic.closure
                == thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM
            && facade.state->acoustic.evidence
                == thermo::AcousticEvidence::CONVERGED_RHOE_ISENTROPE
            && facade.state->acoustic.valid()
            && facade.diagnostics.backendAcousticStatus == "SUCCESS"
            && facade.diagnostics.acousticDerivativeConvergenceComparisons>=2
            && facade.diagnostics.acousticDerivativeRelativeChange<=2.0e-4,
        "facade SUCCESS includes internally converged HEM acoustic evidence");
    const PRPhaseSoundSpeed directAcoustic(
        critical, kij, species, flashOptions(reference));
    const PhaseSoundSpeedResult directAcousticState =
        directAcoustic.evaluateRhoE(
            reference.rho, reference.e, reference.massFractions);
    check(directAcousticState.success()
            && close(
                facade.state->acoustic.pressureDensityDerivative,
                directAcousticState.soundSpeedSquared, 5.0e-8),
        "facade Q acoustics match the authoritative rho-e HEM closure");
    check(facade.diagnostics.modelId == model->modelId(),
        "facade diagnostics carry thermo model provenance");
    check(facade.success()
            && model->stateValidForModel(*facade.state),
        "facade Q uses the configured species mass/mole/MW basis");

    thermo::EquilibriumState wrongMoleBasis = *facade.state;
    wrongMoleBasis.liquid->moleFractions = {0.5, 0.5};
    check(wrongMoleBasis.valid(),
        "wrong mole basis remains structurally well formed");
    check(!model->stateValidForModel(wrongMoleBasis),
        "model-aware validation rejects inconsistent mass/mole fractions");

    thermo::EquilibriumState wrongMolecularWeight = *facade.state;
    constexpr double molecularWeightScale = 1.001;
    wrongMolecularWeight.liquid->molecularWeight *= molecularWeightScale;
    wrongMolecularWeight.liquid->compressibility *= molecularWeightScale;
    wrongMolecularWeight.vapor->molecularWeight *= molecularWeightScale;
    wrongMolecularWeight.vapor->compressibility *= molecularWeightScale;
    check(wrongMolecularWeight.valid(),
        "synchronously altered MW/Z remains structurally well formed");
    check(!model->stateValidForModel(wrongMolecularWeight),
        "model-aware validation rejects an inconsistent molecular weight");

    const ConservedVars conserved = makeConserved(reference);
    const std::vector<ConservedVars> field(2, conserved);
    const std::vector<double> temperatureSeeds(2, reference.temperature);

    Domain serial(2, 0.02, 2, model);
#ifdef RFAM_USE_OPENMP
    const int originalThreadCount = omp_get_max_threads();
    omp_set_dynamic(0);
    omp_set_num_threads(1);
#endif
    check(serial.replaceConservedField(field, temperatureSeeds),
        "serial Domain recovers the PR conserved field");

    Domain parallel(2, 0.02, 2, model);
#ifdef RFAM_USE_OPENMP
    omp_set_num_threads(std::min(4, omp_get_num_procs()));
#endif
    check(parallel.replaceConservedField(field, temperatureSeeds),
        "OpenMP Domain recovers the PR conserved field");
#ifdef RFAM_USE_OPENMP
    omp_set_num_threads(originalThreadCount);
#endif
    check(sameDomain(parallel, capture(serial)),
        "serial and OpenMP U/W/Q recovery agree");
    check(parallel.Q(0).region == thermo::EquilibriumRegion::TWO_PHASE
            && parallel.Q(0).valid(),
        "Domain commits the two-phase Q cache");

    {
        Domain basisGuard(2, 0.02, 2, model);
        check(basisGuard.replaceConservedField(field, temperatureSeeds),
            "species-basis rejection Domain is initialized");
        const DomainSnapshot beforeBasisFailure = capture(basisGuard);
        std::vector<PrimitiveVars> primitives = {
            basisGuard.W(0), basisGuard.W(1)};
        std::vector<thermo::EquilibriumState> equilibrium = {
            basisGuard.Q(0), wrongMoleBasis};
        check(!basisGuard.replaceCompleteField(primitives, equilibrium),
            "Domain rejects a structurally valid Q on the wrong species basis");
        check(sameDomain(basisGuard, beforeBasisFailure),
            "species-basis rejection preserves complete U/Unew/W/Q state");
        check(basisGuard.failure().cell == 1
                && basisGuard.failure().equilibriumStatus
                    == thermo::EquilibriumStatus::RESIDUAL_CONTRACT_FAILED,
            "species-basis rejection reports the first invalid cell");
    }

    const std::string checkpointPath =
        "test_pr_equilibrium_domain_checkpoint.bin";
    bool checkpointWritten = false;
    try
    {
        SectorCheckpointMetadata metadata;
        metadata.time = 2.5e-6;
        metadata.iteration = 17;
        metadata.threadCount = 2;
        metadata.mode = "pr-equilibrium-domain-test";
        metadata.gitCommit = "test-tree";
        metadata.profileHash = 1234;
        SectorCheckpoint::write(checkpointPath, parallel, metadata);
        checkpointWritten = true;
        const SectorCheckpointInfo info =
            SectorCheckpoint::inspect(checkpointPath);
        check(info.metadata.thermoModelId == model->modelId(),
            "checkpoint header carries PR model provenance");

        Domain incompatible(2, 0.02, 2, changedOptionsModel);
        check(incompatible.replaceConservedField(field, temperatureSeeds),
            "provenance mismatch test Domain is initialized");
        const DomainSnapshot beforeMismatch = capture(incompatible);
        bool mismatchRejected = false;
        try
        {
            SectorCheckpoint::read(checkpointPath, incompatible);
        }
        catch (const std::runtime_error&)
        {
            mismatchRejected = true;
        }
        check(mismatchRejected,
            "checkpoint rejects a different PR option fingerprint");
        check(sameDomain(incompatible, beforeMismatch),
            "provenance rejection occurs before Domain mutation");

        Domain caloricMismatch(
            2, 0.02, 2, changedCaloricIdentityModel);
        check(caloricMismatch.replaceConservedField(field, temperatureSeeds),
            "caloric-data mismatch test Domain is initialized");
        const DomainSnapshot beforeCaloricMismatch = capture(caloricMismatch);
        bool caloricMismatchRejected = false;
        try
        {
            SectorCheckpoint::read(checkpointPath, caloricMismatch);
        }
        catch (const std::runtime_error&)
        {
            caloricMismatchRejected = true;
        }
        check(caloricMismatchRejected,
            "checkpoint rejects changed caloric data provenance");
        check(sameDomain(caloricMismatch, beforeCaloricMismatch),
            "caloric provenance rejection precedes Domain mutation");

        Domain restored(2, 0.02, 2, model);
        const SectorCheckpointMetadata loaded =
            SectorCheckpoint::read(checkpointPath, restored);
        check(loaded.iteration == metadata.iteration
                && close(loaded.time, metadata.time),
            "checkpoint metadata round trips");
        check(sameDomain(restored, capture(parallel)),
            "checkpoint reload reconstructs the complete PR U/W/Q field");

        const DomainSnapshot beforeRejectedReplacement = capture(restored);
        std::vector<ConservedVars> invalidField = field;
        invalidField[1].rhoE = std::numeric_limits<double>::quiet_NaN();
        check(!restored.replaceConservedField(
                invalidField, temperatureSeeds),
            "one invalid PR replacement cell rejects the whole field");
        check(sameDomain(restored, beforeRejectedReplacement),
            "failed PR field recovery preserves committed U/W/Q");
        check(restored.failure().cell == 1
                && restored.failure().equilibriumStatus
                    == thermo::EquilibriumStatus::INVALID_INPUT,
            "parallel recovery reports the deterministic first bad cell");
    }
    catch (const std::exception& error)
    {
        check(false, std::string("checkpoint integration threw: ") + error.what());
    }
    if (checkpointWritten) std::remove(checkpointPath.c_str());

    const DomainSnapshot beforeRK2 = capture(parallel);
    std::unique_ptr<StageTwoRejectOperator> rejector(
        new StageTwoRejectOperator());
    StageTwoRejectOperator* rejectorPtr = rejector.get();
    SolverBase solver(parallel);
    solver.setIntegrator(std::unique_ptr<TimeIntegrator>(new RK2Integrator()));
    solver.addOperator(std::move(rejector));
    check(!solver.step(1.0e-7),
        "RK2 rejects an invalid second-stage PR state");
    check(rejectorPtr->calls() == 2,
        "RK2 rejection occurs after a successful PR predictor recovery");
    check(solver.time() == 0.0 && solver.iteration() == 0,
        "rejected PR RK2 step does not advance the clock");
    check(parallel.acceptedSteps() == 0 && parallel.rejectedSteps() == 1,
        "PR RK2 rejection updates the transaction ledger once");
    check(sameDomain(parallel, beforeRK2),
        "PR RK2 rollback restores the exact pre-step U/W/Q field");

    if (failures != 0)
    {
        std::cerr << "PR_EQUILIBRIUM_DOMAIN_FAIL failures="
                  << failures << '\n';
        return 1;
    }
    std::cout << "PR_EQUILIBRIUM_DOMAIN_PASS\n";
    return 0;
}
