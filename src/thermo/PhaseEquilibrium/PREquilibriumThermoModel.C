#include "PREquilibriumThermoModel.H"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <locale>
#include <optional>
#include <sstream>
#include <utility>

#include "NASA7Thermo.H"
#include "SpeciesRegistry.H"
#include "TabulatedIdealCaloric.H"
#include "ThermoDataProvenance.H"

namespace Foam
{
namespace thermo
{
namespace
{

using PRRegion = phaseEquilibrium::EquilibriumRegion;
using PRStatus = phaseEquilibrium::RhoEFlashStatus;

void appendString(
    std::ostringstream& manifest,
    const char* key,
    const std::string& value)
{
    manifest << key << '=' << value.size() << ':' << value << ';';
}

void appendSize(
    std::ostringstream& manifest,
    const char* key,
    std::size_t value)
{
    manifest << key << '=' << value << ';';
}

void appendDouble(
    std::ostringstream& manifest,
    const char* key,
    double value)
{
    manifest << key << '=';
    if (std::isnan(value))
    {
        manifest << "nan";
    }
    else if (std::isinf(value))
    {
        manifest << (value > 0.0 ? "+inf" : "-inf");
    }
    else
    {
        static_assert(sizeof(double) == sizeof(std::uint64_t),
            "RFAM provenance requires 64-bit double precision");
        std::uint64_t bits = 0;
        std::memcpy(&bits, &value, sizeof(bits));
        manifest << std::hex << std::setfill('0') << std::setw(16) << bits
                 << std::dec;
    }
    manifest << ';';
}

std::vector<double> speciesMolecularWeights(
    const std::vector<std::string>& speciesNames)
{
    SpeciesRegistry registry;
    std::vector<double> result;
    result.reserve(speciesNames.size());
    for (const std::string& name : speciesNames)
    {
        if (!registry.has(name)) return {};
        const double molecularWeight =
            registry.molecularWeight(registry.index(name));
        if (!std::isfinite(molecularWeight) || molecularWeight <= 0.0)
            return {};
        result.push_back(molecularWeight);
    }
    return result;
}

std::string makeModelIdImpl(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    const phaseEquilibrium::RhoEFlashOptions& options,
    const phaseEquilibrium::PhaseSoundSpeedOptions& acousticOptions,
    const PREquilibriumThermoDataFingerprints& dataFingerprints)
{
    std::ostringstream manifest;
    manifest.imbue(std::locale::classic());
    manifest << "schema=pr-nasa7-coolprop72-rhoe-acoustic-v6;"
             << "energy=sensible-internal-explicit-chemical-heat;";
    appendString(
        manifest, "data.nasa7.sha256", dataFingerprints.nasa7Sha256);
    appendString(
        manifest,
        "data.tabulatedIdealCaloric.sha256",
        dataFingerprints.tabulatedIdealCaloricSha256);
    appendString(
        manifest,
        "data.speciesRegistry.sha256",
        dataFingerprints.speciesRegistrySha256);
    appendString(
        manifest,
        "implementation.schema",
        dataFingerprints.implementationSchema);
    appendSize(manifest, "species.count", speciesNames.size());
    for (std::size_t i = 0; i < speciesNames.size(); ++i)
    {
        const std::string prefix = "species[" + std::to_string(i) + "].";
        appendString(manifest, (prefix + "name").c_str(), speciesNames[i]);
        const bool available = criticalData.has(speciesNames[i]);
        appendSize(manifest, (prefix + "critical.available").c_str(),
            available ? 1U : 0U);
        if (available)
        {
            const CriticalProperties& properties =
                criticalData.get(speciesNames[i]);
            appendDouble(manifest, (prefix + "Tc").c_str(), properties.Tc);
            appendDouble(manifest, (prefix + "Pc").c_str(), properties.Pc);
            appendDouble(
                manifest, (prefix + "omega").c_str(), properties.omega);
            appendDouble(manifest, (prefix + "MW").c_str(), properties.MW);
        }
    }
    for (std::size_t i = 0; i < speciesNames.size(); ++i)
    {
        for (std::size_t j = i; j < speciesNames.size(); ++j)
        {
            const std::string prefix = "kij[" + std::to_string(i) + ","
                + std::to_string(j) + "].";
            const bool available = kij.has(speciesNames[i], speciesNames[j]);
            appendSize(manifest, (prefix + "available").c_str(),
                available ? 1U : 0U);
            if (available)
            {
                const BinaryInteractionEntry* entry = kij.find(
                    speciesNames[i], speciesNames[j]);
                appendDouble(
                    manifest, (prefix + "value").c_str(), entry->value);
                appendString(
                    manifest,
                    (prefix + "evidence").c_str(),
                    binaryInteractionEvidenceClassName(
                        entry->metadata.evidence));
                appendString(
                    manifest,
                    (prefix + "zeroClass").c_str(),
                    binaryInteractionZeroClassName(
                        entry->metadata.zeroClass));
                appendString(
                    manifest,
                    (prefix + "sourceId").c_str(),
                    entry->metadata.sourceId);
                appendString(
                    manifest,
                    (prefix + "datasetId").c_str(),
                    entry->metadata.datasetId);
                appendDouble(
                    manifest,
                    (prefix + "minimumTemperature").c_str(),
                    entry->metadata.applicability.minimumTemperature);
                appendDouble(
                    manifest,
                    (prefix + "maximumTemperature").c_str(),
                    entry->metadata.applicability.maximumTemperature);
                appendDouble(
                    manifest,
                    (prefix + "minimumPressure").c_str(),
                    entry->metadata.applicability.minimumPressure);
                appendDouble(
                    manifest,
                    (prefix + "maximumPressure").c_str(),
                    entry->metadata.applicability.maximumPressure);
                appendDouble(
                    manifest,
                    (prefix + "absoluteUncertainty").c_str(),
                    entry->metadata.absoluteUncertainty);
            }
        }
    }

    appendDouble(manifest, "rhoe.minimumPressure", options.minimumPressure);
    appendDouble(manifest, "rhoe.maximumPressure", options.maximumPressure);
    appendDouble(
        manifest, "rhoe.minimumTemperature", options.minimumTemperature);
    appendDouble(
        manifest, "rhoe.maximumTemperature", options.maximumTemperature);
    appendDouble(manifest, "rhoe.initialPressure", options.initialPressure);
    appendDouble(
        manifest, "rhoe.initialTemperature", options.initialTemperature);
    appendSize(manifest, "rhoe.pressureScanIntervals",
        options.pressureScanIntervals);
    appendSize(manifest, "rhoe.maximumIterations", options.maximumIterations);
    appendDouble(
        manifest, "rhoe.finiteDifferenceStep", options.finiteDifferenceStep);
    appendDouble(manifest, "rhoe.saturationLogPressureTolerance",
        options.saturationLogPressureTolerance);
    appendDouble(manifest, "rhoe.logStepLimit", options.logStepLimit);
    appendDouble(manifest, "rhoe.densityAbsoluteTolerance",
        options.densityAbsoluteTolerance);
    appendDouble(manifest, "rhoe.densityRelativeTolerance",
        options.densityRelativeTolerance);
    appendDouble(manifest, "rhoe.energyAbsoluteTolerance",
        options.energyAbsoluteTolerance);
    appendDouble(manifest, "rhoe.energyRelativeTolerance",
        options.energyRelativeTolerance);
    appendDouble(manifest, "rhoe.tpRootRelativeTolerance",
        options.tpRootRelativeTolerance);

    const phaseEquilibrium::TPFlashOptions& tp = options.tpOptions;
    appendSize(manifest, "tp.maxFlashIterations", tp.maxFlashIterations);
    appendSize(
        manifest, "tp.maxStabilityIterations", tp.maxStabilityIterations);
    appendSize(manifest, "tp.maxRachfordRiceIterations",
        tp.maxRachfordRiceIterations);
    appendSize(manifest, "tp.saturationScanIntervals",
        tp.saturationScanIntervals);
    appendSize(manifest, "tp.maxSaturationIterations",
        tp.maxSaturationIterations);
    appendDouble(manifest, "tp.fugacityTolerance", tp.fugacityTolerance);
    appendDouble(manifest, "tp.materialTolerance", tp.materialTolerance);
    appendDouble(manifest, "tp.stabilityTolerance", tp.stabilityTolerance);
    appendDouble(manifest, "tp.iterationTolerance", tp.iterationTolerance);
    appendDouble(manifest, "tp.saturationTemperatureTolerance",
        tp.saturationTemperatureTolerance);
    appendDouble(manifest, "tp.successiveSubstitutionDamping",
        tp.successiveSubstitutionDamping);
    appendDouble(manifest, "mixture.compositionAbsolute",
        options.conservedTolerances.compositionAbsolute);
    appendDouble(manifest, "mixture.energyRelative",
        options.conservedTolerances.energyRelative);

    appendDouble(manifest, "acoustic.relativeDensityPerturbation",
        acousticOptions.relativeDensityPerturbation);
    appendDouble(manifest, "acoustic.minimumRelativeDensityPerturbation",
        acousticOptions.minimumRelativeDensityPerturbation);
    appendSize(manifest, "acoustic.entropyScanIntervals",
        acousticOptions.entropyScanIntervals);
    appendSize(manifest, "acoustic.maximumIterations",
        acousticOptions.maximumIterations);
    appendSize(manifest, "acoustic.maximumTopologyRefinements",
        acousticOptions.maximumTopologyRefinements);
    appendSize(manifest, "acoustic.minimumDerivativeConvergenceComparisons",
        acousticOptions.minimumDerivativeConvergenceComparisons);
    appendDouble(manifest, "acoustic.energyCorrectionSpanFactor",
        acousticOptions.energyCorrectionSpanFactor);
    appendDouble(manifest, "acoustic.entropyAbsoluteTolerance",
        acousticOptions.entropyAbsoluteTolerance);
    appendDouble(manifest, "acoustic.energyIntervalRelativeTolerance",
        acousticOptions.energyIntervalRelativeTolerance);
    appendDouble(manifest, "acoustic.maximumFlashDensityErrorFraction",
        acousticOptions.maximumFlashDensityErrorFraction);
    appendDouble(manifest, "acoustic.derivativeRelativeTolerance",
        acousticOptions.derivativeRelativeTolerance);

    return "peng-robinson-nasa7-coolprop72-rho-e-acoustic-v6:sha256="
        + thermoProvenance::sha256Hex(manifest.str());
}

EquilibriumStatus mapPhaseFailure(
    phaseEquilibrium::PhaseEquilibriumStatus status)
{
    using Status = phaseEquilibrium::PhaseEquilibriumStatus;
    switch (status)
    {
        case Status::INVALID_CONFIGURATION:
        case Status::UNKNOWN_SPECIES:
        case Status::INVALID_CRITICAL_PROPERTIES:
        case Status::MISSING_BINARY_INTERACTION:
            return EquilibriumStatus::INVALID_CONFIGURATION;
        case Status::INVALID_COMPOSITION_SIZE:
            return EquilibriumStatus::SPECIES_MISMATCH;
        case Status::NONFINITE_COMPOSITION:
        case Status::NEGATIVE_COMPOSITION:
        case Status::COMPOSITION_NOT_NORMALIZED:
            return EquilibriumStatus::INVALID_INPUT;
        case Status::RACHFORD_RICE_NO_BRACKET:
        case Status::PURE_SATURATION_NO_BRACKET:
            return EquilibriumStatus::NOT_BRACKETED;
        case Status::RACHFORD_RICE_NOT_CONVERGED:
        case Status::STABILITY_NOT_CONVERGED:
        case Status::FLASH_NOT_CONVERGED:
        case Status::PURE_SATURATION_NOT_CONVERGED:
            return EquilibriumStatus::NOT_CONVERGED;
        case Status::NONFINITE_PARAMETER:
        case Status::INVALID_LOG_ARGUMENT:
            return EquilibriumStatus::NONFINITE_PROPERTY;
        case Status::SUCCESS:
        case Status::INVALID_PRESSURE:
        case Status::INVALID_TEMPERATURE:
        case Status::NO_PHYSICAL_ROOT:
        case Status::WILSON_K_FAILURE:
            return EquilibriumStatus::NO_EQUILIBRIUM_STATE;
    }
    return EquilibriumStatus::NO_EQUILIBRIUM_STATE;
}

EquilibriumStatus mapCaloricFailure(
    phaseEquilibrium::PRPhaseCaloricStatus status)
{
    using Status = phaseEquilibrium::PRPhaseCaloricStatus;
    switch (status)
    {
        case Status::INVALID_CONFIGURATION:
        case Status::MISSING_CALORIC_DATA:
        case Status::MOLECULAR_WEIGHT_MISMATCH:
        case Status::MISSING_BINARY_INTERACTION:
            return EquilibriumStatus::INVALID_CONFIGURATION;
        case Status::UNSUPPORTED_CALORIC_RANGE:
            return EquilibriumStatus::UNSUPPORTED_PROPERTY_RANGE;
        case Status::NONFINITE_PROPERTY:
            return EquilibriumStatus::NONFINITE_PROPERTY;
        case Status::EOS_CLOSURE_FAILURE:
            return EquilibriumStatus::RESIDUAL_CONTRACT_FAILED;
        case Status::SUCCESS:
        case Status::INVALID_PRESSURE:
        case Status::INVALID_TEMPERATURE:
        case Status::INVALID_COMPOSITION:
        case Status::INVALID_PHASE_ROOT:
        case Status::PHASE_ROOT_FAILURE:
            return EquilibriumStatus::NO_EQUILIBRIUM_STATE;
    }
    return EquilibriumStatus::NO_EQUILIBRIUM_STATE;
}

EquilibriumStatus mapStatus(
    const phaseEquilibrium::RhoEFlashResult& backend,
    bool speciesMismatch)
{
    switch (backend.status)
    {
        case PRStatus::SUCCESS:
            return EquilibriumStatus::SUCCESS;
        case PRStatus::INVALID_CONFIGURATION:
        case PRStatus::INVALID_OPTIONS:
            return EquilibriumStatus::INVALID_CONFIGURATION;
        case PRStatus::INVALID_DENSITY:
        case PRStatus::NONFINITE_TARGET_ENERGY:
            return EquilibriumStatus::INVALID_INPUT;
        case PRStatus::INVALID_COMPOSITION:
            return speciesMismatch
                ? EquilibriumStatus::SPECIES_MISMATCH
                : EquilibriumStatus::INVALID_INPUT;
        case PRStatus::OUT_OF_CONFIGURED_RANGE:
            return EquilibriumStatus::UNSUPPORTED_PROPERTY_RANGE;
        case PRStatus::REQUIRED_REGION_MISMATCH:
            return EquilibriumStatus::NO_EQUILIBRIUM_STATE;
        case PRStatus::TP_FLASH_FAILURE:
            return mapPhaseFailure(backend.phaseStatus);
        case PRStatus::TP_CLOSURE_FAILURE:
            return backend.tpFlash.region ==
                    PRRegion::PURE_COMPONENT_COEXISTENCE
                ? EquilibriumStatus::UNDERDETERMINED_PHASE_FRACTION
                : EquilibriumStatus::RESIDUAL_CONTRACT_FAILED;
        case PRStatus::CALORIC_FAILURE:
            return mapCaloricFailure(backend.caloricStatus);
        case PRStatus::UNSUPPORTED_CALORIC_RANGE:
            return EquilibriumStatus::UNSUPPORTED_PROPERTY_RANGE;
        case PRStatus::CONSERVED_MIXTURE_FAILURE:
        case PRStatus::RESIDUAL_FAILURE:
            return EquilibriumStatus::RESIDUAL_CONTRACT_FAILED;
        case PRStatus::PURE_SATURATION_NO_BRACKET:
        case PRStatus::NO_VALID_SATURATION_INTERVAL:
        case PRStatus::NO_RHOE_BRACKET:
            return EquilibriumStatus::NOT_BRACKETED;
        case PRStatus::NOT_CONVERGED:
            return EquilibriumStatus::NOT_CONVERGED;
    }
    return EquilibriumStatus::NO_EQUILIBRIUM_STATE;
}

EquilibriumStatus mapAcousticFailure(
    phaseEquilibrium::PhaseSoundSpeedStatus status)
{
    using Status = phaseEquilibrium::PhaseSoundSpeedStatus;
    switch (status)
    {
        case Status::SUCCESS:
            return EquilibriumStatus::SUCCESS;
        case Status::INVALID_CONFIGURATION:
            return EquilibriumStatus::INVALID_CONFIGURATION;
        case Status::INVALID_INPUT:
            return EquilibriumStatus::INVALID_INPUT;
        case Status::BASE_STATE_INVALID:
        case Status::BASE_STATE_MISMATCH:
            return EquilibriumStatus::RESIDUAL_CONTRACT_FAILED;
        case Status::NONFINITE_ENTROPY:
            return EquilibriumStatus::NONFINITE_PROPERTY;
        case Status::NO_ENTROPY_BRACKET:
        case Status::FLASH_FAILURE:
        case Status::NOT_CONVERGED:
        case Status::PHASE_TOPOLOGY_CHANGE:
        case Status::INSUFFICIENT_FLASH_ACCURACY:
        case Status::PERTURBATION_NOT_CONVERGED:
        case Status::NONPOSITIVE_DERIVATIVE:
            return EquilibriumStatus::ACOUSTIC_CLOSURE_FAILED;
    }
    return EquilibriumStatus::ACOUSTIC_CLOSURE_FAILED;
}

EquilibriumPhaseProperties translatePhase(
    const phaseEquilibrium::PRPhaseCaloricProperties& phase)
{
    EquilibriumPhaseProperties translated;
    translated.massFractions = phase.massFractions;
    translated.moleFractions = phase.moleFractions;
    translated.molecularWeight = phase.molecularWeight;
    translated.rho = phase.rho;
    translated.h = phase.h;
    translated.e = phase.e;
    translated.s = phase.s;
    translated.compressibility = phase.compressibility;
    translated.entropyAvailable = std::isfinite(phase.s);
    return translated;
}

std::optional<EquilibriumRegion> translateSinglePhaseRegion(PRRegion region)
{
    switch (region)
    {
        case PRRegion::SINGLE_LIQUID:
            return EquilibriumRegion::SINGLE_LIQUID;
        case PRRegion::SINGLE_VAPOR:
            return EquilibriumRegion::SINGLE_VAPOR;
        case PRRegion::SINGLE_SUPERCRITICAL:
            return EquilibriumRegion::SUPERCRITICAL;
        case PRRegion::UNKNOWN:
        case PRRegion::TWO_PHASE:
        case PRRegion::PURE_COMPONENT_COEXISTENCE:
            return std::nullopt;
    }
    return std::nullopt;
}

std::optional<EquilibriumState> translateState(
    const phaseEquilibrium::RhoEFlashResult& backend,
    const phaseEquilibrium::PhaseSoundSpeedResult& acoustic,
    const std::vector<double>& overallMassFractions)
{
    if (!backend.success() || !acoustic.success()) return std::nullopt;

    EquilibriumState state;
    state.overallMassFractions = overallMassFractions;
    state.p = backend.pressure;
    state.T = backend.temperature;
    state.acoustic.closure = AcousticClosure::HOMOGENEOUS_EQUILIBRIUM;
    state.acoustic.evidence = AcousticEvidence::CONVERGED_RHOE_ISENTROPE;
    state.acoustic.pressureDensityDerivative = acoustic.soundSpeedSquared;

    if (backend.tpFlash.region == PRRegion::TWO_PHASE
        || backend.tpFlash.region == PRRegion::PURE_COMPONENT_COEXISTENCE)
    {
        if (!backend.liquidPhase.has_value()
            || !backend.vaporPhase.has_value()
            || !backend.conservedMixture.has_value())
            return std::nullopt;

        const phaseEquilibrium::ConservedPhaseMixture& mixture =
            *backend.conservedMixture;
        state.region = EquilibriumRegion::TWO_PHASE;
        state.rho = mixture.rho;
        state.h = mixture.h;
        state.e = mixture.e;
        state.s = mixture.s;
        state.vaporMoleFraction = mixture.vaporMoleFraction;
        state.vaporMassFraction = mixture.vaporMassFraction;
        state.vaporVolumeFraction = mixture.vaporVolumeFraction;
        state.liquid = translatePhase(*backend.liquidPhase);
        state.vapor = translatePhase(*backend.vaporPhase);
    }
    else
    {
        const std::optional<EquilibriumRegion> region =
            translateSinglePhaseRegion(backend.tpFlash.region);
        if (!region.has_value() || !backend.singlePhase.has_value())
            return std::nullopt;

        const phaseEquilibrium::PRPhaseCaloricProperties& phase =
            *backend.singlePhase;
        state.region = *region;
        state.rho = phase.rho;
        state.h = phase.h;
        state.e = phase.e;
        state.s = phase.s;
        state.singlePhase = translatePhase(phase);
    }

    state.entropyAvailable = std::isfinite(state.s);
    if (!state.valid()) return std::nullopt;
    return state;
}

double stabilityViolation(
    const phaseEquilibrium::RhoEFlashResult& backend)
{
    if (!backend.success())
        return std::numeric_limits<double>::infinity();

    const auto violation = [](double minimumTpd)
    {
        return std::isfinite(minimumTpd)
            ? std::max(0.0, -minimumTpd)
            : std::numeric_limits<double>::infinity();
    };

    if (backend.tpFlash.region == PRRegion::TWO_PHASE)
    {
        const double liquid = violation(backend.tpFlash.liquidMinimumTpd);
        const double vapor = violation(backend.tpFlash.vaporMinimumTpd);
        if (!std::isfinite(liquid) || !std::isfinite(vapor))
            return std::numeric_limits<double>::infinity();
        return std::max(liquid, vapor);
    }
    if (backend.tpFlash.region == PRRegion::SINGLE_LIQUID
        || backend.tpFlash.region == PRRegion::SINGLE_VAPOR
        || backend.tpFlash.region == PRRegion::SINGLE_SUPERCRITICAL)
        return violation(backend.tpFlash.feedMinimumTpd);

    return std::numeric_limits<double>::infinity();
}

EquilibriumDiagnostics translateDiagnostics(
    const phaseEquilibrium::RhoEFlashResult& backend,
    const std::string& modelId)
{
    EquilibriumDiagnostics diagnostics;
    diagnostics.modelId = modelId;
    diagnostics.backendStatus =
        phaseEquilibrium::rhoEFlashStatusName(backend.status);
    diagnostics.backendPhaseStatus =
        phaseEquilibrium::statusName(backend.phaseStatus);
    diagnostics.backendCaloricStatus =
        phaseEquilibrium::phaseCaloricStatusName(backend.caloricStatus);
    diagnostics.backendAcousticStatus = "NOT_EVALUATED";
    diagnostics.densityResidual = backend.densityResidual;
    diagnostics.energyResidual = backend.energyResidual;
    diagnostics.materialResidual = std::max(
        backend.tpMaterialResidual, backend.mixtureMaterialResidual);
    diagnostics.fugacityResidual = backend.tpFugacityResidual;
    diagnostics.stabilityResidual = stabilityViolation(backend);
    diagnostics.energyIdentityResidual = std::max(
        backend.phaseEnergyIdentityResidual,
        backend.mixtureEnergyIdentityResidual);
    diagnostics.iterations = backend.iterations;
    diagnostics.validSamples = backend.validSamples;
    diagnostics.invalidSamples = backend.invalidSamples;
    return diagnostics;
}

void appendAcousticDiagnostics(
    const phaseEquilibrium::PhaseSoundSpeedResult& acoustic,
    EquilibriumDiagnostics& diagnostics)
{
    diagnostics.backendAcousticStatus =
        phaseEquilibrium::phaseSoundSpeedStatusName(acoustic.status);
    diagnostics.acousticLowerEntropyResidual =
        acoustic.lowerEntropyResidual;
    diagnostics.acousticUpperEntropyResidual =
        acoustic.upperEntropyResidual;
    diagnostics.acousticUsedRelativeDensityPerturbation =
        acoustic.usedRelativeDensityPerturbation;
    diagnostics.acousticDerivativeRelativeChange =
        acoustic.derivativeRelativeChange;
    diagnostics.acousticDerivativeConvergenceComparisons =
        acoustic.derivativeConvergenceComparisons;
}

} // End anonymous namespace

PREquilibriumThermoDataFingerprints
authoritativePREquilibriumThermoDataFingerprints()
{
    const NASA7Thermo& nasa = nasa7Thermo();
    const phaseEquilibrium::TabulatedIdealCaloric& tabulated =
        phaseEquilibrium::tabulatedIdealCaloric();
    const SpeciesRegistry registry;
    PREquilibriumThermoDataFingerprints result;
    result.nasa7Sha256 = nasa.provenanceSha256();
    result.tabulatedIdealCaloricSha256 = tabulated.provenanceSha256();
    result.speciesRegistrySha256 = registry.provenanceSha256();
    result.implementationSchema =
        "rfam.pr-phase-caloric-rhoe-hem-acoustic/v2";
    return result;
}

std::string makePREquilibriumThermoModelId(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    const phaseEquilibrium::RhoEFlashOptions& options,
    const phaseEquilibrium::PhaseSoundSpeedOptions& acousticOptions,
    const PREquilibriumThermoDataFingerprints& dataFingerprints)
{
    return makeModelIdImpl(
        criticalData,
        kij,
        speciesNames,
        options,
        acousticOptions,
        dataFingerprints);
}

PREquilibriumThermoModel::PREquilibriumThermoModel(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    const phaseEquilibrium::RhoEFlashOptions& options,
    const phaseEquilibrium::PhaseSoundSpeedOptions& acousticOptions)
:
    speciesNames_(speciesNames),
    speciesMolecularWeights_(speciesMolecularWeights(speciesNames_)),
    rhoEFlash_(criticalData, kij, speciesNames_, options),
    phaseSoundSpeed_(
        criticalData, kij, speciesNames_, options, acousticOptions),
    modelId_(makePREquilibriumThermoModelId(
        criticalData,
        kij,
        speciesNames,
        options,
        acousticOptions,
        authoritativePREquilibriumThermoDataFingerprints())),
    physicalModelAuthorityFingerprint_(
        phaseEquilibrium::makePRPhasePhysicalModelAuthorityFingerprint(
            criticalData, kij, speciesNames))
{}

bool PREquilibriumThermoModel::stateValidForModel(
    const EquilibriumState& state) const
{
    return validateEquilibriumSpeciesBasis(
            state, speciesMolecularWeights_)
        && state.acoustic.closure
            == AcousticClosure::HOMOGENEOUS_EQUILIBRIUM
        && state.acoustic.evidence
            == AcousticEvidence::CONVERGED_RHOE_ISENTROPE;
}

EquilibriumResult PREquilibriumThermoModel::initializePTY(
    double pressure,
    double temperature,
    const std::vector<double>& overallMassFractions) const
{
    const phaseEquilibrium::RhoEFlashResult backend =
        rhoEFlash_.initializePTYMassFractions(
            pressure, temperature, overallMassFractions);

    EquilibriumResult result;
    result.diagnostics = translateDiagnostics(backend, modelId_);
    result.status = mapStatus(
        backend, overallMassFractions.size() != speciesNames_.size());
    if (result.status != EquilibriumStatus::SUCCESS)
        return result;

    const phaseEquilibrium::PhaseSoundSpeedResult acoustic =
        phaseSoundSpeed_.evaluate(backend, overallMassFractions);
    appendAcousticDiagnostics(acoustic, result.diagnostics);
    if (!acoustic.success())
    {
        result.status = mapAcousticFailure(acoustic.status);
        return result;
    }

    result.state = translateState(backend, acoustic, overallMassFractions);
    if (!result.state.has_value()
        || !stateValidForModel(*result.state))
    {
        result.state.reset();
        result.status = EquilibriumStatus::RESIDUAL_CONTRACT_FAILED;
        result.diagnostics.backendStatus =
            "INVALID_COMMITTED_SPECIES_BASIS";
        return result;
    }
    return result;
}

EquilibriumResult PREquilibriumThermoModel::recoverRhoE(
    double density,
    double internalEnergy,
    const std::vector<double>& overallMassFractions,
    const EquilibriumSeed* seed) const
{
    phaseEquilibrium::RhoEFlashResult backend;
    if (seed)
    {
        phaseEquilibrium::RhoEFlashSeed backendSeed;
        backendSeed.pressure = seed->pressure;
        backendSeed.temperature = seed->temperature;
        backendSeed.localOnly = false;
        backend = rhoEFlash_.flashMassFractions(
            density, internalEnergy, overallMassFractions, backendSeed);
    }
    else
    {
        backend = rhoEFlash_.flashMassFractions(
            density, internalEnergy, overallMassFractions);
    }

    EquilibriumResult result;
    result.diagnostics = translateDiagnostics(backend, modelId_);
    result.status = mapStatus(
        backend, overallMassFractions.size() != speciesNames_.size());
    if (result.status != EquilibriumStatus::SUCCESS)
        return result;

    const phaseEquilibrium::PhaseSoundSpeedResult acoustic =
        phaseSoundSpeed_.evaluate(backend, overallMassFractions);
    appendAcousticDiagnostics(acoustic, result.diagnostics);
    if (!acoustic.success())
    {
        result.status = mapAcousticFailure(acoustic.status);
        return result;
    }

    result.state = translateState(backend, acoustic, overallMassFractions);
    if (!result.state.has_value()
        || !stateValidForModel(*result.state))
    {
        result.state.reset();
        result.status = EquilibriumStatus::RESIDUAL_CONTRACT_FAILED;
        result.diagnostics.backendStatus =
            "INVALID_COMMITTED_SPECIES_BASIS";
        return result;
    }
    return result;
}

} // End namespace thermo
} // End namespace Foam
