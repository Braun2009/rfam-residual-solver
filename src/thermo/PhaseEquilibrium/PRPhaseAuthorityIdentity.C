#include "PRPhaseAuthorityIdentity.H"

#include <cmath>

#include "NASA7Thermo.H"
#include "PRPhaseDisplacementAuthority.H"
#include "PRPhaseRhoERecovery.H"
#include "SpeciesRegistry.H"
#include "TabulatedIdealCaloric.H"
#include "ThermoDataProvenance.H"

#include <cstddef>

namespace Foam
{
namespace phaseEquilibrium
{
namespace
{

void appendPairAuthority(
    std::string& manifest,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    std::size_t first,
    std::size_t second)
{
    const std::string prefix = "pair[" + std::to_string(first) + ","
        + std::to_string(second) + "].";
    if (first == second)
    {
        thermoProvenance::appendSize(manifest, prefix + "available", 1U);
        thermoProvenance::appendDouble(manifest, prefix + "value", 0.0);
        thermoProvenance::appendString(
            manifest, prefix + "authority", "pr-pure-diagonal");
        return;
    }

    const BinaryInteractionEntry* entry = kij.find(
        speciesNames[first], speciesNames[second]);
    thermoProvenance::appendSize(
        manifest, prefix + "available", entry ? 1U : 0U);
    if (!entry) return;

    thermoProvenance::appendDouble(manifest, prefix + "value", entry->value);
    thermoProvenance::appendString(
        manifest,
        prefix + "evidence",
        binaryInteractionEvidenceClassName(entry->metadata.evidence));
    thermoProvenance::appendString(
        manifest,
        prefix + "zeroClass",
        binaryInteractionZeroClassName(entry->metadata.zeroClass));
    thermoProvenance::appendString(
        manifest, prefix + "sourceId", entry->metadata.sourceId);
    thermoProvenance::appendString(
        manifest, prefix + "datasetId", entry->metadata.datasetId);
    thermoProvenance::appendDouble(
        manifest,
        prefix + "minimumTemperature",
        entry->metadata.applicability.minimumTemperature);
    thermoProvenance::appendDouble(
        manifest,
        prefix + "maximumTemperature",
        entry->metadata.applicability.maximumTemperature);
    thermoProvenance::appendDouble(
        manifest,
        prefix + "minimumPressure",
        entry->metadata.applicability.minimumPressure);
    thermoProvenance::appendDouble(
        manifest,
        prefix + "maximumPressure",
        entry->metadata.applicability.maximumPressure);
    thermoProvenance::appendDouble(
        manifest,
        prefix + "absoluteUncertainty",
        entry->metadata.absoluteUncertainty);
}

std::string physicalManifest(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames)
{
    std::string manifest =
        "schema=rfam.pr-phase-physical-authority/v1;"
        "energy=sensible-internal-explicit-chemical-heat;"
        "alpha=peng-robinson-current;";
    thermoProvenance::appendString(
        manifest, "nasa7.sha256", nasa7Thermo().provenanceSha256());
    thermoProvenance::appendString(
        manifest,
        "tabulatedIdealCaloric.sha256",
        tabulatedIdealCaloric().provenanceSha256());
    const SpeciesRegistry registry;
    thermoProvenance::appendString(
        manifest,
        "speciesRegistry.sha256",
        registry.provenanceSha256());
    thermoProvenance::appendSize(
        manifest, "species.count", speciesNames.size());
    for (std::size_t i = 0; i < speciesNames.size(); ++i)
    {
        const std::string prefix = "species[" + std::to_string(i) + "].";
        thermoProvenance::appendString(
            manifest, prefix + "name", speciesNames[i]);
        const bool available = criticalData.has(speciesNames[i]);
        thermoProvenance::appendSize(
            manifest, prefix + "critical.available", available ? 1U : 0U);
        if (!available) continue;
        const CriticalProperties& properties = criticalData.get(
            speciesNames[i]);
        thermoProvenance::appendDouble(
            manifest, prefix + "Tc", properties.Tc);
        thermoProvenance::appendDouble(
            manifest, prefix + "Pc", properties.Pc);
        thermoProvenance::appendDouble(
            manifest, prefix + "omega", properties.omega);
        thermoProvenance::appendDouble(
            manifest, prefix + "MW", properties.MW);
    }
    for (std::size_t i = 0; i < speciesNames.size(); ++i)
    {
        for (std::size_t j = i; j < speciesNames.size(); ++j)
            appendPairAuthority(manifest, kij, speciesNames, i, j);
    }
    return manifest;
}

} // End anonymous namespace

std::string makePRPhasePhysicalModelAuthorityFingerprint(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames)
{
    return thermoProvenance::sha256Hex(
        physicalManifest(criticalData, kij, speciesNames));
}

std::string makePRPhaseRhoEAuthorityFingerprint(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    const PRPhaseRhoEOptions& options)
{
    std::string manifest = "schema=rfam.pr-phase-rhoe-authority/v1;";
    thermoProvenance::appendString(
        manifest,
        "physicalModel",
        makePRPhasePhysicalModelAuthorityFingerprint(
            criticalData, kij, speciesNames));
    thermoProvenance::appendDouble(
        manifest, "minimumTemperature", options.minimumTemperature);
    thermoProvenance::appendDouble(
        manifest, "maximumTemperature", options.maximumTemperature);
    thermoProvenance::appendSigned(
        manifest,
        "temperatureScanIntervals",
        options.temperatureScanIntervals);
    thermoProvenance::appendSigned(
        manifest, "maximumIterations", options.maximumIterations);
    thermoProvenance::appendDouble(
        manifest, "derivativeRelativeStep", options.derivativeRelativeStep);
    thermoProvenance::appendDouble(
        manifest, "energyAbsoluteTolerance", options.energyAbsoluteTolerance);
    thermoProvenance::appendDouble(
        manifest, "energyRelativeTolerance", options.energyRelativeTolerance);
    thermoProvenance::appendDouble(
        manifest,
        "densityRelativeTolerance",
        options.densityRelativeTolerance);
    thermoProvenance::appendDouble(
        manifest, "rootRelativeTolerance", options.rootRelativeTolerance);
    thermoProvenance::appendDouble(
        manifest,
        "derivativeIdentityRelativeTolerance",
        options.derivativeIdentityRelativeTolerance);
    thermoProvenance::appendDouble(
        manifest,
        "derivativeIdentityAbsoluteTolerance",
        options.derivativeIdentityAbsoluteTolerance);
    return thermoProvenance::sha256Hex(manifest);
}

std::string makePRPhaseDisplacementAuthorityFingerprint(
    const std::string& recoveryAuthorityFingerprint,
    const PRPhaseDisplacementBounds& bounds)
{
    std::string manifest =
        "schema=rfam.pr-phase-displacement-authority/v1;";
    thermoProvenance::appendString(
        manifest, "recoveryAuthority", recoveryAuthorityFingerprint);
    thermoProvenance::appendDouble(
        manifest,
        "maximumDensityRelativeDisplacement",
        bounds.maximumDensityRelativeDisplacement);
    thermoProvenance::appendDouble(
        manifest,
        "maximumEnergyRelativeDisplacement",
        bounds.maximumEnergyRelativeDisplacement);
    thermoProvenance::appendDouble(
        manifest,
        "maximumCompositionL1Displacement",
        bounds.maximumCompositionL1Displacement);
    thermoProvenance::appendSigned(
        manifest, "dyadicDepth", bounds.dyadicDepth);
    thermoProvenance::appendDouble(
        manifest, "minimumAdjacentRootGap", bounds.minimumAdjacentRootGap);
    thermoProvenance::appendDouble(
        manifest,
        "minimumDimensionlessDpDrhoMargin",
        bounds.minimumDimensionlessDpDrhoMargin);
    thermoProvenance::appendDouble(
        manifest,
        "maximumCompressibilityRelativeStep",
        bounds.maximumCompressibilityRelativeStep);
    thermoProvenance::appendDouble(
        manifest,
        "maximumPressureRelativeStep",
        bounds.maximumPressureRelativeStep);
    thermoProvenance::appendDouble(
        manifest,
        "maximumTemperatureRelativeStep",
        bounds.maximumTemperatureRelativeStep);
    thermoProvenance::appendDouble(
        manifest,
        "densityClosureRelativeTolerance",
        bounds.densityClosureRelativeTolerance);
    thermoProvenance::appendDouble(
        manifest,
        "energyClosureAbsoluteTolerance",
        bounds.energyClosureAbsoluteTolerance);
    thermoProvenance::appendDouble(
        manifest,
        "energyClosureRelativeTolerance",
        bounds.energyClosureRelativeTolerance);
    thermoProvenance::appendDouble(
        manifest,
        "compositionClosureAbsoluteTolerance",
        bounds.compositionClosureAbsoluteTolerance);
    return thermoProvenance::sha256Hex(manifest);
}

bool validPRPhaseAuthorityFingerprint(const std::string& fingerprint)
{
    if (fingerprint.size() != 64U
        || fingerprint == std::string(64U, '0'))
    {
        return false;
    }
    for (char character : fingerprint)
    {
        if (!((character >= '0' && character <= '9')
            || (character >= 'a' && character <= 'f')))
        {
            return false;
        }
    }
    return true;
}

} // End namespace phaseEquilibrium
} // End namespace Foam
