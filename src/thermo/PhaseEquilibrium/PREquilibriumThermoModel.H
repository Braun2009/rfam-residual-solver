#ifndef PREquilibriumThermoModel_H
#define PREquilibriumThermoModel_H

#include "EquilibriumThermoModel.H"
#include "PRPhaseAuthorityIdentity.H"
#include "PRPhaseSoundSpeed.H"

#include <string>
#include <vector>

namespace Foam
{
namespace thermo
{

struct PREquilibriumThermoDataFingerprints
{
    std::string nasa7Sha256;
    std::string tabulatedIdealCaloricSha256;
    std::string speciesRegistrySha256;
    std::string implementationSchema;
};

PREquilibriumThermoDataFingerprints
authoritativePREquilibriumThermoDataFingerprints();

std::string makePREquilibriumThermoModelId(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    const phaseEquilibrium::RhoEFlashOptions& options,
    const phaseEquilibrium::PhaseSoundSpeedOptions& acousticOptions,
    const PREquilibriumThermoDataFingerprints& dataFingerprints);

class PREquilibriumThermoModel final : public EquilibriumThermoModel
{
    const std::vector<std::string> speciesNames_;
    const std::vector<double> speciesMolecularWeights_;
    const phaseEquilibrium::PRRhoEFlash rhoEFlash_;
    const phaseEquilibrium::PRPhaseSoundSpeed phaseSoundSpeed_;
    const std::string modelId_;
    const std::string physicalModelAuthorityFingerprint_;

public:
    PREquilibriumThermoModel(
        const SpeciesCriticalDB& criticalData,
        const BinaryInteractionDB& kij,
        const std::vector<std::string>& speciesNames,
        const phaseEquilibrium::RhoEFlashOptions& options,
        const phaseEquilibrium::PhaseSoundSpeedOptions& acousticOptions =
            phaseEquilibrium::PhaseSoundSpeedOptions());

    EquilibriumResult recoverRhoE(
        double density,
        double internalEnergy,
        const std::vector<double>& overallMassFractions,
        const EquilibriumSeed* seed = nullptr) const override;

    EquilibriumResult initializePTY(
        double pressure,
        double temperature,
        const std::vector<double>& overallMassFractions) const override;

    const std::string& modelId() const override { return modelId_; }
    std::string phaseLocalPhysicalModelFingerprint() const override
    {
        return physicalModelAuthorityFingerprint_;
    }
    const std::vector<std::string>& speciesNames() const override
    {
        return speciesNames_;
    }
    EnergyConvention energyConvention() const override
    {
        return EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT;
    }
    bool stateValidForModel(const EquilibriumState& state) const override;
    bool supportsPrimitivePTYInitialization() const override { return true; }
};

} // End namespace thermo
} // End namespace Foam

#endif
