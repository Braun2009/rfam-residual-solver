#ifndef CondensedCarbonModel_H
#define CondensedCarbonModel_H

#include "ConservedVars.H"
#include "EquilibriumThermoModel.H"
#include "SpeciesRegistry.H"

#include <array>
#include <limits>
#include <string>
#include <vector>

namespace Foam
{

using SpeciesVector = std::vector<double>;

enum class CondensedCarbonMomentModel
{
    ONE_MOMENT,
    TWO_MOMENT
};

const char* condensedCarbonMomentModelName(CondensedCarbonMomentModel model);

enum class CondensedCarbonApplicabilityStatus
{
    SUPPORTED,
    WARNING_RESEARCH_ONLY,
    OUTSIDE_DILUTE_LIMIT,
    SLIP_MODEL_REQUIRED,
    SOOT_ENERGY_MODEL_REQUIRED,
    FAILED
};

const char* condensedCarbonApplicabilityStatusName(
    CondensedCarbonApplicabilityStatus status);

enum class CondensedCarbonApplicabilityPolicy
{
    DIAGNOSTIC_ONLY,
    FAIL_CLOSED
};

struct CondensedCarbonState
{
    double rhoSoot = 0.0;       // kg C(s) / m3 control volume
    double numberDensity = 0.0; // particles / m3 control volume, two-moment only

    bool valid(CondensedCarbonMomentModel model) const noexcept;
};

struct WallCarbonInventory
{
    double massPerArea = 0.0; // kg C(s) / m2 wall

    bool valid() const noexcept;
};

struct CondensedCarbonExchange
{
    SpeciesVector fluidSpeciesMassRate; // kg / (m3 s), Domain species order
    double sootMassRate = 0.0;           // kg / (m3 s)
    double sootNumberRate = 0.0;         // particles / (m3 s), two-moment only
    double fluidEnergyRate = 0.0;        // W / m3, complete fluid rhoE source
};

struct WallCarbonTransfer
{
    double sootMassFluxToWall = 0.0;   // kg / (m2 s), positive toward wall
    double sootNumberFluxToWall = 0.0; // particles / (m2 s), two-moment only
    double wallHeatFlux = 0.0;         // W / m2; v1 requires zero
};

struct CondensedCarbonApplicabilityLimits
{
    double maximumMassLoading = std::numeric_limits<double>::quiet_NaN();
    double maximumVolumeFraction = std::numeric_limits<double>::quiet_NaN();
    double maximumHeatCapacityRatio = std::numeric_limits<double>::quiet_NaN();
    double maximumStokesNumber = std::numeric_limits<double>::quiet_NaN();

    bool valid() const noexcept;
};

struct CondensedCarbonOptions
{
    CondensedCarbonMomentModel momentModel =
        CondensedCarbonMomentModel::ONE_MOMENT;
    CondensedCarbonApplicabilityPolicy applicabilityPolicy =
        CondensedCarbonApplicabilityPolicy::DIAGNOSTIC_ONLY;

    double carbonMolarMass = 0.0120107;          // kg / mol C atoms
    double carbonFormationEnthalpy = 0.0;        // J / mol at 298.15 K
    double materialDensity = std::numeric_limits<double>::quiet_NaN(); // kg / m3
    double diagnosticSpecificHeat =
        std::numeric_limits<double>::quiet_NaN(); // J / (kg K), audit only
    double oneMomentParticleDiameter =
        std::numeric_limits<double>::quiet_NaN(); // m

    CondensedCarbonState initialState;
    CondensedCarbonState inletState;
    CondensedCarbonApplicabilityLimits limits;
    bool productionValidated = false;

    // These are explicit unsupported switches in v1. A future implementation
    // must add the corresponding conservative flux before enabling one.
    bool molecularDiffusionEnabled = false;
    bool turbulentDiffusionEnabled = false;
    bool thermophoresisEnabled = false;

    void validate() const;
};

struct CondensedCarbonApplicabilityInputs
{
    double rhoFluid = std::numeric_limits<double>::quiet_NaN();
    double rhoSoot = std::numeric_limits<double>::quiet_NaN();
    double fluidCp = std::numeric_limits<double>::quiet_NaN();
    double dynamicViscosity = std::numeric_limits<double>::quiet_NaN();
    double fluidTimeScale = std::numeric_limits<double>::quiet_NaN();
    double particleDiameter = std::numeric_limits<double>::quiet_NaN();
};

struct SootApplicabilityReport
{
    double massLoading = std::numeric_limits<double>::infinity();
    double volumeFraction = std::numeric_limits<double>::infinity();
    double heatCapacityRatio = std::numeric_limits<double>::infinity();
    double stokesNumber = std::numeric_limits<double>::infinity();
    CondensedCarbonApplicabilityStatus status =
        CondensedCarbonApplicabilityStatus::FAILED;

    bool productionSupported() const noexcept
    {
        return status == CondensedCarbonApplicabilityStatus::SUPPORTED;
    }
};

SootApplicabilityReport evaluateCondensedCarbonApplicability(
    const CondensedCarbonOptions& options,
    const CondensedCarbonApplicabilityInputs& inputs);

struct CondensedCarbonExchangeDiagnostics
{
    double fluidMassRate = std::numeric_limits<double>::quiet_NaN();
    double massResidual = std::numeric_limits<double>::infinity();
    std::array<double,4> elementResidual = {
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    double expectedFluidEnergyRate = std::numeric_limits<double>::quiet_NaN();
    double energyResidual = std::numeric_limits<double>::infinity();
};

struct CondensedCarbonExchangeResult
{
    bool accepted = false;
    std::string reason;
    ConservedVars fluidCandidate;
    CondensedCarbonState carbonCandidate;
    CondensedCarbonExchangeDiagnostics diagnostics;
};

class CondensedCarbonExchangeContract
{
public:
    static double expectedFluidEnergyRate(
        const PrimitiveVars& fluid,
        const std::vector<std::string>& speciesNames,
        const SpeciesRegistry& registry,
        const CondensedCarbonOptions& options,
        const CondensedCarbonExchange& exchange);

    static CondensedCarbonExchangeResult evaluate(
        const ConservedVars& fluidBase,
        const PrimitiveVars& fluid,
        const CondensedCarbonState& carbonBase,
        const std::vector<std::string>& speciesNames,
        const SpeciesRegistry& registry,
        const CondensedCarbonOptions& options,
        const CondensedCarbonExchange& exchange,
        double dt);
};

struct CondensedCarbonCellContext
{
    int i;
    int j;
    int k;
    const ConservedVars& fluidConserved;
    const PrimitiveVars& fluidPrimitive;
    const thermo::EquilibriumState& fluidThermo;
    const CondensedCarbonState& carbon;
};

class CondensedCarbonKinetics
{
public:
    virtual ~CondensedCarbonKinetics() = default;
    virtual std::string modelId() const = 0;
    virtual CondensedCarbonExchange evaluate(
        const CondensedCarbonCellContext& context) const = 0;
    virtual double maxStableDt(
        const CondensedCarbonCellContext&) const
    {
        return std::numeric_limits<double>::infinity();
    }
};

struct WallCarbonCellContext
{
    int i;
    int k;
    const ConservedVars& fluidConserved;
    const PrimitiveVars& fluidPrimitive;
    const thermo::EquilibriumState& fluidThermo;
    const CondensedCarbonState& carbon;
    const WallCarbonInventory& wall;
};

class WallCarbonDepositionModel
{
public:
    virtual ~WallCarbonDepositionModel() = default;
    virtual std::string modelId() const = 0;
    virtual WallCarbonTransfer evaluate(
        const WallCarbonCellContext& context) const = 0;
};

std::string condensedCarbonModelManifest(
    const CondensedCarbonOptions& options,
    const SpeciesRegistry& registry,
    const std::string& kineticsModelId,
    const std::string& wallModelId);

std::string condensedCarbonModelFingerprint(
    const CondensedCarbonOptions& options,
    const SpeciesRegistry& registry,
    const std::string& kineticsModelId,
    const std::string& wallModelId);

} // End namespace Foam

#endif
