#ifndef EquilibriumThermoModel_H
#define EquilibriumThermoModel_H

#include "StateContract.H"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

namespace Foam
{
class GERGResidualTable;
class GERGSaturationTable;
namespace thermo
{

enum class EquilibriumStatus
{
    SUCCESS,
    INVALID_CONFIGURATION,
    INVALID_INPUT,
    SPECIES_MISMATCH,
    PRIMITIVE_INITIALIZATION_UNSUPPORTED,
    UNSUPPORTED_PROPERTY_RANGE,
    UNDERDETERMINED_PHASE_FRACTION,
    NO_EQUILIBRIUM_STATE,
    NOT_BRACKETED,
    NOT_CONVERGED,
    ACOUSTIC_CLOSURE_FAILED,
    NONFINITE_PROPERTY,
    RESIDUAL_CONTRACT_FAILED
};

const char* equilibriumStatusName(EquilibriumStatus status);

enum class EquilibriumRegion
{
    UNKNOWN,
    UNCLASSIFIED_SINGLE_PHASE,
    SINGLE_LIQUID,
    SINGLE_VAPOR,
    SUPERCRITICAL,
    TWO_PHASE,
    CONDENSED_SOLID
};

const char* equilibriumRegionName(EquilibriumRegion region);

enum class AcousticClosure
{
    UNKNOWN,
    FROZEN_COMPOSITION_SINGLE_PHASE,
    HOMOGENEOUS_EQUILIBRIUM
};

const char* acousticClosureName(AcousticClosure closure);

enum class AcousticEvidence
{
    UNVERIFIED,
    ANALYTIC_FROZEN_SINGLE_PHASE,
    CONVERGED_RHOE_ISENTROPE,
    // Equilibrium sound speed interpolated from a verified RFAMTAB2 property
    // table (explicitly enabled with RFAM_TETRA_TABLE). The provenance is the
    // table path plus the verification report; see docs/tetra_cold_flow_status.md.
    TABLE_INTERPOLATED
};

const char* acousticEvidenceName(AcousticEvidence evidence);

struct EquilibriumAcousticProperties
{
    AcousticClosure closure = AcousticClosure::UNKNOWN;
    AcousticEvidence evidence = AcousticEvidence::UNVERIFIED;
    double pressureDensityDerivative =
        std::numeric_limits<double>::quiet_NaN();
    // (d a / d p)_s at the state, used by the pressure outlet's second-order
    // Riemann-invariant estimate. Optional: backends that cannot supply it
    // leave it NaN and the outlet falls back to the exact isentrope integral.
    double soundSpeedPressureDerivative =
        std::numeric_limits<double>::quiet_NaN();

    bool valid() const;
    double soundSpeed() const;
};

struct EquilibriumThermodynamicPayload
{
    static constexpr std::uint32_t currentVersion = 1;

    std::uint32_t version = currentVersion;
    AcousticClosure closure = AcousticClosure::UNKNOWN;
    double Cp = std::numeric_limits<double>::quiet_NaN();
    double Cv = std::numeric_limits<double>::quiet_NaN();
    double dpdrho_T = std::numeric_limits<double>::quiet_NaN();
    double dpdT_rho = std::numeric_limits<double>::quiet_NaN();
    double dedrho_T = std::numeric_limits<double>::quiet_NaN();
    double dedT_rho = std::numeric_limits<double>::quiet_NaN();

    bool valid() const;
};

static_assert(
    std::is_trivially_copyable<EquilibriumThermodynamicPayload>::value,
    "per-cell thermodynamic payload must remain fixed-size");
static_assert(
    sizeof(EquilibriumThermodynamicPayload) <= 64,
    "per-cell thermodynamic payload exceeded its memory budget");

struct EquilibriumPhaseProperties
{
    std::vector<double> massFractions;
    std::vector<double> moleFractions;
    double molecularWeight = std::numeric_limits<double>::quiet_NaN();
    double rho = std::numeric_limits<double>::quiet_NaN();
    double h = std::numeric_limits<double>::quiet_NaN();
    double e = std::numeric_limits<double>::quiet_NaN();
    double s = std::numeric_limits<double>::quiet_NaN();
    double compressibility = std::numeric_limits<double>::quiet_NaN();
    bool entropyAvailable = false;

    bool valid() const;
};

// CFD-facing equilibrium state. This is intentionally independent of a PR
// result type and contains only a successfully committed thermodynamic state.
struct EquilibriumState
{
    EquilibriumRegion region = EquilibriumRegion::UNKNOWN;
    std::vector<double> overallMassFractions;
    double p = std::numeric_limits<double>::quiet_NaN();
    double T = std::numeric_limits<double>::quiet_NaN();
    double rho = std::numeric_limits<double>::quiet_NaN();
    double h = std::numeric_limits<double>::quiet_NaN();
    double e = std::numeric_limits<double>::quiet_NaN();
    double s = std::numeric_limits<double>::quiet_NaN();
    double vaporMoleFraction = std::numeric_limits<double>::quiet_NaN();
    double vaporMassFraction = std::numeric_limits<double>::quiet_NaN();
    double vaporVolumeFraction = std::numeric_limits<double>::quiet_NaN();
    // Condensed-phase mass fraction relative to the TOTAL mass (gas is
    // 1 - condensedMassFraction).  Only meaningful for CONDENSED_SOLID.
    double condensedMassFraction = std::numeric_limits<double>::quiet_NaN();
    bool entropyAvailable = false;

    // Diagnostic projection metadata.  A projected Q is structurally valid
    // but does not lie on the authoritative conserved rho-e state. Production
    // flux/characteristic consumers must reject it; explicit diagnostic
    // solvers may inspect the displacement.
    bool conservedStateProjected = false;
    double densityProjectionDelta = 0.0;
    double energyProjectionDelta = 0.0;

    EquilibriumAcousticProperties acoustic;
    std::optional<EquilibriumThermodynamicPayload> thermodynamics;
    std::optional<EquilibriumPhaseProperties> liquid;
    std::optional<EquilibriumPhaseProperties> vapor;
    std::optional<EquilibriumPhaseProperties> singlePhase;
    // Condensed (solid) phase: CO2/H2O frozen out of the gas at cold
    // temperatures (e.g. 185 K LOX near-field with H2O/CO2 dilution).
    // overallMassFractions keeps the TOTAL composition (conservation
    // basis); the condensed phase carries only the condensed species.
    std::optional<EquilibriumPhaseProperties> condensed;

    bool valid() const;
};

struct EquilibriumSpeciesBasisTolerances
{
    double compositionAbsolute = 1.0e-10;
    double molecularWeightRelative = 1.0e-10;
};

// Structural validation cannot infer the ordered species molecular weights.
// Model boundaries must additionally prove that every phase uses one
// consistent mass, mole, and molecular-weight basis.
bool validateEquilibriumSpeciesBasis(
    const EquilibriumState& state,
    const std::vector<double>& orderedSpeciesMolecularWeights,
    const EquilibriumSpeciesBasisTolerances& tolerances = {});

// Diagnostics are returned with every model call but are not stored in every
// Domain cell. A failed result never carries an EquilibriumState.
struct EquilibriumDiagnostics
{
    std::string modelId;
    std::string backendStatus;
    std::string backendPhaseStatus;
    std::string backendCaloricStatus;
    std::string backendAcousticStatus;
    double densityResidual = std::numeric_limits<double>::infinity();
    double energyResidual = std::numeric_limits<double>::infinity();
    double enthalpyResidual = std::numeric_limits<double>::infinity();
    double materialResidual = std::numeric_limits<double>::infinity();
    double fugacityResidual = std::numeric_limits<double>::infinity();
    double stabilityResidual = std::numeric_limits<double>::infinity();
    double energyIdentityResidual = std::numeric_limits<double>::infinity();
    double acousticLowerEntropyResidual =
        std::numeric_limits<double>::infinity();
    double acousticUpperEntropyResidual =
        std::numeric_limits<double>::infinity();
    double acousticUsedRelativeDensityPerturbation =
        std::numeric_limits<double>::quiet_NaN();
    double acousticDerivativeRelativeChange =
        std::numeric_limits<double>::infinity();
    std::size_t acousticDerivativeConvergenceComparisons = 0;
    std::size_t iterations = 0;
    std::size_t validSamples = 0;
    std::size_t invalidSamples = 0;
};

struct EquilibriumResult
{
    EquilibriumStatus status = EquilibriumStatus::INVALID_CONFIGURATION;
    std::optional<EquilibriumState> state;
    EquilibriumDiagnostics diagnostics;

    bool success() const
    {
        return status == EquilibriumStatus::SUCCESS
            && state.has_value() && state->valid();
    }
};

struct EquilibriumSeed
{
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double temperature = std::numeric_limits<double>::quiet_NaN();

    // Optional continuation data from the previously committed cell state.
    // These values are seeds only; the recovery path must still validate the
    // candidate against the current rho-e-Y state and fall back to a full scan.
    EquilibriumRegion previousRegion = EquilibriumRegion::UNKNOWN;
    double previousTemperature =
        std::numeric_limits<double>::quiet_NaN();
    double previousVaporMassFraction =
        std::numeric_limits<double>::quiet_NaN();
};

// A configured model is immutable and every recovery method must be re-entrant.
// Rho-e remains the authoritative conservative recovery operation. Concrete
// models may additionally expose PTY initialization only when it returns the
// same complete caloric, phase, acoustic, and species-basis state contract.
class EquilibriumThermoModel
{
public:
    virtual ~EquilibriumThermoModel() = default;

    virtual EquilibriumResult recoverRhoE(
        double density,
        double internalEnergy,
        const std::vector<double>& overallMassFractions,
        const EquilibriumSeed* seed = nullptr) const = 0;

    virtual EquilibriumResult initializePTY(
        double pressure,
        double temperature,
        const std::vector<double>& overallMassFractions) const;

    // Optional boundary flash. Conservative cell updates must still use
    // recoverRhoE. A model without entropy support rejects explicitly.
    virtual EquilibriumResult initializePSY(
        double /*pressure*/, double /*entropy*/,
        const std::vector<double>& /*overallMassFractions*/,
        const EquilibriumSeed* /*seed*/ = nullptr) const
    {
        EquilibriumResult result;
        result.status = EquilibriumStatus::UNSUPPORTED_PROPERTY_RANGE;
        result.diagnostics.backendStatus = "PRESSURE_ENTROPY_FLASH_UNSUPPORTED";
        return result;
    }

    virtual const std::string& modelId() const = 0;
    virtual const std::vector<std::string>& speciesNames() const = 0;
    virtual EnergyConvention energyConvention() const = 0;

    // Optional immutable identity for phase-local PR authorities. Generic
    // models return an empty value; a displacement owner must reject such a
    // model rather than treating modelId() as a substitute data binding.
    virtual std::string phaseLocalPhysicalModelFingerprint() const
    {
        return {};
    }

    // Register the GERG-2008 residual table (dense-cryogenic production
    // backend).  Models that do not consume it ignore the registration.
    virtual void setGERGTable(const GERGResidualTable* /*table*/) const {}
    virtual void setSaturationTable(const GERGSaturationTable* /*table*/) const {}

    // Validate a state against this model's ordered species basis. Concrete
    // production models cache their molecular weights and override this path;
    // the default keeps test/external implementations fail-closed.
    virtual bool stateValidForModel(const EquilibriumState& state) const;

    // Transitional legacy initialization may start from p,T,Y. Phase-aware
    // models must use a complete EquilibriumState instead.
    virtual bool supportsPrimitivePTYInitialization() const { return false; }

    // Model-local declaration for the transitional committed-Q gas transport
    // bridge. This is not a production ModelCapabilities readiness flag.
    virtual bool supportsCommittedLegacyGasTransport() const { return false; }
};

class LegacySinglePhaseThermoModel final : public EquilibriumThermoModel
{
    std::vector<std::string> speciesNames_;
    std::vector<double> speciesMolecularWeights_;
    bool configurationValid_ = false;
    std::string modelId_;
    mutable const GERGResidualTable* gergTable_ = nullptr;
    mutable const GERGSaturationTable* saturationTable_ = nullptr;

public:
    explicit LegacySinglePhaseThermoModel(
        std::vector<std::string> speciesNames = {});

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
    const std::vector<std::string>& speciesNames() const override
    { return speciesNames_; }
    EnergyConvention energyConvention() const override
    {
        return EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT;
    }
    bool stateValidForModel(const EquilibriumState& state) const override;
    bool supportsPrimitivePTYInitialization() const override { return true; }
    bool supportsCommittedLegacyGasTransport() const override
    { return configurationValid_; }

    void setGERGTable(const GERGResidualTable* table) const override
    { gergTable_ = table; }

    void setSaturationTable(const GERGSaturationTable* table) const override
    { saturationTable_ = table; }
};

std::shared_ptr<const EquilibriumThermoModel>
makeLegacySinglePhaseThermoModel(
    const std::vector<std::string>& speciesNames = {});

} // End namespace thermo
} // End namespace Foam

#endif
