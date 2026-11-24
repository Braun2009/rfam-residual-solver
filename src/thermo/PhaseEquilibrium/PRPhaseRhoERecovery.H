#ifndef PRPhaseRhoERecovery_H
#define PRPhaseRhoERecovery_H

#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "BinaryInteractionDB.H"
#include "PengRobinsonEOS.H"
#include "PRMixtureKernel.H"
#include "PRPhaseCaloric.H"
#include "SpeciesCriticalDB.H"

namespace Foam
{
namespace phaseEquilibrium
{

enum class PRPhaseIdentity
{
    VAPOR,
    SUPERCRITICAL
};

enum class PRPhaseRootIdentity
{
    PR_VAPOR_ROOT,
    PR_UNIQUE_SUPERCRITICAL_ROOT
};

enum class PRPhaseRhoEStatus
{
    SUCCESS,
    INVALID_CONFIGURATION,
    INVALID_PHASE_IDENTITY,
    INVALID_DENSITY,
    INVALID_ENERGY,
    INVALID_TEMPERATURE_SEED,
    INVALID_COMPOSITION,
    UNSUPPORTED_CALORIC_RANGE,
    NONPOSITIVE_PRESSURE,
    ROOT_LOSS,
    ROOT_IDENTITY_MISMATCH,
    NOT_BRACKETED,
    MULTIPLE_ENERGY_ROOTS,
    DERIVATIVE_FAILURE,
    NOT_CONVERGED,
    NONFINITE_PROPERTY,
    RESIDUAL_CONTRACT_FAILED
};

const char* phaseRhoEStatusName(PRPhaseRhoEStatus status);

struct PRPhaseRhoEOptions
{
    double minimumTemperature = 250.0;
    double maximumTemperature = 1500.0;
    int temperatureScanIntervals = 128;
    int maximumIterations = 80;
    double derivativeRelativeStep = 2.0e-7;
    double energyAbsoluteTolerance = 1.0e-5;
    double energyRelativeTolerance = 2.0e-10;
    double densityRelativeTolerance = 2.0e-9;
    double rootRelativeTolerance = 2.0e-12;
    double derivativeIdentityRelativeTolerance = 2.0e-6;
    double derivativeIdentityAbsoluteTolerance = 1.0e-4;
};

struct PRPhaseRhoEDerivatives
{
    double cv = std::numeric_limits<double>::quiet_NaN();
    double dedT_rho = std::numeric_limits<double>::quiet_NaN();
    double dpdT_rho = std::numeric_limits<double>::quiet_NaN();
    double dpdrho_T = std::numeric_limits<double>::quiet_NaN();
    double dedrho_T = std::numeric_limits<double>::quiet_NaN();

    bool valid() const;
};

struct PRPhaseRhoEDiagnostics
{
    int iterations = 0;
    int validScanSamples = 0;
    int invalidScanSamples = 0;
    double densityRelativeResidual =
        std::numeric_limits<double>::infinity();
    double energyRelativeResidual =
        std::numeric_limits<double>::infinity();
    double eosRootRelativeResidual =
        std::numeric_limits<double>::infinity();
    double departureCaloricRelativeResidual =
        std::numeric_limits<double>::infinity();
    double cvIdentityRelativeResidual =
        std::numeric_limits<double>::infinity();
    double energyDensityIdentityRelativeResidual =
        std::numeric_limits<double>::infinity();
};

struct PRPhaseRhoEResult
{
    PRPhaseRhoEStatus status = PRPhaseRhoEStatus::INVALID_CONFIGURATION;
    PRPhaseIdentity phaseIdentity = PRPhaseIdentity::VAPOR;
    PRPhaseRootIdentity rootIdentity =
        PRPhaseRootIdentity::PR_VAPOR_ROOT;
    std::optional<PRPhaseCaloricProperties> state;
    std::optional<PRPhaseRhoEDerivatives> derivatives;
    PRPhaseRhoEDiagnostics diagnostics;

    bool success() const
    {
        return status == PRPhaseRhoEStatus::SUCCESS
            && state.has_value() && derivatives.has_value();
    }
};

// Non-owning thermodynamic kernel for a parent-bound vapor or supercritical
// PR branch. It recovers p and T from phase-local rho-e-Y at fixed density,
// but it does not own parent generations, tokens, quality, phase fractions,
// acoustics, final reflash, Domain publication, or solver admission.
class PRPhaseRhoERecovery
{
private:
    struct FixedDensityEvaluation
    {
        PRPhaseRhoEStatus status = PRPhaseRhoEStatus::ROOT_LOSS;
        std::optional<PRPhaseCaloricProperties> state;
    };

    PRMixtureKernel kernel_;
    PRPhaseCaloric caloric_;
    PengRobinsonEOS eos_;
    std::vector<std::string> speciesNames_;
    PRPhaseRhoEOptions options_;
    std::string physicalModelAuthorityFingerprint_;
    std::string authorityFingerprint_;
    PRPhaseRhoEStatus configurationStatus_;

    static bool compatibleIdentity(
        PRPhaseIdentity phaseIdentity,
        PRPhaseRootIdentity rootIdentity);

    FixedDensityEvaluation evaluateFixedDensity(
        double density,
        double temperature,
        const std::vector<double>& massFractions,
        PRPhaseIdentity phaseIdentity,
        PRPhaseRootIdentity rootIdentity) const;

    bool idealCpMass(
        double temperature,
        const std::vector<double>& massFractions,
        double& cp) const;

    std::optional<PRPhaseRhoEDerivatives> derivatives(
        const PRPhaseCaloricProperties& base,
        PRPhaseIdentity phaseIdentity,
        PRPhaseRootIdentity rootIdentity,
        PRPhaseRhoEDiagnostics& diagnostics) const;

    PRPhaseRhoEResult finalize(
        const PRPhaseCaloricProperties& state,
        double targetDensity,
        double targetEnergy,
        PRPhaseIdentity phaseIdentity,
        PRPhaseRootIdentity rootIdentity,
        PRPhaseRhoEDiagnostics diagnostics) const;

public:
    PRPhaseRhoERecovery(
        const SpeciesCriticalDB& criticalData,
        const BinaryInteractionDB& kij,
        const std::vector<std::string>& speciesNames,
        const PRPhaseRhoEOptions& options = {});

    bool configured() const
    {
        return configurationStatus_ == PRPhaseRhoEStatus::SUCCESS;
    }

    PRPhaseRhoEStatus configurationStatus() const
    {
        return configurationStatus_;
    }

    const std::vector<std::string>& speciesNames() const
    {
        return speciesNames_;
    }

    const PRPhaseRhoEOptions& options() const
    {
        return options_;
    }

    // These identities are immutable construction authorities. They let a
    // token-bound displacement layer prove that it uses the same frozen PR
    // data and phase-local recovery options as this kernel.
    const std::string& physicalModelAuthorityFingerprint() const
    {
        return physicalModelAuthorityFingerprint_;
    }

    const std::string& authorityFingerprint() const
    {
        return authorityFingerprint_;
    }

    PRPhaseRhoEResult recover(
        double density,
        double specificInternalEnergy,
        const std::vector<double>& massFractions,
        PRPhaseIdentity phaseIdentity,
        PRPhaseRootIdentity rootIdentity,
        double temperatureSeed) const;
};

} // End namespace phaseEquilibrium
} // End namespace Foam

#endif
