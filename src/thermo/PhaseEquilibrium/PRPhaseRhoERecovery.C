#include "PRPhaseRhoERecovery.H"

#include "PRPhaseAuthorityIdentity.H"

#include <algorithm>
#include <cmath>
#include <utility>

namespace Foam
{
namespace phaseEquilibrium
{

namespace
{

constexpr int maximumTemperatureScanIntervals = 4096;
constexpr int maximumRecoveryIterations = 512;
constexpr double maximumAbsoluteTolerance = 1.0;
constexpr double maximumRelativeTolerance = 1.0e-2;

bool finitePositive(double value)
{
    return std::isfinite(value) && value > 0.0;
}

double relativeResidual(double lhs, double rhs)
{
    return std::abs(lhs - rhs)
        /std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

bool withinTolerance(
    double lhs,
    double rhs,
    double relativeTolerance,
    double absoluteTolerance)
{
    return std::isfinite(lhs) && std::isfinite(rhs)
        && std::abs(lhs - rhs) <= absoluteTolerance
            + relativeTolerance
                *std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

bool validOptions(const PRPhaseRhoEOptions& options)
{
    return finitePositive(options.minimumTemperature)
        && std::isfinite(options.maximumTemperature)
        && options.maximumTemperature > options.minimumTemperature
        && options.temperatureScanIntervals >= 2
        && options.temperatureScanIntervals
            <= maximumTemperatureScanIntervals
        && options.maximumIterations >= 1
        && options.maximumIterations <= maximumRecoveryIterations
        && finitePositive(options.derivativeRelativeStep)
        && options.derivativeRelativeStep < 0.25
        && std::isfinite(options.energyAbsoluteTolerance)
        && options.energyAbsoluteTolerance >= 0.0
        && options.energyAbsoluteTolerance <= maximumAbsoluteTolerance
        && std::isfinite(options.energyRelativeTolerance)
        && options.energyRelativeTolerance >= 0.0
        && options.energyRelativeTolerance <= maximumRelativeTolerance
        && options.energyAbsoluteTolerance
            + options.energyRelativeTolerance > 0.0
        && finitePositive(options.densityRelativeTolerance)
        && options.densityRelativeTolerance <= maximumRelativeTolerance
        && finitePositive(options.rootRelativeTolerance)
        && options.rootRelativeTolerance <= maximumRelativeTolerance
        && std::isfinite(options.derivativeIdentityRelativeTolerance)
        && options.derivativeIdentityRelativeTolerance >= 0.0
        && options.derivativeIdentityRelativeTolerance
            <= maximumRelativeTolerance
        && std::isfinite(options.derivativeIdentityAbsoluteTolerance)
        && options.derivativeIdentityAbsoluteTolerance >= 0.0
        && options.derivativeIdentityAbsoluteTolerance
            <= maximumAbsoluteTolerance
        && options.derivativeIdentityRelativeTolerance
            + options.derivativeIdentityAbsoluteTolerance > 0.0;
}

PRPhaseRhoEStatus caloricFailureStatus(PRPhaseCaloricStatus status)
{
    switch (status)
    {
        case PRPhaseCaloricStatus::INVALID_COMPOSITION:
            return PRPhaseRhoEStatus::INVALID_COMPOSITION;
        case PRPhaseCaloricStatus::UNSUPPORTED_CALORIC_RANGE:
            return PRPhaseRhoEStatus::UNSUPPORTED_CALORIC_RANGE;
        case PRPhaseCaloricStatus::NONFINITE_PROPERTY:
            return PRPhaseRhoEStatus::NONFINITE_PROPERTY;
        case PRPhaseCaloricStatus::EOS_CLOSURE_FAILURE:
            return PRPhaseRhoEStatus::RESIDUAL_CONTRACT_FAILED;
        case PRPhaseCaloricStatus::INVALID_PRESSURE:
            return PRPhaseRhoEStatus::NONPOSITIVE_PRESSURE;
        case PRPhaseCaloricStatus::PHASE_ROOT_FAILURE:
        case PRPhaseCaloricStatus::MISSING_BINARY_INTERACTION:
            return PRPhaseRhoEStatus::ROOT_LOSS;
        case PRPhaseCaloricStatus::SUCCESS:
        case PRPhaseCaloricStatus::INVALID_CONFIGURATION:
        case PRPhaseCaloricStatus::MISSING_CALORIC_DATA:
        case PRPhaseCaloricStatus::MOLECULAR_WEIGHT_MISMATCH:
        case PRPhaseCaloricStatus::INVALID_TEMPERATURE:
        case PRPhaseCaloricStatus::INVALID_PHASE_ROOT:
            return PRPhaseRhoEStatus::NONFINITE_PROPERTY;
    }
    return PRPhaseRhoEStatus::NONFINITE_PROPERTY;
}

} // End anonymous namespace

const char* phaseRhoEStatusName(PRPhaseRhoEStatus status)
{
    switch (status)
    {
        case PRPhaseRhoEStatus::SUCCESS: return "SUCCESS";
        case PRPhaseRhoEStatus::INVALID_CONFIGURATION:
            return "INVALID_CONFIGURATION";
        case PRPhaseRhoEStatus::INVALID_PHASE_IDENTITY:
            return "INVALID_PHASE_IDENTITY";
        case PRPhaseRhoEStatus::INVALID_DENSITY: return "INVALID_DENSITY";
        case PRPhaseRhoEStatus::INVALID_ENERGY: return "INVALID_ENERGY";
        case PRPhaseRhoEStatus::INVALID_TEMPERATURE_SEED:
            return "INVALID_TEMPERATURE_SEED";
        case PRPhaseRhoEStatus::INVALID_COMPOSITION:
            return "INVALID_COMPOSITION";
        case PRPhaseRhoEStatus::UNSUPPORTED_CALORIC_RANGE:
            return "UNSUPPORTED_CALORIC_RANGE";
        case PRPhaseRhoEStatus::NONPOSITIVE_PRESSURE:
            return "NONPOSITIVE_PRESSURE";
        case PRPhaseRhoEStatus::ROOT_LOSS: return "ROOT_LOSS";
        case PRPhaseRhoEStatus::ROOT_IDENTITY_MISMATCH:
            return "ROOT_IDENTITY_MISMATCH";
        case PRPhaseRhoEStatus::NOT_BRACKETED: return "NOT_BRACKETED";
        case PRPhaseRhoEStatus::MULTIPLE_ENERGY_ROOTS:
            return "MULTIPLE_ENERGY_ROOTS";
        case PRPhaseRhoEStatus::DERIVATIVE_FAILURE:
            return "DERIVATIVE_FAILURE";
        case PRPhaseRhoEStatus::NOT_CONVERGED: return "NOT_CONVERGED";
        case PRPhaseRhoEStatus::NONFINITE_PROPERTY:
            return "NONFINITE_PROPERTY";
        case PRPhaseRhoEStatus::RESIDUAL_CONTRACT_FAILED:
            return "RESIDUAL_CONTRACT_FAILED";
    }
    return "UNKNOWN_STATUS";
}

bool PRPhaseRhoEDerivatives::valid() const
{
    return finitePositive(cv) && finitePositive(dedT_rho)
        && finitePositive(dpdrho_T) && std::isfinite(dpdT_rho)
        && std::isfinite(dedrho_T);
}

bool PRPhaseRhoERecovery::compatibleIdentity(
    PRPhaseIdentity phaseIdentity,
    PRPhaseRootIdentity rootIdentity)
{
    return (phaseIdentity == PRPhaseIdentity::VAPOR
            && rootIdentity == PRPhaseRootIdentity::PR_VAPOR_ROOT)
        || (phaseIdentity == PRPhaseIdentity::SUPERCRITICAL
            && rootIdentity
                == PRPhaseRootIdentity::PR_UNIQUE_SUPERCRITICAL_ROOT);
}

PRPhaseRhoERecovery::PRPhaseRhoERecovery(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    const PRPhaseRhoEOptions& options)
:
    kernel_(criticalData, kij, speciesNames),
    caloric_(criticalData, kij, speciesNames),
    eos_(criticalData, kij),
    speciesNames_(speciesNames),
    options_(options),
    physicalModelAuthorityFingerprint_(
        makePRPhasePhysicalModelAuthorityFingerprint(
            criticalData, kij, speciesNames)),
    authorityFingerprint_(makePRPhaseRhoEAuthorityFingerprint(
        criticalData, kij, speciesNames, options)),
    configurationStatus_(PRPhaseRhoEStatus::INVALID_CONFIGURATION)
{
    if (!kernel_.configured() || !caloric_.configured()
        || !validOptions(options_)
        || !validPRPhaseAuthorityFingerprint(
            physicalModelAuthorityFingerprint_)
        || !validPRPhaseAuthorityFingerprint(authorityFingerprint_))
    {
        return;
    }

    // A configured phase kernel may defer missing-pair rejection until a pair
    // becomes active. This recovery object freezes the whole potential set.
    for (std::size_t i = 0; i < speciesNames_.size(); ++i)
    {
        for (std::size_t j = i + 1; j < speciesNames_.size(); ++j)
        {
            if (!kij.has(speciesNames_[i], speciesNames_[j])) return;
        }
    }

    try
    {
        eos_.setSpeciesNames(speciesNames_);
    }
    catch (...)
    {
        return;
    }
    configurationStatus_ = PRPhaseRhoEStatus::SUCCESS;
}

PRPhaseRhoERecovery::FixedDensityEvaluation
PRPhaseRhoERecovery::evaluateFixedDensity(
    double density,
    double temperature,
    const std::vector<double>& massFractions,
    PRPhaseIdentity phaseIdentity,
    PRPhaseRootIdentity rootIdentity) const
{
    FixedDensityEvaluation result;
    if (!finitePositive(density))
    {
        result.status = PRPhaseRhoEStatus::INVALID_DENSITY;
        return result;
    }
    if (!finitePositive(temperature))
    {
        result.status = PRPhaseRhoEStatus::INVALID_TEMPERATURE_SEED;
        return result;
    }
    if (!compatibleIdentity(phaseIdentity, rootIdentity))
    {
        result.status = PRPhaseRhoEStatus::INVALID_PHASE_IDENTITY;
        return result;
    }

    std::vector<double> moleFractions;
    if (kernel_.massToMoleFractions(massFractions, moleFractions)
        != PhaseEquilibriumStatus::SUCCESS)
    {
        result.status = PRPhaseRhoEStatus::INVALID_COMPOSITION;
        return result;
    }

    const double pressure = eos_.pressure(
        density, temperature, massFractions);
    if (!finitePositive(pressure))
    {
        result.status = PRPhaseRhoEStatus::NONPOSITIVE_PRESSURE;
        return result;
    }

    const PRRootSet roots = kernel_.roots(
        pressure, temperature, moleFractions);
    if (!roots.valid() || roots.physicalRoots.empty())
    {
        result.status = PRPhaseRhoEStatus::ROOT_LOSS;
        return result;
    }
    if (rootIdentity == PRPhaseRootIdentity::PR_UNIQUE_SUPERCRITICAL_ROOT
        && roots.physicalRoots.size() != 1U)
    {
        result.status = PRPhaseRhoEStatus::ROOT_IDENTITY_MISMATCH;
        return result;
    }

    const double retainedRoot = roots.physicalRoots.back();
    const PRPhaseCaloricResult evaluated = caloric_.evaluate(
        pressure, temperature, moleFractions, PhaseRoot::VAPOR);
    if (!evaluated.success())
    {
        result.status = caloricFailureStatus(evaluated.status);
        return result;
    }
    if (evaluated.state->root != PhaseRoot::VAPOR
        || relativeResidual(
                evaluated.state->compressibility, retainedRoot)
            > options_.rootRelativeTolerance
        || relativeResidual(evaluated.state->rho, density)
            > options_.densityRelativeTolerance)
    {
        result.status = PRPhaseRhoEStatus::ROOT_IDENTITY_MISMATCH;
        return result;
    }

    result.status = PRPhaseRhoEStatus::SUCCESS;
    result.state = *evaluated.state;
    return result;
}

bool PRPhaseRhoERecovery::idealCpMass(
    double temperature,
    const std::vector<double>& massFractions,
    double& cp) const
{
    const std::optional<double> corrected = caloric_.idealCpMass(
        temperature, massFractions);
    cp = corrected.value_or(std::numeric_limits<double>::quiet_NaN());
    return corrected.has_value();
}

std::optional<PRPhaseRhoEDerivatives> PRPhaseRhoERecovery::derivatives(
    const PRPhaseCaloricProperties& base,
    PRPhaseIdentity phaseIdentity,
    PRPhaseRootIdentity rootIdentity,
    PRPhaseRhoEDiagnostics& diagnostics) const
{
    const double densityStep = options_.derivativeRelativeStep*base.rho;
    const double temperatureStep =
        options_.derivativeRelativeStep*base.temperature;
    if (!finitePositive(densityStep) || !finitePositive(temperatureStep)
        || base.rho <= densityStep || base.temperature <= temperatureStep)
    {
        return std::nullopt;
    }

    const FixedDensityEvaluation rhoMinus = evaluateFixedDensity(
        base.rho - densityStep, base.temperature, base.massFractions,
        phaseIdentity, rootIdentity);
    const FixedDensityEvaluation rhoPlus = evaluateFixedDensity(
        base.rho + densityStep, base.temperature, base.massFractions,
        phaseIdentity, rootIdentity);
    const FixedDensityEvaluation temperatureMinus = evaluateFixedDensity(
        base.rho, base.temperature - temperatureStep, base.massFractions,
        phaseIdentity, rootIdentity);
    const FixedDensityEvaluation temperaturePlus = evaluateFixedDensity(
        base.rho, base.temperature + temperatureStep, base.massFractions,
        phaseIdentity, rootIdentity);
    if (!rhoMinus.state || !rhoPlus.state
        || !temperatureMinus.state || !temperaturePlus.state)
    {
        return std::nullopt;
    }

    PRPhaseRhoEDerivatives result;
    result.dpdrho_T =
        (rhoPlus.state->pressure - rhoMinus.state->pressure)
        /(2.0*densityStep);
    result.dedrho_T = (rhoPlus.state->e - rhoMinus.state->e)
        /(2.0*densityStep);
    result.dpdT_rho =
        (temperaturePlus.state->pressure - temperatureMinus.state->pressure)
        /(2.0*temperatureStep);
    result.dedT_rho =
        (temperaturePlus.state->e - temperatureMinus.state->e)
        /(2.0*temperatureStep);

    double idealCp = std::numeric_limits<double>::quiet_NaN();
    if (!idealCpMass(base.temperature, base.massFractions, idealCp))
        return std::nullopt;
    const double molecularWeight = eos_.MW(base.massFractions);
    const double departureCv = eos_.cvDeparture(
        base.pressure, base.temperature, base.compressibility,
        base.massFractions);
    if (!finitePositive(molecularWeight) || !std::isfinite(departureCv))
        return std::nullopt;

    const double idealCv = idealCp
        - PRMixtureKernel::universalGasConstant/molecularWeight;
    result.cv = idealCv + departureCv;
    if (!result.valid()) return std::nullopt;

    diagnostics.departureCaloricRelativeResidual = relativeResidual(
        result.dedT_rho - idealCv, departureCv);
    diagnostics.cvIdentityRelativeResidual = relativeResidual(
        result.dedT_rho, result.cv);
    const double expectedEnergyDensityDerivative =
        (base.pressure - base.temperature*result.dpdT_rho)
        /(base.rho*base.rho);
    diagnostics.energyDensityIdentityRelativeResidual = relativeResidual(
        result.dedrho_T, expectedEnergyDensityDerivative);

    return result;
}

PRPhaseRhoEResult PRPhaseRhoERecovery::finalize(
    const PRPhaseCaloricProperties& state,
    double targetDensity,
    double targetEnergy,
    PRPhaseIdentity phaseIdentity,
    PRPhaseRootIdentity rootIdentity,
    PRPhaseRhoEDiagnostics diagnostics) const
{
    PRPhaseRhoEResult result;
    result.phaseIdentity = phaseIdentity;
    result.rootIdentity = rootIdentity;
    result.diagnostics = diagnostics;

    result.diagnostics.densityRelativeResidual = relativeResidual(
        state.rho, targetDensity);
    result.diagnostics.energyRelativeResidual = relativeResidual(
        state.e, targetEnergy);

    const PRRootSet roots = kernel_.roots(
        state.pressure, state.temperature, state.moleFractions);
    if (!roots.valid() || roots.physicalRoots.empty())
    {
        result.status = PRPhaseRhoEStatus::ROOT_LOSS;
        return result;
    }
    if (rootIdentity == PRPhaseRootIdentity::PR_UNIQUE_SUPERCRITICAL_ROOT
        && roots.physicalRoots.size() != 1U)
    {
        result.status = PRPhaseRhoEStatus::ROOT_IDENTITY_MISMATCH;
        return result;
    }
    result.diagnostics.eosRootRelativeResidual = relativeResidual(
        state.compressibility, roots.physicalRoots.back());

    const std::optional<PRPhaseRhoEDerivatives> phaseDerivatives =
        derivatives(state, phaseIdentity, rootIdentity, result.diagnostics);
    if (!phaseDerivatives)
    {
        result.status = PRPhaseRhoEStatus::DERIVATIVE_FAILURE;
        return result;
    }

    double idealCp = std::numeric_limits<double>::quiet_NaN();
    const double molecularWeight = eos_.MW(state.massFractions);
    const double departureCv = eos_.cvDeparture(
        state.pressure, state.temperature, state.compressibility,
        state.massFractions);
    if (!idealCpMass(state.temperature, state.massFractions, idealCp)
        || !finitePositive(molecularWeight) || !std::isfinite(departureCv))
    {
        result.status = PRPhaseRhoEStatus::DERIVATIVE_FAILURE;
        return result;
    }
    const double idealCv = idealCp
        - PRMixtureKernel::universalGasConstant/molecularWeight;
    const double energyDensityIdentity =
        (state.pressure
            - state.temperature*phaseDerivatives->dpdT_rho)
        /(state.rho*state.rho);

    const double energyTolerance = options_.energyAbsoluteTolerance
        + options_.energyRelativeTolerance
            *std::max(std::abs(targetEnergy), 1.0);
    if (std::abs(state.e - targetEnergy) > energyTolerance
        || result.diagnostics.densityRelativeResidual
            > options_.densityRelativeTolerance
        || result.diagnostics.eosRootRelativeResidual
            > options_.rootRelativeTolerance
        || !withinTolerance(
            phaseDerivatives->dedT_rho,
            phaseDerivatives->cv,
            options_.derivativeIdentityRelativeTolerance,
            options_.derivativeIdentityAbsoluteTolerance)
        || !withinTolerance(
            phaseDerivatives->dedT_rho - idealCv,
            departureCv,
            options_.derivativeIdentityRelativeTolerance,
            options_.derivativeIdentityAbsoluteTolerance)
        || !withinTolerance(
            phaseDerivatives->dedrho_T,
            energyDensityIdentity,
            options_.derivativeIdentityRelativeTolerance,
            options_.derivativeIdentityAbsoluteTolerance))
    {
        result.status = PRPhaseRhoEStatus::RESIDUAL_CONTRACT_FAILED;
        return result;
    }

    result.status = PRPhaseRhoEStatus::SUCCESS;
    result.state = state;
    result.derivatives = *phaseDerivatives;
    return result;
}

PRPhaseRhoEResult PRPhaseRhoERecovery::recover(
    double density,
    double specificInternalEnergy,
    const std::vector<double>& massFractions,
    PRPhaseIdentity phaseIdentity,
    PRPhaseRootIdentity rootIdentity,
    double temperatureSeed) const
{
    PRPhaseRhoEResult result;
    result.phaseIdentity = phaseIdentity;
    result.rootIdentity = rootIdentity;
    result.status = configurationStatus_;
    if (!configured()) return result;
    if (!compatibleIdentity(phaseIdentity, rootIdentity))
    {
        result.status = PRPhaseRhoEStatus::INVALID_PHASE_IDENTITY;
        return result;
    }
    if (!finitePositive(density))
    {
        result.status = PRPhaseRhoEStatus::INVALID_DENSITY;
        return result;
    }
    if (!std::isfinite(specificInternalEnergy))
    {
        result.status = PRPhaseRhoEStatus::INVALID_ENERGY;
        return result;
    }
    if (!finitePositive(temperatureSeed))
    {
        result.status = PRPhaseRhoEStatus::INVALID_TEMPERATURE_SEED;
        return result;
    }

    std::vector<double> moleFractions;
    if (kernel_.massToMoleFractions(massFractions, moleFractions)
        != PhaseEquilibriumStatus::SUCCESS)
    {
        result.status = PRPhaseRhoEStatus::INVALID_COMPOSITION;
        return result;
    }
    const PRPhaseCaloricTemperatureRange supported =
        caloric_.supportedTemperatureRange(moleFractions);
    if (!supported.success())
    {
        result.status = PRPhaseRhoEStatus::UNSUPPORTED_CALORIC_RANGE;
        return result;
    }

    const double minimumTemperature = std::max(
        options_.minimumTemperature, supported.minimum);
    const double maximumTemperature = std::min(
        options_.maximumTemperature, supported.maximum);
    if (!(maximumTemperature > minimumTemperature))
    {
        result.status = PRPhaseRhoEStatus::UNSUPPORTED_CALORIC_RANGE;
        return result;
    }

    struct Bracket
    {
        FixedDensityEvaluation low;
        FixedDensityEvaluation high;
    };
    std::vector<Bracket> brackets;
    std::vector<FixedDensityEvaluation> exactRoots;
    std::optional<FixedDensityEvaluation> previous;
    bool insideExactRootBand = false;
    const double energyTolerance = options_.energyAbsoluteTolerance
        + options_.energyRelativeTolerance
            *std::max(std::abs(specificInternalEnergy), 1.0);

    for (int i = 0; i <= options_.temperatureScanIntervals; ++i)
    {
        const double fraction = static_cast<double>(i)
            /static_cast<double>(options_.temperatureScanIntervals);
        const double temperature = minimumTemperature
            + fraction*(maximumTemperature - minimumTemperature);
        FixedDensityEvaluation current = evaluateFixedDensity(
            density, temperature, massFractions, phaseIdentity, rootIdentity);
        if (!current.state)
        {
            ++result.diagnostics.invalidScanSamples;
            previous.reset();
            insideExactRootBand = false;
            continue;
        }

        ++result.diagnostics.validScanSamples;
        const double residual = current.state->e - specificInternalEnergy;
        if (std::abs(residual) <= energyTolerance)
        {
            if (!insideExactRootBand)
            {
                exactRoots.push_back(current);
            }
            else
            {
                const double retainedResidual =
                    exactRoots.back().state->e - specificInternalEnergy;
                if (std::abs(residual) < std::abs(retainedResidual))
                    exactRoots.back() = current;
            }
            insideExactRootBand = true;
            previous.reset();
            continue;
        }
        insideExactRootBand = false;
        if (previous && previous->state)
        {
            const double previousResidual =
                previous->state->e - specificInternalEnergy;
            if ((previousResidual < 0.0 && residual > 0.0)
                || (previousResidual > 0.0 && residual < 0.0))
            {
                brackets.push_back({*previous, current});
            }
        }
        previous = std::move(current);
    }

    const std::size_t rootCount = exactRoots.size() + brackets.size();
    if (rootCount == 0U)
    {
        result.status = PRPhaseRhoEStatus::NOT_BRACKETED;
        return result;
    }
    if (rootCount != 1U)
    {
        result.status = PRPhaseRhoEStatus::MULTIPLE_ENERGY_ROOTS;
        return result;
    }
    if (!exactRoots.empty())
    {
        return finalize(
            *exactRoots.front().state, density, specificInternalEnergy,
            phaseIdentity, rootIdentity, result.diagnostics);
    }

    FixedDensityEvaluation low = brackets.front().low;
    FixedDensityEvaluation high = brackets.front().high;
    double lowTemperature = low.state->temperature;
    double highTemperature = high.state->temperature;
    double lowResidual = low.state->e - specificInternalEnergy;
    double highResidual = high.state->e - specificInternalEnergy;
    if (lowResidual > 0.0)
    {
        std::swap(lowTemperature, highTemperature);
        std::swap(lowResidual, highResidual);
    }
    if (!(lowResidual < 0.0 && highResidual > 0.0))
    {
        result.status = PRPhaseRhoEStatus::NOT_BRACKETED;
        return result;
    }

    double trialTemperature = temperatureSeed;
    if (!(trialTemperature > lowTemperature
        && trialTemperature < highTemperature))
    {
        trialTemperature = 0.5*(lowTemperature + highTemperature);
    }

    for (int iteration = 0; iteration < options_.maximumIterations;
        ++iteration)
    {
        const FixedDensityEvaluation trial = evaluateFixedDensity(
            density, trialTemperature, massFractions,
            phaseIdentity, rootIdentity);
        result.diagnostics.iterations = iteration + 1;
        if (!trial.state)
        {
            result.status = trial.status;
            return result;
        }

        const double residual = trial.state->e - specificInternalEnergy;
        if (std::abs(residual) <= energyTolerance)
        {
            return finalize(
                *trial.state, density, specificInternalEnergy,
                phaseIdentity, rootIdentity, result.diagnostics);
        }
        if (residual < 0.0)
        {
            lowTemperature = trialTemperature;
            lowResidual = residual;
        }
        else
        {
            highTemperature = trialTemperature;
            highResidual = residual;
        }

        const double span = highTemperature - lowTemperature;
        if (!(span > 0.0))
        {
            result.status = PRPhaseRhoEStatus::NOT_CONVERGED;
            return result;
        }

        double nextTemperature = 0.5*(lowTemperature + highTemperature);
        const double derivativeStep = options_.derivativeRelativeStep
            *std::max(std::abs(trialTemperature), 1.0);
        const double derivativeLowTemperature = std::max(
            lowTemperature, trialTemperature - derivativeStep);
        const double derivativeHighTemperature = std::min(
            highTemperature, trialTemperature + derivativeStep);
        if (derivativeHighTemperature > derivativeLowTemperature)
        {
            const FixedDensityEvaluation derivativeLow = evaluateFixedDensity(
                density, derivativeLowTemperature, massFractions,
                phaseIdentity, rootIdentity);
            const FixedDensityEvaluation derivativeHigh = evaluateFixedDensity(
                density, derivativeHighTemperature, massFractions,
                phaseIdentity, rootIdentity);
            if (derivativeLow.state && derivativeHigh.state)
            {
                const double derivative =
                    (derivativeHigh.state->e - derivativeLow.state->e)
                    /(derivativeHighTemperature
                        - derivativeLowTemperature);
                const double newtonTemperature =
                    trialTemperature - residual/derivative;
                const double guard = 0.1*span;
                if (finitePositive(derivative)
                    && newtonTemperature > lowTemperature + guard
                    && newtonTemperature < highTemperature - guard)
                {
                    nextTemperature = newtonTemperature;
                }
            }
        }
        trialTemperature = nextTemperature;
    }

    result.status = PRPhaseRhoEStatus::NOT_CONVERGED;
    return result;
}

} // End namespace phaseEquilibrium
} // End namespace Foam
