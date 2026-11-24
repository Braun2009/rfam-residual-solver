#include "MixingRules/BinaryInteractionDB.H"
#include "MixingRules/SpeciesCriticalDB.H"
#include "PhaseEquilibrium/PREquilibriumThermoModel.H"
#include "PhaseEquilibrium/PRMixtureKernel.H"
#include "PhaseEquilibrium/PRPhaseCaloric.H"
#include "PhaseEquilibrium/PRTPFlash.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

using Foam::BinaryInteractionDB;
using Foam::SpeciesCriticalDB;
using Foam::phaseEquilibrium::PhaseEquilibriumStatus;
using Foam::phaseEquilibrium::PhaseRoot;
using Foam::phaseEquilibrium::PRMixtureKernel;
using Foam::phaseEquilibrium::PRPhaseCaloric;
using Foam::phaseEquilibrium::PRPhaseCaloricProperties;
using Foam::phaseEquilibrium::PRPhaseCaloricResult;
using Foam::phaseEquilibrium::PRTPFlash;
using Foam::phaseEquilibrium::RhoEFlashOptions;
using Foam::phaseEquilibrium::TPFlashResult;
using Foam::thermo::EquilibriumRegion;
using Foam::thermo::EquilibriumResult;
using Foam::thermo::PREquilibriumThermoModel;

int failures = 0;

void check(bool condition, const std::string& label)
{
    std::cout << (condition ? "  PASS " : "  FAIL ") << label << '\n';
    if (!condition) ++failures;
}

bool finitePositive(double value)
{
    return std::isfinite(value) && value > 0.0;
}

bool close(double lhs, double rhs, double relative, double absolute = 0.0)
{
    return std::abs(lhs - rhs)
        <= absolute + relative*std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

bool closeVector(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs,
    double tolerance)
{
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        if (!close(lhs[i], rhs[i], 0.0, tolerance)) return false;
    }
    return true;
}

RhoEFlashOptions fixtureOptions()
{
    RhoEFlashOptions options;
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

std::vector<double> twoPhaseComposition(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& interactions)
{
    PRTPFlash flash(critical, interactions, {"O2", "C3H8"});
    const TPFlashResult split =
        flash.flashMoleFractions(2.0e6, 300.0, {0.2, 0.8});
    if (!split.valid()
        || split.region != Foam::phaseEquilibrium::EquilibriumRegion::TWO_PHASE)
    {
        throw std::runtime_error("two-phase fixture flash failed");
    }

    return split.feedMassFractions;
}

struct Evaluation
{
    bool valid = false;
    PRPhaseCaloricProperties state;
};

Evaluation evaluateVapor(
    const PRMixtureKernel& kernel,
    const PRPhaseCaloric& caloric,
    double pressure,
    double temperature,
    const std::vector<double>& massFractions)
{
    Evaluation result;
    std::vector<double> moleFractions;
    if (kernel.massToMoleFractions(massFractions, moleFractions)
        != PhaseEquilibriumStatus::SUCCESS)
    {
        return result;
    }
    const PRPhaseCaloricResult evaluated = caloric.evaluate(
        pressure, temperature, moleFractions, PhaseRoot::VAPOR);
    if (!evaluated.success() || evaluated.state->root != PhaseRoot::VAPOR)
        return result;
    result.valid = true;
    result.state = *evaluated.state;
    return result;
}

struct RecoveryResult
{
    bool success = false;
    std::optional<PRPhaseCaloricProperties> state;
    int iterations = 0;
    double densityResidual = std::numeric_limits<double>::infinity();
    double energyResidual = std::numeric_limits<double>::infinity();
};

std::array<double, 2> residuals(
    const PRPhaseCaloricProperties& state,
    double targetDensity,
    double targetEnergy)
{
    const double energyScale = std::max(std::abs(targetEnergy), 1.0e5);
    return {
        std::log(state.rho/targetDensity),
        (state.e - targetEnergy)/energyScale};
}

double normSquared(const std::array<double, 2>& residual)
{
    return residual[0]*residual[0] + residual[1]*residual[1];
}

RecoveryResult recoverVaporRhoE(
    const PRMixtureKernel& kernel,
    const PRPhaseCaloric& caloric,
    double targetDensity,
    double targetEnergy,
    const std::vector<double>& massFractions,
    double pressureSeed,
    double temperatureSeed)
{
    RecoveryResult result;
    if (!finitePositive(targetDensity) || !std::isfinite(targetEnergy)
        || !finitePositive(pressureSeed) || !finitePositive(temperatureSeed))
    {
        return result;
    }

    constexpr double minimumPressure = 1.0e5;
    constexpr double maximumPressure = 1.0e7;
    constexpr double minimumTemperature = 250.0;
    constexpr double maximumTemperature = 1500.0;
    constexpr double perturbation = 1.0e-5;
    constexpr double tolerance = 2.0e-10;

    double logPressure = std::log(pressureSeed);
    double logTemperature = std::log(temperatureSeed);
    for (int iteration = 0; iteration < 30; ++iteration)
    {
        const double pressure = std::exp(logPressure);
        const double temperature = std::exp(logTemperature);
        const Evaluation base = evaluateVapor(
            kernel, caloric, pressure, temperature, massFractions);
        if (!base.valid) return result;

        const std::array<double, 2> baseResidual = residuals(
            base.state, targetDensity, targetEnergy);
        result.iterations = iteration + 1;
        result.densityResidual = std::abs(baseResidual[0]);
        result.energyResidual = std::abs(baseResidual[1]);
        if (result.densityResidual <= tolerance
            && result.energyResidual <= tolerance)
        {
            result.success = true;
            result.state = base.state;
            return result;
        }

        const Evaluation pressurePerturbed = evaluateVapor(
            kernel, caloric, std::exp(logPressure + perturbation),
            temperature, massFractions);
        const Evaluation temperaturePerturbed = evaluateVapor(
            kernel, caloric, pressure,
            std::exp(logTemperature + perturbation), massFractions);
        if (!pressurePerturbed.valid || !temperaturePerturbed.valid)
            return result;

        const std::array<double, 2> pressureResidual = residuals(
            pressurePerturbed.state, targetDensity, targetEnergy);
        const std::array<double, 2> temperatureResidual = residuals(
            temperaturePerturbed.state, targetDensity, targetEnergy);
        const double j00 = (pressureResidual[0] - baseResidual[0])
            /perturbation;
        const double j10 = (pressureResidual[1] - baseResidual[1])
            /perturbation;
        const double j01 = (temperatureResidual[0] - baseResidual[0])
            /perturbation;
        const double j11 = (temperatureResidual[1] - baseResidual[1])
            /perturbation;
        const double determinant = j00*j11 - j01*j10;
        if (!std::isfinite(determinant) || std::abs(determinant) < 1.0e-12)
            return result;

        double pressureStep =
            (-baseResidual[0]*j11 + j01*baseResidual[1])/determinant;
        double temperatureStep =
            (-j00*baseResidual[1] + j10*baseResidual[0])/determinant;
        pressureStep = std::max(-0.25, std::min(0.25, pressureStep));
        temperatureStep = std::max(-0.25, std::min(0.25, temperatureStep));

        bool accepted = false;
        const double baseNorm = normSquared(baseResidual);
        for (int lineSearch = 0; lineSearch < 10; ++lineSearch)
        {
            const double scale = std::ldexp(1.0, -lineSearch);
            const double trialLogPressure = std::max(
                std::log(minimumPressure),
                std::min(
                    std::log(maximumPressure),
                    logPressure + scale*pressureStep));
            const double trialLogTemperature = std::max(
                std::log(minimumTemperature),
                std::min(
                    std::log(maximumTemperature),
                    logTemperature + scale*temperatureStep));
            const Evaluation trial = evaluateVapor(
                kernel, caloric, std::exp(trialLogPressure),
                std::exp(trialLogTemperature), massFractions);
            if (!trial.valid) continue;
            if (normSquared(residuals(
                    trial.state, targetDensity, targetEnergy)) < baseNorm)
            {
                logPressure = trialLogPressure;
                logTemperature = trialLogTemperature;
                accepted = true;
                break;
            }
        }
        if (!accepted) return result;
    }
    return result;
}

std::optional<PRPhaseCaloricProperties> stateAtFixedDensity(
    const PRMixtureKernel& kernel,
    const PRPhaseCaloric& caloric,
    double targetDensity,
    double temperature,
    const std::vector<double>& massFractions,
    double pressureSeed)
{
    if (!finitePositive(targetDensity) || !finitePositive(temperature)
        || !finitePositive(pressureSeed))
    {
        return std::nullopt;
    }

    double logPressure = std::log(pressureSeed);
    constexpr double perturbation = 1.0e-5;
    for (int iteration = 0; iteration < 30; ++iteration)
    {
        const Evaluation base = evaluateVapor(
            kernel, caloric, std::exp(logPressure),
            temperature, massFractions);
        if (!base.valid) return std::nullopt;
        const double residual = std::log(base.state.rho/targetDensity);
        if (std::abs(residual) <= 2.0e-12) return base.state;

        const Evaluation perturbed = evaluateVapor(
            kernel, caloric, std::exp(logPressure + perturbation),
            temperature, massFractions);
        if (!perturbed.valid) return std::nullopt;
        const double derivative =
            (std::log(perturbed.state.rho/targetDensity) - residual)
            /perturbation;
        if (!std::isfinite(derivative) || std::abs(derivative) < 1.0e-12)
            return std::nullopt;
        const double step = std::max(-0.2, std::min(0.2, -residual/derivative));
        logPressure = std::max(
            std::log(1.0e5),
            std::min(std::log(1.0e7), logPressure + step));
    }
    return std::nullopt;
}

std::vector<double> compositionPerturbation(
    const std::vector<double>& initial)
{
    std::vector<double> result = initial;
    constexpr double displacement = 1.0e-7;
    result[0] += displacement;
    result[1] -= displacement;
    return result;
}

} // End anonymous namespace

int main()
{
    using namespace Foam;

    std::cout << "=== Offline vapor-root rho-e-Y feasibility ===\n";

    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB interactions = BinaryInteractionDB::createDefault();
    const std::vector<std::string> species = {"O2", "C3H8"};
    PREquilibriumThermoModel owner(
        critical, interactions, species, fixtureOptions());
    const std::vector<double> overallY =
        twoPhaseComposition(critical, interactions);
    const EquilibriumResult parent = owner.initializePTY(2.0e6, 300.0, overallY);
    check(parent.success(), "authoritative two-phase parent initializes");
    check(
        parent.success() && parent.state->region == EquilibriumRegion::TWO_PHASE
            && parent.state->vapor.has_value(),
        "authoritative parent carries a committed vapor payload");
    if (!parent.success() || parent.state->region != EquilibriumRegion::TWO_PHASE
        || !parent.state->vapor)
    {
        return 1;
    }

    const auto& vapor = *parent.state->vapor;
    PRMixtureKernel kernel(
        critical, interactions, species);
    PRPhaseCaloric caloric(
        critical, interactions, species);
    check(kernel.configured() && caloric.configured(),
        "frozen PR/kij/caloric phase backend configures");

    const Evaluation parentEvaluation = evaluateVapor(
        kernel, caloric, parent.state->p, parent.state->T,
        vapor.massFractions);
    check(parentEvaluation.valid,
        "explicit vapor root evaluates without whole-mixture acoustics");
    if (parentEvaluation.valid)
    {
        check(close(parentEvaluation.state.rho, vapor.rho, 2.0e-10),
            "parent vapor density is reproduced");
        check(close(parentEvaluation.state.e, vapor.e, 2.0e-10),
            "parent vapor internal energy is reproduced");
        check(close(parentEvaluation.state.h, vapor.h, 2.0e-10),
            "parent vapor enthalpy is reproduced");
        check(close(parentEvaluation.state.compressibility,
                vapor.compressibility, 2.0e-10),
            "parent vapor PR root identity is reproduced");
        check(closeVector(parentEvaluation.state.massFractions,
                vapor.massFractions, 2.0e-12),
            "parent vapor ordered composition is reproduced");
    }

    const EquilibriumResult wholeMixtureRecovery = owner.recoverRhoE(
        vapor.rho, vapor.e, vapor.massFractions);
    check(
        wholeMixtureRecovery.status
                == thermo::EquilibriumStatus::ACOUSTIC_CLOSURE_FAILED
            && !wholeMixtureRecovery.state.has_value()
            && wholeMixtureRecovery.diagnostics.backendStatus == "SUCCESS"
            && wholeMixtureRecovery.diagnostics.backendPhaseStatus == "SUCCESS"
            && wholeMixtureRecovery.diagnostics.backendAcousticStatus
                == "PHASE_TOPOLOGY_CHANGE",
        "whole-mixture recovery remains correctly blocked by acoustic topology");

    const RecoveryResult parentRecovered = recoverVaporRhoE(
        kernel, caloric, vapor.rho, vapor.e, vapor.massFractions,
        1.04*parent.state->p, 1.02*parent.state->T);
    check(parentRecovered.success,
        "biased seed recovers the committed vapor rho-e-Y payload");
    if (parentRecovered.state)
    {
        check(close(parentRecovered.state->pressure, parent.state->p, 2.0e-9)
                && close(parentRecovered.state->temperature,
                    parent.state->T, 2.0e-9),
            "parent recovery returns the committed p-T branch");
    }

    const RecoveryResult alternateSeedRecovered = recoverVaporRhoE(
        kernel, caloric, vapor.rho, vapor.e, vapor.massFractions,
        0.96*parent.state->p, 0.98*parent.state->T);
    check(alternateSeedRecovered.success,
        "an alternate numerical seed recovers the parent vapor payload");
    check(
        alternateSeedRecovered.state.has_value()
            && parentRecovered.state.has_value()
            && close(alternateSeedRecovered.state->pressure,
                parentRecovered.state->pressure, 2.0e-9)
            && close(alternateSeedRecovered.state->temperature,
                parentRecovered.state->temperature, 2.0e-9)
            && close(alternateSeedRecovered.state->compressibility,
                parentRecovered.state->compressibility, 2.0e-9),
        "reasonable seed changes do not change the recovered vapor root");

    const std::vector<double> displacedY =
        compositionPerturbation(vapor.massFractions);
    double displacedSum = 0.0;
    for (double fraction : displacedY) displacedSum += fraction;
    check(close(displacedSum, 1.0, 0.0, 2.0e-12)
            && *std::min_element(displacedY.begin(), displacedY.end()) >= 0.0,
        "normalized composition perturbation preserves mass and positivity");

    const auto target = stateAtFixedDensity(
        kernel, caloric, vapor.rho, parent.state->T + 2.0,
        displacedY, parent.state->p);
    check(target.has_value(),
        "a nearby displaced vapor target exists at fixed phase density");
    if (target)
    {
        const RecoveryResult displacedRecovered = recoverVaporRhoE(
            kernel, caloric, vapor.rho, target->e, displacedY,
            parent.state->p, parent.state->T);
        check(displacedRecovered.success,
            "fixed-density displaced rho-e-Y target recovers");
        if (displacedRecovered.state)
        {
            check(close(displacedRecovered.state->pressure,
                    target->pressure, 2.0e-9)
                    && close(displacedRecovered.state->temperature,
                        target->temperature, 2.0e-9)
                    && close(displacedRecovered.state->rho,
                        vapor.rho, 2.0e-10)
                    && close(displacedRecovered.state->e,
                        target->e, 2.0e-10),
                "displaced recovery closes p-T-rho-e on the vapor root");
        }
    }

    std::vector<double> invalidY = vapor.massFractions;
    invalidY[0] = -1.0e-8;
    const RecoveryResult invalidComposition = recoverVaporRhoE(
        kernel, caloric, vapor.rho, vapor.e, invalidY,
        parent.state->p, parent.state->T);
    check(!invalidComposition.success && !invalidComposition.state.has_value(),
        "invalid composition fails without exposing a phase state");

    const RecoveryResult invalidDensity = recoverVaporRhoE(
        kernel, caloric, 0.0, vapor.e, vapor.massFractions,
        parent.state->p, parent.state->T);
    check(!invalidDensity.success && !invalidDensity.state.has_value(),
        "invalid density fails without exposing a phase state");

    const RecoveryResult unreachableEnergy = recoverVaporRhoE(
        kernel, caloric, vapor.rho, 1.0e12, vapor.massFractions,
        parent.state->p, parent.state->T);
    check(!unreachableEnergy.success && !unreachableEnergy.state.has_value(),
        "out-of-envelope energy fails without exposing a phase state");

    std::cout << std::setprecision(17)
              << "parent_p=" << parent.state->p << '\n'
              << "parent_T=" << parent.state->T << '\n'
              << "parent_vapor_rho=" << vapor.rho << '\n'
              << "parent_vapor_e=" << vapor.e << '\n'
              << "parent_recovery_iterations="
              << parentRecovered.iterations << '\n';
    if (target)
    {
        std::cout << "displaced_target_p=" << target->pressure << '\n'
                  << "displaced_target_T=" << target->temperature << '\n'
                  << "displaced_target_e=" << target->e << '\n';
    }

    if (failures != 0)
    {
        std::cerr << "PR_PHASE_CALORIC_RHO_E_FEASIBILITY_FAIL failures="
                  << failures << '\n';
        return 1;
    }

    std::cout
        << "PR_PHASE_CALORIC_RHO_E_FEASIBILITY_PASS checks=21\n"
        << "SCOPE offline-diagnostic-only runtime_api=0 "
        << "phaseAwareChemistry=0 production_ready=0\n";
    return 0;
}
