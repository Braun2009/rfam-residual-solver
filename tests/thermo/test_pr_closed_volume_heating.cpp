#include "PhaseEquilibrium/ConservedPhaseMixture.H"
#include "PhaseEquilibrium/PRPhaseCaloric.H"
#include "PhaseEquilibrium/PRRhoEFlash.H"
#include "PhaseEquilibrium/PRTPFlash.H"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

using namespace Foam;
using namespace Foam::phaseEquilibrium;

namespace
{

int checks = 0;
int failures = 0;

void check(bool condition, const std::string& label)
{
    ++checks;
    std::cout << "  " << (condition ? "PASS " : "FAIL ") << label << '\n';
    if (!condition) ++failures;
}

bool close(
    double lhs,
    double rhs,
    double relativeTolerance = 1.0e-9,
    double absoluteTolerance = 1.0e-9)
{
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) return false;
    return std::abs(lhs - rhs) <= absoluteTolerance
        + relativeTolerance*std::max(std::abs(lhs), std::abs(rhs));
}

ConservedPhaseProperties conserved(
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

RhoEFlashOptions oxygenOptions()
{
    RhoEFlashOptions result;
    result.minimumPressure = 1.0e6;
    result.maximumPressure = 4.9e6;
    result.minimumTemperature = 100.0;
    result.maximumTemperature = 180.0;
    result.initialPressure = 2.0e6;
    result.initialTemperature = 132.0;
    result.pressureScanIntervals = 96;
    result.maximumIterations = 120;
    result.finiteDifferenceStep = 1.0e-4;
    result.saturationLogPressureTolerance = 1.0e-10;
    result.logStepLimit = 0.75;
    result.densityAbsoluteTolerance = 1.0e-7;
    result.densityRelativeTolerance = 2.0e-8;
    result.energyAbsoluteTolerance = 1.0e-4;
    result.energyRelativeTolerance = 2.0e-10;
    return result;
}

struct SaturatedReference
{
    bool valid = false;
    double density = std::numeric_limits<double>::quiet_NaN();
    double energy = std::numeric_limits<double>::quiet_NaN();
};

SaturatedReference saturatedOxygen(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& kij,
    double pressure,
    double vaporMassFraction)
{
    SaturatedReference result;
    PRTPFlash tp(critical, kij, {"O2"});
    const PureSaturationResult saturation = tp.solvePureSaturationTemperature(
        pressure, 100.0, 180.0);
    if (!saturation.valid()) return result;

    PRPhaseCaloric caloric(critical, kij, {"O2"});
    const PRPhaseCaloricResult liquid = caloric.evaluate(
        pressure, saturation.temperature, {1.0}, PhaseRoot::LIQUID);
    const PRPhaseCaloricResult vapor = caloric.evaluate(
        pressure, saturation.temperature, {1.0}, PhaseRoot::VAPOR);
    if (!liquid.success() || !vapor.success()) return result;

    const ConservedMixtureResult mixture = assembleConservedPhaseMixture(
        pressure,
        vaporMassFraction,
        conserved(*liquid.state),
        conserved(*vapor.state),
        {1.0});
    if (!mixture.success()) return result;

    result.valid = true;
    result.density = mixture.state->rho;
    result.energy = mixture.state->e;
    return result;
}

enum class HeatingStepStatus
{
    SUCCESS,
    INVALID_INPUT,
    STAGE_FLASH_FAILED,
    FINAL_FLASH_FAILED,
    ENERGY_LEDGER_FAILED
};

struct HeatingBoxSnapshot
{
    double specificEnergy = std::numeric_limits<double>::quiet_NaN();
    double heatLedger = 0.0;
    double time = 0.0;
    std::size_t acceptedSteps = 0;
    std::size_t rejectedSteps = 0;
    RhoEFlashResult equilibrium;
};

class ClosedVolumeHeatingBox
{
private:
    const PRRhoEFlash& flash_;
    double mass_;
    double volume_;
    double density_;
    double initialRecoveredEnergy_;
    HeatingBoxSnapshot committed_;

    RhoEFlashResult reflash(
        double specificEnergy,
        const RhoEFlashResult& seedState) const
    {
        RhoEFlashSeed seed;
        seed.pressure = seedState.pressure;
        seed.temperature = seedState.temperature;
        seed.localOnly = true;
        seed.requiredRegion = seedState.tpFlash.region;
        return flash_.flashMassFractions(
            density_, specificEnergy, {1.0}, seed);
    }

    HeatingStepStatus reject(HeatingStepStatus status)
    {
        ++committed_.rejectedSteps;
        return status;
    }

public:
    ClosedVolumeHeatingBox(
        const PRRhoEFlash& flash,
        double mass,
        double volume,
        double initialSpecificEnergy,
        RhoEFlashResult initialEquilibrium)
    :
        flash_(flash),
        mass_(mass),
        volume_(volume),
        density_(mass/volume),
        initialRecoveredEnergy_(initialEquilibrium.energy)
    {
        committed_.specificEnergy = initialSpecificEnergy;
        committed_.equilibrium = std::move(initialEquilibrium);
    }

    bool valid() const
    {
        return std::isfinite(mass_) && mass_ > 0.0
            && std::isfinite(volume_) && volume_ > 0.0
            && std::isfinite(density_) && density_ > 0.0
            && committed_.equilibrium.success()
            && committed_.equilibrium.tpFlash.region
                == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE
            && committed_.equilibrium.conservedMixture.has_value()
            && close(committed_.equilibrium.targetDensity, density_, 2.0e-8)
            && close(
                committed_.equilibrium.targetEnergy,
                committed_.specificEnergy,
                2.0e-10,
                1.0e-4);
    }

    HeatingStepStatus advanceRK2(
        double dt,
        double wallConductance,
        double wallTemperature)
    {
        if (!valid() || !std::isfinite(dt) || !(dt > 0.0)
            || !std::isfinite(wallConductance) || wallConductance < 0.0
            || !std::isfinite(wallTemperature) || !(wallTemperature > 0.0))
            return reject(HeatingStepStatus::INVALID_INPUT);

        const double stageOnePower = wallConductance
            *(wallTemperature - committed_.equilibrium.temperature);
        const double stageSpecificEnergy = committed_.specificEnergy
            + dt*stageOnePower/mass_;
        if (!std::isfinite(stageSpecificEnergy))
            return reject(HeatingStepStatus::STAGE_FLASH_FAILED);

        RhoEFlashResult stage = reflash(
            stageSpecificEnergy, committed_.equilibrium);
        if (!stage.success())
            return reject(HeatingStepStatus::STAGE_FLASH_FAILED);

        const double stageTwoPower = wallConductance
            *(wallTemperature - stage.temperature);
        const double heatIncrement = 0.5*dt
            *(stageOnePower + stageTwoPower);
        const double candidateSpecificEnergy = committed_.specificEnergy
            + heatIncrement/mass_;
        if (!std::isfinite(heatIncrement)
            || !std::isfinite(candidateSpecificEnergy))
            return reject(HeatingStepStatus::FINAL_FLASH_FAILED);

        RhoEFlashResult candidate = reflash(
            candidateSpecificEnergy, stage);
        if (!candidate.success())
            return reject(HeatingStepStatus::FINAL_FLASH_FAILED);

        const double candidateLedger = committed_.heatLedger + heatIncrement;
        const double incrementalEnergyResidual = mass_
            *(candidate.energy - committed_.equilibrium.energy)
            - heatIncrement;
        const double cumulativeEnergyResidual = mass_
            *(candidate.energy - initialRecoveredEnergy_)
            - candidateLedger;
        if (!std::isfinite(candidateLedger)
            || !close(incrementalEnergyResidual, 0.0, 0.0, 1.0e-3)
            || !close(cumulativeEnergyResidual, 0.0, 0.0, 1.0e-3))
            return reject(HeatingStepStatus::ENERGY_LEDGER_FAILED);

        committed_.specificEnergy = candidateSpecificEnergy;
        committed_.heatLedger = candidateLedger;
        committed_.time += dt;
        ++committed_.acceptedSteps;
        committed_.equilibrium = std::move(candidate);
        return HeatingStepStatus::SUCCESS;
    }

    double mass() const { return mass_; }
    double volume() const { return volume_; }
    double density() const { return density_; }
    const HeatingBoxSnapshot& state() const { return committed_; }
};

double vaporMassFraction(const HeatingBoxSnapshot& state)
{
    if (!state.equilibrium.conservedMixture)
        return std::numeric_limits<double>::quiet_NaN();
    return state.equilibrium.conservedMixture->vaporMassFraction;
}

bool sameValue(double lhs, double rhs)
{
    return lhs == rhs || (std::isnan(lhs) && std::isnan(rhs));
}

bool sameVector(const std::vector<double>& lhs, const std::vector<double>& rhs)
{
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t i = 0; i < lhs.size(); ++i)
        if (!sameValue(lhs[i], rhs[i])) return false;
    return true;
}

bool samePhaseState(const PRPhaseState& lhs, const PRPhaseState& rhs)
{
    return lhs.status == rhs.status
        && lhs.requestedRoot == rhs.requestedRoot
        && sameVector(lhs.moleFractions, rhs.moleFractions)
        && sameVector(
            lhs.lnFugacityCoefficients, rhs.lnFugacityCoefficients)
        && sameValue(lhs.compressibility, rhs.compressibility)
        && sameValue(lhs.residualGibbsOverRT, rhs.residualGibbsOverRT)
        && lhs.physicalRootCount == rhs.physicalRootCount;
}

bool sameTPFlash(const TPFlashResult& lhs, const TPFlashResult& rhs)
{
    return lhs.status == rhs.status && lhs.region == rhs.region
        && lhs.phaseFractionDetermined == rhs.phaseFractionDetermined
        && sameValue(lhs.pressure, rhs.pressure)
        && sameValue(lhs.temperature, rhs.temperature)
        && sameValue(lhs.vaporMoleFraction, rhs.vaporMoleFraction)
        && sameVector(lhs.feedMoleFractions, rhs.feedMoleFractions)
        && sameVector(lhs.feedMassFractions, rhs.feedMassFractions)
        && sameVector(lhs.liquidMoleFractions, rhs.liquidMoleFractions)
        && sameVector(lhs.vaporMoleFractions, rhs.vaporMoleFractions)
        && sameVector(
            lhs.singlePhaseMoleFractions, rhs.singlePhaseMoleFractions)
        && sameVector(lhs.liquidMassFractions, rhs.liquidMassFractions)
        && sameVector(lhs.vaporMassFractions, rhs.vaporMassFractions)
        && sameVector(
            lhs.singlePhaseMassFractions, rhs.singlePhaseMassFractions)
        && samePhaseState(lhs.liquidState, rhs.liquidState)
        && samePhaseState(lhs.vaporState, rhs.vaporState)
        && samePhaseState(lhs.singlePhaseState, rhs.singlePhaseState)
        && sameValue(lhs.feedMinimumTpd, rhs.feedMinimumTpd)
        && sameValue(lhs.liquidMinimumTpd, rhs.liquidMinimumTpd)
        && sameValue(lhs.vaporMinimumTpd, rhs.vaporMinimumTpd)
        && sameValue(
            lhs.maximumLogFugacityResidual,
            rhs.maximumLogFugacityResidual)
        && sameValue(
            lhs.maximumMoleMaterialResidual,
            rhs.maximumMoleMaterialResidual)
        && sameValue(lhs.rachfordRiceResidual, rhs.rachfordRiceResidual)
        && lhs.flashIterations == rhs.flashIterations;
}

bool sameCaloricState(
    const PRPhaseCaloricProperties& lhs,
    const PRPhaseCaloricProperties& rhs)
{
    return lhs.root == rhs.root
        && sameVector(lhs.moleFractions, rhs.moleFractions)
        && sameVector(lhs.massFractions, rhs.massFractions)
        && sameValue(lhs.pressure, rhs.pressure)
        && sameValue(lhs.temperature, rhs.temperature)
        && sameValue(lhs.molecularWeight, rhs.molecularWeight)
        && sameValue(lhs.compressibility, rhs.compressibility)
        && sameValue(lhs.rho, rhs.rho) && sameValue(lhs.h, rhs.h)
        && sameValue(lhs.e, rhs.e) && sameValue(lhs.s, rhs.s);
}

bool sameOptionalCaloric(
    const std::optional<PRPhaseCaloricProperties>& lhs,
    const std::optional<PRPhaseCaloricProperties>& rhs)
{
    return lhs.has_value() == rhs.has_value()
        && (!lhs || sameCaloricState(*lhs, *rhs));
}

bool sameOptionalMixture(
    const std::optional<ConservedPhaseMixture>& lhs,
    const std::optional<ConservedPhaseMixture>& rhs)
{
    return lhs.has_value() == rhs.has_value()
        && (!lhs
            || (sameValue(lhs->vaporMoleFraction, rhs->vaporMoleFraction)
                && sameValue(
                    lhs->vaporMassFraction, rhs->vaporMassFraction)
                && sameValue(
                    lhs->vaporVolumeFraction, rhs->vaporVolumeFraction)
                && sameValue(lhs->rho, rhs->rho)
                && sameValue(lhs->h, rhs->h)
                && sameValue(lhs->e, rhs->e)
                && sameValue(lhs->s, rhs->s)
                && sameVector(
                    lhs->overallMassFractions,
                    rhs->overallMassFractions)));
}

bool sameFlashResult(const RhoEFlashResult& lhs, const RhoEFlashResult& rhs)
{
    return lhs.status == rhs.status && lhs.phaseStatus == rhs.phaseStatus
        && lhs.caloricStatus == rhs.caloricStatus
        && lhs.conservedStatus == rhs.conservedStatus
        && sameValue(lhs.targetDensity, rhs.targetDensity)
        && sameValue(lhs.targetEnergy, rhs.targetEnergy)
        && sameValue(lhs.pressure, rhs.pressure)
        && sameValue(lhs.temperature, rhs.temperature)
        && sameValue(lhs.density, rhs.density)
        && sameValue(lhs.energy, rhs.energy)
        && sameValue(lhs.entropy, rhs.entropy)
        && sameValue(lhs.densityResidual, rhs.densityResidual)
        && sameValue(
            lhs.densityRelativeResidual, rhs.densityRelativeResidual)
        && sameValue(lhs.energyResidual, rhs.energyResidual)
        && sameValue(lhs.energyRelativeResidual, rhs.energyRelativeResidual)
        && sameValue(lhs.tpFugacityResidual, rhs.tpFugacityResidual)
        && sameValue(lhs.tpMaterialResidual, rhs.tpMaterialResidual)
        && sameValue(
            lhs.phaseEnergyIdentityResidual,
            rhs.phaseEnergyIdentityResidual)
        && sameValue(
            lhs.mixtureMaterialResidual, rhs.mixtureMaterialResidual)
        && sameValue(
            lhs.mixtureEnergyIdentityResidual,
            rhs.mixtureEnergyIdentityResidual)
        && sameValue(lhs.vaporMoleFraction, rhs.vaporMoleFraction)
        && lhs.pureVolumeEnergyLeverRule == rhs.pureVolumeEnergyLeverRule
        && lhs.iterations == rhs.iterations
        && lhs.validSamples == rhs.validSamples
        && lhs.invalidSamples == rhs.invalidSamples
        && sameTPFlash(lhs.tpFlash, rhs.tpFlash)
        && sameOptionalCaloric(lhs.liquidPhase, rhs.liquidPhase)
        && sameOptionalCaloric(lhs.vaporPhase, rhs.vaporPhase)
        && sameOptionalCaloric(lhs.singlePhase, rhs.singlePhase)
        && sameOptionalMixture(lhs.conservedMixture, rhs.conservedMixture);
}

bool sameCommittedState(
    const HeatingBoxSnapshot& lhs,
    const HeatingBoxSnapshot& rhs)
{
    return sameValue(lhs.specificEnergy, rhs.specificEnergy)
        && sameValue(lhs.heatLedger, rhs.heatLedger)
        && sameValue(lhs.time, rhs.time)
        && lhs.acceptedSteps == rhs.acceptedSteps
        && sameFlashResult(lhs.equilibrium, rhs.equilibrium);
}

struct Trajectory
{
    bool success = false;
    bool monotone = false;
    double maximumFlashEnergyResidual = 0.0;
    HeatingBoxSnapshot finalState;
};

Trajectory runTrajectory(
    const PRRhoEFlash& flash,
    const RhoEFlashResult& initial,
    double initialSpecificEnergy,
    std::size_t stepCount,
    double totalTime,
    double wallConductance,
    double wallTemperature)
{
    Trajectory result;
    constexpr double mass = 1.0;
    const double volume = mass/initial.targetDensity;
    ClosedVolumeHeatingBox box(
        flash, mass, volume, initialSpecificEnergy, initial);
    if (!box.valid() || stepCount == 0) return result;

    double previousQuality = vaporMassFraction(box.state());
    double previousPressure = box.state().equilibrium.pressure;
    double previousTemperature = box.state().equilibrium.temperature;
    result.monotone = std::isfinite(previousQuality);
    const double dt = totalTime/static_cast<double>(stepCount);
    for (std::size_t step = 0; step < stepCount; ++step)
    {
        if (box.advanceRK2(dt, wallConductance, wallTemperature)
            != HeatingStepStatus::SUCCESS)
            return result;
        const HeatingBoxSnapshot& state = box.state();
        const double quality = vaporMassFraction(state);
        result.monotone = result.monotone
            && state.equilibrium.tpFlash.region
                == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE
            && std::isfinite(quality) && quality >= previousQuality
            && quality >= 0.0 && quality <= 1.0
            && state.equilibrium.pressure >= previousPressure
            && state.equilibrium.temperature >= previousTemperature;
        result.maximumFlashEnergyResidual = std::max(
            result.maximumFlashEnergyResidual,
            std::abs(state.equilibrium.energyResidual));
        previousQuality = quality;
        previousPressure = state.equilibrium.pressure;
        previousTemperature = state.equilibrium.temperature;
    }

    result.finalState = box.state();
    result.success = box.valid()
        && result.finalState.acceptedSteps == stepCount
        && result.finalState.rejectedSteps == 0
        && close(result.finalState.time, totalTime, 1.0e-12, 1.0e-14)
        && close(
            result.finalState.heatLedger,
            mass*(result.finalState.equilibrium.energy - initial.energy),
            0.0,
            1.0e-3);
    return result;
}

double convergenceRatio(double coarse, double medium, double fine)
{
    const double numerator = std::abs(coarse - medium);
    const double denominator = std::abs(medium - fine);
    return denominator > 0.0
        ? numerator/denominator
        : std::numeric_limits<double>::infinity();
}

bool secondOrderRatios(
    double coarse,
    double medium,
    double fine,
    double reference)
{
    const double first = convergenceRatio(coarse, medium, fine);
    const double second = convergenceRatio(medium, fine, reference);
    return first >= 3.2 && first <= 5.0
        && second >= 3.2 && second <= 5.0;
}

} // End anonymous namespace

int main()
{
    std::cout << "=== Pure-O2 closed-volume HEM heating Box A1 ===\n";

    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB kij = BinaryInteractionDB::createDefault();
    const SaturatedReference reference = saturatedOxygen(
        critical, kij, 2.0e6, 0.10);
    check(reference.valid,
        "initial pure-O2 saturation state is constructible");
    if (!reference.valid) return 1;

    const RhoEFlashOptions options = oxygenOptions();
    PRRhoEFlash flash(critical, kij, {"O2"}, options);
    const RhoEFlashResult initial = flash.flashMassFractions(
        reference.density, reference.energy, {1.0});
    check(initial.success()
            && initial.tpFlash.region
                == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE,
        "rho-e flash owns the initial two-phase state");
    if (!initial.success()) return 1;

    constexpr double totalTime = 0.4;
    constexpr double wallConductance = 5000.0;
    constexpr double wallTemperature = 150.0;
    const Trajectory coarse = runTrajectory(
        flash, initial, reference.energy, 4, totalTime,
        wallConductance, wallTemperature);
    const Trajectory medium = runTrajectory(
        flash, initial, reference.energy, 8, totalTime,
        wallConductance, wallTemperature);
    const Trajectory fine = runTrajectory(
        flash, initial, reference.energy, 16, totalTime,
        wallConductance, wallTemperature);
    const Trajectory referenceTrajectory = runTrajectory(
        flash, initial, reference.energy, 32, totalTime,
        wallConductance, wallTemperature);

    check(coarse.success && medium.success && fine.success
            && referenceTrajectory.success,
        "all timestep refinements complete without a hidden fallback");
    check(referenceTrajectory.monotone,
        "heating raises saturation pressure, temperature, and vapor quality monotonically");
    const double initialQuality = initial.conservedMixture
        ? initial.conservedMixture->vaporMassFraction
        : std::numeric_limits<double>::quiet_NaN();
    check(referenceTrajectory.finalState.heatLedger > 0.0
            && vaporMassFraction(referenceTrajectory.finalState)
                > initialQuality,
        "wall heat produces evaporation rather than an independent quality update");
    check(referenceTrajectory.maximumFlashEnergyResidual < 2.0e-4,
        "every accepted reflash closes its target specific energy");
    check(close(
            referenceTrajectory.finalState.equilibrium.conservedMixture->e,
            referenceTrajectory.finalState.equilibrium.conservedMixture->h
                - referenceTrajectory.finalState.equilibrium.pressure
                    /referenceTrajectory.finalState.equilibrium.conservedMixture->rho,
            2.0e-10,
            1.0e-5),
        "final phase mixture retains e=h-p/rho on the shared energy basis");
    check(referenceTrajectory.finalState.equilibrium.conservedMixture
                ->overallMassFractions.size() == 1
            && close(
                referenceTrajectory.finalState.equilibrium.conservedMixture
                    ->overallMassFractions[0],
                1.0,
                0.0,
                2.0e-12),
        "fixed-mass Box A preserves the pure-O2 species inventory");

    const double energyRatioOne = convergenceRatio(
        coarse.finalState.specificEnergy,
        medium.finalState.specificEnergy,
        fine.finalState.specificEnergy);
    const double energyRatioTwo = convergenceRatio(
        medium.finalState.specificEnergy,
        fine.finalState.specificEnergy,
        referenceTrajectory.finalState.specificEnergy);
    std::cout << "    dt Richardson ratios e=" << energyRatioOne
              << ',' << energyRatioTwo
              << " p=" << convergenceRatio(
                    coarse.finalState.equilibrium.pressure,
                    medium.finalState.equilibrium.pressure,
                    fine.finalState.equilibrium.pressure)
              << ',' << convergenceRatio(
                    medium.finalState.equilibrium.pressure,
                    fine.finalState.equilibrium.pressure,
                    referenceTrajectory.finalState.equilibrium.pressure)
              << " T=" << convergenceRatio(
                    coarse.finalState.equilibrium.temperature,
                    medium.finalState.equilibrium.temperature,
                    fine.finalState.equilibrium.temperature)
              << ',' << convergenceRatio(
                    medium.finalState.equilibrium.temperature,
                    fine.finalState.equilibrium.temperature,
                    referenceTrajectory.finalState.equilibrium.temperature)
              << " quality=" << convergenceRatio(
                    vaporMassFraction(coarse.finalState),
                    vaporMassFraction(medium.finalState),
                    vaporMassFraction(fine.finalState))
              << ',' << convergenceRatio(
                    vaporMassFraction(medium.finalState),
                    vaporMassFraction(fine.finalState),
                    vaporMassFraction(referenceTrajectory.finalState))
              << '\n';
    check(secondOrderRatios(
                coarse.finalState.specificEnergy,
                medium.finalState.specificEnergy,
                fine.finalState.specificEnergy,
                referenceTrajectory.finalState.specificEnergy)
            && secondOrderRatios(
                coarse.finalState.equilibrium.pressure,
                medium.finalState.equilibrium.pressure,
                fine.finalState.equilibrium.pressure,
                referenceTrajectory.finalState.equilibrium.pressure)
            && secondOrderRatios(
                coarse.finalState.equilibrium.temperature,
                medium.finalState.equilibrium.temperature,
                fine.finalState.equilibrium.temperature,
                referenceTrajectory.finalState.equilibrium.temperature)
            && secondOrderRatios(
                vaporMassFraction(coarse.finalState),
                vaporMassFraction(medium.finalState),
                vaporMassFraction(fine.finalState),
                vaporMassFraction(referenceTrajectory.finalState))
            && secondOrderRatios(
                coarse.finalState.heatLedger,
                medium.finalState.heatLedger,
                fine.finalState.heatLedger,
                referenceTrajectory.finalState.heatLedger),
        "RK2 wall-heating e/p/T/quality/heat converge at second-order scale");

    constexpr double mass = 1.0;
    ClosedVolumeHeatingBox rollbackBox(
        flash, mass, mass/reference.density, reference.energy, initial);
    const HeatingBoxSnapshot beforeFailure = rollbackBox.state();
    const HeatingStepStatus failure = rollbackBox.advanceRK2(
        100.0, wallConductance, 1000.0);
    const HeatingBoxSnapshot afterFailure = rollbackBox.state();
    check(failure == HeatingStepStatus::STAGE_FLASH_FAILED
            || failure == HeatingStepStatus::FINAL_FLASH_FAILED,
        "an out-of-envelope heating candidate fails closed");
    check(sameCommittedState(beforeFailure, afterFailure)
            && afterFailure.rejectedSteps
                == beforeFailure.rejectedSteps + 1,
        "failed reflash rolls back energy, phase state, time, and heat ledger");

    const HeatingBoxSnapshot beforeInvalid = rollbackBox.state();
    const HeatingStepStatus invalid = rollbackBox.advanceRK2(
        std::numeric_limits<double>::quiet_NaN(),
        wallConductance,
        wallTemperature);
    check(invalid == HeatingStepStatus::INVALID_INPUT
            && sameCommittedState(beforeInvalid, rollbackBox.state())
            && rollbackBox.state().rejectedSteps
                == beforeInvalid.rejectedSteps + 1,
        "invalid timestep is rejected without publishing a candidate");

    if (failures != 0)
    {
        std::cerr << "PR_CLOSED_VOLUME_HEATING_FAIL failures="
                  << failures << " checks=" << checks << '\n';
        return 1;
    }

    std::cout << "PR_CLOSED_VOLUME_HEATING_PASS checks=" << checks << '\n'
              << "NOTE: initially two-phase rigid Box A1 only; external caloric "
              << "validation, transport, chemistry, and production capability "
              << "flags remain closed.\n";
    return 0;
}
