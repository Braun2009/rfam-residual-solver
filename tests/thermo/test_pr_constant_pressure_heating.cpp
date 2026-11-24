#include "PhaseEquilibrium/PRPHFlash.H"
#include "PhaseEquilibrium/PRPhaseCaloric.H"
#include "PhaseEquilibrium/PRRhoEFlash.H"
#include "PhaseEquilibrium/PRTPFlash.H"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <optional>
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
    double absoluteTolerance = 1.0e-8)
{
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) return false;
    return std::abs(lhs - rhs) <= absoluteTolerance
        + relativeTolerance*std::max(std::abs(lhs), std::abs(rhs));
}

PHFlashOptions phOptions()
{
    PHFlashOptions result;
    result.minimumTemperature = 100.0;
    result.maximumTemperature = 180.0;
    result.scanIntervals = 64;
    result.maximumIterations = 160;
    result.enthalpyAbsoluteTolerance = 1.0e-5;
    result.enthalpyRelativeTolerance = 1.0e-10;
    return result;
}

RhoEFlashOptions rhoEOptions()
{
    RhoEFlashOptions result;
    result.minimumPressure = 2.0e6;
    result.maximumPressure = 4.9e6;
    result.minimumTemperature = 100.0;
    result.maximumTemperature = 180.0;
    result.initialPressure = 4.0e6;
    result.initialTemperature = 130.0;
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

struct PhasePoint
{
    bool valid = false;
    EquilibriumRegion region = EquilibriumRegion::UNKNOWN;
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double rho = std::numeric_limits<double>::quiet_NaN();
    double h = std::numeric_limits<double>::quiet_NaN();
    double e = std::numeric_limits<double>::quiet_NaN();
    double quality = std::numeric_limits<double>::quiet_NaN();
};

PhasePoint pointFromPH(const PHFlashResult& result)
{
    PhasePoint point;
    if (!result.success()) return point;
    point.region = result.tpFlash.region;
    point.temperature = result.temperature;
    if (result.conservedMixture)
    {
        point.rho = result.conservedMixture->rho;
        point.h = result.conservedMixture->h;
        point.e = result.conservedMixture->e;
        point.quality = result.conservedMixture->vaporMassFraction;
    }
    else if (result.singlePhase)
    {
        point.rho = result.singlePhase->rho;
        point.h = result.singlePhase->h;
        point.e = result.singlePhase->e;
        if (point.region == EquilibriumRegion::SINGLE_LIQUID)
            point.quality = 0.0;
        else if (point.region == EquilibriumRegion::SINGLE_VAPOR)
            point.quality = 1.0;
    }
    point.valid = std::isfinite(point.temperature)
        && std::isfinite(point.rho) && point.rho > 0.0
        && std::isfinite(point.h) && std::isfinite(point.e)
        && std::isfinite(point.quality)
        && point.quality >= 0.0 && point.quality <= 1.0;
    return point;
}

PhasePoint pointFromRhoE(const RhoEFlashResult& result)
{
    PhasePoint point;
    if (!result.success()) return point;
    point.region = result.tpFlash.region;
    point.temperature = result.temperature;
    point.rho = result.density;
    point.e = result.energy;
    if (result.conservedMixture)
    {
        point.h = result.conservedMixture->h;
        point.quality = result.conservedMixture->vaporMassFraction;
    }
    else if (result.singlePhase)
    {
        point.h = result.singlePhase->h;
        if (point.region == EquilibriumRegion::SINGLE_LIQUID)
            point.quality = 0.0;
        else if (point.region == EquilibriumRegion::SINGLE_VAPOR)
            point.quality = 1.0;
    }
    point.valid = std::isfinite(point.temperature)
        && std::isfinite(point.rho) && point.rho > 0.0
        && std::isfinite(point.h) && std::isfinite(point.e)
        && std::isfinite(point.quality)
        && point.quality >= 0.0 && point.quality <= 1.0;
    return point;
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

bool samePHFlashResult(const PHFlashResult& lhs, const PHFlashResult& rhs)
{
    return lhs.status == rhs.status && lhs.phaseStatus == rhs.phaseStatus
        && lhs.caloricStatus == rhs.caloricStatus
        && lhs.conservedStatus == rhs.conservedStatus
        && sameValue(lhs.pressure, rhs.pressure)
        && sameValue(lhs.targetEnthalpy, rhs.targetEnthalpy)
        && sameValue(lhs.temperature, rhs.temperature)
        && sameValue(lhs.vaporMoleFraction, rhs.vaporMoleFraction)
        && sameValue(lhs.computedEnthalpy, rhs.computedEnthalpy)
        && sameValue(lhs.enthalpyResidual, rhs.enthalpyResidual)
        && sameValue(
            lhs.enthalpyRelativeResidual, rhs.enthalpyRelativeResidual)
        && sameValue(lhs.entropy, rhs.entropy)
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
        && lhs.pureComponentEnthalpyLeverRule
            == rhs.pureComponentEnthalpyLeverRule
        && lhs.iterations == rhs.iterations
        && lhs.validTemperatureSamples == rhs.validTemperatureSamples
        && lhs.invalidTemperatureSamples == rhs.invalidTemperatureSamples
        && sameTPFlash(lhs.tpFlash, rhs.tpFlash)
        && sameOptionalCaloric(lhs.liquidPhase, rhs.liquidPhase)
        && sameOptionalCaloric(lhs.vaporPhase, rhs.vaporPhase)
        && sameOptionalCaloric(lhs.singlePhase, rhs.singlePhase)
        && sameOptionalMixture(lhs.conservedMixture, rhs.conservedMixture);
}

bool sameRhoEFlashResult(
    const RhoEFlashResult& lhs,
    const RhoEFlashResult& rhs)
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

enum class PistonStepStatus
{
    SUCCESS,
    INVALID_HEAT,
    PH_FLASH_FAILED,
    RHO_E_REFLASH_FAILED,
    CROSS_FLASH_MISMATCH,
    FIRST_LAW_FAILED
};

struct PistonSnapshot
{
    double targetEnthalpy = std::numeric_limits<double>::quiet_NaN();
    double heatLedger = 0.0;
    double boundaryWorkLedger = 0.0;
    double time = 0.0;
    std::size_t acceptedSteps = 0;
    std::size_t rejectedSteps = 0;
    PHFlashResult ph;
    RhoEFlashResult equilibrium;
};

class ConstantPressureHeatingBox
{
private:
    const PRPHFlash& phFlash_;
    const PRRhoEFlash& rhoEFlash_;
    double pressure_;
    double mass_;
    double initialEnthalpy_;
    double initialEnergy_;
    double initialVolume_;
    PistonSnapshot committed_;

    PistonStepStatus reject(PistonStepStatus status)
    {
        ++committed_.rejectedSteps;
        return status;
    }

public:
    ConstantPressureHeatingBox(
        const PRPHFlash& phFlash,
        const PRRhoEFlash& rhoEFlash,
        double pressure,
        double mass,
        double targetEnthalpy,
        PHFlashResult initialPH,
        RhoEFlashResult initialEquilibrium)
    :
        phFlash_(phFlash),
        rhoEFlash_(rhoEFlash),
        pressure_(pressure),
        mass_(mass),
        initialEnthalpy_(pointFromRhoE(initialEquilibrium).h),
        initialEnergy_(initialEquilibrium.energy),
        initialVolume_(mass/initialEquilibrium.density)
    {
        committed_.targetEnthalpy = targetEnthalpy;
        committed_.ph = std::move(initialPH);
        committed_.equilibrium = std::move(initialEquilibrium);
    }

    bool valid() const
    {
        const PhasePoint phPoint = pointFromPH(committed_.ph);
        const PhasePoint rhoEPoint = pointFromRhoE(committed_.equilibrium);
        return std::isfinite(pressure_) && pressure_ > 0.0
            && std::isfinite(mass_) && mass_ > 0.0
            && phPoint.valid && rhoEPoint.valid
            && phPoint.region == rhoEPoint.region
            && close(phPoint.temperature, rhoEPoint.temperature, 2.0e-7)
            && close(phPoint.rho, rhoEPoint.rho, 2.0e-8)
            && close(phPoint.h, rhoEPoint.h, 2.0e-9, 2.0e-4)
            && close(phPoint.e, rhoEPoint.e, 2.0e-9, 2.0e-4)
            && close(phPoint.quality, rhoEPoint.quality, 2.0e-8, 1.0e-9)
            && close(committed_.ph.pressure, pressure_, 0.0, 0.0)
            && close(
                committed_.equilibrium.pressure,
                pressure_,
                2.0e-7,
                1.0e-2);
    }

    PistonStepStatus addHeat(double heatIncrement, double elapsedTime)
    {
        if (!valid() || !std::isfinite(heatIncrement)
            || !(heatIncrement > 0.0) || !std::isfinite(elapsedTime)
            || !(elapsedTime > 0.0))
            return reject(PistonStepStatus::INVALID_HEAT);

        const double candidateTargetEnthalpy = committed_.targetEnthalpy
            + heatIncrement/mass_;
        PHFlashResult candidatePH = phFlash_.flashMassFractions(
            pressure_, candidateTargetEnthalpy, {1.0});
        const PhasePoint phPoint = pointFromPH(candidatePH);
        if (!phPoint.valid)
            return reject(PistonStepStatus::PH_FLASH_FAILED);

        RhoEFlashResult candidateEquilibrium =
            rhoEFlash_.flashMassFractions(phPoint.rho, phPoint.e, {1.0});
        const PhasePoint candidate = pointFromRhoE(candidateEquilibrium);
        if (!candidate.valid)
            return reject(PistonStepStatus::RHO_E_REFLASH_FAILED);
        if (candidate.region != phPoint.region
            || !close(candidate.temperature, phPoint.temperature, 2.0e-7)
            || !close(candidate.rho, phPoint.rho, 2.0e-8)
            || !close(candidate.h, phPoint.h, 2.0e-9, 2.0e-4)
            || !close(candidate.e, phPoint.e, 2.0e-9, 2.0e-4)
            || !close(candidate.quality, phPoint.quality, 2.0e-8, 1.0e-9)
            || !close(candidateEquilibrium.pressure, pressure_, 2.0e-7, 1.0e-2))
            return reject(PistonStepStatus::CROSS_FLASH_MISMATCH);

        const PhasePoint previous = pointFromRhoE(committed_.equilibrium);
        const double previousVolume = mass_/previous.rho;
        const double candidateVolume = mass_/candidate.rho;
        const double workIncrement = pressure_
            *(candidateVolume - previousVolume);
        const double incrementalResidual = heatIncrement
            - (mass_*(candidate.e - previous.e) + workIncrement);
        const double candidateHeatLedger = committed_.heatLedger
            + heatIncrement;
        const double candidateWorkLedger = committed_.boundaryWorkLedger
            + workIncrement;
        const double cumulativeResidual = candidateHeatLedger
            - (mass_*(candidate.e - initialEnergy_)
                + candidateWorkLedger);
        const double enthalpyResidual = candidateHeatLedger
            - mass_*(candidate.h - initialEnthalpy_);
        if (!close(incrementalResidual, 0.0, 0.0, 2.0e-3)
            || !close(cumulativeResidual, 0.0, 0.0, 2.0e-3)
            || !close(enthalpyResidual, 0.0, 0.0, 2.0e-3))
            return reject(PistonStepStatus::FIRST_LAW_FAILED);

        committed_.targetEnthalpy = candidateTargetEnthalpy;
        committed_.heatLedger = candidateHeatLedger;
        committed_.boundaryWorkLedger = candidateWorkLedger;
        committed_.time += elapsedTime;
        ++committed_.acceptedSteps;
        committed_.ph = std::move(candidatePH);
        committed_.equilibrium = std::move(candidateEquilibrium);
        return PistonStepStatus::SUCCESS;
    }

    double pressure() const { return pressure_; }
    double mass() const { return mass_; }
    double initialVolume() const { return initialVolume_; }
    const PistonSnapshot& state() const { return committed_; }
};

bool sameCommittedState(
    const PistonSnapshot& lhs,
    const PistonSnapshot& rhs)
{
    return sameValue(lhs.targetEnthalpy, rhs.targetEnthalpy)
        && sameValue(lhs.heatLedger, rhs.heatLedger)
        && sameValue(lhs.boundaryWorkLedger, rhs.boundaryWorkLedger)
        && sameValue(lhs.time, rhs.time)
        && lhs.acceptedSteps == rhs.acceptedSteps
        && samePHFlashResult(lhs.ph, rhs.ph)
        && sameRhoEFlashResult(lhs.equilibrium, rhs.equilibrium);
}

struct RunResult
{
    bool success = false;
    bool crossedBoundary = false;
    bool monotoneQuality = false;
    PistonSnapshot finalState;
};

RunResult runHeatingPath(
    const PRPHFlash& phFlash,
    const PRRhoEFlash& rhoEFlash,
    double pressure,
    double mass,
    double initialEnthalpy,
    const PHFlashResult& initialPH,
    const RhoEFlashResult& initialEquilibrium,
    double totalHeat,
    std::size_t stepCount)
{
    RunResult result;
    if (stepCount == 0) return result;
    ConstantPressureHeatingBox box(
        phFlash,
        rhoEFlash,
        pressure,
        mass,
        initialEnthalpy,
        initialPH,
        initialEquilibrium);
    if (!box.valid()) return result;

    double previousQuality = pointFromRhoE(box.state().equilibrium).quality;
    result.monotoneQuality = previousQuality == 0.0;
    for (std::size_t step = 0; step < stepCount; ++step)
    {
        if (box.addHeat(
                totalHeat/static_cast<double>(stepCount),
                1.0/static_cast<double>(stepCount))
            != PistonStepStatus::SUCCESS)
            return result;
        const PhasePoint point = pointFromRhoE(box.state().equilibrium);
        result.crossedBoundary = result.crossedBoundary
            || point.region == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE;
        result.monotoneQuality = result.monotoneQuality
            && point.quality >= previousQuality
            && point.quality >= 0.0 && point.quality <= 1.0;
        previousQuality = point.quality;
    }

    result.finalState = box.state();
    const PhasePoint finalPoint = pointFromRhoE(result.finalState.equilibrium);
    const double finalVolume = mass/finalPoint.rho;
    result.success = box.valid()
        && result.finalState.acceptedSteps == stepCount
        && result.finalState.rejectedSteps == 0
        && close(result.finalState.heatLedger, totalHeat, 1.0e-12, 1.0e-7)
        && close(
            result.finalState.heatLedger,
            mass*(finalPoint.e - initialEquilibrium.energy)
                + pressure*(finalVolume - box.initialVolume()),
            0.0,
            2.0e-3);
    return result;
}

} // End anonymous namespace

int main()
{
    std::cout << "=== Pure-O2 constant-pressure HEM heating Box A2 ===\n";

    constexpr double pressure = 4.0e6;
    constexpr double mass = 1.0;
    constexpr double initialTemperature = 120.0;
    constexpr double targetQuality = 0.5;

    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB kij = BinaryInteractionDB::createDefault();
    PRPhaseCaloric caloric(critical, kij, {"O2"});
    const PRPhaseCaloricResult initialCaloric = caloric.evaluate(
        pressure, initialTemperature, {1.0}, PhaseRoot::LIQUID);
    check(initialCaloric.success(),
        "4 MPa, 120 K subcooled O2 liquid is constructible");
    if (!initialCaloric.success()) return 1;

    PRTPFlash tp(critical, kij, {"O2"});
    const PureSaturationResult saturation = tp.solvePureSaturationTemperature(
        pressure, 100.0, 180.0);
    check(saturation.valid() && saturation.temperature > initialTemperature,
        "4 MPa saturation lies above the subcooled initial temperature");
    if (!saturation.valid()) return 1;
    const PRPhaseCaloricResult saturatedLiquid = caloric.evaluate(
        pressure, saturation.temperature, {1.0}, PhaseRoot::LIQUID);
    const PRPhaseCaloricResult saturatedVapor = caloric.evaluate(
        pressure, saturation.temperature, {1.0}, PhaseRoot::VAPOR);
    check(saturatedLiquid.success() && saturatedVapor.success()
            && saturatedVapor.state->h > saturatedLiquid.state->h,
        "saturation endpoints provide a positive latent enthalpy");
    if (!saturatedLiquid.success() || !saturatedVapor.success()) return 1;

    PRPHFlash phFlash(critical, kij, {"O2"}, phOptions());
    PRRhoEFlash rhoEFlash(critical, kij, {"O2"}, rhoEOptions());
    const double initialEnthalpy = initialCaloric.state->h;
    const PHFlashResult initialPH = phFlash.flashMassFractions(
        pressure, initialEnthalpy, {1.0});
    const PhasePoint initialPHPoint = pointFromPH(initialPH);
    check(initialPHPoint.valid
            && initialPHPoint.region == EquilibriumRegion::SINGLE_LIQUID
            && initialPHPoint.quality == 0.0,
        "PH initialization identifies a stable single liquid");
    if (!initialPHPoint.valid) return 1;

    const RhoEFlashResult initialEquilibrium =
        rhoEFlash.flashMassFractions(
            initialPHPoint.rho, initialPHPoint.e, {1.0});
    const PhasePoint initialRhoEPoint = pointFromRhoE(initialEquilibrium);
    check(initialRhoEPoint.valid
            && initialRhoEPoint.region == EquilibriumRegion::SINGLE_LIQUID
            && close(initialEquilibrium.pressure, pressure, 2.0e-7, 1.0e-2),
        "authoritative rho-e flash reproduces the subcooled PH state");
    if (!initialRhoEPoint.valid) return 1;

    const double targetEnthalpy = (1.0 - targetQuality)
            *saturatedLiquid.state->h
        + targetQuality*saturatedVapor.state->h;
    const double totalHeat = mass*(targetEnthalpy - initialEnthalpy);
    check(std::isfinite(totalHeat) && totalHeat > 0.0,
        "prescribed heat spans sensible heating plus half the latent interval");

    const RunResult coarse = runHeatingPath(
        phFlash, rhoEFlash, pressure, mass, initialEnthalpy,
        initialPH, initialEquilibrium, totalHeat, 8);
    const RunResult fine = runHeatingPath(
        phFlash, rhoEFlash, pressure, mass, initialEnthalpy,
        initialPH, initialEquilibrium, totalHeat, 32);
    check(coarse.success && fine.success,
        "coarse and fine heat partitions close every PH/rho-e transaction");
    check(fine.crossedBoundary && fine.monotoneQuality,
        "subcooled O2 crosses into coexistence with monotone derived quality");

    const PhasePoint coarseFinal = pointFromRhoE(coarse.finalState.equilibrium);
    const PhasePoint fineFinal = pointFromRhoE(fine.finalState.equilibrium);
    std::cout << "    final p=" << fine.finalState.equilibrium.pressure
              << " T=" << fineFinal.temperature
              << " quality=" << fineFinal.quality
              << " Q=" << fine.finalState.heatLedger
              << " W=" << fine.finalState.boundaryWorkLedger << '\n';
    check(fineFinal.region == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE
            && close(fineFinal.quality, targetQuality, 2.0e-8, 1.0e-8),
        "final rho-e state recovers the prescribed 0.5 mass quality");
    check(close(coarseFinal.temperature, fineFinal.temperature, 2.0e-7)
            && close(coarseFinal.quality, fineFinal.quality, 2.0e-8, 1.0e-8)
            && close(
                coarse.finalState.heatLedger,
                fine.finalState.heatLedger,
                1.0e-12,
                1.0e-7)
            && close(
                coarse.finalState.boundaryWorkLedger,
                fine.finalState.boundaryWorkLedger,
                2.0e-8,
                2.0e-3),
        "final state and first-law ledgers are invariant to heat partition");
    check(fine.finalState.boundaryWorkLedger > 0.0,
        "positive expansion work is separated from internal-energy change");

    bool qualityMatrixClosed = true;
    double previousVolume = mass/initialRhoEPoint.rho;
    for (double quality : {1.0e-6, 0.01, 0.5, 0.99, 1.0 - 1.0e-6})
    {
        const double requestedEnthalpy = (1.0 - quality)
                *saturatedLiquid.state->h
            + quality*saturatedVapor.state->h;
        const RunResult path = runHeatingPath(
            phFlash,
            rhoEFlash,
            pressure,
            mass,
            initialEnthalpy,
            initialPH,
            initialEquilibrium,
            mass*(requestedEnthalpy - initialEnthalpy),
            8);
        const PhasePoint final = pointFromRhoE(path.finalState.equilibrium);
        const double finalVolume = mass/final.rho;
        std::cout << "    boundary quality=" << quality
                  << " recovered=" << final.quality
                  << " V=" << finalVolume << '\n';
        qualityMatrixClosed = qualityMatrixClosed
            && path.success && path.crossedBoundary && path.monotoneQuality
            && final.region == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE
            && close(final.quality, quality, 2.0e-8, 2.0e-8)
            && finalVolume > previousVolume;
        previousVolume = finalVolume;
    }
    check(qualityMatrixClosed,
        "PH-to-rho-e path closes from incipient boiling to the vapor endpoint");

    const PRPhaseCaloricResult superheatedVapor = caloric.evaluate(
        pressure, 160.0, {1.0}, PhaseRoot::VAPOR);
    const RunResult vaporPath = superheatedVapor.success()
        ? runHeatingPath(
            phFlash,
            rhoEFlash,
            pressure,
            mass,
            initialEnthalpy,
            initialPH,
            initialEquilibrium,
            mass*(superheatedVapor.state->h - initialEnthalpy),
            16)
        : RunResult{};
    const PhasePoint vaporFinal = pointFromRhoE(
        vaporPath.finalState.equilibrium);
    check(superheatedVapor.success() && vaporPath.success
            && vaporFinal.region == EquilibriumRegion::SINGLE_VAPOR
            && vaporFinal.quality == 1.0
            && vaporPath.monotoneQuality,
        "continued heating exits coexistence into a single vapor without clipping quality");

    ConstantPressureHeatingBox rollbackBox(
        phFlash,
        rhoEFlash,
        pressure,
        mass,
        initialEnthalpy,
        initialPH,
        initialEquilibrium);
    const PistonSnapshot beforeFailure = rollbackBox.state();
    const PistonStepStatus failed = rollbackBox.addHeat(1.0e12, 1.0);
    const PistonSnapshot afterFailure = rollbackBox.state();
    check(failed == PistonStepStatus::PH_FLASH_FAILED,
        "out-of-envelope enthalpy fails before rho-e publication");
    check(sameCommittedState(beforeFailure, afterFailure)
            && afterFailure.rejectedSteps
                == beforeFailure.rejectedSteps + 1,
        "failed PH/rho-e transaction rolls back state and both energy ledgers");

    RhoEFlashOptions limitedOptions = rhoEOptions();
    limitedOptions.minimumPressure = 3.5e6;
    limitedOptions.maximumPressure = 4.5e6;
    limitedOptions.minimumTemperature = 110.0;
    limitedOptions.maximumTemperature = 130.0;
    limitedOptions.initialTemperature = initialTemperature;
    PRRhoEFlash limitedRhoE(critical, kij, {"O2"}, limitedOptions);
    const RhoEFlashResult limitedInitial = limitedRhoE.flashMassFractions(
        initialPHPoint.rho, initialPHPoint.e, {1.0});
    ConstantPressureHeatingBox rhoERollbackBox(
        phFlash,
        limitedRhoE,
        pressure,
        mass,
        initialEnthalpy,
        initialPH,
        limitedInitial);
    const PistonSnapshot beforeRhoEFailure = rhoERollbackBox.state();
    const PistonStepStatus rhoEFailure = rhoERollbackBox.addHeat(
        totalHeat, 1.0);
    const PistonSnapshot afterRhoEFailure = rhoERollbackBox.state();
    check(limitedInitial.success()
            && rhoEFailure == PistonStepStatus::RHO_E_REFLASH_FAILED,
        "a valid PH candidate is rejected when authoritative rho-e is out of range");
    check(sameCommittedState(beforeRhoEFailure, afterRhoEFailure)
            && afterRhoEFailure.rejectedSteps
                == beforeRhoEFailure.rejectedSteps + 1,
        "rho-e failure rolls back the PH candidate and both energy ledgers");

    if (failures != 0)
    {
        std::cerr << "PR_CONSTANT_PRESSURE_HEATING_FAIL failures="
                  << failures << " checks=" << checks << '\n';
        return 1;
    }

    std::cout << "PR_CONSTANT_PRESSURE_HEATING_PASS checks=" << checks << '\n'
              << "NOTE: pure-O2 piston Box A2 only; external LOX caloric "
              << "validation, phase-aware transport/chemistry, and production "
              << "capability flags remain closed.\n";
    return 0;
}
