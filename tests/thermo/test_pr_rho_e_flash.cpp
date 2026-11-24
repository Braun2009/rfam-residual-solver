#include "PhaseEquilibrium/PRRhoEFlash.H"
#include "PhaseEquilibrium/PRPhaseSoundSpeed.H"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace Foam;
using namespace Foam::phaseEquilibrium;

namespace
{

int failures = 0;

void check(bool condition, const std::string& label)
{
    std::cout << "  " << (condition ? "PASS " : "FAIL ") << label << '\n';
    if (!condition) ++failures;
}

bool close(double lhs, double rhs, double relative = 1.0e-7)
{
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) return false;
    const double scale = std::max({std::abs(lhs), std::abs(rhs), 1.0});
    return std::abs(lhs - rhs) <= relative*scale;
}

ConservedPhaseProperties conserved(const PRPhaseCaloricProperties& phase)
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

RhoEFlashOptions options(
    double minimumPressure,
    double maximumPressure,
    double minimumTemperature,
    double maximumTemperature)
{
    RhoEFlashOptions result;
    result.minimumPressure = minimumPressure;
    result.maximumPressure = maximumPressure;
    result.minimumTemperature = minimumTemperature;
    result.maximumTemperature = maximumTemperature;
    result.pressureScanIntervals = 96;
    result.maximumIterations = 120;
    result.finiteDifferenceStep = 1.0e-4;
    result.logStepLimit = 0.75;
    result.densityAbsoluteTolerance = 1.0e-7;
    result.densityRelativeTolerance = 2.0e-8;
    result.energyAbsoluteTolerance = 1.0e-4;
    result.energyRelativeTolerance = 2.0e-10;
    return result;
}

struct ReferenceState
{
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double density = std::numeric_limits<double>::quiet_NaN();
    double energy = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> massFractions;
};

struct PureIsentropePoint
{
    bool valid = false;
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double density = std::numeric_limits<double>::quiet_NaN();
    double vaporMassFraction = std::numeric_limits<double>::quiet_NaN();
};

struct IndependentAcousticReference
{
    bool valid = false;
    double soundSpeed = std::numeric_limits<double>::quiet_NaN();
    double soundSpeedSquared = std::numeric_limits<double>::quiet_NaN();
    double cv = std::numeric_limits<double>::quiet_NaN();
    double pressureDensityDerivative =
        std::numeric_limits<double>::quiet_NaN();
    double pressureTemperatureDerivative =
        std::numeric_limits<double>::quiet_NaN();
};

struct EquilibriumPTPoint
{
    bool valid = false;
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double temperature = std::numeric_limits<double>::quiet_NaN();
    double density = std::numeric_limits<double>::quiet_NaN();
    double energy = std::numeric_limits<double>::quiet_NaN();
    double entropy = std::numeric_limits<double>::quiet_NaN();
};

ReferenceState pureSaturatedState(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& kij,
    const std::string& species,
    double pressure,
    double beta)
{
    ReferenceState result;
    result.pressure = pressure;
    result.massFractions = {1.0};

    PRTPFlash tp(critical, kij, {species});
    const PureSaturationResult saturation = tp.solvePureSaturationTemperature(
        pressure, 80.0, 500.0);
    if (!saturation.valid()) return result;
    result.temperature = saturation.temperature;

    PRPhaseCaloric caloric(critical, kij, {species});
    const PRPhaseCaloricResult liquid = caloric.evaluate(
        pressure, result.temperature, result.massFractions, PhaseRoot::LIQUID);
    const PRPhaseCaloricResult vapor = caloric.evaluate(
        pressure, result.temperature, result.massFractions, PhaseRoot::VAPOR);
    if (!liquid.success() || !vapor.success()) return ReferenceState{};

    const ConservedMixtureResult mixture = assembleConservedPhaseMixture(
        pressure, beta, conserved(*liquid.state), conserved(*vapor.state),
        result.massFractions);
    if (!mixture.success()) return ReferenceState{};
    result.density = mixture.state->rho;
    result.energy = mixture.state->e;
    return result;
}

ReferenceState pureSinglePhaseState(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& kij,
    const std::string& species,
    double pressure,
    double temperature,
    PhaseRoot root)
{
    ReferenceState result;
    result.pressure = pressure;
    result.temperature = temperature;
    result.massFractions = {1.0};

    PRPhaseCaloric caloric(critical, kij, {species});
    const PRPhaseCaloricResult phase = caloric.evaluate(
        pressure, temperature, result.massFractions, root);
    if (!phase.success()) return ReferenceState{};
    result.density = phase.state->rho;
    result.energy = phase.state->e;
    return result;
}

ReferenceState mixtureTwoPhaseState(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& kij)
{
    ReferenceState result;
    result.pressure = 2.0e6;
    result.temperature = 300.0;
    PRTPFlash tp(critical, kij, {"O2", "C3H8"});
    const TPFlashResult split = tp.flashMoleFractions(
        result.pressure, result.temperature, {0.2, 0.8});
    if (!split.valid() || split.region != EquilibriumRegion::TWO_PHASE)
        return result;

    PRPhaseCaloric caloric(critical, kij, {"O2", "C3H8"});
    const PRPhaseCaloricResult liquid = caloric.evaluate(
        result.pressure, result.temperature, split.liquidMoleFractions,
        PhaseRoot::LIQUID);
    const PRPhaseCaloricResult vapor = caloric.evaluate(
        result.pressure, result.temperature, split.vaporMoleFractions,
        PhaseRoot::VAPOR);
    if (!liquid.success() || !vapor.success()) return ReferenceState{};

    result.massFractions = split.feedMassFractions;
    const ConservedMixtureResult mixture = assembleConservedPhaseMixture(
        result.pressure, split.vaporMoleFraction,
        conserved(*liquid.state), conserved(*vapor.state), result.massFractions);
    if (!mixture.success()) return ReferenceState{};
    result.density = mixture.state->rho;
    result.energy = mixture.state->e;
    return result;
}

PhaseSoundSpeedResult evaluateAcousticAtDelta(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& species,
    const RhoEFlashOptions& rhoEOptions,
    const RhoEFlashResult& base,
    const std::vector<double>& massFractions,
    double relativeDensityPerturbation)
{
    PhaseSoundSpeedOptions acousticOptions;
    acousticOptions.relativeDensityPerturbation = relativeDensityPerturbation;
    PRPhaseSoundSpeed acoustic(
        critical, kij, species, rhoEOptions, acousticOptions);
    return acoustic.evaluate(base, massFractions);
}

std::vector<PhaseSoundSpeedResult> acousticConvergenceSequence(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& species,
    const RhoEFlashOptions& rhoEOptions,
    const RhoEFlashResult& base,
    const std::vector<double>& massFractions,
    const std::vector<double>& densityPerturbations)
{
    std::vector<PhaseSoundSpeedResult> result;
    result.reserve(densityPerturbations.size());
    for (double perturbation : densityPerturbations)
    {
        result.push_back(evaluateAcousticAtDelta(
            critical, kij, species, rhoEOptions, base, massFractions,
            perturbation));
    }
    return result;
}

double relativeDifference(double lhs, double rhs)
{
    if (!std::isfinite(lhs) || !std::isfinite(rhs))
        return std::numeric_limits<double>::infinity();
    return std::abs(lhs - rhs)
        /std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

bool acousticSequenceStable(
    const std::vector<PhaseSoundSpeedResult>& sequence,
    double relativeTolerance)
{
    if (sequence.empty()) return false;
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    for (const PhaseSoundSpeedResult& state : sequence)
    {
        if (!state.success()) return false;
        minimum = std::min(minimum, state.soundSpeed);
        maximum = std::max(maximum, state.soundSpeed);
    }
    return relativeDifference(minimum, maximum) <= relativeTolerance;
}

void printAcousticSequence(
    const std::string& label,
    const std::vector<double>& densityPerturbations,
    const std::vector<PhaseSoundSpeedResult>& sequence)
{
    std::cout << "    " << label << " acoustic convergence";
    for (std::size_t i = 0; i < sequence.size(); ++i)
    {
        std::cout << " delta=" << densityPerturbations[i]
                  << " status=" << phaseSoundSpeedStatusName(sequence[i].status)
                  << " a=" << sequence[i].soundSpeed;
    }
    std::cout << '\n';
}

PureIsentropePoint pureSaturationIsentropePoint(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& kij,
    const std::string& species,
    double pressure,
    double targetEntropy,
    double minimumTemperature,
    double maximumTemperature)
{
    PureIsentropePoint result;
    result.pressure = pressure;

    PRTPFlash tp(critical, kij, {species});
    const PureSaturationResult saturation = tp.solvePureSaturationTemperature(
        pressure, minimumTemperature, maximumTemperature);
    if (!saturation.valid()) return result;

    PRPhaseCaloric caloric(critical, kij, {species});
    const PRPhaseCaloricResult liquid = caloric.evaluate(
        pressure, saturation.temperature, {1.0}, PhaseRoot::LIQUID);
    const PRPhaseCaloricResult vapor = caloric.evaluate(
        pressure, saturation.temperature, {1.0}, PhaseRoot::VAPOR);
    if (!liquid.success() || !vapor.success()) return result;

    const double entropySpan = vapor.state->s - liquid.state->s;
    if (!std::isfinite(entropySpan) || entropySpan <= 0.0) return result;
    const double quality = (targetEntropy - liquid.state->s)/entropySpan;
    if (!std::isfinite(quality) || quality < 0.0 || quality > 1.0)
        return result;

    const ConservedMixtureResult mixture = assembleConservedPhaseMixture(
        pressure, quality, conserved(*liquid.state), conserved(*vapor.state),
        {1.0});
    if (!mixture.success() || !std::isfinite(mixture.state->s)) return result;

    result.valid = true;
    result.density = mixture.state->rho;
    result.vaporMassFraction = mixture.state->vaporMassFraction;
    return result;
}

double independentPureIsentropicSoundSpeed(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& kij,
    const std::string& species,
    double basePressure,
    double targetEntropy,
    double minimumTemperature,
    double maximumTemperature,
    double relativePressurePerturbation)
{
    const PureIsentropePoint lower = pureSaturationIsentropePoint(
        critical, kij, species,
        basePressure*(1.0 - relativePressurePerturbation), targetEntropy,
        minimumTemperature, maximumTemperature);
    const PureIsentropePoint upper = pureSaturationIsentropePoint(
        critical, kij, species,
        basePressure*(1.0 + relativePressurePerturbation), targetEntropy,
        minimumTemperature, maximumTemperature);
    if (!lower.valid || !upper.valid || !(upper.density > lower.density))
        return std::numeric_limits<double>::quiet_NaN();
    const double derivative = (upper.pressure - lower.pressure)
        /(upper.density - lower.density);
    return derivative > 0.0 && std::isfinite(derivative)
        ? std::sqrt(derivative)
        : std::numeric_limits<double>::quiet_NaN();
}

double fivePointDerivative(
    double minusTwo,
    double minusOne,
    double plusOne,
    double plusTwo,
    double step)
{
    if (!std::isfinite(minusTwo) || !std::isfinite(minusOne)
        || !std::isfinite(plusOne) || !std::isfinite(plusTwo)
        || !std::isfinite(step) || !(step > 0.0))
        return std::numeric_limits<double>::quiet_NaN();
    return (minusTwo - 8.0*minusOne + 8.0*plusOne - plusTwo)
        /(12.0*step);
}

IndependentAcousticReference independentSinglePhaseSoundSpeed(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& kij,
    const std::string& species,
    double pressure,
    double temperature,
    PhaseRoot root,
    double relativeStep)
{
    IndependentAcousticReference result;
    PRPhaseCaloric caloric(critical, kij, {species});
    const PRPhaseCaloricResult base = caloric.evaluate(
        pressure, temperature, {1.0}, root);
    if (!base.success()) return result;

    PengRobinsonEOS eos(critical, kij);
    eos.setSpeciesNames({species});
    const std::vector<double> massFractions = {1.0};
    const double density = base.state->rho;
    const double densityStep = relativeStep*density;
    const double temperatureStep = relativeStep*temperature;
    if (!(density - 2.0*densityStep > 0.0)
        || !(temperature - 2.0*temperatureStep > 0.0))
        return result;

    const double pRhoMinusTwo = eos.pressure(
        density - 2.0*densityStep, temperature, massFractions);
    const double pRhoMinusOne = eos.pressure(
        density - densityStep, temperature, massFractions);
    const double pRhoPlusOne = eos.pressure(
        density + densityStep, temperature, massFractions);
    const double pRhoPlusTwo = eos.pressure(
        density + 2.0*densityStep, temperature, massFractions);
    result.pressureDensityDerivative = fivePointDerivative(
        pRhoMinusTwo, pRhoMinusOne, pRhoPlusOne, pRhoPlusTwo,
        densityStep);

    const double pTMinusTwo = eos.pressure(
        density, temperature - 2.0*temperatureStep, massFractions);
    const double pTMinusOne = eos.pressure(
        density, temperature - temperatureStep, massFractions);
    const double pTPlusOne = eos.pressure(
        density, temperature + temperatureStep, massFractions);
    const double pTPlusTwo = eos.pressure(
        density, temperature + 2.0*temperatureStep, massFractions);
    result.pressureTemperatureDerivative = fivePointDerivative(
        pTMinusTwo, pTMinusOne, pTPlusOne, pTPlusTwo,
        temperatureStep);

    auto energyAtFixedDensity = [&](double sampleTemperature)
    {
        const double samplePressure = eos.pressure(
            density, sampleTemperature, massFractions);
        const PRPhaseCaloricResult phase = caloric.evaluate(
            samplePressure, sampleTemperature, {1.0}, root);
        if (!phase.success()
            || relativeDifference(phase.state->rho, density) > 1.0e-8)
            return std::numeric_limits<double>::quiet_NaN();
        return phase.state->e;
    };
    result.cv = fivePointDerivative(
        energyAtFixedDensity(temperature - 2.0*temperatureStep),
        energyAtFixedDensity(temperature - temperatureStep),
        energyAtFixedDensity(temperature + temperatureStep),
        energyAtFixedDensity(temperature + 2.0*temperatureStep),
        temperatureStep);

    result.soundSpeedSquared = result.pressureDensityDerivative
        + temperature*result.pressureTemperatureDerivative
            *result.pressureTemperatureDerivative
            /(density*density*result.cv);
    if (!std::isfinite(result.pressureDensityDerivative)
        || !std::isfinite(result.pressureTemperatureDerivative)
        || !std::isfinite(result.cv) || result.cv <= 0.0
        || !std::isfinite(result.soundSpeedSquared)
        || result.soundSpeedSquared <= 0.0)
        return result;
    result.soundSpeed = std::sqrt(result.soundSpeedSquared);
    result.valid = true;
    return result;
}

EquilibriumPTPoint equilibriumMixtureAtPT(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& kij,
    double pressure,
    double temperature,
    const std::vector<double>& massFractions)
{
    EquilibriumPTPoint result;
    result.pressure = pressure;
    result.temperature = temperature;
    PRTPFlash tp(critical, kij, {"O2", "C3H8"});
    const TPFlashResult split = tp.flashMassFractions(
        pressure, temperature, massFractions);
    if (!split.valid() || split.region != EquilibriumRegion::TWO_PHASE
        || !split.phaseSplitReady())
        return result;

    PRPhaseCaloric caloric(critical, kij, {"O2", "C3H8"});
    const PRPhaseCaloricResult liquid = caloric.evaluate(
        pressure, temperature, split.liquidMoleFractions,
        PhaseRoot::LIQUID);
    const PRPhaseCaloricResult vapor = caloric.evaluate(
        pressure, temperature, split.vaporMoleFractions,
        PhaseRoot::VAPOR);
    if (!liquid.success() || !vapor.success()) return result;
    const ConservedMixtureResult mixture = assembleConservedPhaseMixture(
        pressure, split.vaporMoleFraction,
        conserved(*liquid.state), conserved(*vapor.state), massFractions);
    if (!mixture.success() || !std::isfinite(mixture.state->s)) return result;

    result.valid = true;
    result.density = mixture.state->rho;
    result.energy = mixture.state->e;
    result.entropy = mixture.state->s;
    return result;
}

EquilibriumPTPoint solveMixtureEntropyAtPressure(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& kij,
    double pressure,
    double targetEntropy,
    const std::vector<double>& massFractions,
    double temperatureSeed)
{
    const double entropyTolerance = 2.0e-7;
    EquilibriumPTPoint center = equilibriumMixtureAtPT(
        critical, kij, pressure, temperatureSeed, massFractions);
    if (!center.valid) return EquilibriumPTPoint{};
    const double centerResidual = center.entropy - targetEntropy;
    if (std::abs(centerResidual) <= entropyTolerance) return center;

    EquilibriumPTPoint lower;
    EquilibriumPTPoint upper;
    double lowerResidual = std::numeric_limits<double>::quiet_NaN();
    double upperResidual = std::numeric_limits<double>::quiet_NaN();
    double temperatureSpan = 0.02;
    bool bracketed = false;
    for (int expansion = 0; expansion < 7; ++expansion)
    {
        lower = equilibriumMixtureAtPT(
            critical, kij, pressure, temperatureSeed - temperatureSpan,
            massFractions);
        upper = equilibriumMixtureAtPT(
            critical, kij, pressure, temperatureSeed + temperatureSpan,
            massFractions);
        if (!lower.valid || !upper.valid) return EquilibriumPTPoint{};
        lowerResidual = lower.entropy - targetEntropy;
        upperResidual = upper.entropy - targetEntropy;
        bracketed = (lowerResidual < 0.0 && upperResidual > 0.0)
            || (lowerResidual > 0.0 && upperResidual < 0.0);
        if (bracketed) break;
        temperatureSpan *= 2.0;
    }
    if (!bracketed) return EquilibriumPTPoint{};

    for (int iteration = 0; iteration < 60; ++iteration)
    {
        const double middleTemperature =
            0.5*(lower.temperature + upper.temperature);
        EquilibriumPTPoint middle = equilibriumMixtureAtPT(
            critical, kij, pressure, middleTemperature, massFractions);
        if (!middle.valid) return EquilibriumPTPoint{};
        const double middleResidual = middle.entropy - targetEntropy;
        if (std::abs(middleResidual) <= entropyTolerance) return middle;
        const bool lowerToMiddle =
            (lowerResidual < 0.0 && middleResidual > 0.0)
            || (lowerResidual > 0.0 && middleResidual < 0.0);
        if (lowerToMiddle)
        {
            upper = middle;
            upperResidual = middleResidual;
        }
        else
        {
            lower = middle;
            lowerResidual = middleResidual;
        }
    }
    return EquilibriumPTPoint{};
}

IndependentAcousticReference independentMixtureSoundSpeed(
    const SpeciesCriticalDB& critical,
    const BinaryInteractionDB& kij,
    double basePressure,
    double baseTemperature,
    double targetEntropy,
    const std::vector<double>& massFractions,
    double relativePressurePerturbation)
{
    IndependentAcousticReference result;
    const EquilibriumPTPoint lower = solveMixtureEntropyAtPressure(
        critical, kij,
        basePressure*(1.0 - relativePressurePerturbation), targetEntropy,
        massFractions, baseTemperature);
    const EquilibriumPTPoint upper = solveMixtureEntropyAtPressure(
        critical, kij,
        basePressure*(1.0 + relativePressurePerturbation), targetEntropy,
        massFractions, baseTemperature);
    if (!lower.valid || !upper.valid || !(upper.density > lower.density))
        return result;
    result.soundSpeedSquared = (upper.pressure - lower.pressure)
        /(upper.density - lower.density);
    if (!std::isfinite(result.soundSpeedSquared)
        || result.soundSpeedSquared <= 0.0)
        return result;
    result.soundSpeed = std::sqrt(result.soundSpeedSquared);
    result.valid = true;
    return result;
}

void printFailure(const RhoEFlashResult& result)
{
    std::cout << "    status=" << rhoEFlashStatusName(result.status)
              << " phase=" << statusName(result.phaseStatus)
              << " caloric=" << phaseCaloricStatusName(result.caloricStatus)
              << " p=" << result.pressure
              << " T=" << result.temperature
              << " rho=" << result.density
              << " e=" << result.energy
              << " beta=" << result.vaporMoleFraction
              << " rhoResidual=" << result.densityResidual
              << " eResidual=" << result.energyResidual
              << " validSamples=" << result.validSamples
              << " invalidSamples=" << result.invalidSamples << '\n';
}

void printAcousticFailure(const PhaseSoundSpeedResult& result)
{
    std::cout << "    acoustic status="
              << phaseSoundSpeedStatusName(result.status)
              << " base=" << rhoEFlashStatusName(result.baseFlashStatus)
              << " lower=" << rhoEFlashStatusName(result.lowerFlashStatus)
              << " upper=" << rhoEFlashStatusName(result.upperFlashStatus)
              << " targetS=" << result.targetEntropy
              << " lowerSres=" << result.lowerEntropyResidual
              << " upperSres=" << result.upperEntropyResidual
              << " lowerSamples=" << result.lowerValidSamples
              << " upperSamples=" << result.upperValidSamples
              << " lowerP=" << result.lowerPressure
              << " upperP=" << result.upperPressure;
    if (result.baseState)
        std::cout << " baseP=" << result.baseState->pressure
                  << " baseT=" << result.baseState->temperature
                  << " baseRho=" << result.baseState->density
                  << " baseE=" << result.baseState->energy;
    std::cout << '\n';
}

} // End anonymous namespace

int main()
{
    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB kij = BinaryInteractionDB::createDefault();
    const std::vector<double> densityPerturbations = {
        1.0e-3, 3.0e-4, 1.0e-4, 3.0e-5};

    std::cout << "=== Strict Peng-Robinson rho-E flash ===\n";

    const ReferenceState propane = pureSaturatedState(
        critical, kij, "C3H8", 1.0e6, 0.37);
    check(std::isfinite(propane.temperature) && std::isfinite(propane.density)
            && std::isfinite(propane.energy),
        "pure C3H8 saturation reference is constructible");
    if (std::isfinite(propane.density))
    {
        PRRhoEFlash solver(critical, kij, {"C3H8"},
            options(5.0e5, 1.5e6, 280.0, 360.0));
        const RhoEFlashResult recovered = solver.flashMassFractions(
            propane.density, propane.energy, propane.massFractions);
        if (!recovered.success()) printFailure(recovered);
        check(recovered.success(),
            "pure C3H8 saturated rho-E inversion succeeds");
        check(recovered.pureVolumeEnergyLeverRule
                && recovered.liquidPhase.has_value()
                && recovered.vaporPhase.has_value()
                && recovered.conservedMixture.has_value(),
            "pure rho-E result exposes both phases and the lever rule");
        check(close(recovered.pressure, propane.pressure, 2.0e-7)
                && close(recovered.temperature, propane.temperature, 2.0e-7),
            "pure rho-E recovers saturation pressure and temperature");
        check(recovered.densityRelativeResidual < 2.0e-8
                && recovered.energyRelativeResidual < 2.0e-9,
            "pure rho-E closes density and energy residuals");
        PRPhaseSoundSpeed acoustic(
            critical, kij, {"C3H8"}, options(5.0e5, 1.5e6, 280.0, 360.0));
        const PhaseSoundSpeedResult acousticState = acoustic.evaluate(
            recovered, propane.massFractions);
        if (!acousticState.success()) printAcousticFailure(acousticState);
        check(acousticState.success(),
            "pure two-phase HEM sound-speed closure succeeds");
        const bool pureAcousticResiduals =
            std::abs(acousticState.lowerEntropyResidual) < 2.0e-6
            && std::abs(acousticState.upperEntropyResidual) < 2.0e-6
            && acousticState.soundSpeedSquared > 0.0;
        if (!pureAcousticResiduals) printAcousticFailure(acousticState);
        check(pureAcousticResiduals,
            "pure HEM acoustic samples close entropy and positivity");
        const std::vector<PhaseSoundSpeedResult> propaneSequence =
            acousticConvergenceSequence(
                critical, kij, {"C3H8"},
                options(5.0e5, 1.5e6, 280.0, 360.0), recovered,
                propane.massFractions, densityPerturbations);
        printAcousticSequence(
            "C3H8", densityPerturbations, propaneSequence);
        check(std::all_of(
                propaneSequence.begin(), propaneSequence.end(),
                [](const PhaseSoundSpeedResult& state)
                {
                    return state.success();
                }),
            "pure C3H8 HEM acoustics close over density-step refinement");
        check(acousticSequenceStable(propaneSequence, 5.0e-5),
            "pure C3H8 HEM sound speed converges with density-step refinement");
        const double propaneIndependentSoundSpeed =
            independentPureIsentropicSoundSpeed(
                critical, kij, "C3H8", recovered.pressure,
                recovered.entropy, 280.0, 360.0, 3.0e-5);
        std::cout << "    C3H8 independent saturation-isentrope a="
                  << propaneIndependentSoundSpeed << '\n';
        check(std::isfinite(propaneIndependentSoundSpeed)
                && propaneIndependentSoundSpeed > 0.0,
            "pure C3H8 independent saturation isentrope is constructible");
        check(relativeDifference(
                    propaneSequence[2].soundSpeed,
                    propaneIndependentSoundSpeed) < 5.0e-5,
            "pure C3H8 HEM sound speed matches the independent isentrope");
        if (recovered.conservedMixture)
        {
            check(close(recovered.conservedMixture->e,
                        recovered.conservedMixture->h
                            - recovered.pressure/recovered.conservedMixture->rho,
                        1.0e-10),
                "pure rho-E conserved state closes e=h-p/rho");
        }
    }

    const ReferenceState oxygen = pureSaturatedState(
        critical, kij, "O2", 4.0e6, 0.41);
    check(std::isfinite(oxygen.temperature) && std::isfinite(oxygen.density)
            && std::isfinite(oxygen.energy),
        "low-temperature O2 saturation reference is constructible");
    if (std::isfinite(oxygen.density))
    {
        PRRhoEFlash solver(critical, kij, {"O2"},
            options(2.0e6, 5.0e6, 100.0, 220.0));
        const RhoEFlashResult recovered = solver.flashMassFractions(
            oxygen.density, oxygen.energy, oxygen.massFractions);
        if (!recovered.success()) printFailure(recovered);
        check(recovered.success(),
            "low-temperature O2 saturated rho-E inversion succeeds");
        check(recovered.densityRelativeResidual < 2.0e-8
                && recovered.energyRelativeResidual < 2.0e-9,
            "O2 rho-E closes density and energy residuals");
        PRPhaseSoundSpeed acoustic(
            critical, kij, {"O2"}, options(2.0e6, 5.0e6, 100.0, 220.0));
        const PhaseSoundSpeedResult acousticState = acoustic.evaluate(
            recovered, oxygen.massFractions);
        if (!acousticState.success()) printAcousticFailure(acousticState);
        check(acousticState.success(),
            "low-temperature O2 HEM sound-speed closure succeeds");
        const std::vector<PhaseSoundSpeedResult> oxygenSequence =
            acousticConvergenceSequence(
                critical, kij, {"O2"},
                options(2.0e6, 5.0e6, 100.0, 220.0), recovered,
                oxygen.massFractions, densityPerturbations);
        printAcousticSequence("O2", densityPerturbations, oxygenSequence);
        check(std::all_of(
                oxygenSequence.begin(), oxygenSequence.end(),
                [](const PhaseSoundSpeedResult& state)
                {
                    return state.success();
                }),
            "low-temperature O2 HEM acoustics close over density-step refinement");
        check(acousticSequenceStable(oxygenSequence, 5.0e-5),
            "low-temperature O2 HEM sound speed converges with density-step refinement");
        const double oxygenIndependentSoundSpeed =
            independentPureIsentropicSoundSpeed(
                critical, kij, "O2", recovered.pressure,
                recovered.entropy, 100.0, 220.0, 3.0e-5);
        std::cout << "    O2 independent saturation-isentrope a="
                  << oxygenIndependentSoundSpeed << '\n';
        check(std::isfinite(oxygenIndependentSoundSpeed)
                && oxygenIndependentSoundSpeed > 0.0,
            "low-temperature O2 independent saturation isentrope is constructible");
        check(relativeDifference(
                    oxygenSequence[2].soundSpeed,
                    oxygenIndependentSoundSpeed) < 5.0e-5,
            "low-temperature O2 HEM sound speed matches the independent isentrope");

        PhaseSoundSpeedOptions boundaryAcousticOptions;
        boundaryAcousticOptions.relativeDensityPerturbation = 1.0e-3;
        boundaryAcousticOptions.minimumRelativeDensityPerturbation = 1.0e-8;
        boundaryAcousticOptions.maximumTopologyRefinements = 20;
        PRPhaseSoundSpeed boundaryAcoustic(
            critical, kij, {"O2"},
            options(2.0e6, 5.0e6, 100.0, 220.0),
            boundaryAcousticOptions);

        const std::vector<double> boundaryQualities = {
            1.0e-3, 1.0e-2, 0.5, 0.99, 0.999};
        bool boundaryQualityClosure = true;
        bool boundaryReferenceClosure = true;
        bool warmSeedClosure = true;
        std::vector<double> boundaryQualitySoundSpeeds;
        for (double quality : boundaryQualities)
        {
            const ReferenceState reference = pureSaturatedState(
                critical, kij, "O2", 4.0e6, quality);
            const RhoEFlashResult boundaryBase = solver.flashMassFractions(
                reference.density, reference.energy,
                reference.massFractions);
            if (!boundaryBase.success()) printFailure(boundaryBase);
            const RhoEFlashResult warmBoundaryBase =
                solver.flashMassFractions(
                    reference.density, reference.energy,
                    reference.massFractions,
                    {reference.pressure, reference.temperature, true,
                     EquilibriumRegion::PURE_COMPONENT_COEXISTENCE});
            warmSeedClosure = warmSeedClosure
                && warmBoundaryBase.success()
                && close(warmBoundaryBase.pressure, boundaryBase.pressure)
                && close(warmBoundaryBase.temperature, boundaryBase.temperature)
                && warmBoundaryBase.validSamples < boundaryBase.validSamples;
            const PhaseSoundSpeedResult boundaryState =
                boundaryAcoustic.evaluate(
                    boundaryBase, reference.massFractions);
            if (!boundaryState.success()) printAcousticFailure(boundaryState);
            std::cout << "    O2 coexistence quality=" << quality
                      << " a=" << boundaryState.soundSpeed
                      << " requestedDelta="
                      << boundaryState.requestedRelativeDensityPerturbation
                      << " usedDelta="
                      << boundaryState.usedRelativeDensityPerturbation
                      << " stencilRefinements="
                      << boundaryState.stencilRefinements
                      << " topologyRefinements="
                      << boundaryState.topologyRefinements << '\n';
            boundaryQualityClosure = boundaryQualityClosure
                && boundaryBase.success() && boundaryState.success()
                && boundaryState.baseRegion
                    == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE
                && boundaryState.lowerRegion == boundaryState.baseRegion
                && boundaryState.upperRegion == boundaryState.baseRegion;
            boundaryQualitySoundSpeeds.push_back(boundaryState.soundSpeed);

            if (quality == 1.0e-2 || quality == 0.5 || quality == 0.99)
            {
                const double independent =
                    independentPureIsentropicSoundSpeed(
                        critical, kij, "O2", boundaryBase.pressure,
                        boundaryBase.entropy, 100.0, 220.0, 3.0e-5);
                std::cout << "      independent a=" << independent << '\n';
                boundaryReferenceClosure = boundaryReferenceClosure
                    && std::isfinite(independent)
                    && relativeDifference(
                        boundaryState.soundSpeed, independent) < 2.0e-4;
            }
        }
        check(boundaryQualityClosure,
            "O2 HEM acoustics retain two-phase topology toward both quality limits");
        check(boundaryReferenceClosure,
            "O2 HEM quality sweep matches independent saturation isentropes");
        check(warmSeedClosure,
            "pure rho-E warm seeds preserve the cold solution with fewer samples");
        check(relativeDifference(
                    boundaryQualitySoundSpeeds[0],
                    boundaryQualitySoundSpeeds[1]) < 3.0e-2
                && relativeDifference(
                    boundaryQualitySoundSpeeds[3],
                    boundaryQualitySoundSpeeds[4]) < 1.0e-2,
            "O2 two-phase acoustics approach finite liquid and vapor quality limits");

        bool exactBoundaryRejected = true;
        for (double quality : {0.0, 1.0})
        {
            const ReferenceState reference = pureSaturatedState(
                critical, kij, "O2", 4.0e6, quality);
            const RhoEFlashResult boundaryBase = solver.flashMassFractions(
                reference.density, reference.energy,
                reference.massFractions);
            const PhaseSoundSpeedResult boundaryState =
                boundaryAcoustic.evaluate(
                    boundaryBase, reference.massFractions);
            std::cout << "    O2 exact boundary quality=" << quality
                      << " acousticStatus="
                      << phaseSoundSpeedStatusName(boundaryState.status)
                      << " usedDelta="
                      << boundaryState.usedRelativeDensityPerturbation
                      << '\n';
            exactBoundaryRejected = exactBoundaryRejected
                && boundaryBase.success() && !boundaryState.success()
                && boundaryState.status
                    == PhaseSoundSpeedStatus::PHASE_TOPOLOGY_CHANGE;
        }
        check(exactBoundaryRejected,
            "HEM acoustics reject a centered derivative at exact phase boundaries");

        bool endpointSeedInvariant = true;
        for (double quality : {1.0 - 1.0e-6})
        {
            const ReferenceState reference = pureSaturatedState(
                critical, kij, "O2", 4.0e6, quality);
            const RhoEFlashResult cold = solver.flashMassFractions(
                reference.density, reference.energy,
                reference.massFractions);
            const RhoEFlashResult distantSeed = solver.flashMassFractions(
                reference.density, reference.energy,
                reference.massFractions, {4.8e6, 110.0});
            endpointSeedInvariant = endpointSeedInvariant
                && cold.success() && distantSeed.success()
                && close(cold.pressure, distantSeed.pressure, 2.0e-7)
                && close(cold.temperature, distantSeed.temperature, 2.0e-7)
                && cold.tpFlash.region == distantSeed.tpFlash.region;
        }
        check(endpointSeedInvariant,
            "pure rho-E near-endpoint solution is invariant to distant seed hints");

        const std::vector<double> saturationOffsets = {
            1.0, 0.3, 0.1, 0.03};
        bool liquidLimitClosure = true;
        bool vaporLimitClosure = true;
        bool boundaryStencilRefined = false;
        std::vector<double> liquidLimitSoundSpeeds;
        std::vector<double> vaporLimitSoundSpeeds;
        for (double offset : saturationOffsets)
        {
            const ReferenceState liquidReference = pureSinglePhaseState(
                critical, kij, "O2", 4.0e6, oxygen.temperature - offset,
                PhaseRoot::LIQUID);
            const RhoEFlashResult liquidBase = solver.flashMassFractions(
                liquidReference.density, liquidReference.energy,
                liquidReference.massFractions);
            if (!liquidBase.success()) printFailure(liquidBase);
            const PhaseSoundSpeedResult liquidAcoustic =
                boundaryAcoustic.evaluate(
                    liquidBase, liquidReference.massFractions);
            if (!liquidAcoustic.success()) printAcousticFailure(liquidAcoustic);
            std::cout << "    O2 liquid-side dT=" << offset
                      << " a=" << liquidAcoustic.soundSpeed
                      << " usedDelta="
                      << liquidAcoustic.usedRelativeDensityPerturbation
                      << " stencilRefinements="
                      << liquidAcoustic.stencilRefinements
                      << " topologyRefinements="
                      << liquidAcoustic.topologyRefinements << '\n';
            liquidLimitClosure = liquidLimitClosure
                && liquidBase.success() && liquidBase.singlePhase.has_value()
                && liquidBase.tpFlash.region == EquilibriumRegion::SINGLE_LIQUID
                && liquidAcoustic.success()
                && liquidAcoustic.lowerRegion == EquilibriumRegion::SINGLE_LIQUID
                && liquidAcoustic.upperRegion == EquilibriumRegion::SINGLE_LIQUID;
            liquidLimitSoundSpeeds.push_back(liquidAcoustic.soundSpeed);
            boundaryStencilRefined = boundaryStencilRefined
                || liquidAcoustic.stencilRefinements > 0;

            const ReferenceState vaporReference = pureSinglePhaseState(
                critical, kij, "O2", 4.0e6, oxygen.temperature + offset,
                PhaseRoot::VAPOR);
            const RhoEFlashResult vaporBase = solver.flashMassFractions(
                vaporReference.density, vaporReference.energy,
                vaporReference.massFractions);
            if (!vaporBase.success()) printFailure(vaporBase);
            const PhaseSoundSpeedResult vaporAcoustic =
                boundaryAcoustic.evaluate(
                    vaporBase, vaporReference.massFractions);
            if (!vaporAcoustic.success()) printAcousticFailure(vaporAcoustic);
            std::cout << "    O2 vapor-side dT=" << offset
                      << " a=" << vaporAcoustic.soundSpeed
                      << " usedDelta="
                      << vaporAcoustic.usedRelativeDensityPerturbation
                      << " stencilRefinements="
                      << vaporAcoustic.stencilRefinements
                      << " topologyRefinements="
                      << vaporAcoustic.topologyRefinements << '\n';
            vaporLimitClosure = vaporLimitClosure
                && vaporBase.success() && vaporBase.singlePhase.has_value()
                && vaporBase.tpFlash.region == EquilibriumRegion::SINGLE_VAPOR
                && vaporAcoustic.success()
                && vaporAcoustic.lowerRegion == EquilibriumRegion::SINGLE_VAPOR
                && vaporAcoustic.upperRegion == EquilibriumRegion::SINGLE_VAPOR;
            vaporLimitSoundSpeeds.push_back(vaporAcoustic.soundSpeed);
            boundaryStencilRefined = boundaryStencilRefined
                || vaporAcoustic.stencilRefinements > 0;
        }
        check(liquidLimitClosure,
            "O2 liquid-side acoustic limit remains on the liquid branch");
        check(vaporLimitClosure,
            "O2 vapor-side acoustic limit remains on the vapor branch");
        check(boundaryStencilRefined,
            "O2 boundary acoustics refine the stencil instead of mixing phase branches");
        check(relativeDifference(
                    liquidLimitSoundSpeeds[2], liquidLimitSoundSpeeds[3])
                    < 1.0e-2
                && relativeDifference(
                    vaporLimitSoundSpeeds[2], vaporLimitSoundSpeeds[3])
                    < 1.0e-2,
            "O2 single-phase acoustics approach finite one-sided limits");
        check(relativeDifference(
                    boundaryQualitySoundSpeeds.front(),
                    liquidLimitSoundSpeeds.back()) > 0.1
                && relativeDifference(
                    boundaryQualitySoundSpeeds.back(),
                    vaporLimitSoundSpeeds.back()) > 0.1,
            "O2 boundary limits remain distinct instead of being averaged");
    }

    PRPhaseCaloric caloric(critical, kij, {"C3H8"});
    const PRPhaseCaloricResult singleReference = caloric.evaluate(
        1.0e6, 433.7, {1.0}, PhaseRoot::VAPOR);
    check(singleReference.success(),
        "single-phase C3H8 rho-E reference is constructible");
    if (singleReference.success())
    {
        PRRhoEFlash solver(critical, kij, {"C3H8"},
            options(5.0e5, 1.5e6, 400.0, 460.0));
        const RhoEFlashResult recovered = solver.flashMassFractions(
            singleReference.state->rho, singleReference.state->e, {1.0});
        if (!recovered.success()) printFailure(recovered);
        check(recovered.success() && recovered.singlePhase.has_value(),
            "single-phase rho-E inversion succeeds without a phase mixture");
        check(!recovered.liquidPhase.has_value()
                && !recovered.vaporPhase.has_value()
                && !recovered.conservedMixture.has_value(),
            "single-phase rho-E result does not invent two phases");
        check(close(recovered.pressure, 1.0e6, 2.0e-7)
                && close(recovered.temperature, 433.7, 2.0e-7),
            "single-phase rho-E recovers pressure and temperature");
    }

    const PRPhaseCaloricResult analyticReference = caloric.evaluate(
        4.0e6, 450.0, {1.0}, PhaseRoot::VAPOR);
    check(analyticReference.success(),
        "single-phase analytic acoustic reference is constructible");
    if (analyticReference.success())
    {
        const RhoEFlashOptions analyticRhoEOptions =
            options(2.0e6, 6.0e6, 400.0, 500.0);
        PRRhoEFlash solver(
            critical, kij, {"C3H8"}, analyticRhoEOptions);
        const RhoEFlashResult recovered = solver.flashMassFractions(
            analyticReference.state->rho, analyticReference.state->e,
            {1.0}, {4.0e6, 450.0, true,
                    EquilibriumRegion::SINGLE_VAPOR});
        if (!recovered.success()) printFailure(recovered);
        check(recovered.success() && recovered.singlePhase.has_value(),
            "single-phase analytic acoustic state recovers through rho-E");
        const RhoEFlashResult wrongRegion = solver.flashMassFractions(
            analyticReference.state->rho, analyticReference.state->e,
            {1.0}, {4.0e6, 450.0, true,
                    EquilibriumRegion::SINGLE_LIQUID});
        check(wrongRegion.status
                == RhoEFlashStatus::REQUIRED_REGION_MISMATCH,
            "strict local rho-E reports a required phase-region mismatch");

        const PhaseSoundSpeedResult acoustic = evaluateAcousticAtDelta(
            critical, kij, {"C3H8"}, analyticRhoEOptions,
            recovered, {1.0}, 1.0e-4);
        if (!acoustic.success()) printAcousticFailure(acoustic);
        const IndependentAcousticReference coarse =
            independentSinglePhaseSoundSpeed(
                critical, kij, "C3H8", 4.0e6, 450.0,
                PhaseRoot::VAPOR, 1.0e-4);
        const IndependentAcousticReference fine =
            independentSinglePhaseSoundSpeed(
                critical, kij, "C3H8", 4.0e6, 450.0,
                PhaseRoot::VAPOR, 3.0e-5);
        std::cout << "    C3H8 single-phase acoustic="
                  << acoustic.soundSpeed
                  << " independentCoarse=" << coarse.soundSpeed
                  << " independentFine=" << fine.soundSpeed
                  << " Cv=" << fine.cv << '\n';
        check(coarse.valid && fine.valid
                && relativeDifference(
                    coarse.soundSpeed, fine.soundSpeed) < 2.0e-5,
            "single-phase analytic sound-speed reference converges");
        check(acoustic.success()
                && relativeDifference(
                    acoustic.soundSpeed, fine.soundSpeed) < 5.0e-5,
            "single-phase HEM sound speed matches the analytic identity");
        // Frozen corrected-PR value from the 2026-07-17 Vortex audit.
        check(std::abs(fine.soundSpeed - 258.281279994) < 2.0e-3,
            "single-phase analytic sound speed matches the frozen PR reference");
    }

    const ReferenceState mixture = mixtureTwoPhaseState(critical, kij);
    check(std::isfinite(mixture.density) && std::isfinite(mixture.energy)
            && mixture.massFractions.size() == 2,
        "O2/C3H8 two-phase rho-E reference is constructible");
    if (std::isfinite(mixture.density))
    {
        PRRhoEFlash solver(critical, kij, {"O2", "C3H8"},
            options(1.0e6, 3.0e6, 280.0, 320.0));
        const RhoEFlashResult recovered = solver.flashMassFractions(
            mixture.density, mixture.energy, mixture.massFractions);
        if (!recovered.success()) printFailure(recovered);
        const RhoEFlashResult distantSeedRecovered =
            solver.flashMassFractions(
                mixture.density, mixture.energy, mixture.massFractions,
                {2.8e6, 315.0});
        check(recovered.success()
                && recovered.liquidPhase.has_value()
                && recovered.vaporPhase.has_value()
                && recovered.conservedMixture.has_value(),
            "O2/C3H8 two-phase rho-E inversion succeeds");
        check(distantSeedRecovered.success()
                && close(
                    recovered.pressure, distantSeedRecovered.pressure,
                    1.0e-10)
                && close(
                    recovered.temperature, distantSeedRecovered.temperature,
                    1.0e-10)
                && recovered.tpFlash.region
                    == distantSeedRecovered.tpFlash.region,
            "multicomponent rho-E root is invariant to nonlocal seed hints");
        check(close(recovered.pressure, mixture.pressure, 2.0e-6)
                && close(recovered.temperature, mixture.temperature, 2.0e-6),
            "mixture rho-E recovers pressure and temperature");
        check(recovered.tpFugacityResidual
                <= solver.tpFlash().options().fugacityTolerance
                && recovered.tpMaterialResidual
                    <= solver.tpFlash().options().materialTolerance,
            "mixture rho-E preserves TP fugacity and material closure");
        PRPhaseSoundSpeed acoustic(
            critical, kij, {"O2", "C3H8"},
            options(1.0e6, 3.0e6, 280.0, 320.0));
        const PhaseSoundSpeedResult acousticState = acoustic.evaluate(
            recovered, mixture.massFractions);
        if (!acousticState.success()) printAcousticFailure(acousticState);
        check(acousticState.success(),
            "O2/C3H8 two-phase HEM sound-speed closure succeeds");
        const bool mixtureAcousticResiduals =
            std::abs(acousticState.lowerEntropyResidual) < 2.0e-6
            && std::abs(acousticState.upperEntropyResidual) < 2.0e-6;
        if (!mixtureAcousticResiduals) printAcousticFailure(acousticState);
        check(mixtureAcousticResiduals,
            "mixture HEM acoustic samples close entropy");
        const std::vector<PhaseSoundSpeedResult> mixtureSequence =
            acousticConvergenceSequence(
                critical, kij, {"O2", "C3H8"},
                options(1.0e6, 3.0e6, 280.0, 320.0), recovered,
                mixture.massFractions, densityPerturbations);
        printAcousticSequence(
            "O2/C3H8", densityPerturbations, mixtureSequence);
        check(std::all_of(
                mixtureSequence.begin(), mixtureSequence.end(),
                [](const PhaseSoundSpeedResult& state)
                {
                    return state.success();
                }),
            "O2/C3H8 HEM acoustics close over density-step refinement");
        check(acousticSequenceStable(mixtureSequence, 5.0e-5),
            "O2/C3H8 HEM sound speed converges with density-step refinement");
        check(std::all_of(
                mixtureSequence.begin(),mixtureSequence.end(),
                [](const PhaseSoundSpeedResult& state)
                {
                    return state.derivativeConvergenceComparisons>=2
                        && std::isfinite(state.derivativeRelativeChange)
                        && state.derivativeRelativeChange<=2.0e-4;
                }),
            "every successful HEM result proves internal perturbation convergence");

        PhaseSoundSpeedOptions impossibleConvergence;
        impossibleConvergence.maximumTopologyRefinements=2;
        impossibleConvergence.minimumDerivativeConvergenceComparisons=2;
        impossibleConvergence.derivativeRelativeTolerance=1.0e-16;
        PRPhaseSoundSpeed convergenceGate(
            critical,kij,{"O2","C3H8"},
            options(1.0e6,3.0e6,280.0,320.0),impossibleConvergence);
        const PhaseSoundSpeedResult rejectedUnconverged=
            convergenceGate.evaluate(recovered,mixture.massFractions);
        check(rejectedUnconverged.status
                ==PhaseSoundSpeedStatus::PERTURBATION_NOT_CONVERGED
                && !rejectedUnconverged.success(),
            "HEM acoustics fail closed when perturbation convergence is unresolved");

        PhaseSoundSpeedOptions invalidConvergence=impossibleConvergence;
        invalidConvergence.minimumDerivativeConvergenceComparisons=3;
        PRPhaseSoundSpeed invalidConvergenceGate(
            critical,kij,{"O2","C3H8"},
            options(1.0e6,3.0e6,280.0,320.0),invalidConvergence);
        check(invalidConvergenceGate.configurationStatus()
                ==PhaseSoundSpeedStatus::INVALID_CONFIGURATION,
            "acoustic options reject an impossible convergence matrix");

        const EquilibriumPTPoint independentBase = equilibriumMixtureAtPT(
            critical, kij, mixture.pressure, mixture.temperature,
            mixture.massFractions);
        const IndependentAcousticReference mixtureReferenceCoarse =
            independentMixtureSoundSpeed(
                critical, kij, mixture.pressure, mixture.temperature,
                independentBase.entropy, mixture.massFractions, 1.0e-4);
        const IndependentAcousticReference mixtureReferenceFine =
            independentMixtureSoundSpeed(
                critical, kij, mixture.pressure, mixture.temperature,
                independentBase.entropy, mixture.massFractions, 3.0e-5);
        std::cout << "    O2/C3H8 independent TP-isentrope coarse a="
                  << mixtureReferenceCoarse.soundSpeed
                  << " fine a=" << mixtureReferenceFine.soundSpeed << '\n';
        check(independentBase.valid
                && relativeDifference(
                    independentBase.entropy, recovered.entropy) < 1.0e-9,
            "mixture rho-E entropy matches independent TP assembly");
        check(mixtureReferenceCoarse.valid && mixtureReferenceFine.valid
                && relativeDifference(
                    mixtureReferenceCoarse.soundSpeed,
                    mixtureReferenceFine.soundSpeed) < 5.0e-5,
            "O2/C3H8 independent TP-isentrope sound speed converges");
        check(mixtureReferenceFine.valid
                && relativeDifference(
                    mixtureSequence[2].soundSpeed,
                    mixtureReferenceFine.soundSpeed) < 2.0e-4,
            "O2/C3H8 HEM sound speed matches the independent TP isentrope");
        const PhaseSoundSpeedResult directAcoustic = acoustic.evaluateRhoE(
            mixture.density, mixture.energy, mixture.massFractions);
        check(directAcoustic.success()
                && relativeDifference(
                    directAcoustic.soundSpeed,
                    mixtureSequence[2].soundSpeed) < 5.0e-5,
            "authoritative rho-E acoustic entry point matches verified-base evaluation");

        RhoEFlashResult tamperedBase = recovered;
        tamperedBase.pressure *= 1.01;
        const PhaseSoundSpeedResult rejectedTamperedBase = acoustic.evaluate(
            tamperedBase, mixture.massFractions);
        check(rejectedTamperedBase.status
                == PhaseSoundSpeedStatus::BASE_STATE_MISMATCH,
            "HEM acoustics reject a base not reproduced by the owned model");
        std::vector<double> mismatchedComposition = mixture.massFractions;
        std::swap(mismatchedComposition[0], mismatchedComposition[1]);
        const PhaseSoundSpeedResult rejectedComposition = acoustic.evaluate(
            recovered, mismatchedComposition);
        check(rejectedComposition.status
                == PhaseSoundSpeedStatus::BASE_STATE_MISMATCH,
            "HEM acoustics reject a composition inconsistent with the base");
        RhoEFlashResult missingEntropy = recovered;
        missingEntropy.entropy = std::numeric_limits<double>::quiet_NaN();
        const PhaseSoundSpeedResult rejected = acoustic.evaluate(
            missingEntropy, mixture.massFractions);
        check(rejected.status == PhaseSoundSpeedStatus::NONFINITE_ENTROPY
                && !rejected.success(),
            "HEM acoustics reject a state without entropy instead of guessing");
    }

    PRRhoEFlash invalidRange(critical, kij, {"O2"},
        options(4.0e6, 4.5e6, 70.0, 75.0));
    const RhoEFlashResult outOfTable = invalidRange.flashMassFractions(
        900.0, -1.0e5, {1.0});
    check(!outOfTable.success()
            && !outOfTable.singlePhase.has_value()
            && !outOfTable.conservedMixture.has_value(),
        "out-of-range rho-E request fails closed without a state");

    PRRhoEFlash invalidInput(critical, kij, {"C3H8"},
        options(5.0e5, 1.5e6, 280.0, 360.0));
    const RhoEFlashResult nanTarget = invalidInput.flashMassFractions(
        std::numeric_limits<double>::quiet_NaN(), 0.0, {1.0});
    check(nanTarget.status == RhoEFlashStatus::INVALID_DENSITY
            && !nanTarget.success(),
        "NaN density is rejected before any thermodynamic evaluation");

    if (failures != 0)
    {
        std::cerr << "PR_RHO_E_FLASH_FAIL failures=" << failures << '\n';
        return 1;
    }

    std::cout << "PR_RHO_E_FLASH_PASS\n"
              << "NOTE: this enables only converged rho-e isentropic acoustic "
              << "evidence; complete flash, boundary, transport, chemistry, "
              << "and production gates remain closed.\n";
    return 0;
}
