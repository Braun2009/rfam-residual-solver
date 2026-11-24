#include "PhaseEquilibrium/PRPhaseCaloric.H"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "NASA7Thermo.H"
#include "PhaseEquilibrium/TabulatedIdealCaloric.H"

namespace
{

int failures = 0;

void check(bool condition, const std::string& label)
{
    std::cout << (condition ? "  PASS " : "  FAIL ") << label << '\n';
    if (!condition) ++failures;
}

bool close(double lhs, double rhs, double relative = 1.0e-10)
{
    const double scale = std::max({std::abs(lhs), std::abs(rhs), 1.0});
    return std::abs(lhs - rhs) <= relative*scale;
}

struct IdealCaloricSample
{
    bool valid = false;
    double h = std::numeric_limits<double>::quiet_NaN();
    double s = std::numeric_limits<double>::quiet_NaN();
};

IdealCaloricSample idealCaloricSample(
    const Foam::phaseEquilibrium::PRPhaseCaloric& caloric,
    const Foam::PengRobinsonEOS& eos,
    double molecularWeight,
    double pressure,
    double temperature)
{
    using namespace Foam;
    using namespace Foam::phaseEquilibrium;
    const std::vector<double> pure = {1.0};
    const PRPhaseCaloricResult result = caloric.evaluate(
        pressure, temperature, pure, PhaseRoot::VAPOR);
    IdealCaloricSample sample;
    if (!result.success()) return sample;
    const double departureH = eos.enthalpyDeparture(
        pressure,
        temperature,
        result.state->compressibility,
        pure);
    const double departureS = eos.entropyDeparture(
        pressure,
        temperature,
        result.state->compressibility,
        pure);
    sample.h = result.state->h - departureH;
    sample.s = result.state->s - departureS
        + PRMixtureKernel::universalGasConstant/molecularWeight
            *std::log(pressure/1.0e5);
    sample.valid = std::isfinite(sample.h) && std::isfinite(sample.s);
    return sample;
}

double expectedHandoffCp(const std::string& species, double temperature)
{
    using namespace Foam;
    using namespace Foam::phaseEquilibrium;
    const NASA7Thermo& nasa = nasa7Thermo();
    const TabulatedIdealCaloric& table = tabulatedIdealCaloric();
    const double molecularWeight = nasa.MW(species);
    const bool oxygen = species == "O2";
    const double minimum = oxygen ? 195.0 : 295.0;
    const double maximum = oxygen
        ? 200.0 : NASA7Thermo::sensibleReferenceTemperature;
    const double nasaMinimum = oxygen ? 200.0 : 300.0;
    const double targetCp =
        nasa.cp_mole(species, nasaMinimum)/molecularWeight;

    if (temperature < maximum)
    {
        const TabulatedIdealCaloricResult state =
            table.evaluate(species, temperature);
        const TabulatedIdealCaloricResult boundary =
            table.evaluate(species, maximum);
        if (!state.success() || !boundary.success())
            return std::numeric_limits<double>::quiet_NaN();
        const double u = std::clamp(
            (temperature - minimum)/(maximum - minimum), 0.0, 1.0);
        const double blend = u*u*(3.0 - 2.0*u);
        return state.state->cp
            + (targetCp - boundary.state->cp)*blend;
    }
    if (!oxygen && temperature < nasaMinimum) return targetCp;
    return nasa.cp_mole(species, temperature)/molecularWeight;
}

void checkDerivativeIdentity(
    const std::string& species,
    const std::string& location,
    const Foam::phaseEquilibrium::PRPhaseCaloric& caloric,
    const Foam::PengRobinsonEOS& eos,
    double molecularWeight,
    double pressure,
    double temperature)
{
    constexpr double step = 1.0e-3;
    const IdealCaloricSample lower = idealCaloricSample(
        caloric, eos, molecularWeight, pressure, temperature - step);
    const IdealCaloricSample upper = idealCaloricSample(
        caloric, eos, molecularWeight, pressure, temperature + step);
    const std::string label = species + " " + location;
    check(lower.valid && upper.valid, label + " derivative samples are valid");
    if (!lower.valid || !upper.valid) return;

    const double expectedCp = expectedHandoffCp(species, temperature);
    const std::optional<double> cp = caloric.idealCpMass(
        temperature, {1.0});
    check(cp.has_value() && close(*cp, expectedCp, 2.0e-12),
        label + " exposes the corrected ideal-Cp authority");
    if (!cp) return;
    const double dhdT = (upper.h - lower.h)/(2.0*step);
    const double dsdT = (upper.s - lower.s)/(2.0*step);
    check(close(dhdT, *cp, 8.0e-7), label + " satisfies dh/dT=Cp");
    check(close(dsdT, *cp/temperature, 8.0e-7),
        label + " satisfies ds/dT=Cp/T");
}

void checkBoundaryC1(
    const std::string& species,
    const std::string& boundaryLabel,
    const Foam::phaseEquilibrium::PRPhaseCaloric& caloric,
    const Foam::PengRobinsonEOS& eos,
    double molecularWeight,
    double pressure,
    double temperature)
{
    constexpr double step = 1.0e-4;
    const IdealCaloricSample center = idealCaloricSample(
        caloric, eos, molecularWeight, pressure, temperature);
    const IdealCaloricSample left = idealCaloricSample(
        caloric, eos, molecularWeight, pressure, temperature - step);
    const IdealCaloricSample right = idealCaloricSample(
        caloric, eos, molecularWeight, pressure, temperature + step);
    const IdealCaloricSample leftLimit = idealCaloricSample(
        caloric,
        eos,
        molecularWeight,
        pressure,
        std::nextafter(temperature, 0.0));
    const IdealCaloricSample rightLimit = idealCaloricSample(
        caloric,
        eos,
        molecularWeight,
        pressure,
        std::nextafter(temperature, std::numeric_limits<double>::infinity()));
    const bool valid = center.valid && left.valid && right.valid
        && leftLimit.valid && rightLimit.valid;
    check(valid, boundaryLabel + " C1 samples are valid");
    if (!valid) return;

    check(std::abs(leftLimit.h - center.h) < 1.0e-7
            && std::abs(rightLimit.h - center.h) < 1.0e-7,
        boundaryLabel + " h is continuous");
    check(std::abs(leftLimit.s - center.s) < 1.0e-9
            && std::abs(rightLimit.s - center.s) < 1.0e-9,
        boundaryLabel + " s is continuous");

    const double expectedCp = expectedHandoffCp(species, temperature);
    const std::optional<double> cp = caloric.idealCpMass(
        temperature, {1.0});
    check(cp.has_value() && close(*cp, expectedCp, 2.0e-12),
        boundaryLabel + " exposes the corrected ideal-Cp authority");
    if (!cp) return;
    const double leftDh = (center.h - left.h)/step;
    const double rightDh = (right.h - center.h)/step;
    const double leftDs = (center.s - left.s)/step;
    const double rightDs = (right.s - center.s)/step;
    check(close(leftDh, *cp, 2.0e-6)
            && close(rightDh, *cp, 2.0e-6)
            && close(leftDh, rightDh, 2.0e-6),
        boundaryLabel + " removes the Cp jump");
    check(close(leftDs, *cp/temperature, 2.0e-6)
            && close(rightDs, *cp/temperature, 2.0e-6),
        boundaryLabel + " keeps ds/dT=Cp/T on both sides");
}

} // End anonymous namespace

int main()
{
    using namespace Foam;
    using namespace Foam::phaseEquilibrium;

    std::cout << "=== Explicit-root PR phase caloric adapter ===\n";

    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB kij = BinaryInteractionDB::createDefault();

    double oxygenMinimum = 0.0;
    double oxygenMaximum = 0.0;
    check(nasa7Thermo().temperatureRange(
            "O2", oxygenMinimum, oxygenMaximum)
            && oxygenMinimum == 200.0 && oxygenMaximum == 3500.0,
        "NASA7 exposes the authoritative O2 caloric range");

    PRPhaseCaloric propane(critical, kij, {"C3H8"});
    check(propane.configured(), "pure C3H8 adapter configures");

    constexpr double propanePressure = 1.0e6;
    constexpr double propaneTemperature = 300.0;
    const std::vector<double> pure = {1.0};
    const auto liquid = propane.evaluate(
        propanePressure, propaneTemperature, pure, PhaseRoot::LIQUID);
    const auto vapor = propane.evaluate(
        propanePressure, propaneTemperature, pure, PhaseRoot::VAPOR);
    check(liquid.success(), "C3H8 liquid root is supported at 300 K");
    check(vapor.success(), "C3H8 vapor root is supported at 300 K");

    PengRobinsonEOS referenceEos(critical, kij);
    referenceEos.setSpeciesNames({"C3H8"});
    const auto referenceLiquid = referenceEos.solveZ(
        propanePressure, propaneTemperature, pure,
        PengRobinsonEOS::PhaseHint::LIQUID);
    const auto referenceVapor = referenceEos.solveZ(
        propanePressure, propaneTemperature, pure,
        PengRobinsonEOS::PhaseHint::GAS);
    check(referenceLiquid.valid && referenceVapor.valid,
        "reference EOS exposes both explicit C3H8 roots");

    if (liquid.state && vapor.state &&
        referenceLiquid.valid && referenceVapor.valid)
    {
        const double MW = critical.get("C3H8").MW;
        const double rhoLiquid = propanePressure*MW/
            (referenceLiquid.Z*PRMixtureKernel::universalGasConstant
                *propaneTemperature);
        const double rhoVapor = propanePressure*MW/
            (referenceVapor.Z*PRMixtureKernel::universalGasConstant
                *propaneTemperature);
        check(close(liquid.state->compressibility, referenceLiquid.Z),
            "LIQUID request selects the smallest physical root");
        check(close(vapor.state->compressibility, referenceVapor.Z),
            "VAPOR request selects the largest physical root");
        check(close(liquid.state->rho, rhoLiquid),
            "liquid density follows its specified root");
        check(close(vapor.state->rho, rhoVapor),
            "vapor density follows its specified root");
        check(liquid.state->rho > vapor.state->rho,
            "liquid and vapor roots remain physically distinct");

        const double expectedLiquidH =
            nasa7Thermo().hsens_mass_mix({"C3H8"}, pure,
                propaneTemperature)
            + referenceEos.enthalpyDeparture(
                propanePressure, propaneTemperature,
                referenceLiquid.Z, pure);
        check(close(liquid.state->h, expectedLiquidH),
            "phase enthalpy uses NASA7 sensible plus corrected PR departure");
        check(close(liquid.state->e,
            liquid.state->h - propanePressure/liquid.state->rho, 1.0e-13),
            "liquid e=h-p/rho identity closes");
        check(close(vapor.state->e,
            vapor.state->h - propanePressure/vapor.state->rho, 1.0e-13),
            "vapor e=h-p/rho identity closes");
        check(std::isfinite(liquid.state->s)
                && std::isfinite(vapor.state->s),
            "both explicit roots expose finite entropy on one reference");

        for (const auto& phase :
             {std::pair<const char*, PhaseRoot>{"liquid", PhaseRoot::LIQUID},
              std::pair<const char*, PhaseRoot>{"vapor", PhaseRoot::VAPOR}})
        {
            constexpr double deltaT = 1.0e-3;
            const auto lower = propane.evaluate(
                propanePressure, propaneTemperature - deltaT, pure,
                phase.second);
            const auto upper = propane.evaluate(
                propanePressure, propaneTemperature + deltaT, pure,
                phase.second);
            const bool derivativeValid = lower.success() && upper.success();
            check(derivativeValid,
                std::string(phase.first)
                    + " entropy/enthalpy derivative samples are valid");
            if (derivativeValid)
            {
                const double dhdT = (upper.state->h - lower.state->h)
                    /(2.0*deltaT);
                const double TdsdT = propaneTemperature
                    *(upper.state->s - lower.state->s)/(2.0*deltaT);
                check(std::abs(dhdT - TdsdT)
                        /std::max({std::abs(dhdT), std::abs(TdsdT), 1.0})
                        < 2.0e-5,
                    std::string(phase.first)
                        + " satisfies dh=T ds along fixed-p root");
            }
        }
    }

    const auto propaneBelow = propane.evaluate(
        propanePressure, 299.999, pure, PhaseRoot::LIQUID);
    check(propaneBelow.success(),
        "C3H8 uses its low-temperature table below the NASA range");
    const auto propaneOutOfTable = propane.evaluate(
        propanePressure, 79.999, pure, PhaseRoot::LIQUID);
    check(propaneOutOfTable.status ==
        PRPhaseCaloricStatus::UNSUPPORTED_CALORIC_RANGE,
        "C3H8 below its frozen table fails closed");
    check(!propaneOutOfTable.state.has_value() &&
        propaneOutOfTable.supportedTemperatureMinimum == 80.0,
        "C3H8 table boundary is reported without clamping");

    PRPhaseCaloric oxygen(critical, kij, {"O2"});
    check(oxygen.configured(), "pure O2 adapter configures");
    const auto oxygenAtLower = oxygen.evaluate(
        4.0e6, 200.0, pure, PhaseRoot::LIQUID);
    const auto oxygenBelow = oxygen.evaluate(
        4.0e6, 199.999, pure, PhaseRoot::LIQUID);
    const auto oxygenAtUpper = oxygen.evaluate(
        4.0e6, 3500.0, pure, PhaseRoot::VAPOR);
    const auto oxygenAbove = oxygen.evaluate(
        4.0e6, 3500.001, pure, PhaseRoot::VAPOR);
    check(oxygenAtLower.success(), "O2 is supported at its 200 K boundary");
    check(oxygenBelow.success(),
        "O2 below 200 K uses the frozen low-temperature ideal caloric table");
    const auto oxygenOutOfTable = oxygen.evaluate(
        4.0e6, 79.999, pure, PhaseRoot::LIQUID);
    check(oxygenOutOfTable.status ==
        PRPhaseCaloricStatus::UNSUPPORTED_CALORIC_RANGE,
        "O2 below its frozen table fails closed");
    check(!oxygenOutOfTable.state.has_value() &&
        oxygenOutOfTable.supportedTemperatureMinimum == 80.0,
        "O2 table boundary is reported without clamping");
    check(oxygenBelow.state && oxygenAtLower.state &&
        std::abs(oxygenBelow.state->h - oxygenAtLower.state->h) < 2.0,
        "O2 low-table/NASA splice is enthalpy-continuous");
    check(oxygenBelow.state && oxygenAtLower.state &&
        std::abs(oxygenBelow.state->s - oxygenAtLower.state->s) < 2.0e-2,
        "O2 low-table/NASA splice is entropy-continuous");
    check(oxygenAtUpper.success(), "O2 is supported at its 3500 K boundary");
    check(oxygenAbove.status ==
        PRPhaseCaloricStatus::UNSUPPORTED_CALORIC_RANGE,
        "O2 above 3500 K is rejected without extrapolation");

    constexpr double handoffPressure = 1.0e5;
    PengRobinsonEOS oxygenReferenceEos(critical, kij);
    oxygenReferenceEos.setSpeciesNames({"O2"});
    const double oxygenMW = critical.get("O2").MW;
    const double propaneMW = critical.get("C3H8").MW;

    checkBoundaryC1(
        "O2", "O2 195 K handoff start", oxygen, oxygenReferenceEos,
        oxygenMW, handoffPressure, 195.0);
    checkBoundaryC1(
        "O2", "O2 200 K table/NASA handoff", oxygen, oxygenReferenceEos,
        oxygenMW, handoffPressure, 200.0);
    checkDerivativeIdentity(
        "O2", "below handoff", oxygen, oxygenReferenceEos,
        oxygenMW, handoffPressure, 194.0);
    checkDerivativeIdentity(
        "O2", "inside handoff", oxygen, oxygenReferenceEos,
        oxygenMW, handoffPressure, 197.5);
    checkDerivativeIdentity(
        "O2", "above handoff", oxygen, oxygenReferenceEos,
        oxygenMW, handoffPressure, 201.0);

    checkBoundaryC1(
        "C3H8", "C3H8 295 K handoff start", propane, referenceEos,
        propaneMW, handoffPressure, 295.0);
    checkBoundaryC1(
        "C3H8", "C3H8 298.15 K table/bridge handoff", propane,
        referenceEos, propaneMW, handoffPressure,
        NASA7Thermo::sensibleReferenceTemperature);
    checkBoundaryC1(
        "C3H8", "C3H8 300 K bridge/NASA handoff", propane, referenceEos,
        propaneMW, handoffPressure, 300.0);
    checkDerivativeIdentity(
        "C3H8", "below handoff", propane, referenceEos,
        propaneMW, handoffPressure, 294.0);
    checkDerivativeIdentity(
        "C3H8", "inside handoff", propane, referenceEos,
        propaneMW, handoffPressure, 296.5);
    checkDerivativeIdentity(
        "C3H8", "above handoff entropy bridge", propane, referenceEos,
        propaneMW, handoffPressure, 299.0);
    checkDerivativeIdentity(
        "C3H8", "NASA side", propane, referenceEos,
        propaneMW, handoffPressure, 301.0);

    const TabulatedIdealCaloric& rawTable = tabulatedIdealCaloric();
    const auto rawOxygen185 = rawTable.evaluate("O2", 185.0);
    const auto rawOxygen190 = rawTable.evaluate("O2", 190.0);
    const auto idealOxygen185 = idealCaloricSample(
        oxygen, oxygenReferenceEos, oxygenMW, handoffPressure, 185.0);
    const auto idealOxygen190 = idealCaloricSample(
        oxygen, oxygenReferenceEos, oxygenMW, handoffPressure, 190.0);
    check(rawOxygen185.success() && rawOxygen190.success()
            && idealOxygen185.valid && idealOxygen190.valid
            && close(
                idealOxygen190.h - idealOxygen185.h,
                rawOxygen190.state->sensibleEnthalpy
                    - rawOxygen185.state->sensibleEnthalpy,
                2.0e-10)
            && close(
                idealOxygen190.s - idealOxygen185.s,
                rawOxygen190.state->sensibleEntropy
                    - rawOxygen185.state->sensibleEntropy,
                2.0e-10),
        "O2 table increments below 195 K are unchanged");

    const auto rawPropane285 = rawTable.evaluate("C3H8", 285.0);
    const auto rawPropane290 = rawTable.evaluate("C3H8", 290.0);
    const auto idealPropane285 = idealCaloricSample(
        propane, referenceEos, propaneMW, handoffPressure, 285.0);
    const auto idealPropane290 = idealCaloricSample(
        propane, referenceEos, propaneMW, handoffPressure, 290.0);
    check(rawPropane285.success() && rawPropane290.success()
            && idealPropane285.valid && idealPropane290.valid
            && close(
                idealPropane290.h - idealPropane285.h,
                rawPropane290.state->sensibleEnthalpy
                    - rawPropane285.state->sensibleEnthalpy,
                2.0e-10)
            && close(
                idealPropane290.s - idealPropane285.s,
                rawPropane290.state->sensibleEntropy
                    - rawPropane285.state->sensibleEntropy,
                2.0e-10),
        "C3H8 table increments below 295 K are unchanged");

    const auto idealOxygen205 = idealCaloricSample(
        oxygen, oxygenReferenceEos, oxygenMW, handoffPressure, 205.0);
    check(idealOxygen205.valid
            && close(
                idealOxygen205.h,
                nasa7Thermo().hsens_mole("O2", 205.0)/oxygenMW)
            && close(
                idealOxygen205.s,
                nasa7Thermo().s_mole("O2", 205.0)/oxygenMW),
        "O2 NASA behavior above 200 K is unchanged");

    const auto idealPropane305 = idealCaloricSample(
        propane, referenceEos, propaneMW, handoffPressure, 305.0);
    check(idealPropane305.valid
            && close(
                idealPropane305.h,
                nasa7Thermo().hsens_mole("C3H8", 305.0)/propaneMW)
            && close(
                idealPropane305.s,
                nasa7Thermo().s_mole("C3H8", 305.0)/propaneMW),
        "C3H8 NASA behavior above 300 K is unchanged");

    const auto idealPropane299 = idealCaloricSample(
        propane, referenceEos, propaneMW, handoffPressure, 299.0);
    const double propaneCp300 =
        nasa7Thermo().cp_mole("C3H8", 300.0)/propaneMW;
    const double expectedPropaneEntropy299 =
        nasa7Thermo().s_mole("C3H8", 300.0)/propaneMW
        + propaneCp300*std::log(299.0/300.0);
    check(idealPropane299.valid
            && close(
                idealPropane299.h,
                nasa7Thermo().hsens_mole("C3H8", 299.0)/propaneMW)
            && close(idealPropane299.s, expectedPropaneEntropy299),
        "C3H8 298.15-300 K bridge uses Cp*ln(T/300) entropy");

    PRPhaseCaloric mixture(critical, kij, {"O2", "C3H8"});
    check(mixture.configured(), "O2/C3H8 adapter has sourced kij");
    const auto inactivePropane = mixture.evaluate(
        4.0e6, 200.0, {1.0, 0.0}, PhaseRoot::VAPOR);
    check(inactivePropane.success(),
        "zero-fraction C3H8 does not narrow the O2 caloric range");
    const auto activePropane = mixture.evaluate(
        4.0e6, 250.0, {0.999, 0.001}, PhaseRoot::VAPOR);
    check(activePropane.success(),
        "present C3H8 uses the low-temperature caloric table at 250 K");
    const auto activeCryogenicPropane = mixture.evaluate(
        4.0e6, 130.0, {0.999, 0.001}, PhaseRoot::LIQUID);
    check(activeCryogenicPropane.success(),
        "present C3H8 retains a sourced ideal caloric term at 130 K");

    PRPhaseCaloric missingPair(critical, kij, {"O2", "H2O"});
    check(missingPair.configured(),
        "missing kij does not invalidate every lower-dimensional subspace");
    const auto missingPairInactive = missingPair.evaluate(
        1.0e5, 300.0, {1.0, 0.0}, PhaseRoot::VAPOR);
    check(missingPairInactive.success(),
        "inactive missing kij leaves the pure O2 subspace usable");
    const auto missingPairActive = missingPair.evaluate(
        1.0e5, 300.0, {0.5, 0.5}, PhaseRoot::VAPOR);
    check(missingPairActive.status ==
            PRPhaseCaloricStatus::MISSING_BINARY_INTERACTION,
        "active missing kij fails the caloric state explicitly");
    check(!missingPairActive.state.has_value(),
        "active missing kij exposes no caloric state");

    const std::vector<double> equimolar = {0.5, 0.5};
    const auto mixedState = mixture.evaluate(
        1.0e6, 450.0, equimolar, PhaseRoot::VAPOR);
    check(mixedState.success(), "supported O2/C3H8 mixture returns a state");
    if (mixedState.state)
    {
        const double oxygenMW = critical.get("O2").MW;
        const double propaneMW = critical.get("C3H8").MW;
        const double expectedMW = 0.5*(oxygenMW + propaneMW);
        check(close(mixedState.state->molecularWeight, expectedMW),
            "mixture molecular weight remains on a mole basis");
        check(close(mixedState.state->massFractions[0],
            0.5*oxygenMW/expectedMW) &&
            close(mixedState.state->massFractions[1],
                0.5*propaneMW/expectedMW),
            "returned phase composition is converted to mass fractions");
    }

    BinaryInteractionDB mutableKij = BinaryInteractionDB::createDefault();
    PRPhaseCaloric frozenCaloric(
        critical, mutableKij, {"O2", "C3H8"});
    const auto frozenBefore = frozenCaloric.evaluate(
        1.0e6, 450.0, equimolar, PhaseRoot::VAPOR);
    mutableKij.add("O2", "C3H8", 0.9);
    const auto frozenAfter = frozenCaloric.evaluate(
        1.0e6, 450.0, equimolar, PhaseRoot::VAPOR);
    check(frozenBefore.success() && frozenAfter.success() &&
            close(frozenBefore.state->rho, frozenAfter.state->rho) &&
            close(frozenBefore.state->h, frozenAfter.state->h),
        "phase caloric freezes EOS and kij data at construction");

    const auto stableRejected = oxygen.evaluate(
        4.0e6, 200.0, pure, PhaseRoot::STABLE);
    check(stableRejected.status == PRPhaseCaloricStatus::INVALID_PHASE_ROOT,
        "adapter requires an explicit LIQUID or VAPOR root");
    check(!stableRejected.state.has_value(),
        "invalid root request exposes no state");

    if (failures != 0)
    {
        std::cerr << "PR_PHASE_CALORIC_FAIL failures=" << failures << '\n';
        return 1;
    }

    std::cout << "PR_PHASE_CALORIC_PASS\n"
              << "NOTE: explicit-root phase calorics are not a phase-stability "
              << "test or a complete flash.\n";
    return 0;
}
