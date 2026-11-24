#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "pr_high_temperature_residual_handoff_prototype.H"

namespace diagnostic = rfamHighTemperatureDiagnostic;

namespace
{

int checks = 0;

void require(bool condition, const std::string& message)
{
    ++checks;
    if (!condition)
    {
        std::cerr << "PR_HIGH_TEMPERATURE_HANDOFF_FAIL " << message << '\n';
        std::exit(2);
    }
}

double relativeError(double actual, double expected)
{
    return std::abs(actual - expected)
        /std::max({std::abs(actual), std::abs(expected), 1.0});
}

diagnostic::Component stableComponent(
    const std::string& name,
    double massFraction)
{
    if (name == "O2")
        return {name, massFraction, 0.0319988, 154.58, 5.043e6, 0.0222, true};
    if (name == "N2")
        return {name, massFraction, 0.0280134, 126.19, 3.396e6, 0.0372, true};
    if (name == "H2")
        return {name, massFraction, 0.00201588, 32.98, 1.293e6, -0.2160, true};
    if (name == "C3H8")
        return {name, massFraction, 0.04409562, 369.89, 4.251e6, 0.1523, true};
    if (name == "H2O")
        return {name, massFraction, 0.01801528, 647.10, 22.064e6, 0.3443, true};
    if (name == "CO")
        return {name, massFraction, 0.0280101, 132.86, 3.494e6, 0.0497, true};
    if (name == "CO2")
        return {name, massFraction, 0.0440095, 304.13, 7.377e6, 0.2236, true};
    return {};
}

diagnostic::InteractionMap legacyDiagnosticInteractions()
{
    diagnostic::InteractionMap values;
    const auto add = [&values](
        const std::string& first,
        const std::string& second,
        double value)
    {
        values[diagnostic::interactionKey(first, second)] = value;
    };
    add("N2", "O2", 0.0);
    add("N2", "C3H8", 0.08);
    add("N2", "CO2", -0.017);
    add("O2", "C3H8", 0.055);
    add("O2", "H2", 0.01);
    add("O2", "CO2", 0.02);
    add("C3H8", "CO2", 0.13);
    add("H2", "CO", 0.12);
    add("CO", "CO2", 0.07);
    add("CO2", "H2O", 0.19);
    return values;
}

struct RegressionCase
{
    std::string name;
    std::vector<diagnostic::Component> components;
};

std::vector<RegressionCase> regressionCases()
{
    const double radicalMassFraction = 1.0e-8;
    const double stableScale = 1.0 - radicalMassFraction;
    return {
        {"pure-CO", {stableComponent("CO", 1.0)}},
        {"pure-N2", {stableComponent("N2", 1.0)}},
        {"CO-N2", {
            stableComponent("CO", 0.5), stableComponent("N2", 0.5)}},
        {"CO2-N2", {
            stableComponent("CO2", 0.5), stableComponent("N2", 0.5)}},
        {"H2O-N2", {
            stableComponent("H2O", 0.5), stableComponent("N2", 0.5)}},
        {"O2-N2", {
            stableComponent("O2", 0.5), stableComponent("N2", 0.5)}},
        {"reference-crossing", {
            stableComponent("O2", 0.70),
            stableComponent("N2", 0.15),
            stableComponent("C3H8", 0.15)}},
        {"full-products", {
            stableComponent("O2", 0.0943289863862),
            stableComponent("H2", 0.1252197516978),
            stableComponent("H2O", 0.254843792564),
            stableComponent("CO", 0.299197804074),
            stableComponent("CO2", 0.226409665278)}},
        {"full-products-trace-OH", {
            stableComponent("O2", stableScale*0.0943289863862),
            stableComponent("H2", stableScale*0.1252197516978),
            stableComponent("H2O", stableScale*0.254843792564),
            stableComponent("CO", stableScale*0.299197804074),
            stableComponent("CO2", stableScale*0.226409665278),
            {"OH-ideal-trace", radicalMassFraction, 0.01700734,
             std::numeric_limits<double>::quiet_NaN(),
             std::numeric_limits<double>::quiet_NaN(),
             std::numeric_limits<double>::quiet_NaN(), false}}}
    };
}

void testSwitchConstraints()
{
    const double coefficients[] = {1.0, 0.0, 0.0, 0.0,
                                   -35.0, 84.0, -70.0, 20.0};
    double valueAtOne = 0.0;
    double firstAtOne = 0.0;
    double secondAtOne = 0.0;
    double thirdAtOne = 0.0;
    for (int power = 0; power <= 7; ++power)
    {
        valueAtOne += coefficients[power];
        if (power >= 1) firstAtOne += power*coefficients[power];
        if (power >= 2)
            secondAtOne += power*(power - 1)*coefficients[power];
        if (power >= 3)
            thirdAtOne += power*(power - 1)*(power - 2)
                *coefficients[power];
    }
    require(coefficients[0] == 1.0, "lower value coefficient");
    require(coefficients[1] == 0.0 && coefficients[2] == 0.0
            && coefficients[3] == 0.0,
            "lower derivative coefficients");
    require(valueAtOne == 0.0, "upper value constraint");
    require(firstAtOne == 0.0, "upper first-derivative constraint");
    require(secondAtOne == 0.0, "upper second-derivative constraint");
    require(thirdAtOne == 0.0, "upper third-derivative constraint");

    double previous = 1.0;
    for (int index = 0; index <= 1000; ++index)
    {
        const double temperature = diagnostic::handoffTemperatureLower
            +(diagnostic::handoffTemperatureUpper
              - diagnostic::handoffTemperatureLower)*index/1000.0;
        const diagnostic::SwitchDerivatives fade =
            diagnostic::residualFade(temperature);
        require(fade.value >= -1.0e-14 && fade.value <= 1.0 + 1.0e-14,
                "fade bounds");
        require(fade.value <= previous + 1.0e-14, "fade monotonicity");
        previous = fade.value;
    }
}

void testZeroOrderingAndThreshold()
{
    const std::vector<diagnostic::Component> stable = {
        stableComponent("O2", 1.0), stableComponent("N2", 1.0),
        stableComponent("H2", 1.0), stableComponent("C3H8", 1.0),
        stableComponent("H2O", 1.0), stableComponent("CO", 1.0),
        stableComponent("CO2", 1.0)};
    double lowest = std::numeric_limits<double>::infinity();
    std::string limiting;
    for (const diagnostic::Component& component : stable)
    {
        const double zero = diagnostic::alphaZeroTemperature(component);
        if (zero < lowest)
        {
            lowest = zero;
            limiting = component.name;
        }
    }
    require(limiting == "CO", "lowest hazardous species");
    require(std::abs(lowest - 1376.8158788294697) < 1.0e-9,
            "CO zero temperature");
    require(diagnostic::handoffTemperatureUpper < lowest,
            "handoff ends below lowest zero");

    const diagnostic::Component nitrogen = stableComponent("N2", 1.0);
    const diagnostic::Component oxygen = stableComponent("O2", 1.0);
    const double zero = diagnostic::alphaZeroTemperature(nitrogen);
    const diagnostic::RootDerivatives nitrogenRoot =
        diagnostic::originalAlphaRoot(nitrogen, zero);
    const diagnostic::RootDerivatives oxygenRoot =
        diagnostic::originalAlphaRoot(oxygen, zero);
    const double diagonalHalfWidth = std::sqrt(
        diagnostic::principalRootDerivativeThreshold
        /diagnostic::alphaScale(nitrogen))/std::abs(nitrogenRoot.first);
    const double crossHalfWidth =
        diagnostic::principalRootDerivativeThreshold
        /(std::sqrt(diagnostic::alphaScale(nitrogen)
                   *diagnostic::alphaScale(oxygen))
          *std::abs(oxygenRoot.value)*std::abs(nitrogenRoot.first));
    const double ulp = std::nextafter(zero,
        std::numeric_limits<double>::infinity()) - zero;
    require(diagonalHalfWidth > ulp,
            "diagonal threshold covers representable temperatures");
    require(crossHalfWidth < ulp,
            "cross threshold is narrower than one representable step");
}

void testPotentialDerivatives()
{
    const diagnostic::InteractionMap interactions =
        legacyDiagnosticInteractions();
    const std::vector<diagnostic::Component> components = {
        stableComponent("O2", 0.70),
        stableComponent("N2", 0.15),
        stableComponent("C3H8", 0.15)};
    const double temperature = 1175.0;
    const double molecularWeight =
        diagnostic::mixtureMolecularWeight(components);
    const double volume = molecularWeight/10.0;
    const diagnostic::PotentialDerivatives center =
        diagnostic::handoffResidualPotential(
            components, interactions, temperature, volume);
    require(center.valid, "center potential valid");

    const double temperatureStep = 0.02;
    const double volumeStep = 1.0e-5*volume;
    const double volumeSecondStep = 1.0e-3*volume;
    const auto potential = [&](double testTemperature, double testVolume)
    {
        return diagnostic::handoffResidualPotential(
            components, interactions, testTemperature, testVolume).value;
    };
    const double temperatureFirst =
        (potential(temperature + temperatureStep, volume)
         - potential(temperature - temperatureStep, volume))
        /(2.0*temperatureStep);
    const double temperatureSecond =
        (potential(temperature + temperatureStep, volume)
         - 2.0*potential(temperature, volume)
         + potential(temperature - temperatureStep, volume))
        /(temperatureStep*temperatureStep);
    const double volumeFirst =
        (potential(temperature, volume + volumeStep)
         - potential(temperature, volume - volumeStep))
        /(2.0*volumeStep);
    const double volumeSecond =
        (potential(temperature, volume + volumeSecondStep)
         - 2.0*potential(temperature, volume)
         + potential(temperature, volume - volumeSecondStep))
        /(volumeSecondStep*volumeSecondStep);
    const double temperatureVolume =
        (potential(temperature + temperatureStep, volume + volumeStep)
         - potential(temperature + temperatureStep, volume - volumeStep)
         - potential(temperature - temperatureStep, volume + volumeStep)
         + potential(temperature - temperatureStep, volume - volumeStep))
        /(4.0*temperatureStep*volumeStep);

    require(relativeError(center.temperature, temperatureFirst) < 2.0e-8,
            "potential temperature derivative");
    require(relativeError(center.temperature2, temperatureSecond) < 2.0e-6,
            "potential second temperature derivative");
    require(relativeError(center.volume, volumeFirst) < 2.0e-9,
            "potential volume derivative");
    require(relativeError(center.volume2, volumeSecond) < 2.0e-6,
            "potential second volume derivative");
    require(relativeError(center.temperatureVolume, temperatureVolume)
            < 2.0e-8,
            "potential mixed derivative");
}

void testRegressionMatrix()
{
    const diagnostic::InteractionMap interactions =
        legacyDiagnosticInteractions();
    const std::vector<RegressionCase> cases = regressionCases();
    for (const RegressionCase& testCase : cases)
    {
        for (const double density : {1.0, 10.0})
        {
            for (int index = 0; index <= 620; ++index)
            {
                const double temperature = 900.0 + 5.0*index;
                const diagnostic::ThermodynamicState state =
                    diagnostic::evaluateHandoff(
                        testCase.components, interactions, density, temperature);
                require(state.valid, testCase.name + " handoff validity");
                require(state.pressure > 0.0,
                        testCase.name + " positive pressure");
                require(state.cv > 0.0, testCase.name + " positive Cv");
                require(state.cp > 0.0, testCase.name + " positive Cp");
                require(state.dpdrhoIsothermal > 0.0,
                        testCase.name + " positive dpdrho");
                require(state.soundSpeedSquared > 0.0,
                        testCase.name + " positive sound speed squared");
            }
        }

        const double molecularWeight =
            diagnostic::mixtureMolecularWeight(testCase.components);
        const double volume = molecularWeight/10.0;
        const diagnostic::PotentialDerivatives lower =
            diagnostic::handoffResidualPotential(
                testCase.components, interactions,
                diagnostic::handoffTemperatureLower, volume);
        const diagnostic::PotentialDerivatives lowerReference =
            diagnostic::directResidualPotential(
                testCase.components, interactions,
                diagnostic::handoffTemperatureLower, volume,
                diagnostic::RootPolicy::PositiveGuarded);
        require(lower.valid && lowerReference.valid,
                testCase.name + " lower potentials valid");
        require(relativeError(lower.value, lowerReference.value) < 1.0e-14
                && relativeError(lower.temperature,
                                 lowerReference.temperature) < 1.0e-14
                && relativeError(lower.temperature2,
                                 lowerReference.temperature2) < 1.0e-14
                && relativeError(lower.volume,
                                 lowerReference.volume) < 1.0e-14
                && relativeError(lower.volume2,
                                 lowerReference.volume2) < 1.0e-14
                && relativeError(lower.temperatureVolume,
                                 lowerReference.temperatureVolume) < 1.0e-14,
                testCase.name + " lower C2 authority match");

        const diagnostic::PotentialDerivatives upper =
            diagnostic::handoffResidualPotential(
                testCase.components, interactions,
                diagnostic::handoffTemperatureUpper, volume);
        require(upper.valid && upper.value == 0.0
                && upper.temperature == 0.0
                && upper.temperature2 == 0.0
                && upper.volume == 0.0 && upper.volume2 == 0.0
                && upper.temperatureVolume == 0.0,
                testCase.name + " upper zero residual authority");

        const double highTemperature = 2000.0;
        const diagnostic::ThermodynamicState high =
            diagnostic::evaluateHandoff(
                testCase.components, interactions, 10.0, highTemperature);
        const double highVolume = molecularWeight/10.0;
        const double idealCv = 2.5*diagnostic::universalGasConstant;
        require(relativeError(
                    high.pressure,
                    diagnostic::universalGasConstant*highTemperature/highVolume)
                < 1.0e-14,
                testCase.name + " ideal pressure limit");
        require(relativeError(high.cv, idealCv) < 1.0e-14,
                testCase.name + " ideal Cv limit");
        require(relativeError(
                    high.cp, idealCv + diagnostic::universalGasConstant)
                < 1.0e-14,
                testCase.name + " ideal Cp limit");
    }
}

void testAlternativeRoutes()
{
    const diagnostic::InteractionMap interactions =
        legacyDiagnosticInteractions();
    const std::vector<diagnostic::Component> reference = {
        stableComponent("O2", 0.70),
        stableComponent("N2", 0.15),
        stableComponent("C3H8", 0.15)};
    const double zero = diagnostic::alphaZeroTemperature(
        stableComponent("N2", 1.0));
    const double epsilon = 1.0e-8;
    const double density = 10.0;
    const double molecularWeight = diagnostic::mixtureMolecularWeight(reference);

    const diagnostic::ThermodynamicState principalBelow =
        diagnostic::evaluateDirect(
            reference, interactions, density, zero - epsilon,
            diagnostic::RootPolicy::PrincipalAbsolute);
    const diagnostic::ThermodynamicState principalAbove =
        diagnostic::evaluateDirect(
            reference, interactions, density, zero + epsilon,
            diagnostic::RootPolicy::PrincipalAbsolute);
    const diagnostic::ThermodynamicState signedBelow =
        diagnostic::evaluateDirect(
            reference, interactions, density, zero - epsilon,
            diagnostic::RootPolicy::SignedDiagnostic);
    const diagnostic::ThermodynamicState signedAbove =
        diagnostic::evaluateDirect(
            reference, interactions, density, zero + epsilon,
            diagnostic::RootPolicy::SignedDiagnostic);
    const diagnostic::ThermodynamicState constrainedBelow =
        diagnostic::evaluateDirect(
            reference, interactions, density, zero - epsilon,
            diagnostic::RootPolicy::ConstrainedPositive);
    const diagnostic::ThermodynamicState constrainedAbove =
        diagnostic::evaluateDirect(
            reference, interactions, density, zero + epsilon,
            diagnostic::RootPolicy::ConstrainedPositive);
    require(principalBelow.valid && principalAbove.valid
            && signedBelow.valid && signedAbove.valid
            && constrainedBelow.valid && constrainedAbove.valid,
            "alternative route states valid");

    const double principalJump =
        (principalAbove.internalEnergy - principalBelow.internalEnergy)
        /molecularWeight;
    const double signedJump =
        (signedAbove.internalEnergy - signedBelow.internalEnergy)
        /molecularWeight;
    const double constrainedJump =
        (constrainedAbove.internalEnergy
         - constrainedBelow.internalEnergy)/molecularWeight;
    require(std::abs(principalJump - 171.55024922639132) < 5.0e-3,
            "principal-root jump reproduction");
    require(std::abs(signedJump) < 1.0e-2,
            "signed-root diagnostic continuity");
    require(std::abs(constrainedJump) < 1.0e-2,
            "constrained-alpha continuity");

    for (const RegressionCase& testCase : regressionCases())
    {
        for (const double densityValue : {1.0, 10.0})
        {
            for (int index = 0; index <= 620; ++index)
            {
                const double temperature = 900.0 + 5.0*index;
                const diagnostic::ThermodynamicState state =
                    diagnostic::evaluateDirect(
                        testCase.components, interactions,
                        densityValue, temperature,
                        diagnostic::RootPolicy::ConstrainedPositive);
                require(state.valid,
                        testCase.name + " constrained route validity");
                require(state.cv > 0.0 && state.cp > 0.0
                        && state.dpdrhoIsothermal > 0.0
                        && state.soundSpeedSquared > 0.0,
                        testCase.name + " constrained route stability");
            }
        }
    }
}

} // namespace

int main()
{
    std::cout << std::setprecision(17);
    testSwitchConstraints();
    testZeroOrderingAndThreshold();
    testPotentialDerivatives();
    testRegressionMatrix();
    testAlternativeRoutes();
    std::cout << "PR_HIGH_TEMPERATURE_HANDOFF_PASS checks="
              << checks << '\n';
    return 0;
}
