#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "PhaseEquilibrium/PRResidualHelmholtzAuthority.H"
#include "SpeciesCriticalDB.H"
#include "BinaryInteractionDB.H"
#include "pr_high_temperature_residual_handoff_prototype.H"

namespace phase = Foam::phaseEquilibrium;
namespace diagnostic = rfamHighTemperatureDiagnostic;

namespace
{

int checks = 0;

void require(bool condition, const std::string& message)
{
    ++checks;
    if (!condition)
    {
        std::cerr << "PR_RESIDUAL_HELMHOLTZ_FAIL " << message << '\n';
        std::exit(2);
    }
}

double relativeError(double actual, double expected)
{
    return std::abs(actual - expected)
        /std::max({std::abs(actual), std::abs(expected), 1.0e-14});
}

std::string scientific(double value)
{
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(17) << value;
    return stream.str();
}

const std::vector<std::string>& names()
{
    static const std::vector<std::string> value = {
        "O2", "N2", "C3H8", "H2", "H2O", "CO", "CO2"};
    return value;
}

double knownKij(const std::string& first, const std::string& second)
{
    const auto key = first < second
        ? std::make_pair(first, second)
        : std::make_pair(second, first);
    static const std::map<std::pair<std::string, std::string>, double> values = {
        {{"N2", "O2"}, 0.0},
        {{"N2", "C3H8"}, 0.08},
        {{"O2", "C3H8"}, 0.055},
        {{"O2", "H2"}, 0.01},
        {{"O2", "CO2"}, 0.02},
        {{"C3H8", "CO2"}, 0.13},
        {{"H2", "CO"}, 0.12},
        {{"CO", "CO2"}, 0.07},
        {{"CO2", "H2O"}, 0.19}};
    const auto iterator = values.find(key);
    return iterator == values.end() ? 0.0 : iterator->second;
}

Foam::BinaryInteractionDB completeInteractions()
{
    Foam::BinaryInteractionDB result;
    for (std::size_t i = 0; i < names().size(); ++i)
        for (std::size_t j = i + 1; j < names().size(); ++j)
            result.add(names()[i], names()[j],
                       knownKij(names()[i], names()[j]));
    return result;
}

phase::PRResidualHelmholtzAuthority authority()
{
    return phase::PRResidualHelmholtzAuthority(
        Foam::SpeciesCriticalDB::createDefault(), completeInteractions(), names());
}

double molecularWeight(const std::vector<double>& moleFractions)
{
    const Foam::SpeciesCriticalDB critical =
        Foam::SpeciesCriticalDB::createDefault();
    double result = 0.0;
    for (std::size_t i = 0; i < moleFractions.size(); ++i)
        result += moleFractions[i]*critical.get(names()[i]).MW;
    return result;
}

phase::PRResidualHelmholtzRequest request(
    const std::vector<double>& moleFractions,
    double temperature,
    double density,
    double idealCv = 3.0*phase::PRResidualHelmholtzAuthority::
        universalGasConstant)
{
    const double mw = molecularWeight(moleFractions);
    const double volume = mw/density;
    phase::PRResidualHelmholtzRequest result;
    result.temperature = temperature;
    result.density = density;
    result.moleFractions = moleFractions;
    result.ideal.cv = idealCv;
    result.ideal.internalEnergy = idealCv*(temperature - 300.0);
    result.ideal.entropy = idealCv*std::log(temperature/300.0)
        + phase::PRResidualHelmholtzAuthority::universalGasConstant
            *std::log(volume/0.01);
    return result;
}

std::vector<double> referenceComposition()
{
    return {0.70, 0.15, 0.15, 0.0, 0.0, 0.0, 0.0};
}

std::vector<double> productComposition()
{
    return {0.0943289863862, 0.0, 0.0, 0.1252197516978,
            0.254843792564, 0.299197804074, 0.226409665278};
}

diagnostic::Component prototypeComponent(
    const std::string& name, double massFraction)
{
    if (name == "O2")
        return {name, massFraction, 0.0319988, 154.58, 5.043e6, 0.0222, true};
    if (name == "N2")
        return {name, massFraction, 0.0280134, 126.19, 3.396e6, 0.0372, true};
    if (name == "C3H8")
        return {name, massFraction, 0.04409562, 369.89, 4.251e6, 0.1523, true};
    if (name == "H2")
        return {name, massFraction, 0.00201588, 32.98, 1.293e6, -0.2160, true};
    if (name == "H2O")
        return {name, massFraction, 0.01801528, 647.10, 22.064e6, 0.3443, true};
    if (name == "CO")
        return {name, massFraction, 0.0280101, 132.86, 3.494e6, 0.0497, true};
    if (name == "CO2")
        return {name, massFraction, 0.0440095, 304.13, 7.377e6, 0.2236, true};
    return {};
}

diagnostic::InteractionMap prototypeInteractions()
{
    diagnostic::InteractionMap result;
    for (std::size_t i = 0; i < names().size(); ++i)
        for (std::size_t j = i + 1; j < names().size(); ++j)
            result[diagnostic::interactionKey(names()[i], names()[j])] =
                knownKij(names()[i], names()[j]);
    return result;
}

void testCapabilityAndConfiguration()
{
    const phase::PRResidualHelmholtzAuthority model = authority();
    static_assert(!phase::PRResidualHelmholtzAuthority::productionCapability,
                  "candidate must remain disabled");
    require(model.configured(), "complete configuration accepted");
    require(model.minimumCertifiedAlphaRoot() >
                phase::PRResidualHelmholtzAuthority::minimumAlphaRootMargin,
            "global alpha-root margin recorded");
    require(model.limitingSpecies() == "CO", "CO is limiting species");

    const auto cap = model.evaluate(request(referenceComposition(), 1175.0, 10.0));
    require(cap.status == phase::PRResidualHelmholtzStatus::CAPABILITY_DISABLED,
            "runtime entry remains capability disabled");
    require(!cap.state.has_value(), "disabled entry publishes no state");

    Foam::BinaryInteractionDB missing;
    missing.add("O2", "N2", 0.0);
    const phase::PRResidualHelmholtzAuthority missingModel(
        Foam::SpeciesCriticalDB::createDefault(), missing,
        {"O2", "N2", "C3H8"});
    require(missingModel.configurationStatus()
                == phase::PRResidualHelmholtzStatus::MISSING_BINARY_INTERACTION,
            "missing pair fails at configuration");

    Foam::BinaryInteractionDB invalid = completeInteractions();
    invalid.add("O2", "N2", 1.0);
    const phase::PRResidualHelmholtzAuthority invalidModel(
        Foam::SpeciesCriticalDB::createDefault(), invalid, {"O2", "N2"});
    require(invalidModel.configurationStatus()
                == phase::PRResidualHelmholtzStatus::INVALID_BINARY_INTERACTION,
            "kij >= 1 fails at configuration");

    const phase::PRResidualHelmholtzAuthority duplicateModel(
        Foam::SpeciesCriticalDB::createDefault(), completeInteractions(),
        {"O2", "O2"});
    require(duplicateModel.configurationStatus()
                == phase::PRResidualHelmholtzStatus::DUPLICATE_SPECIES,
            "duplicate species fails at configuration");

    Foam::SpeciesCriticalDB badCritical =
        Foam::SpeciesCriticalDB::createDefault();
    badCritical.add("X", 100.0, 1.0e6, 0.5, 0.03);
    const phase::PRResidualHelmholtzAuthority rootExit(
        badCritical, Foam::BinaryInteractionDB(), {"X"});
    require(rootExit.configurationStatus()
                == phase::PRResidualHelmholtzStatus::ALPHA_ROOT_INTERVAL_EXIT,
            "alpha root outside fixed interval fails closed");
}

void testHandoffAndHighBranch()
{
    const phase::PRResidualHelmholtzAuthority model = authority();
    const std::vector<double> x = referenceComposition();
    const auto lower = model.evaluateDiagnostic(request(x, 1000.0, 10.0));
    const auto interior = model.evaluateDiagnostic(request(x, 1175.0, 10.0));
    const auto upper = model.evaluateDiagnostic(request(x, 1325.0, 10.0));
    const auto high = model.evaluateDiagnostic(request(x, 4000.0, 10.0));
    require(lower.success() && interior.success() && upper.success()
                && high.success(), "handoff states valid");
    require(lower.state->token.alphaEvaluated,
            "lower side evaluates PR alpha");
    require(interior.state->token.alphaEvaluated,
            "interior side evaluates PR alpha");
    require(!upper.state->token.alphaEvaluated
                && !high.state->token.alphaEvaluated,
            "upper branch never evaluates PR alpha");
    require(lower.state->token.fade == 1.0
                && lower.state->token.dFadeDT == 0.0
                && lower.state->token.d2FadeDT2 == 0.0
                && lower.state->token.d3FadeDT3 == 0.0,
            "lower C3 fade endpoint");
    require(upper.state->token.fade == 0.0
                && upper.state->token.dFadeDT == 0.0
                && upper.state->token.d2FadeDT2 == 0.0
                && upper.state->token.d3FadeDT3 == 0.0,
            "upper C3 fade endpoint");
    require(upper.state->token.residualHelmholtz == 0.0
                && upper.state->token.dResidualHelmholtzDT == 0.0
                && upper.state->token.d2ResidualHelmholtzDT2 == 0.0
                && upper.state->token.dResidualHelmholtzDv == 0.0
                && upper.state->token.d2ResidualHelmholtzDv2 == 0.0
                && upper.state->token.d2ResidualHelmholtzDTDv == 0.0,
            "upper residual is exactly zero");
    require(high.state->pressure > 0.0 && high.state->cv > 0.0
                && high.state->cp > 0.0
                && high.state->soundSpeedSquared > 0.0,
            "zero-residual high branch remains thermodynamically finite");
}

void testPotentialDerivatives()
{
    const phase::PRResidualHelmholtzAuthority model = authority();
    const std::vector<double> x = referenceComposition();
    const double temperature = 1175.0;
    const double density = 10.0;
    const double mw = molecularWeight(x);
    const double volume = mw/density;
    const double dT = 0.02;
    const double dV = 1.0e-5*volume;
    const double volumeSecondStep = 1.0e-3*volume;
    const auto evaluateAt = [&](double T, double V)
    {
        return model.evaluateDiagnostic(request(x, T, mw/V));
    };
    const auto center = evaluateAt(temperature, volume);
    const auto plusT = evaluateAt(temperature + dT, volume);
    const auto minusT = evaluateAt(temperature - dT, volume);
    const auto plusV = evaluateAt(temperature, volume + dV);
    const auto minusV = evaluateAt(temperature, volume - dV);
    const auto plusV2 = evaluateAt(temperature, volume + volumeSecondStep);
    const auto minusV2 = evaluateAt(temperature, volume - volumeSecondStep);
    const auto plusTV = evaluateAt(temperature + dT, volume + dV);
    const auto plusTminusV = evaluateAt(temperature + dT, volume - dV);
    const auto minusTplusV = evaluateAt(temperature - dT, volume + dV);
    const auto minusTminusV = evaluateAt(temperature - dT, volume - dV);
    require(center.success() && plusT.success() && minusT.success()
                && plusV.success() && minusV.success() && plusV2.success()
                && minusV2.success() && plusTV.success()
                && plusTminusV.success() && minusTplusV.success()
                && minusTminusV.success(), "finite-difference states valid");

    const auto value = [](const phase::PRResidualHelmholtzResult& result)
    {
        return result.state->token.residualHelmholtz;
    };
    const double dT1 = (value(plusT) - value(minusT))/(2.0*dT);
    const double dT2 = (value(plusT) - 2.0*value(center) + value(minusT))
        /(dT*dT);
    const double dV1 = (value(plusV) - value(minusV))/(2.0*dV);
    const double secondVolumeDerivative =
        (value(plusV2) - 2.0*value(center) + value(minusV2))
        /(volumeSecondStep*volumeSecondStep);
    const double dTV = (value(plusTV) - value(plusTminusV)
        - value(minusTplusV) + value(minusTminusV))/(4.0*dT*dV);
    require(relativeError(center.state->token.dResidualHelmholtzDT, dT1)
                < 2.0e-8, "A_T derivative");
    require(relativeError(center.state->token.d2ResidualHelmholtzDT2, dT2)
                < 2.0e-6, "A_TT derivative");
    require(relativeError(center.state->token.dResidualHelmholtzDv, dV1)
                < 2.0e-9, "A_v derivative");
    require(relativeError(
                center.state->token.d2ResidualHelmholtzDv2,
                secondVolumeDerivative)
                < 2.0e-6, "A_vv derivative");
    const double mixedRelativeError = relativeError(
        center.state->token.d2ResidualHelmholtzDTDv, dTV);
    require(mixedRelativeError < 5.0e-8,
            "A_Tv derivative relative-error="
                + scientific(mixedRelativeError)
                + " analytic="
                + scientific(center.state->token.d2ResidualHelmholtzDTDv)
                + " finite-difference=" + scientific(dTV));
}

void testThermodynamicIdentities()
{
    const phase::PRResidualHelmholtzAuthority model = authority();
    const std::vector<std::vector<double>> compositions = {
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        referenceComposition(), productComposition()};
    for (const std::vector<double>& x : compositions)
    {
        for (double density : {1.0, 10.0, 25.0})
        {
            for (int index = 0; index <= 500; ++index)
            {
                const double temperature = 1000.0 + 10.0*index;
                const auto result = model.evaluateDiagnostic(
                    request(x, temperature, density, 3.0*phase::
                            PRResidualHelmholtzAuthority::universalGasConstant));
                require(result.success(), "interval state succeeds");
                const auto& state = *result.state;
                require(state.pressure > 0.0 && state.cv > 0.0
                            && state.cp > 0.0
                            && state.isothermalPressureDensityDerivative > 0.0
                            && state.soundSpeedSquared > 0.0,
                        "interval stability margins positive");
                require(relativeError(
                            state.enthalpy,
                            state.internalEnergy + state.pressure/state.density)
                            < 2.0e-14,
                        "h=e+p/rho identity");
                require(std::isfinite(state.pressureTemperatureDerivative)
                            && std::isfinite(state.pressureMolarVolumeDerivative),
                        "mechanical derivatives finite");
            }
        }
    }
}

phase::PRResidualHelmholtzResult evaluateAtTemperatureVolume(
    const phase::PRResidualHelmholtzAuthority& model,
    const std::vector<double>& moleFractions,
    double temperature,
    double molarVolume)
{
    return model.evaluateDiagnostic(request(
        moleFractions, temperature,
        molecularWeight(moleFractions)/molarVolume));
}

phase::PRResidualHelmholtzResult solveVolumeAtPressure(
    const phase::PRResidualHelmholtzAuthority& model,
    const std::vector<double>& moleFractions,
    double temperature,
    double targetPressure,
    double volumeSeed)
{
    double volume = volumeSeed;
    phase::PRResidualHelmholtzResult state;
    for (int iteration = 0; iteration < 12; ++iteration)
    {
        state = evaluateAtTemperatureVolume(
            model, moleFractions, temperature, volume);
        if (!state.success()) return state;
        const double correction =
            (state.state->pressure - targetPressure)
            /state.state->pressureMolarVolumeDerivative;
        volume -= correction;
        if (std::abs(correction) <= 2.0e-14*volume) break;
    }
    return evaluateAtTemperatureVolume(
        model, moleFractions, temperature, volume);
}

phase::PRResidualHelmholtzResult solveTemperatureAtEntropy(
    const phase::PRResidualHelmholtzAuthority& model,
    const std::vector<double>& moleFractions,
    double molarVolume,
    double targetEntropy,
    double temperatureSeed)
{
    double temperature = temperatureSeed;
    phase::PRResidualHelmholtzResult state;
    for (int iteration = 0; iteration < 12; ++iteration)
    {
        state = evaluateAtTemperatureVolume(
            model, moleFractions, temperature, molarVolume);
        if (!state.success()) return state;
        const double entropyTemperatureDerivative =
            state.state->cv/temperature;
        const double correction =
            (state.state->entropy - targetEntropy)
            /entropyTemperatureDerivative;
        temperature -= correction;
        if (std::abs(correction) <= 2.0e-14*temperature) break;
    }
    return evaluateAtTemperatureVolume(
        model, moleFractions, temperature, molarVolume);
}

void testDerivedDerivativeConsistency()
{
    const phase::PRResidualHelmholtzAuthority model = authority();
    const std::vector<double> x = referenceComposition();
    const double temperature = 1175.0;
    const double density = 10.0;
    const double mw = molecularWeight(x);
    const double volume = mw/density;
    const double temperatureStep = 0.02;
    const double volumeStep = 1.0e-5*volume;
    const auto center = evaluateAtTemperatureVolume(
        model, x, temperature, volume);
    const auto plusT = evaluateAtTemperatureVolume(
        model, x, temperature + temperatureStep, volume);
    const auto minusT = evaluateAtTemperatureVolume(
        model, x, temperature - temperatureStep, volume);
    const auto plusV = evaluateAtTemperatureVolume(
        model, x, temperature, volume + volumeStep);
    const auto minusV = evaluateAtTemperatureVolume(
        model, x, temperature, volume - volumeStep);
    require(center.success() && plusT.success() && minusT.success()
                && plusV.success() && minusV.success(),
            "derived finite-difference states valid");

    const double pressureTemperature =
        (plusT.state->pressure - minusT.state->pressure)
        /(2.0*temperatureStep);
    const double pressureVolume =
        (plusV.state->pressure - minusV.state->pressure)
        /(2.0*volumeStep);
    const double cv =
        (plusT.state->internalEnergy - minusT.state->internalEnergy)
        /(2.0*temperatureStep);
    const double entropyTemperature =
        (plusT.state->entropy - minusT.state->entropy)
        /(2.0*temperatureStep);
    const double entropyVolume =
        (plusV.state->entropy - minusV.state->entropy)
        /(2.0*volumeStep);
    require(relativeError(
                center.state->pressureTemperatureDerivative,
                pressureTemperature) < 2.0e-8,
            "p_T finite difference");
    require(relativeError(
                center.state->pressureMolarVolumeDerivative,
                pressureVolume) < 2.0e-9,
            "p_v finite difference");
    require(relativeError(center.state->cv, cv) < 2.0e-8,
            "Cv from constant-volume energy derivative");
    require(relativeError(
                center.state->cv/temperature,
                entropyTemperature) < 2.0e-8,
            "Cv from constant-volume entropy derivative");
    require(relativeError(
                center.state->pressureTemperatureDerivative/mw,
                entropyVolume) < 2.0e-8,
            "Maxwell entropy-volume derivative");

    const auto constantPressurePlus = solveVolumeAtPressure(
        model, x, temperature + temperatureStep,
        center.state->pressure, volume);
    const auto constantPressureMinus = solveVolumeAtPressure(
        model, x, temperature - temperatureStep,
        center.state->pressure, volume);
    require(constantPressurePlus.success() && constantPressureMinus.success(),
            "constant-pressure solves succeed");
    require(relativeError(constantPressurePlus.state->pressure,
                          center.state->pressure) < 2.0e-13
                && relativeError(constantPressureMinus.state->pressure,
                                 center.state->pressure) < 2.0e-13,
            "constant-pressure solves close");
    const double cp =
        (constantPressurePlus.state->enthalpy
         - constantPressureMinus.state->enthalpy)
        /(2.0*temperatureStep);
    require(relativeError(center.state->cp, cp) < 3.0e-8,
            "Cp from constant-pressure enthalpy derivative");

    const double densityStep = 1.0e-5*density;
    const double lowerDensity = density - densityStep;
    const double upperDensity = density + densityStep;
    const auto isentropicLower = solveTemperatureAtEntropy(
        model, x, mw/lowerDensity, center.state->entropy, temperature);
    const auto isentropicUpper = solveTemperatureAtEntropy(
        model, x, mw/upperDensity, center.state->entropy, temperature);
    require(isentropicLower.success() && isentropicUpper.success(),
            "isentropic solves succeed");
    require(relativeError(isentropicLower.state->entropy,
                          center.state->entropy) < 2.0e-13
                && relativeError(isentropicUpper.state->entropy,
                                 center.state->entropy) < 2.0e-13,
            "isentropic solves close");
    const double soundSpeedSquared =
        (isentropicUpper.state->pressure - isentropicLower.state->pressure)
        /(upperDensity - lowerDensity);
    require(relativeError(center.state->soundSpeedSquared, soundSpeedSquared)
                < 3.0e-8,
            "sound speed from isentropic pressure-density derivative");
}

void testFailClosedInputs()
{
    const phase::PRResidualHelmholtzAuthority model = authority();
    const std::vector<double> x = referenceComposition();
    auto assertFailure = [&](const phase::PRResidualHelmholtzRequest& value,
                             phase::PRResidualHelmholtzStatus status,
                             const std::string& message)
    {
        const auto result = model.evaluateDiagnostic(value);
        require(result.status == status, message + " status");
        require(!result.state.has_value(), message + " no publication");
    };

    assertFailure(request(x, 999.0, 10.0),
                  phase::PRResidualHelmholtzStatus::
                      OUTSIDE_CERTIFIED_TEMPERATURE_INTERVAL,
                  "low temperature");
    assertFailure(request(x, 6000.1, 10.0),
                  phase::PRResidualHelmholtzStatus::
                      OUTSIDE_CERTIFIED_TEMPERATURE_INTERVAL,
                  "high temperature");
    assertFailure(request(x, 1175.0, 0.0),
                  phase::PRResidualHelmholtzStatus::INVALID_DENSITY,
                  "zero density");
    auto nonnormalized = request(x, 1175.0, 10.0);
    nonnormalized.moleFractions[0] += 1.0e-5;
    assertFailure(nonnormalized,
                  phase::PRResidualHelmholtzStatus::COMPOSITION_NOT_NORMALIZED,
                  "non-normalized composition");
    auto negative = request(x, 1175.0, 10.0);
    negative.moleFractions[0] = -1.0e-4;
    assertFailure(negative,
                  phase::PRResidualHelmholtzStatus::NEGATIVE_COMPOSITION,
                  "negative composition");
    auto badIdeal = request(x, 1175.0, 10.0);
    badIdeal.ideal.cv = -1.0;
    assertFailure(badIdeal,
                  phase::PRResidualHelmholtzStatus::INVALID_IDEAL_CALORIC_STATE,
                  "invalid ideal caloric state");
    auto tooDense = request(x, 1000.0, 1.0e6);
    const auto denseResult = model.evaluateDiagnostic(tooDense);
    require(!denseResult.success() && !denseResult.state.has_value(),
            "close-packed state fails atomically");
    auto tooSmallCv = request(x, 1175.0, 10.0, 1.0e-12);
    const auto cvResult = model.evaluateDiagnostic(tooSmallCv);
    require(!cvResult.success() && !cvResult.state.has_value(),
            "nonpositive total Cv fails atomically");
}

void testPrototypeParity()
{
    const phase::PRResidualHelmholtzAuthority model = authority();
    const std::vector<double> x = referenceComposition();
    const double mw = molecularWeight(x);
    std::vector<diagnostic::Component> components;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        if (x[i] == 0.0) continue;
        const Foam::SpeciesCriticalDB critical =
            Foam::SpeciesCriticalDB::createDefault();
        const double massFraction = x[i]
            *critical.get(names()[i]).MW/mw;
        components.push_back(prototypeComponent(names()[i], massFraction));
    }
    const auto diagnosticState = diagnostic::handoffResidualPotential(
        components, prototypeInteractions(), 1175.0, mw/10.0);
    const auto candidate = model.evaluateDiagnostic(
        request(x, 1175.0, 10.0));
    require(diagnosticState.valid && candidate.success(),
            "prototype parity states valid");
    const auto& token = candidate.state->token;
    require(relativeError(token.residualHelmholtz, diagnosticState.value)
                < 1.0e-12, "residual potential parity");
    require(relativeError(token.dResidualHelmholtzDT, diagnosticState.temperature)
                < 1.0e-12, "temperature derivative parity");
    require(relativeError(token.d2ResidualHelmholtzDT2,
                          diagnosticState.temperature2) < 1.0e-12,
            "second temperature derivative parity");
    require(relativeError(token.dResidualHelmholtzDv, diagnosticState.volume)
                < 1.0e-12, "volume derivative parity");
    require(relativeError(token.d2ResidualHelmholtzDv2,
                          diagnosticState.volume2) < 1.0e-12,
            "second volume derivative parity");
    require(relativeError(token.d2ResidualHelmholtzDTDv,
                          diagnosticState.temperatureVolume) < 1.0e-12,
            "mixed derivative parity");
}

} // namespace

int main()
{
    std::cout << std::setprecision(17);
    testCapabilityAndConfiguration();
    testHandoffAndHighBranch();
    testPotentialDerivatives();
    testThermodynamicIdentities();
    testDerivedDerivativeConsistency();
    testFailClosedInputs();
    testPrototypeParity();
    std::cout << "PR_RESIDUAL_HELMHOLTZ_PASS checks=" << checks << '\n';
    return 0;
}
