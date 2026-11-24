#ifndef RFAM_DIAGNOSTIC_PR_HIGH_TEMPERATURE_RESIDUAL_HANDOFF_PROTOTYPE_H
#define RFAM_DIAGNOSTIC_PR_HIGH_TEMPERATURE_RESIDUAL_HANDOFF_PROTOTYPE_H

// Offline diagnostic prototype only. This header is intentionally outside the
// production thermo tree and is not registered in CMake.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace rfamHighTemperatureDiagnostic
{

constexpr double universalGasConstant = 8.314462618;
constexpr double principalRootDerivativeThreshold = 1.0e-30;
constexpr double handoffTemperatureLower = 1000.0;
constexpr double handoffTemperatureUpper = 1325.0;

struct Component
{
    std::string name;
    double massFraction = 0.0;
    double molecularWeight = std::numeric_limits<double>::quiet_NaN();
    double criticalTemperature = std::numeric_limits<double>::quiet_NaN();
    double criticalPressure = std::numeric_limits<double>::quiet_NaN();
    double acentricFactor = std::numeric_limits<double>::quiet_NaN();
    bool residualActive = true;
};

using InteractionMap =
    std::map<std::pair<std::string, std::string>, double>;

inline std::pair<std::string, std::string> interactionKey(
    const std::string& first,
    const std::string& second)
{
    return first < second
        ? std::make_pair(first, second)
        : std::make_pair(second, first);
}

inline double interactionValue(
    const InteractionMap& interactions,
    const std::string& first,
    const std::string& second)
{
    const auto iterator = interactions.find(interactionKey(first, second));
    return iterator == interactions.end() ? 0.0 : iterator->second;
}

struct SwitchDerivatives
{
    double value = 0.0;
    double first = 0.0;
    double second = 0.0;
    double third = 0.0;
};

// Unique degree-seven fade satisfying w(0)=1, w(1)=0 and zero first,
// second, and third derivatives at both ends.
inline SwitchDerivatives residualFade(double temperature)
{
    if (temperature <= handoffTemperatureLower)
        return {1.0, 0.0, 0.0, 0.0};
    if (temperature >= handoffTemperatureUpper)
        return {0.0, 0.0, 0.0, 0.0};

    const double width =
        handoffTemperatureUpper - handoffTemperatureLower;
    const double s = (temperature - handoffTemperatureLower)/width;
    const double s2 = s*s;
    const double s3 = s2*s;
    const double s4 = s3*s;
    const double s5 = s4*s;
    const double s6 = s5*s;
    const double s7 = s6*s;

    SwitchDerivatives result;
    result.value = 1.0 - 35.0*s4 + 84.0*s5 - 70.0*s6 + 20.0*s7;
    result.first =
        (-140.0*s3 + 420.0*s4 - 420.0*s5 + 140.0*s6)/width;
    result.second =
        (-420.0*s2 + 1680.0*s3 - 2100.0*s4 + 840.0*s5)
        /(width*width);
    result.third =
        (-840.0*s + 5040.0*s2 - 8400.0*s3 + 4200.0*s4)
        /(width*width*width);
    return result;
}

struct RootDerivatives
{
    double value = std::numeric_limits<double>::quiet_NaN();
    double first = std::numeric_limits<double>::quiet_NaN();
    double second = std::numeric_limits<double>::quiet_NaN();
};

inline double kappa(const Component& component)
{
    const double omega = component.acentricFactor;
    return 0.37464 + 1.54226*omega - 0.26992*omega*omega;
}

inline double alphaZeroTemperature(const Component& component)
{
    if (!component.residualActive) return std::numeric_limits<double>::infinity();
    const double value = kappa(component);
    if (!(std::isfinite(value) && value > 0.0
          && std::isfinite(component.criticalTemperature)
          && component.criticalTemperature > 0.0))
        return std::numeric_limits<double>::infinity();
    const double ratio = (1.0 + value)/value;
    return component.criticalTemperature*ratio*ratio;
}

inline double alphaScale(const Component& component)
{
    return 0.45724*universalGasConstant*universalGasConstant
        *component.criticalTemperature*component.criticalTemperature
        /component.criticalPressure;
}

inline double coVolume(const Component& component)
{
    return 0.07780*universalGasConstant*component.criticalTemperature
        /component.criticalPressure;
}

inline RootDerivatives originalAlphaRoot(
    const Component& component,
    double temperature)
{
    const double value = kappa(component);
    const double rootTemperature = std::sqrt(temperature);
    const double rootCritical = std::sqrt(component.criticalTemperature);
    RootDerivatives result;
    result.value = 1.0 + value*(1.0 - rootTemperature/rootCritical);
    result.first = -value/(2.0*rootTemperature*rootCritical);
    result.second = value/(4.0*rootCritical*temperature*rootTemperature);
    return result;
}

// Secondary-route prototype: a positive log-quadratic continuation of the
// alpha root. It matches g, dg/dT, and d2g/dT2 at 1000 K, hence is C2 there,
// but has no high-temperature physical authority.
inline RootDerivatives constrainedAlphaRoot(
    const Component& component,
    double temperature)
{
    if (temperature <= handoffTemperatureLower)
        return originalAlphaRoot(component, temperature);

    const RootDerivatives anchor =
        originalAlphaRoot(component, handoffTemperatureLower);
    if (!(std::isfinite(anchor.value) && anchor.value > 0.0))
        return {};

    const double logFirst = anchor.first/anchor.value;
    const double logSecond = anchor.second/anchor.value
        - logFirst*logFirst;
    const double delta = temperature - handoffTemperatureLower;
    const double exponent = std::log(anchor.value) + logFirst*delta
        + 0.5*logSecond*delta*delta;
    const double logDerivative = logFirst + logSecond*delta;

    RootDerivatives result;
    result.value = std::exp(exponent);
    result.first = result.value*logDerivative;
    result.second = result.value
        *(logDerivative*logDerivative + logSecond);
    return result;
}

enum class RootPolicy
{
    PositiveGuarded,
    PrincipalAbsolute,
    SignedDiagnostic,
    ConstrainedPositive
};

struct MixtureParameters
{
    bool valid = false;
    std::string failure;
    double molecularWeight = std::numeric_limits<double>::quiet_NaN();
    double a = std::numeric_limits<double>::quiet_NaN();
    double aTemperature = std::numeric_limits<double>::quiet_NaN();
    double aTemperature2 = std::numeric_limits<double>::quiet_NaN();
    double b = std::numeric_limits<double>::quiet_NaN();
    double lowestActiveZero = std::numeric_limits<double>::infinity();
    double minimumActiveRoot = std::numeric_limits<double>::infinity();
};

inline bool finitePositive(double value)
{
    return std::isfinite(value) && value > 0.0;
}

inline MixtureParameters mixtureParameters(
    const std::vector<Component>& components,
    const InteractionMap& interactions,
    double temperature,
    RootPolicy policy)
{
    MixtureParameters result;
    if (!finitePositive(temperature) || components.empty())
    {
        result.failure = "invalid-request";
        return result;
    }

    double massSum = 0.0;
    double moleSum = 0.0;
    for (const Component& component : components)
    {
        if (!std::isfinite(component.massFraction)
            || component.massFraction < 0.0
            || !finitePositive(component.molecularWeight))
        {
            result.failure = "invalid-composition";
            return result;
        }
        massSum += component.massFraction;
        moleSum += component.massFraction/component.molecularWeight;
        if (component.massFraction > 0.0 && component.residualActive)
        {
            if (!finitePositive(component.criticalTemperature)
                || !finitePositive(component.criticalPressure)
                || !std::isfinite(component.acentricFactor))
            {
                result.failure = "missing-critical-data";
                return result;
            }
            result.lowestActiveZero = std::min(
                result.lowestActiveZero,
                alphaZeroTemperature(component));
        }
    }
    if (!std::isfinite(massSum) || std::abs(massSum - 1.0) > 1.0e-12
        || !finitePositive(moleSum))
    {
        result.failure = "invalid-composition";
        return result;
    }

    result.molecularWeight = 1.0/moleSum;
    std::vector<double> moleFractions(components.size(), 0.0);
    std::vector<RootDerivatives> roots(components.size());
    std::vector<double> rootScales(components.size(), 0.0);
    for (std::size_t i = 0; i < components.size(); ++i)
    {
        const Component& component = components[i];
        moleFractions[i] = component.massFraction
            /component.molecularWeight/moleSum;
        if (!component.residualActive || component.massFraction == 0.0)
            continue;

        roots[i] = policy == RootPolicy::ConstrainedPositive
            ? constrainedAlphaRoot(component, temperature)
            : originalAlphaRoot(component, temperature);
        if (!std::isfinite(roots[i].value)
            || !std::isfinite(roots[i].first)
            || !std::isfinite(roots[i].second))
        {
            result.failure = "invalid-alpha-root";
            return result;
        }
        if (policy == RootPolicy::PositiveGuarded && roots[i].value <= 0.0)
        {
            result.failure = "alpha-domain-exit";
            return result;
        }
        result.minimumActiveRoot = std::min(
            result.minimumActiveRoot, roots[i].value);
        rootScales[i] = std::sqrt(alphaScale(component));
    }

    result.b = 0.0;
    result.a = 0.0;
    result.aTemperature = 0.0;
    result.aTemperature2 = 0.0;
    for (std::size_t i = 0; i < components.size(); ++i)
    {
        if (components[i].residualActive)
            result.b += moleFractions[i]*coVolume(components[i]);
        for (std::size_t j = 0; j < components.size(); ++j)
        {
            if (!components[i].residualActive
                || !components[j].residualActive)
                continue;
            RootDerivatives first = roots[i];
            RootDerivatives second = roots[j];
            if (policy == RootPolicy::PrincipalAbsolute)
            {
                const double signFirst = first.value < 0.0 ? -1.0 : 1.0;
                const double signSecond = second.value < 0.0 ? -1.0 : 1.0;
                first.value = std::abs(first.value);
                first.first *= signFirst;
                first.second *= signFirst;
                second.value = std::abs(second.value);
                second.first *= signSecond;
                second.second *= signSecond;
            }
            const double scale = rootScales[i]*rootScales[j]
                *(1.0 - interactionValue(
                    interactions, components[i].name, components[j].name))
                *moleFractions[i]*moleFractions[j];
            result.a += scale*first.value*second.value;
            result.aTemperature += scale
                *(first.first*second.value
                  + first.value*second.first);
            result.aTemperature2 += scale
                *(first.second*second.value
                  + 2.0*first.first*second.first
                  + first.value*second.second);
        }
    }

    if (!finitePositive(result.molecularWeight)
        || !finitePositive(result.b)
        || !std::isfinite(result.a)
        || !std::isfinite(result.aTemperature)
        || !std::isfinite(result.aTemperature2))
    {
        result.failure = "invalid-mixture-parameters";
        return result;
    }
    result.valid = true;
    result.failure.clear();
    return result;
}

struct PotentialDerivatives
{
    bool valid = false;
    std::string failure;
    double value = 0.0;
    double temperature = 0.0;
    double temperature2 = 0.0;
    double volume = 0.0;
    double volume2 = 0.0;
    double temperatureVolume = 0.0;
    double molecularWeight = std::numeric_limits<double>::quiet_NaN();
};

inline PotentialDerivatives prResidualPotential(
    const MixtureParameters& mixture,
    double temperature,
    double molarVolume)
{
    PotentialDerivatives result;
    result.molecularWeight = mixture.molecularWeight;
    if (!mixture.valid || !finitePositive(temperature)
        || !finitePositive(molarVolume) || molarVolume <= mixture.b)
    {
        result.failure = mixture.valid ? "invalid-volume" : mixture.failure;
        return result;
    }

    const double rootTwo = std::sqrt(2.0);
    const double denominator = molarVolume*molarVolume
        + 2.0*mixture.b*molarVolume - mixture.b*mixture.b;
    const double logVolume = std::log(1.0 - mixture.b/molarVolume);
    const double logAttraction = std::log(
        (molarVolume + (1.0 + rootTwo)*mixture.b)
        /(molarVolume + (1.0 - rootTwo)*mixture.b));
    const double attractionScale = 1.0/(2.0*rootTwo*mixture.b);

    result.value = -universalGasConstant*temperature*logVolume
        - mixture.a*attractionScale*logAttraction;
    result.temperature = -universalGasConstant*logVolume
        - mixture.aTemperature*attractionScale*logAttraction;
    result.temperature2 = -mixture.aTemperature2
        *attractionScale*logAttraction;
    result.volume = -universalGasConstant*temperature*mixture.b
        /(molarVolume*(molarVolume - mixture.b))
        + mixture.a/denominator;
    result.volume2 = universalGasConstant*temperature*mixture.b
        *(2.0*molarVolume - mixture.b)
        /(molarVolume*molarVolume
          *(molarVolume - mixture.b)*(molarVolume - mixture.b))
        - 2.0*mixture.a*(molarVolume + mixture.b)
          /(denominator*denominator);
    result.temperatureVolume = -universalGasConstant*mixture.b
        /(molarVolume*(molarVolume - mixture.b))
        + mixture.aTemperature/denominator;
    result.valid = std::isfinite(result.value)
        && std::isfinite(result.temperature)
        && std::isfinite(result.temperature2)
        && std::isfinite(result.volume)
        && std::isfinite(result.volume2)
        && std::isfinite(result.temperatureVolume);
    if (!result.valid) result.failure = "nonfinite-potential";
    return result;
}

inline PotentialDerivatives zeroResidualPotential(double molecularWeight)
{
    PotentialDerivatives result;
    result.valid = finitePositive(molecularWeight);
    result.molecularWeight = molecularWeight;
    if (!result.valid) result.failure = "invalid-molecular-weight";
    return result;
}

inline double mixtureMolecularWeight(
    const std::vector<Component>& components)
{
    double mass = 0.0;
    double moles = 0.0;
    for (const Component& component : components)
    {
        if (!std::isfinite(component.massFraction)
            || component.massFraction < 0.0
            || !finitePositive(component.molecularWeight))
            return std::numeric_limits<double>::quiet_NaN();
        mass += component.massFraction;
        moles += component.massFraction/component.molecularWeight;
    }
    return std::abs(mass - 1.0) <= 1.0e-12 && finitePositive(moles)
        ? 1.0/moles : std::numeric_limits<double>::quiet_NaN();
}

inline PotentialDerivatives handoffResidualPotential(
    const std::vector<Component>& components,
    const InteractionMap& interactions,
    double temperature,
    double molarVolume)
{
    const double molecularWeight = mixtureMolecularWeight(components);
    if (!finitePositive(molecularWeight))
        return zeroResidualPotential(molecularWeight);

    // The high side does not evaluate PR alpha at all. This is the critical
    // structural distinction from evaluating through a zero and multiplying
    // the resulting property by zero afterward.
    if (temperature >= handoffTemperatureUpper)
        return zeroResidualPotential(molecularWeight);

    const MixtureParameters mixture = mixtureParameters(
        components, interactions, temperature, RootPolicy::PositiveGuarded);
    PotentialDerivatives base =
        prResidualPotential(mixture, temperature, molarVolume);
    if (!base.valid) return base;

    const SwitchDerivatives fade = residualFade(temperature);
    PotentialDerivatives result;
    result.valid = true;
    result.molecularWeight = molecularWeight;
    result.value = fade.value*base.value;
    result.temperature = fade.value*base.temperature
        + fade.first*base.value;
    result.temperature2 = fade.value*base.temperature2
        + 2.0*fade.first*base.temperature
        + fade.second*base.value;
    result.volume = fade.value*base.volume;
    result.volume2 = fade.value*base.volume2;
    result.temperatureVolume = fade.value*base.temperatureVolume
        + fade.first*base.volume;
    return result;
}

inline PotentialDerivatives directResidualPotential(
    const std::vector<Component>& components,
    const InteractionMap& interactions,
    double temperature,
    double molarVolume,
    RootPolicy policy)
{
    const MixtureParameters mixture =
        mixtureParameters(components, interactions, temperature, policy);
    return prResidualPotential(mixture, temperature, molarVolume);
}

struct ThermodynamicState
{
    bool valid = false;
    std::string failure;
    double molecularWeight = std::numeric_limits<double>::quiet_NaN();
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double internalEnergy = std::numeric_limits<double>::quiet_NaN();
    double enthalpy = std::numeric_limits<double>::quiet_NaN();
    double entropy = std::numeric_limits<double>::quiet_NaN();
    double cv = std::numeric_limits<double>::quiet_NaN();
    double cp = std::numeric_limits<double>::quiet_NaN();
    double pressureTemperature = std::numeric_limits<double>::quiet_NaN();
    double pressureVolume = std::numeric_limits<double>::quiet_NaN();
    double dpdrhoIsothermal = std::numeric_limits<double>::quiet_NaN();
    double soundSpeedSquared = std::numeric_limits<double>::quiet_NaN();
    double residualPotential = std::numeric_limits<double>::quiet_NaN();
};

inline ThermodynamicState deriveState(
    const PotentialDerivatives& residual,
    double temperature,
    double molarVolume,
    double idealCv = 2.5*universalGasConstant)
{
    ThermodynamicState result;
    result.molecularWeight = residual.molecularWeight;
    if (!residual.valid || !finitePositive(temperature)
        || !finitePositive(molarVolume) || !finitePositive(idealCv)
        || !finitePositive(residual.molecularWeight))
    {
        result.failure = residual.valid
            ? "invalid-state-request" : residual.failure;
        return result;
    }

    result.residualPotential = residual.value;
    result.pressure = universalGasConstant*temperature/molarVolume
        - residual.volume;
    result.internalEnergy = idealCv*temperature
        + residual.value - temperature*residual.temperature;
    result.enthalpy = result.internalEnergy + result.pressure*molarVolume;
    result.entropy = idealCv*std::log(temperature)
        + universalGasConstant*std::log(molarVolume)
        - residual.temperature;
    result.cv = idealCv - temperature*residual.temperature2;
    result.pressureTemperature = universalGasConstant/molarVolume
        - residual.temperatureVolume;
    result.pressureVolume = -universalGasConstant*temperature
        /(molarVolume*molarVolume) - residual.volume2;
    if (!std::isfinite(result.pressureVolume)
        || result.pressureVolume == 0.0)
    {
        result.failure = "singular-mechanical-derivative";
        return result;
    }
    result.cp = result.cv - temperature*result.pressureTemperature
        *result.pressureTemperature/result.pressureVolume;
    result.dpdrhoIsothermal = -result.pressureVolume*molarVolume*molarVolume
        /residual.molecularWeight;
    if (result.cv == 0.0)
    {
        result.failure = "singular-caloric-derivative";
        return result;
    }
    const double pressureVolumeEntropy = result.pressureVolume
        - temperature*result.pressureTemperature
          *result.pressureTemperature/result.cv;
    result.soundSpeedSquared = -molarVolume*molarVolume
        /residual.molecularWeight*pressureVolumeEntropy;

    result.valid = std::isfinite(result.pressure)
        && std::isfinite(result.internalEnergy)
        && std::isfinite(result.enthalpy)
        && std::isfinite(result.entropy)
        && std::isfinite(result.cv) && std::isfinite(result.cp)
        && std::isfinite(result.dpdrhoIsothermal)
        && std::isfinite(result.soundSpeedSquared);
    if (!result.valid) result.failure = "nonfinite-derived-state";
    return result;
}

inline ThermodynamicState evaluateHandoff(
    const std::vector<Component>& components,
    const InteractionMap& interactions,
    double density,
    double temperature)
{
    const double molecularWeight = mixtureMolecularWeight(components);
    const double molarVolume = molecularWeight/density;
    return deriveState(
        handoffResidualPotential(
            components, interactions, temperature, molarVolume),
        temperature,
        molarVolume);
}

inline ThermodynamicState evaluateDirect(
    const std::vector<Component>& components,
    const InteractionMap& interactions,
    double density,
    double temperature,
    RootPolicy policy)
{
    const double molecularWeight = mixtureMolecularWeight(components);
    const double molarVolume = molecularWeight/density;
    return deriveState(
        directResidualPotential(
            components, interactions, temperature, molarVolume, policy),
        temperature,
        molarVolume);
}

} // namespace rfamHighTemperatureDiagnostic

#endif
