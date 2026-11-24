#include "PRResidualHelmholtzAuthority.H"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <set>
#include <utility>

namespace Foam
{
namespace phaseEquilibrium
{

namespace
{

double quietNaN()
{
    return std::numeric_limits<double>::quiet_NaN();
}

} // End anonymous namespace

const char* residualHelmholtzStatusName(PRResidualHelmholtzStatus status)
{
    switch (status)
    {
        case PRResidualHelmholtzStatus::SUCCESS: return "SUCCESS";
        case PRResidualHelmholtzStatus::CAPABILITY_DISABLED:
            return "CAPABILITY_DISABLED";
        case PRResidualHelmholtzStatus::INVALID_CONFIGURATION:
            return "INVALID_CONFIGURATION";
        case PRResidualHelmholtzStatus::EMPTY_SPECIES_SET:
            return "EMPTY_SPECIES_SET";
        case PRResidualHelmholtzStatus::DUPLICATE_SPECIES:
            return "DUPLICATE_SPECIES";
        case PRResidualHelmholtzStatus::UNKNOWN_SPECIES:
            return "UNKNOWN_SPECIES";
        case PRResidualHelmholtzStatus::INVALID_CRITICAL_PROPERTIES:
            return "INVALID_CRITICAL_PROPERTIES";
        case PRResidualHelmholtzStatus::MISSING_BINARY_INTERACTION:
            return "MISSING_BINARY_INTERACTION";
        case PRResidualHelmholtzStatus::INVALID_BINARY_INTERACTION:
            return "INVALID_BINARY_INTERACTION";
        case PRResidualHelmholtzStatus::ALPHA_ROOT_INTERVAL_EXIT:
            return "ALPHA_ROOT_INTERVAL_EXIT";
        case PRResidualHelmholtzStatus::INVALID_TEMPERATURE:
            return "INVALID_TEMPERATURE";
        case PRResidualHelmholtzStatus::OUTSIDE_CERTIFIED_TEMPERATURE_INTERVAL:
            return "OUTSIDE_CERTIFIED_TEMPERATURE_INTERVAL";
        case PRResidualHelmholtzStatus::INVALID_DENSITY:
            return "INVALID_DENSITY";
        case PRResidualHelmholtzStatus::INVALID_COMPOSITION_SIZE:
            return "INVALID_COMPOSITION_SIZE";
        case PRResidualHelmholtzStatus::NONFINITE_COMPOSITION:
            return "NONFINITE_COMPOSITION";
        case PRResidualHelmholtzStatus::NEGATIVE_COMPOSITION:
            return "NEGATIVE_COMPOSITION";
        case PRResidualHelmholtzStatus::COMPOSITION_NOT_NORMALIZED:
            return "COMPOSITION_NOT_NORMALIZED";
        case PRResidualHelmholtzStatus::INVALID_IDEAL_CALORIC_STATE:
            return "INVALID_IDEAL_CALORIC_STATE";
        case PRResidualHelmholtzStatus::INVALID_MOLAR_VOLUME:
            return "INVALID_MOLAR_VOLUME";
        case PRResidualHelmholtzStatus::NONFINITE_POTENTIAL:
            return "NONFINITE_POTENTIAL";
        case PRResidualHelmholtzStatus::NONPOSITIVE_PRESSURE:
            return "NONPOSITIVE_PRESSURE";
        case PRResidualHelmholtzStatus::NONPOSITIVE_MOLAR_CV:
            return "NONPOSITIVE_MOLAR_CV";
        case PRResidualHelmholtzStatus::MECHANICAL_INSTABILITY:
            return "MECHANICAL_INSTABILITY";
        case PRResidualHelmholtzStatus::NONPOSITIVE_MOLAR_CP:
            return "NONPOSITIVE_MOLAR_CP";
        case PRResidualHelmholtzStatus::NONPOSITIVE_SOUND_SPEED_SQUARED:
            return "NONPOSITIVE_SOUND_SPEED_SQUARED";
    }
    return "UNKNOWN";
}

bool PRResidualHelmholtzAuthority::finitePositive(double value)
{
    return std::isfinite(value) && value > 0.0;
}

double PRResidualHelmholtzAuthority::kappa(
    const CriticalProperties& properties)
{
    const double omega = properties.omega;
    return 0.37464 + 1.54226*omega - 0.26992*omega*omega;
}

PRResidualHelmholtzAuthority::RootDerivatives
PRResidualHelmholtzAuthority::alphaRoot(
    const CriticalProperties& properties,
    double temperature)
{
    RootDerivatives result;
    if (!finitePositive(temperature)
        || !finitePositive(properties.Tc)
        || !std::isfinite(properties.omega))
        return result;

    const double coefficient = kappa(properties);
    const double rootTemperature = std::sqrt(temperature);
    const double rootCriticalTemperature = std::sqrt(properties.Tc);
    result.value = 1.0 + coefficient
        *(1.0 - rootTemperature/rootCriticalTemperature);
    result.first = -coefficient
        /(2.0*rootTemperature*rootCriticalTemperature);
    result.second = coefficient
        /(4.0*rootCriticalTemperature*temperature*rootTemperature);
    return result;
}

PRResidualHelmholtzAuthority::FadeDerivatives
PRResidualHelmholtzAuthority::fade(double temperature)
{
    if (temperature <= handoffTemperatureLower)
        return {1.0, 0.0, 0.0, 0.0};
    if (temperature >= handoffTemperatureUpper)
        return {0.0, 0.0, 0.0, 0.0};

    const double width =
        handoffTemperatureUpper - handoffTemperatureLower;
    const double s = (temperature - handoffTemperatureLower)/width;
    const double q = 1.0 - s;
    const double s2 = s*s;
    const double q2 = q*q;

    FadeDerivatives result;
    result.value = q2*q2*(35.0 - 84.0*q + 70.0*q2 - 20.0*q2*q);
    result.first = -140.0*s2*s*q2*q/width;
    result.second = -420.0*s2*q2*(1.0 - 2.0*s)/(width*width);
    result.third = -840.0*s*q*(1.0 - 5.0*s + 5.0*s2)
        /(width*width*width);
    return result;
}

PRResidualHelmholtzAuthority::PRResidualHelmholtzAuthority(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& binaryInteraction,
    const std::vector<std::string>& speciesNames)
:
    speciesNames_(speciesNames),
    configurationStatus_(PRResidualHelmholtzStatus::INVALID_CONFIGURATION),
    minimumCertifiedAlphaRoot_(quietNaN())
{
    if (speciesNames_.empty())
    {
        configurationStatus_ = PRResidualHelmholtzStatus::EMPTY_SPECIES_SET;
        return;
    }

    std::set<std::string> uniqueNames;
    properties_.reserve(speciesNames_.size());
    minimumCertifiedAlphaRoot_ = std::numeric_limits<double>::infinity();
    for (const std::string& name : speciesNames_)
    {
        if (name.empty() || !uniqueNames.insert(name).second)
        {
            configurationStatus_ = PRResidualHelmholtzStatus::DUPLICATE_SPECIES;
            limitingSpecies_ = name;
            return;
        }
        if (!criticalData.has(name))
        {
            configurationStatus_ = PRResidualHelmholtzStatus::UNKNOWN_SPECIES;
            limitingSpecies_ = name;
            return;
        }

        const CriticalProperties properties = criticalData.get(name);
        if (!finitePositive(properties.Tc)
            || !finitePositive(properties.Pc)
            || !finitePositive(properties.MW)
            || !std::isfinite(properties.omega))
        {
            configurationStatus_ =
                PRResidualHelmholtzStatus::INVALID_CRITICAL_PROPERTIES;
            limitingSpecies_ = name;
            return;
        }

        const RootDerivatives lower =
            alphaRoot(properties, handoffTemperatureLower);
        const RootDerivatives upper =
            alphaRoot(properties, handoffTemperatureUpper);
        const double minimumRoot = std::min(lower.value, upper.value);
        if (!std::isfinite(lower.value) || !std::isfinite(upper.value)
            || minimumRoot < minimumAlphaRootMargin)
        {
            configurationStatus_ =
                PRResidualHelmholtzStatus::ALPHA_ROOT_INTERVAL_EXIT;
            limitingSpecies_ = name;
            minimumCertifiedAlphaRoot_ = minimumRoot;
            return;
        }
        if (minimumRoot < minimumCertifiedAlphaRoot_)
        {
            minimumCertifiedAlphaRoot_ = minimumRoot;
            limitingSpecies_ = name;
        }
        properties_.push_back(properties);
    }

    const std::size_t count = speciesNames_.size();
    binaryInteraction_.assign(count*count, 0.0);
    for (std::size_t i = 0; i < count; ++i)
    {
        for (std::size_t j = i + 1; j < count; ++j)
        {
            const BinaryInteractionEntry* entry = binaryInteraction.find(
                speciesNames_[i], speciesNames_[j]);
            if (!entry)
            {
                configurationStatus_ =
                    PRResidualHelmholtzStatus::MISSING_BINARY_INTERACTION;
                limitingSpecies_ = speciesNames_[i] + "|" + speciesNames_[j];
                return;
            }
            if (!std::isfinite(entry->value) || !(1.0 - entry->value > 0.0))
            {
                configurationStatus_ =
                    PRResidualHelmholtzStatus::INVALID_BINARY_INTERACTION;
                limitingSpecies_ = speciesNames_[i] + "|" + speciesNames_[j];
                return;
            }
            binaryInteraction_[i*count + j] = entry->value;
            binaryInteraction_[j*count + i] = entry->value;
        }
    }

    configurationStatus_ = PRResidualHelmholtzStatus::SUCCESS;
}

PRResidualHelmholtzStatus
PRResidualHelmholtzAuthority::validateComposition(
    const std::vector<double>& moleFractions) const
{
    if (moleFractions.size() != speciesNames_.size())
        return PRResidualHelmholtzStatus::INVALID_COMPOSITION_SIZE;

    double sum = 0.0;
    for (double fraction : moleFractions)
    {
        if (!std::isfinite(fraction))
            return PRResidualHelmholtzStatus::NONFINITE_COMPOSITION;
        if (fraction < 0.0)
            return PRResidualHelmholtzStatus::NEGATIVE_COMPOSITION;
        sum += fraction;
    }
    if (!std::isfinite(sum))
        return PRResidualHelmholtzStatus::NONFINITE_COMPOSITION;
    if (std::abs(sum - 1.0) > compositionTolerance)
        return PRResidualHelmholtzStatus::COMPOSITION_NOT_NORMALIZED;
    return PRResidualHelmholtzStatus::SUCCESS;
}

PRResidualHelmholtzResult PRResidualHelmholtzAuthority::evaluate(
    const PRResidualHelmholtzRequest&) const
{
    PRResidualHelmholtzResult result;
    result.status = PRResidualHelmholtzStatus::CAPABILITY_DISABLED;
    return result;
}

PRResidualHelmholtzResult
PRResidualHelmholtzAuthority::evaluateDiagnostic(
    const PRResidualHelmholtzRequest& request) const
{
    return evaluateImpl(request);
}

PRResidualHelmholtzResult PRResidualHelmholtzAuthority::evaluateImpl(
    const PRResidualHelmholtzRequest& request) const
{
    PRResidualHelmholtzResult result;
    result.status = configurationStatus_;
    result.limitingSpecies = limitingSpecies_;
    if (!configured()) return result;

    if (!finitePositive(request.temperature))
    {
        result.status = PRResidualHelmholtzStatus::INVALID_TEMPERATURE;
        return result;
    }
    if (request.temperature < handoffTemperatureLower
        || request.temperature > certifiedTemperatureMaximum)
    {
        result.status = PRResidualHelmholtzStatus::
            OUTSIDE_CERTIFIED_TEMPERATURE_INTERVAL;
        return result;
    }
    if (!finitePositive(request.density))
    {
        result.status = PRResidualHelmholtzStatus::INVALID_DENSITY;
        return result;
    }
    result.status = validateComposition(request.moleFractions);
    if (result.status != PRResidualHelmholtzStatus::SUCCESS) return result;
    if (!std::isfinite(request.ideal.internalEnergy)
        || !std::isfinite(request.ideal.entropy)
        || !finitePositive(request.ideal.cv))
    {
        result.status =
            PRResidualHelmholtzStatus::INVALID_IDEAL_CALORIC_STATE;
        return result;
    }

    const std::size_t count = speciesNames_.size();
    double molecularWeight = 0.0;
    for (std::size_t i = 0; i < count; ++i)
        molecularWeight += request.moleFractions[i]*properties_[i].MW;
    const double molarVolume = molecularWeight/request.density;
    if (!finitePositive(molecularWeight) || !finitePositive(molarVolume))
    {
        result.status = PRResidualHelmholtzStatus::INVALID_MOLAR_VOLUME;
        return result;
    }

    PRResidualHelmholtzToken token;
    const FadeDerivatives switchValue = fade(request.temperature);
    token.fade = switchValue.value;
    token.dFadeDT = switchValue.first;
    token.d2FadeDT2 = switchValue.second;
    token.d3FadeDT3 = switchValue.third;

    if (request.temperature < handoffTemperatureUpper)
    {
        token.alphaEvaluated = true;
        std::vector<RootDerivatives> roots(count);
        std::vector<double> rootScales(count, 0.0);
        for (std::size_t i = 0; i < count; ++i)
        {
            roots[i] = alphaRoot(properties_[i], request.temperature);
            if (!std::isfinite(roots[i].value)
                || !std::isfinite(roots[i].first)
                || !std::isfinite(roots[i].second)
                || roots[i].value < minimumAlphaRootMargin)
            {
                result.status =
                    PRResidualHelmholtzStatus::ALPHA_ROOT_INTERVAL_EXIT;
                result.limitingSpecies = speciesNames_[i];
                return result;
            }
            rootScales[i] = std::sqrt(
                0.45724*universalGasConstant*universalGasConstant
                *properties_[i].Tc*properties_[i].Tc/properties_[i].Pc);
            token.coVolume += request.moleFractions[i]
                *0.07780*universalGasConstant
                *properties_[i].Tc/properties_[i].Pc;
        }

        for (std::size_t i = 0; i < count; ++i)
        {
            for (std::size_t j = 0; j < count; ++j)
            {
                const double scale = request.moleFractions[i]
                    *request.moleFractions[j]*rootScales[i]*rootScales[j]
                    *(1.0 - binaryInteraction_[i*count + j]);
                token.attraction += scale*roots[i].value*roots[j].value;
                token.dAttractionDT += scale
                    *(roots[i].first*roots[j].value
                      + roots[i].value*roots[j].first);
                token.d2AttractionDT2 += scale
                    *(roots[i].second*roots[j].value
                      + 2.0*roots[i].first*roots[j].first
                      + roots[i].value*roots[j].second);
            }
        }

        if (!finitePositive(token.coVolume)
            || !finitePositive(token.attraction)
            || !std::isfinite(token.dAttractionDT)
            || !std::isfinite(token.d2AttractionDT2)
            || !(molarVolume > token.coVolume))
        {
            result.status = PRResidualHelmholtzStatus::INVALID_MOLAR_VOLUME;
            return result;
        }

        const double rootTwo = std::sqrt(2.0);
        const double denominator = molarVolume*molarVolume
            + 2.0*token.coVolume*molarVolume
            - token.coVolume*token.coVolume;
        const double repulsiveArgument = 1.0 - token.coVolume/molarVolume;
        const double attractionNumerator = molarVolume
            + (1.0 + rootTwo)*token.coVolume;
        const double attractionDenominator = molarVolume
            + (1.0 - rootTwo)*token.coVolume;
        if (!finitePositive(denominator)
            || !finitePositive(repulsiveArgument)
            || !finitePositive(attractionNumerator)
            || !finitePositive(attractionDenominator))
        {
            result.status = PRResidualHelmholtzStatus::INVALID_MOLAR_VOLUME;
            return result;
        }

        const double logVolume = std::log(repulsiveArgument);
        const double logAttraction = std::log(
            attractionNumerator/attractionDenominator);
        const double attractionScale =
            1.0/(2.0*rootTwo*token.coVolume);

        const double base = -universalGasConstant*request.temperature*logVolume
            - token.attraction*attractionScale*logAttraction;
        const double baseTemperature = -universalGasConstant*logVolume
            - token.dAttractionDT*attractionScale*logAttraction;
        const double baseTemperature2 = -token.d2AttractionDT2
            *attractionScale*logAttraction;
        const double baseVolume = -universalGasConstant*request.temperature
            *token.coVolume/(molarVolume*(molarVolume - token.coVolume))
            + token.attraction/denominator;
        const double baseVolume2 = universalGasConstant*request.temperature
            *token.coVolume*(2.0*molarVolume - token.coVolume)
            /(molarVolume*molarVolume
              *(molarVolume - token.coVolume)*(molarVolume - token.coVolume))
            - 2.0*token.attraction*(molarVolume + token.coVolume)
              /(denominator*denominator);
        const double baseTemperatureVolume = -universalGasConstant
            *token.coVolume/(molarVolume*(molarVolume - token.coVolume))
            + token.dAttractionDT/denominator;

        token.residualHelmholtz = switchValue.value*base;
        token.dResidualHelmholtzDT = switchValue.value*baseTemperature
            + switchValue.first*base;
        token.d2ResidualHelmholtzDT2 = switchValue.value*baseTemperature2
            + 2.0*switchValue.first*baseTemperature
            + switchValue.second*base;
        token.dResidualHelmholtzDv = switchValue.value*baseVolume;
        token.d2ResidualHelmholtzDv2 = switchValue.value*baseVolume2;
        token.d2ResidualHelmholtzDTDv =
            switchValue.value*baseTemperatureVolume
            + switchValue.first*baseVolume;
    }

    const double tokenValues[] = {
        token.residualHelmholtz,
        token.dResidualHelmholtzDT,
        token.d2ResidualHelmholtzDT2,
        token.dResidualHelmholtzDv,
        token.d2ResidualHelmholtzDv2,
        token.d2ResidualHelmholtzDTDv,
        token.attraction,
        token.dAttractionDT,
        token.d2AttractionDT2,
        token.coVolume,
        token.fade,
        token.dFadeDT,
        token.d2FadeDT2,
        token.d3FadeDT3};
    if (std::any_of(
            std::begin(tokenValues), std::end(tokenValues),
            [](double value) { return !std::isfinite(value); }))
    {
        result.status = PRResidualHelmholtzStatus::NONFINITE_POTENTIAL;
        return result;
    }

    PRResidualHelmholtzState candidate;
    candidate.moleFractions = request.moleFractions;
    candidate.temperature = request.temperature;
    candidate.density = request.density;
    candidate.molecularWeight = molecularWeight;
    candidate.molarVolume = molarVolume;
    candidate.token = token;

    candidate.pressure = universalGasConstant*request.temperature/molarVolume
        - token.dResidualHelmholtzDv;
    if (!finitePositive(candidate.pressure))
    {
        result.status = PRResidualHelmholtzStatus::NONPOSITIVE_PRESSURE;
        return result;
    }

    const double internalEnergyMolar = request.ideal.internalEnergy
        + token.residualHelmholtz
        - request.temperature*token.dResidualHelmholtzDT;
    const double entropyMolar = request.ideal.entropy
        - token.dResidualHelmholtzDT;
    const double cvMolar = request.ideal.cv
        - request.temperature*token.d2ResidualHelmholtzDT2;
    if (!std::isfinite(internalEnergyMolar)
        || !std::isfinite(entropyMolar)
        || !std::isfinite(cvMolar))
    {
        result.status = PRResidualHelmholtzStatus::NONFINITE_POTENTIAL;
        return result;
    }
    if (!(cvMolar >= minimumMolarCv))
    {
        result.status = PRResidualHelmholtzStatus::NONPOSITIVE_MOLAR_CV;
        return result;
    }

    candidate.pressureTemperatureDerivative =
        universalGasConstant/molarVolume
        - token.d2ResidualHelmholtzDTDv;
    candidate.pressureMolarVolumeDerivative =
        -universalGasConstant*request.temperature/(molarVolume*molarVolume)
        - token.d2ResidualHelmholtzDv2;
    if (!std::isfinite(candidate.pressureTemperatureDerivative)
        || !std::isfinite(candidate.pressureMolarVolumeDerivative)
        || !(candidate.pressureMolarVolumeDerivative < 0.0))
    {
        result.status = PRResidualHelmholtzStatus::MECHANICAL_INSTABILITY;
        return result;
    }

    const double cpMolar = cvMolar
        - request.temperature*candidate.pressureTemperatureDerivative
          *candidate.pressureTemperatureDerivative
          /candidate.pressureMolarVolumeDerivative;
    if (!finitePositive(cpMolar))
    {
        result.status = PRResidualHelmholtzStatus::NONPOSITIVE_MOLAR_CP;
        return result;
    }

    candidate.isothermalPressureDensityDerivative =
        -candidate.pressureMolarVolumeDerivative*molarVolume*molarVolume
        /molecularWeight;
    const double isentropicPressureMolarVolumeDerivative =
        candidate.pressureMolarVolumeDerivative
        - request.temperature*candidate.pressureTemperatureDerivative
          *candidate.pressureTemperatureDerivative/cvMolar;
    candidate.soundSpeedSquared = -molarVolume*molarVolume/molecularWeight
        *isentropicPressureMolarVolumeDerivative;
    if (!std::isfinite(candidate.isothermalPressureDensityDerivative)
        || !std::isfinite(candidate.soundSpeedSquared)
        || !(candidate.soundSpeedSquared >= minimumSoundSpeedSquared))
    {
        result.status =
            PRResidualHelmholtzStatus::NONPOSITIVE_SOUND_SPEED_SQUARED;
        return result;
    }

    candidate.internalEnergy = internalEnergyMolar/molecularWeight;
    candidate.enthalpy = (internalEnergyMolar
        + candidate.pressure*molarVolume)/molecularWeight;
    candidate.entropy = entropyMolar/molecularWeight;
    candidate.cv = cvMolar/molecularWeight;
    candidate.cp = cpMolar/molecularWeight;
    if (!std::isfinite(candidate.internalEnergy)
        || !std::isfinite(candidate.enthalpy)
        || !std::isfinite(candidate.entropy)
        || !finitePositive(candidate.cv)
        || !finitePositive(candidate.cp))
    {
        result.status = PRResidualHelmholtzStatus::NONFINITE_POTENTIAL;
        return result;
    }

    result.status = PRResidualHelmholtzStatus::SUCCESS;
    result.limitingSpecies.clear();
    result.state = std::move(candidate);
    return result;
}

} // End namespace phaseEquilibrium
} // End namespace Foam
