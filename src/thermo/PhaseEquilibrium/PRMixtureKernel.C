#include "PRMixtureKernel.H"

#include "CubicRoots.H"
#include "PengRobinsonAlpha.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

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

bool finitePositive(double value)
{
    return std::isfinite(value) && value > 0.0;
}

bool sameRoot(double lhs, double rhs)
{
    return std::abs(lhs - rhs)
        <= 1.0e-11*std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

} // End anonymous namespace

const char* statusName(PhaseEquilibriumStatus status)
{
    switch (status)
    {
        case PhaseEquilibriumStatus::SUCCESS: return "SUCCESS";
        case PhaseEquilibriumStatus::INVALID_CONFIGURATION:
            return "INVALID_CONFIGURATION";
        case PhaseEquilibriumStatus::UNKNOWN_SPECIES: return "UNKNOWN_SPECIES";
        case PhaseEquilibriumStatus::INVALID_CRITICAL_PROPERTIES:
            return "INVALID_CRITICAL_PROPERTIES";
        case PhaseEquilibriumStatus::MISSING_BINARY_INTERACTION:
            return "MISSING_BINARY_INTERACTION";
        case PhaseEquilibriumStatus::INVALID_PRESSURE: return "INVALID_PRESSURE";
        case PhaseEquilibriumStatus::INVALID_TEMPERATURE:
            return "INVALID_TEMPERATURE";
        case PhaseEquilibriumStatus::INVALID_COMPOSITION_SIZE:
            return "INVALID_COMPOSITION_SIZE";
        case PhaseEquilibriumStatus::NONFINITE_COMPOSITION:
            return "NONFINITE_COMPOSITION";
        case PhaseEquilibriumStatus::NEGATIVE_COMPOSITION:
            return "NEGATIVE_COMPOSITION";
        case PhaseEquilibriumStatus::COMPOSITION_NOT_NORMALIZED:
            return "COMPOSITION_NOT_NORMALIZED";
        case PhaseEquilibriumStatus::NONFINITE_PARAMETER:
            return "NONFINITE_PARAMETER";
        case PhaseEquilibriumStatus::NO_PHYSICAL_ROOT: return "NO_PHYSICAL_ROOT";
        case PhaseEquilibriumStatus::INVALID_LOG_ARGUMENT:
            return "INVALID_LOG_ARGUMENT";
        case PhaseEquilibriumStatus::WILSON_K_FAILURE: return "WILSON_K_FAILURE";
        case PhaseEquilibriumStatus::RACHFORD_RICE_NO_BRACKET:
            return "RACHFORD_RICE_NO_BRACKET";
        case PhaseEquilibriumStatus::RACHFORD_RICE_NOT_CONVERGED:
            return "RACHFORD_RICE_NOT_CONVERGED";
        case PhaseEquilibriumStatus::STABILITY_NOT_CONVERGED:
            return "STABILITY_NOT_CONVERGED";
        case PhaseEquilibriumStatus::FLASH_NOT_CONVERGED:
            return "FLASH_NOT_CONVERGED";
        case PhaseEquilibriumStatus::PURE_SATURATION_NO_BRACKET:
            return "PURE_SATURATION_NO_BRACKET";
        case PhaseEquilibriumStatus::PURE_SATURATION_NOT_CONVERGED:
            return "PURE_SATURATION_NOT_CONVERGED";
    }
    return "UNKNOWN_STATUS";
}

PRMixtureKernel::PRMixtureKernel(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    double compositionTolerance)
:
    speciesNames_(speciesNames),
    configurationStatus_(PhaseEquilibriumStatus::SUCCESS),
    compositionTolerance_(compositionTolerance)
{
    if (speciesNames_.empty() || !finitePositive(compositionTolerance_))
    {
        configurationStatus_ = PhaseEquilibriumStatus::INVALID_CONFIGURATION;
        return;
    }

    std::set<std::string> uniqueNames;
    properties_.reserve(speciesNames_.size());
    for (const std::string& name : speciesNames_)
    {
        if (name.empty() || !uniqueNames.insert(name).second)
        {
            configurationStatus_ = PhaseEquilibriumStatus::INVALID_CONFIGURATION;
            properties_.clear();
            return;
        }
        if (!criticalData.has(name))
        {
            configurationStatus_ = PhaseEquilibriumStatus::UNKNOWN_SPECIES;
            properties_.clear();
            return;
        }

        const CriticalProperties properties = criticalData.get(name);
        if (!finitePositive(properties.Tc) || !finitePositive(properties.Pc)
            || !finitePositive(properties.MW) || !std::isfinite(properties.omega))
        {
            configurationStatus_ =
                PhaseEquilibriumStatus::INVALID_CRITICAL_PROPERTIES;
            properties_.clear();
            return;
        }
        properties_.push_back(properties);
    }

    const std::size_t count = speciesNames_.size();
    binaryInteraction_.assign(
        count*count, std::numeric_limits<double>::quiet_NaN());
    hasBinaryInteraction_.assign(count*count, false);
    for (std::size_t i = 0; i < count; ++i)
    {
        binaryInteraction_[i*count + i] = 0.0;
        hasBinaryInteraction_[i*count + i] = true;
        for (std::size_t j = i + 1; j < count; ++j)
        {
            const bool available = kij.has(
                speciesNames_[i], speciesNames_[j]);
            hasBinaryInteraction_[i*count + j] = available;
            hasBinaryInteraction_[j*count + i] = available;
            if (!available) continue;

            const double value = kij.get(
                speciesNames_[i], speciesNames_[j]);
            if (!std::isfinite(value))
            {
                configurationStatus_ =
                    PhaseEquilibriumStatus::NONFINITE_PARAMETER;
                return;
            }
            binaryInteraction_[i*count + j] = value;
            binaryInteraction_[j*count + i] = value;
        }
    }

}

PhaseEquilibriumStatus PRMixtureKernel::validateMoleFractions(
    const std::vector<double>& moleFractions) const
{
    if (!configured()) return configurationStatus_;
    if (moleFractions.size() != properties_.size())
        return PhaseEquilibriumStatus::INVALID_COMPOSITION_SIZE;

    double sum = 0.0;
    for (double fraction : moleFractions)
    {
        if (!std::isfinite(fraction))
            return PhaseEquilibriumStatus::NONFINITE_COMPOSITION;
        if (fraction < 0.0)
            return PhaseEquilibriumStatus::NEGATIVE_COMPOSITION;
        sum += fraction;
    }
    if (!std::isfinite(sum))
        return PhaseEquilibriumStatus::NONFINITE_COMPOSITION;
    if (std::abs(sum - 1.0) > compositionTolerance_)
        return PhaseEquilibriumStatus::COMPOSITION_NOT_NORMALIZED;
    return PhaseEquilibriumStatus::SUCCESS;
}

PhaseEquilibriumStatus PRMixtureKernel::massToMoleFractions(
    const std::vector<double>& massFractions,
    std::vector<double>& moleFractions) const
{
    moleFractions.clear();
    const PhaseEquilibriumStatus inputStatus =
        validateMoleFractions(massFractions);
    if (inputStatus != PhaseEquilibriumStatus::SUCCESS) return inputStatus;

    double sumMoles = 0.0;
    std::vector<double> unnormalised(properties_.size(), 0.0);
    for (std::size_t i = 0; i < properties_.size(); ++i)
    {
        unnormalised[i] = massFractions[i]/properties_[i].MW;
        sumMoles += unnormalised[i];
    }
    if (!finitePositive(sumMoles))
        return PhaseEquilibriumStatus::NONFINITE_PARAMETER;

    moleFractions.resize(properties_.size());
    for (std::size_t i = 0; i < properties_.size(); ++i)
        moleFractions[i] = unnormalised[i]/sumMoles;
    return PhaseEquilibriumStatus::SUCCESS;
}

PhaseEquilibriumStatus PRMixtureKernel::moleToMassFractions(
    const std::vector<double>& moleFractions,
    std::vector<double>& massFractions) const
{
    massFractions.clear();
    const PhaseEquilibriumStatus inputStatus =
        validateMoleFractions(moleFractions);
    if (inputStatus != PhaseEquilibriumStatus::SUCCESS) return inputStatus;

    double mixtureMW = 0.0;
    for (std::size_t i = 0; i < properties_.size(); ++i)
        mixtureMW += moleFractions[i]*properties_[i].MW;
    if (!finitePositive(mixtureMW))
        return PhaseEquilibriumStatus::NONFINITE_PARAMETER;

    massFractions.resize(properties_.size());
    for (std::size_t i = 0; i < properties_.size(); ++i)
        massFractions[i] = moleFractions[i]*properties_[i].MW/mixtureMW;
    return PhaseEquilibriumStatus::SUCCESS;
}

double PRMixtureKernel::mixtureMolecularWeight(
    const std::vector<double>& moleFractions) const
{
    if (validateMoleFractions(moleFractions)
        != PhaseEquilibriumStatus::SUCCESS)
        return quietNaN();

    double mixtureMW = 0.0;
    for (std::size_t i = 0; i < properties_.size(); ++i)
        mixtureMW += moleFractions[i]*properties_[i].MW;
    return finitePositive(mixtureMW) ? mixtureMW : quietNaN();
}

PhaseEquilibriumStatus PRMixtureKernel::wilsonKValues(
    double pressure,
    double temperature,
    std::vector<double>& K) const
{
    K.clear();
    if (!configured()) return configurationStatus_;
    if (!finitePositive(pressure))
        return PhaseEquilibriumStatus::INVALID_PRESSURE;
    if (!finitePositive(temperature))
        return PhaseEquilibriumStatus::INVALID_TEMPERATURE;

    const double maximumLog = std::log(std::numeric_limits<double>::max());
    const double minimumLog = std::log(std::numeric_limits<double>::min());
    K.resize(properties_.size());
    for (std::size_t i = 0; i < properties_.size(); ++i)
    {
        const CriticalProperties& properties = properties_[i];
        const double logK = std::log(properties.Pc/pressure)
            + 5.373*(1.0 + properties.omega)
                *(1.0 - properties.Tc/temperature);
        if (!std::isfinite(logK) || logK >= maximumLog || logK <= minimumLog)
        {
            K.clear();
            return PhaseEquilibriumStatus::WILSON_K_FAILURE;
        }
        K[i] = std::exp(logK);
        if (!finitePositive(K[i]))
        {
            K.clear();
            return PhaseEquilibriumStatus::WILSON_K_FAILURE;
        }
    }
    return PhaseEquilibriumStatus::SUCCESS;
}

PRMixtureParameters PRMixtureKernel::mixtureParameters(
    double pressure,
    double temperature,
    const std::vector<double>& moleFractions) const
{
    PRMixtureParameters result;
    if (!configured())
    {
        result.status = configurationStatus_;
        return result;
    }
    if (!finitePositive(pressure))
    {
        result.status = PhaseEquilibriumStatus::INVALID_PRESSURE;
        return result;
    }
    if (!finitePositive(temperature))
    {
        result.status = PhaseEquilibriumStatus::INVALID_TEMPERATURE;
        return result;
    }
    result.status = validateMoleFractions(moleFractions);
    if (result.status != PhaseEquilibriumStatus::SUCCESS) return result;

    const std::size_t count = properties_.size();
    result.moleFractions = moleFractions;
    result.speciesA.resize(count);
    result.speciesB.resize(count);
    result.pairA.resize(count*count);
    result.sumPairA.assign(count, 0.0);

    for (std::size_t i = 0; i < count; ++i)
    {
        const PRAlphaDerivatives alpha = evaluatePRAlpha(
            properties_[i].Tc,
            properties_[i].Pc,
            properties_[i].omega,
            temperature);
        result.speciesA[i] = alpha.a;
        result.speciesB[i] = 0.07780*universalGasConstant
            *properties_[i].Tc/properties_[i].Pc;
        if (!finitePositive(result.speciesA[i])
            || !finitePositive(result.speciesB[i]))
        {
            result.status = PhaseEquilibriumStatus::NONFINITE_PARAMETER;
            return result;
        }
    }

    result.a = 0.0;
    result.b = 0.0;
    for (std::size_t i = 0; i < count; ++i)
    {
        result.b += moleFractions[i]*result.speciesB[i];
        for (std::size_t j = 0; j < count; ++j)
        {
            const bool unlikeActivePair = i != j
                && moleFractions[i] > 0.0 && moleFractions[j] > 0.0;
            if (unlikeActivePair
                && !hasBinaryInteraction_[i*count + j])
            {
                result.status =
                    PhaseEquilibriumStatus::MISSING_BINARY_INTERACTION;
                return result;
            }
            // If either component is exactly absent, this pair contributes
            // exactly zero to a_mix and every active-component sumPairA.
            // Do not query BinaryInteractionDB::get(), whose legacy fallback
            // would otherwise turn a missing pair into an unaudited zero kij.
            if (i != j && !unlikeActivePair)
            {
                result.pairA[i*count + j] = 0.0;
                continue;
            }
            const double oneMinusKij = i == j
                ? 1.0
                : 1.0 - binaryInteraction_[i*count + j];
            const double pair = std::sqrt(
                result.speciesA[i]*result.speciesA[j])*oneMinusKij;
            if (!std::isfinite(pair))
            {
                result.status = PhaseEquilibriumStatus::NONFINITE_PARAMETER;
                return result;
            }
            result.pairA[i*count + j] = pair;
            result.sumPairA[i] += moleFractions[j]*pair;
            result.a += moleFractions[i]*moleFractions[j]*pair;
        }
    }

    if (!finitePositive(result.a) || !finitePositive(result.b))
    {
        result.status = PhaseEquilibriumStatus::NONFINITE_PARAMETER;
        return result;
    }

    result.A = result.a*pressure
        /(universalGasConstant*universalGasConstant
            *temperature*temperature);
    result.B = result.b*pressure/(universalGasConstant*temperature);
    if (!finitePositive(result.A) || !finitePositive(result.B))
    {
        result.status = PhaseEquilibriumStatus::NONFINITE_PARAMETER;
        return result;
    }
    result.status = PhaseEquilibriumStatus::SUCCESS;
    return result;
}

PRRootSet PRMixtureKernel::roots(
    double pressure,
    double temperature,
    const std::vector<double>& moleFractions) const
{
    PRRootSet result;
    result.mixture = mixtureParameters(pressure, temperature, moleFractions);
    result.status = result.mixture.status;
    if (!result.mixture.valid()) return result;

    const double A = result.mixture.A;
    const double B = result.mixture.B;
    const double cubicA2 = B - 1.0;
    const double cubicA1 = A - 2.0*B - 3.0*B*B;
    const double cubicA0 = -A*B + B*B + B*B*B;
    double root1 = quietNaN();
    double root2 = quietNaN();
    double root3 = quietNaN();
    const int rootCount = CubicRoots::solve(
        cubicA2, cubicA1, cubicA0, root1, root2, root3);
    const double candidates[3] = {root1, root2, root3};
    const double rootMargin = 64.0*std::numeric_limits<double>::epsilon()
        *std::max(1.0, std::abs(B));

    for (int i = 0; i < rootCount && i < 3; ++i)
    {
        const double root = candidates[i];
        if (std::isfinite(root) && root > B + rootMargin)
            result.physicalRoots.push_back(root);
    }
    std::sort(result.physicalRoots.begin(), result.physicalRoots.end());
    result.physicalRoots.erase(
        std::unique(
            result.physicalRoots.begin(), result.physicalRoots.end(), sameRoot),
        result.physicalRoots.end());

    if (result.physicalRoots.empty())
    {
        result.status = PhaseEquilibriumStatus::NO_PHYSICAL_ROOT;
        return result;
    }
    result.status = PhaseEquilibriumStatus::SUCCESS;
    return result;
}

PRPhaseState PRMixtureKernel::evaluateAtRoot(
    const PRMixtureParameters& mixture,
    double compressibility,
    PhaseRoot requestedRoot,
    std::size_t rootCount) const
{
    PRPhaseState result;
    result.requestedRoot = requestedRoot;
    result.moleFractions = mixture.moleFractions;
    result.compressibility = compressibility;
    result.physicalRootCount = rootCount;
    result.status = mixture.status;
    if (!mixture.valid()) return result;

    const double sqrtTwo = std::sqrt(2.0);
    const double zMinusB = compressibility - mixture.B;
    const double logarithmNumerator =
        compressibility + (1.0 + sqrtTwo)*mixture.B;
    const double logarithmDenominator =
        compressibility + (1.0 - sqrtTwo)*mixture.B;
    if (!finitePositive(zMinusB) || !finitePositive(logarithmNumerator)
        || !finitePositive(logarithmDenominator))
    {
        result.status = PhaseEquilibriumStatus::INVALID_LOG_ARGUMENT;
        return result;
    }

    const double logarithm =
        std::log(logarithmNumerator/logarithmDenominator);
    const double attractionCoefficient =
        mixture.A/(2.0*sqrtTwo*mixture.B);
    if (!std::isfinite(logarithm) || !std::isfinite(attractionCoefficient))
    {
        result.status = PhaseEquilibriumStatus::INVALID_LOG_ARGUMENT;
        return result;
    }

    result.lnFugacityCoefficients.resize(properties_.size());
    result.residualGibbsOverRT = 0.0;
    for (std::size_t i = 0; i < properties_.size(); ++i)
    {
        const double bRatio = mixture.speciesB[i]/mixture.b;
        const double attractionRatio =
            2.0*mixture.sumPairA[i]/mixture.a - bRatio;
        const double lnPhi = bRatio*(compressibility - 1.0)
            - std::log(zMinusB)
            - attractionCoefficient*attractionRatio*logarithm;
        if (!std::isfinite(lnPhi))
        {
            result.status = PhaseEquilibriumStatus::NONFINITE_PARAMETER;
            result.lnFugacityCoefficients.clear();
            return result;
        }
        result.lnFugacityCoefficients[i] = lnPhi;
        result.residualGibbsOverRT += mixture.moleFractions[i]*lnPhi;
    }
    if (!std::isfinite(result.residualGibbsOverRT))
    {
        result.status = PhaseEquilibriumStatus::NONFINITE_PARAMETER;
        return result;
    }
    result.status = PhaseEquilibriumStatus::SUCCESS;
    return result;
}

PRPhaseState PRMixtureKernel::phaseState(
    double pressure,
    double temperature,
    const std::vector<double>& moleFractions,
    PhaseRoot root) const
{
    const PRRootSet rootSet = roots(pressure, temperature, moleFractions);
    if (!rootSet.valid())
    {
        PRPhaseState result;
        result.status = rootSet.status;
        result.requestedRoot = root;
        result.moleFractions = moleFractions;
        return result;
    }

    if (root == PhaseRoot::LIQUID)
        return evaluateAtRoot(
            rootSet.mixture, rootSet.physicalRoots.front(), root,
            rootSet.physicalRoots.size());
    if (root == PhaseRoot::VAPOR)
        return evaluateAtRoot(
            rootSet.mixture, rootSet.physicalRoots.back(), root,
            rootSet.physicalRoots.size());

    PRPhaseState stable;
    stable.status = PhaseEquilibriumStatus::NO_PHYSICAL_ROOT;
    stable.requestedRoot = PhaseRoot::STABLE;
    double minimumGibbs = std::numeric_limits<double>::infinity();
    for (double physicalRoot : rootSet.physicalRoots)
    {
        PRPhaseState candidate = evaluateAtRoot(
            rootSet.mixture, physicalRoot, PhaseRoot::STABLE,
            rootSet.physicalRoots.size());
        if (candidate.valid() && candidate.residualGibbsOverRT < minimumGibbs)
        {
            minimumGibbs = candidate.residualGibbsOverRT;
            stable = candidate;
        }
    }
    return stable;
}

} // End namespace phaseEquilibrium
} // End namespace Foam
