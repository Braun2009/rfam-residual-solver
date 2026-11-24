#include "PRPhaseRootSignature.H"

#include <algorithm>
#include <cmath>

namespace Foam
{
namespace phaseEquilibrium
{

namespace
{

bool finitePositive(double value)
{
    return std::isfinite(value) && value > 0.0;
}

double rootDistance(double lhs, double rhs)
{
    return std::abs(lhs - rhs)
        /std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

} // End anonymous namespace

const char* phaseRootSignatureStatusName(
    PRPhaseRootSignatureStatus status)
{
    switch (status)
    {
        case PRPhaseRootSignatureStatus::SUCCESS: return "SUCCESS";
        case PRPhaseRootSignatureStatus::INVALID_CONFIGURATION:
            return "INVALID_CONFIGURATION";
        case PRPhaseRootSignatureStatus::INVALID_RECOVERY_RESULT:
            return "INVALID_RECOVERY_RESULT";
        case PRPhaseRootSignatureStatus::ROOT_SET_FAILURE:
            return "ROOT_SET_FAILURE";
        case PRPhaseRootSignatureStatus::RETAINED_ROOT_NOT_FOUND:
            return "RETAINED_ROOT_NOT_FOUND";
        case PRPhaseRootSignatureStatus::RETAINED_ROOT_AMBIGUOUS:
            return "RETAINED_ROOT_AMBIGUOUS";
        case PRPhaseRootSignatureStatus::NONFINITE_MARGIN:
            return "NONFINITE_MARGIN";
    }
    return "UNKNOWN_STATUS";
}

bool PRPhaseRootSignature::valid() const
{
    return status == PRPhaseRootSignatureStatus::SUCCESS
        && physicalRootCount > 0
        && retainedOrdinal < physicalRootCount
        && finitePositive(retainedCompressibility)
        && ((physicalRootCount == 1 && !hasAdjacentRoot
                && std::isinf(adjacentRootGap) && adjacentRootGap > 0.0)
            || (physicalRootCount > 1 && hasAdjacentRoot
                && finitePositive(adjacentRootGap)))
        && finitePositive(dimensionlessDpDrhoMargin);
}

PRPhaseRootSignatureEvaluator::PRPhaseRootSignatureEvaluator(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    double rootRelativeTolerance)
:
    kernel_(criticalData, kij, speciesNames),
    rootRelativeTolerance_(rootRelativeTolerance),
    configured_(
        kernel_.configured()
        && finitePositive(rootRelativeTolerance_)
        && rootRelativeTolerance_ <= 1.0e-2)
{}

PRPhaseRootSignature PRPhaseRootSignatureEvaluator::evaluate(
    const PRPhaseRhoEResult& recovered) const
{
    PRPhaseRootSignature result;
    if (!configured_) return result;
    if (!recovered.success() || !recovered.state || !recovered.derivatives)
    {
        result.status =
            PRPhaseRootSignatureStatus::INVALID_RECOVERY_RESULT;
        return result;
    }

    const PRPhaseCaloricProperties& state = *recovered.state;
    const PRRootSet roots = kernel_.roots(
        state.pressure, state.temperature, state.moleFractions);
    if (!roots.valid() || roots.physicalRoots.empty())
    {
        result.status = PRPhaseRootSignatureStatus::ROOT_SET_FAILURE;
        return result;
    }

    std::vector<std::size_t> matches;
    for (std::size_t index = 0; index < roots.physicalRoots.size(); ++index)
    {
        if (rootDistance(
                state.compressibility, roots.physicalRoots[index])
            <= rootRelativeTolerance_)
        {
            matches.push_back(index);
        }
    }
    if (matches.empty())
    {
        result.status =
            PRPhaseRootSignatureStatus::RETAINED_ROOT_NOT_FOUND;
        return result;
    }
    if (matches.size() != 1U)
    {
        result.status =
            PRPhaseRootSignatureStatus::RETAINED_ROOT_AMBIGUOUS;
        return result;
    }

    result.physicalRootCount = roots.physicalRoots.size();
    result.retainedOrdinal = matches.front();
    result.retainedCompressibility = state.compressibility;
    if (result.physicalRootCount > 1U)
    {
        result.hasAdjacentRoot = true;
        result.adjacentRootGap = std::numeric_limits<double>::infinity();
        if (result.retainedOrdinal > 0U)
        {
            result.adjacentRootGap = std::min(
                result.adjacentRootGap,
                state.compressibility
                    - roots.physicalRoots[result.retainedOrdinal - 1U]);
        }
        if (result.retainedOrdinal + 1U < result.physicalRootCount)
        {
            result.adjacentRootGap = std::min(
                result.adjacentRootGap,
                roots.physicalRoots[result.retainedOrdinal + 1U]
                    - state.compressibility);
        }
    }
    result.dimensionlessDpDrhoMargin =
        state.rho/state.pressure*recovered.derivatives->dpdrho_T;
    if (!finitePositive(result.retainedCompressibility)
        || (result.hasAdjacentRoot
            && !finitePositive(result.adjacentRootGap))
        || !finitePositive(result.dimensionlessDpDrhoMargin))
    {
        result.status = PRPhaseRootSignatureStatus::NONFINITE_MARGIN;
        return result;
    }

    result.status = PRPhaseRootSignatureStatus::SUCCESS;
    return result;
}

} // End namespace phaseEquilibrium
} // End namespace Foam
