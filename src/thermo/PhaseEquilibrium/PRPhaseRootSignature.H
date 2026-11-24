#ifndef PRPhaseRootSignature_H
#define PRPhaseRootSignature_H

#include <cstddef>
#include <limits>
#include <vector>

#include "PRMixtureKernel.H"
#include "PRPhaseRhoERecovery.H"

namespace Foam
{
namespace phaseEquilibrium
{

enum class PRPhaseRootSignatureStatus
{
    SUCCESS,
    INVALID_CONFIGURATION,
    INVALID_RECOVERY_RESULT,
    ROOT_SET_FAILURE,
    RETAINED_ROOT_NOT_FOUND,
    RETAINED_ROOT_AMBIGUOUS,
    NONFINITE_MARGIN
};

const char* phaseRootSignatureStatusName(
    PRPhaseRootSignatureStatus status);

struct PRPhaseRootSignature
{
    PRPhaseRootSignatureStatus status =
        PRPhaseRootSignatureStatus::INVALID_CONFIGURATION;
    std::size_t physicalRootCount = 0;
    std::size_t retainedOrdinal = std::numeric_limits<std::size_t>::max();
    double retainedCompressibility =
        std::numeric_limits<double>::quiet_NaN();
    bool hasAdjacentRoot = false;
    double adjacentRootGap = std::numeric_limits<double>::infinity();
    double dimensionlessDpDrhoMargin =
        std::numeric_limits<double>::quiet_NaN();

    bool valid() const;
};

// Computes a branch signature from a successful phase-local rho-e-Y result.
// The retained ordinal is zero-based in the kernel's ascending physical-root
// order. A unique root has no adjacent root and reports an infinite gap.
class PRPhaseRootSignatureEvaluator
{
private:
    PRMixtureKernel kernel_;
    double rootRelativeTolerance_;
    bool configured_;

public:
    PRPhaseRootSignatureEvaluator(
        const SpeciesCriticalDB& criticalData,
        const BinaryInteractionDB& kij,
        const std::vector<std::string>& speciesNames,
        double rootRelativeTolerance);

    bool configured() const
    {
        return configured_;
    }

    PRPhaseRootSignature evaluate(
        const PRPhaseRhoEResult& recovered) const;
};

} // End namespace phaseEquilibrium
} // End namespace Foam

#endif
