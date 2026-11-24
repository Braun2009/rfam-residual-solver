#ifndef PRPhaseAuthorityIdentity_H
#define PRPhaseAuthorityIdentity_H

#include <string>
#include <vector>

#include "BinaryInteractionDB.H"
#include "SpeciesCriticalDB.H"

namespace Foam
{
namespace phaseEquilibrium
{

struct PRPhaseRhoEOptions;
struct PRPhaseDisplacementBounds;

// Physical-model identity shared by the whole-mixture PR owner, the
// phase-local rho-e recovery, and any bounded displacement authority. The
// recovery identity additionally binds every numerical recovery option.
std::string makePRPhasePhysicalModelAuthorityFingerprint(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames);

std::string makePRPhaseRhoEAuthorityFingerprint(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij,
    const std::vector<std::string>& speciesNames,
    const PRPhaseRhoEOptions& options);

std::string makePRPhaseDisplacementAuthorityFingerprint(
    const std::string& recoveryAuthorityFingerprint,
    const PRPhaseDisplacementBounds& bounds);

bool validPRPhaseAuthorityFingerprint(const std::string& fingerprint);

} // End namespace phaseEquilibrium
} // End namespace Foam

#endif
