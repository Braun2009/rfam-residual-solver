#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

#include "EquationOfState/PengRobinsonEOS.H"
#include "MixingRules/BinaryInteractionDB.H"
#include "MixingRules/SpeciesCriticalDB.H"
#include "PhaseEquilibrium/PRMixtureKernel.H"
#include "PhaseEquilibrium/PRPhaseCaloric.H"
#include "PhaseEquilibrium/PRPhaseDisplacementAuthority.H"
#include "PhaseEquilibrium/PRPhaseRootSignature.H"
#include "PhaseEquilibrium/PRTPFlash.H"

namespace
{

using Foam::BinaryInteractionDB;
using Foam::PengRobinsonEOS;
using Foam::SpeciesCriticalDB;
using Foam::phaseEquilibrium::EquilibriumRegion;
using Foam::phaseEquilibrium::PhaseEquilibriumStatus;
using Foam::phaseEquilibrium::PhaseRoot;
using Foam::phaseEquilibrium::PRMixtureKernel;
using Foam::phaseEquilibrium::PRPhaseCaloric;
using Foam::phaseEquilibrium::PRPhaseCaloricProperties;
using Foam::phaseEquilibrium::PRPhaseDisplacementAuthority;
using Foam::phaseEquilibrium::PRPhaseDisplacementBounds;
using Foam::phaseEquilibrium::PRPhaseDisplacementDiagnostics;
using Foam::phaseEquilibrium::PRPhaseDisplacementResult;
using Foam::phaseEquilibrium::PRPhaseDisplacementStatus;
using Foam::phaseEquilibrium::PRPhaseIdentity;
using Foam::phaseEquilibrium::PRPhaseRhoEOptions;
using Foam::phaseEquilibrium::PRPhaseRhoERecovery;
using Foam::phaseEquilibrium::PRPhaseRhoERequest;
using Foam::phaseEquilibrium::PRPhaseRhoEResult;
using Foam::phaseEquilibrium::PRPhaseRootIdentity;
using Foam::phaseEquilibrium::PRPhaseRootSignature;
using Foam::phaseEquilibrium::PRPhaseRootSignatureEvaluator;
using Foam::phaseEquilibrium::PRPhaseRootSignatureStatus;
using Foam::phaseEquilibrium::PRTPFlash;

static_assert(
    !std::is_default_constructible<PRPhaseDisplacementBounds>::value,
    "displacement bounds must never be implicit");

int checks = 0;
int failures = 0;

void check(bool condition, const std::string& label)
{
    ++checks;
    std::cout << (condition ? "  PASS " : "  FAIL ") << label << '\n';
    if (!condition) ++failures;
}

bool close(
    double lhs,
    double rhs,
    double relativeTolerance,
    double absoluteTolerance = 0.0)
{
    return std::isfinite(lhs) && std::isfinite(rhs)
        && std::abs(lhs - rhs) <= absoluteTolerance
            + relativeTolerance
                *std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

bool sameVector(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs)
{
    return lhs == rhs;
}

bool samePublishedResult(
    const PRPhaseRhoEResult& lhs,
    const PRPhaseRhoEResult& rhs)
{
    if (lhs.status != rhs.status
        || lhs.phaseIdentity != rhs.phaseIdentity
        || lhs.rootIdentity != rhs.rootIdentity
        || lhs.state.has_value() != rhs.state.has_value()
        || lhs.derivatives.has_value() != rhs.derivatives.has_value())
    {
        return false;
    }
    if (lhs.state)
    {
        const auto& a = *lhs.state;
        const auto& b = *rhs.state;
        if (a.root != b.root
            || !sameVector(a.moleFractions, b.moleFractions)
            || !sameVector(a.massFractions, b.massFractions)
            || a.pressure != b.pressure
            || a.temperature != b.temperature
            || a.molecularWeight != b.molecularWeight
            || a.compressibility != b.compressibility
            || a.rho != b.rho || a.h != b.h || a.e != b.e || a.s != b.s)
        {
            return false;
        }
    }
    if (lhs.derivatives)
    {
        const auto& a = *lhs.derivatives;
        const auto& b = *rhs.derivatives;
        if (a.cv != b.cv || a.dedT_rho != b.dedT_rho
            || a.dpdT_rho != b.dpdT_rho
            || a.dpdrho_T != b.dpdrho_T
            || a.dedrho_T != b.dedrho_T)
        {
            return false;
        }
    }
    return lhs.diagnostics.iterations == rhs.diagnostics.iterations
        && lhs.diagnostics.validScanSamples
            == rhs.diagnostics.validScanSamples
        && lhs.diagnostics.invalidScanSamples
            == rhs.diagnostics.invalidScanSamples
        && lhs.diagnostics.densityRelativeResidual
            == rhs.diagnostics.densityRelativeResidual
        && lhs.diagnostics.energyRelativeResidual
            == rhs.diagnostics.energyRelativeResidual
        && lhs.diagnostics.eosRootRelativeResidual
            == rhs.diagnostics.eosRootRelativeResidual
        && lhs.diagnostics.departureCaloricRelativeResidual
            == rhs.diagnostics.departureCaloricRelativeResidual
        && lhs.diagnostics.cvIdentityRelativeResidual
            == rhs.diagnostics.cvIdentityRelativeResidual
        && lhs.diagnostics.energyDensityIdentityRelativeResidual
            == rhs.diagnostics.energyDensityIdentityRelativeResidual;
}

PRPhaseRhoEOptions recoveryOptions(double minimumTemperature = 250.0)
{
    PRPhaseRhoEOptions options;
    options.minimumTemperature = minimumTemperature;
    options.maximumTemperature = 500.0;
    options.temperatureScanIntervals = 128;
    options.maximumIterations = 80;
    options.derivativeRelativeStep = 2.0e-7;
    options.energyAbsoluteTolerance = 1.0e-5;
    options.energyRelativeTolerance = 2.0e-10;
    options.densityRelativeTolerance = 2.0e-9;
    options.rootRelativeTolerance = 2.0e-12;
    options.derivativeIdentityRelativeTolerance = 2.0e-6;
    options.derivativeIdentityAbsoluteTolerance = 1.0e-4;
    return options;
}

PRPhaseDisplacementBounds normalBounds()
{
    return PRPhaseDisplacementBounds(
        0.02,
        0.05,
        1.0e-3,
        3,
        1.0e-5,
        1.0e-4,
        0.05,
        0.05,
        0.05,
        2.0e-9,
        1.0e-5,
        2.0e-10,
        2.0e-12);
}

PRPhaseDisplacementBounds strictStepBounds()
{
    return PRPhaseDisplacementBounds(
        0.02,
        0.05,
        1.0e-3,
        3,
        1.0e-5,
        1.0e-4,
        1.0e-14,
        1.0e-14,
        1.0e-14,
        2.0e-9,
        1.0e-5,
        2.0e-10,
        2.0e-12);
}

std::optional<PRPhaseCaloricProperties> fixedDensityVaporState(
    const PRMixtureKernel& kernel,
    const PRPhaseCaloric& caloric,
    const PengRobinsonEOS& eos,
    double density,
    double temperature,
    const std::vector<double>& massFractions)
{
    std::vector<double> moleFractions;
    if (kernel.massToMoleFractions(massFractions, moleFractions)
        != PhaseEquilibriumStatus::SUCCESS)
    {
        return std::nullopt;
    }
    const double pressure = eos.pressure(
        density, temperature, massFractions);
    if (!(std::isfinite(pressure) && pressure > 0.0)) return std::nullopt;
    const auto roots = kernel.roots(pressure, temperature, moleFractions);
    const auto evaluated = caloric.evaluate(
        pressure, temperature, moleFractions, PhaseRoot::VAPOR);
    if (!roots.valid() || roots.physicalRoots.empty()
        || !evaluated.success()
        || !close(evaluated.state->rho, density, 2.0e-9)
        || !close(
            evaluated.state->compressibility,
            roots.physicalRoots.back(), 2.0e-12))
    {
        return std::nullopt;
    }
    return *evaluated.state;
}

bool sameSuccessfulDisplacement(
    const PRPhaseDisplacementResult& lhs,
    const PRPhaseDisplacementResult& rhs)
{
    return lhs.success() && rhs.success()
        && samePublishedResult(*lhs.requested, *rhs.requested)
        && lhs.rootSignature->physicalRootCount
            == rhs.rootSignature->physicalRootCount
        && lhs.rootSignature->retainedOrdinal
            == rhs.rootSignature->retainedOrdinal
        && lhs.rootSignature->retainedCompressibility
            == rhs.rootSignature->retainedCompressibility
        && lhs.rootSignature->hasAdjacentRoot
            == rhs.rootSignature->hasAdjacentRoot
        && lhs.rootSignature->adjacentRootGap
            == rhs.rootSignature->adjacentRootGap
        && lhs.rootSignature->dimensionlessDpDrhoMargin
            == rhs.rootSignature->dimensionlessDpDrhoMargin
        && lhs.diagnostics.attemptedWaypoints
            == rhs.diagnostics.attemptedWaypoints
        && lhs.diagnostics.acceptedWaypoints
            == rhs.diagnostics.acceptedWaypoints
        && lhs.diagnostics.maximumObservedCompressibilityRelativeStep
            == rhs.diagnostics.maximumObservedCompressibilityRelativeStep
        && lhs.diagnostics.maximumObservedPressureRelativeStep
            == rhs.diagnostics.maximumObservedPressureRelativeStep
        && lhs.diagnostics.maximumObservedTemperatureRelativeStep
            == rhs.diagnostics.maximumObservedTemperatureRelativeStep;
}

} // End anonymous namespace

int main()
{
    std::cout << "=== Offline PR phase displacement authority ===\n";

    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB kij = BinaryInteractionDB::createDefault();

    const std::vector<std::string> propaneSpecies = {"C3H8"};
    PRPhaseCaloric propaneCaloric(critical, kij, propaneSpecies);
    PRPhaseRhoERecovery propaneRecovery(
        critical, kij, propaneSpecies, recoveryOptions(240.0));
    const auto propaneTarget = propaneCaloric.evaluate(
        1.0e6, 300.0, {1.0}, PhaseRoot::VAPOR);
    PRPhaseRhoEResult propaneRecovered;
    if (propaneTarget.success())
    {
        propaneRecovered = propaneRecovery.recover(
            propaneTarget.state->rho,
            propaneTarget.state->e,
            propaneTarget.state->massFractions,
            PRPhaseIdentity::VAPOR,
            PRPhaseRootIdentity::PR_VAPOR_ROOT,
            305.0);
    }
    PRPhaseRootSignatureEvaluator propaneSignatures(
        critical, kij, propaneSpecies, 2.0e-12);
    const PRPhaseRootSignature propaneSignature =
        propaneSignatures.evaluate(propaneRecovered);
    check(
        propaneSignature.valid()
            && propaneSignature.physicalRootCount == 3U
            && propaneSignature.retainedOrdinal == 2U
            && propaneSignature.hasAdjacentRoot
            && propaneSignature.adjacentRootGap > 0.0
            && propaneSignature.dimensionlessDpDrhoMargin > 0.0,
        "three-root propane retains ordinal two with positive gap and margin");

    const std::vector<std::string> oxygenSpecies = {"O2"};
    PRPhaseCaloric oxygenCaloric(critical, kij, oxygenSpecies);
    PRPhaseRhoERecovery oxygenRecovery(
        critical, kij, oxygenSpecies, recoveryOptions(160.0));
    const auto oxygenTarget = oxygenCaloric.evaluate(
        6.0e6, 180.0, {1.0}, PhaseRoot::VAPOR);
    PRPhaseRhoEResult oxygenRecovered;
    if (oxygenTarget.success())
    {
        oxygenRecovered = oxygenRecovery.recover(
            oxygenTarget.state->rho,
            oxygenTarget.state->e,
            oxygenTarget.state->massFractions,
            PRPhaseIdentity::SUPERCRITICAL,
            PRPhaseRootIdentity::PR_UNIQUE_SUPERCRITICAL_ROOT,
            180.0);
    }
    PRPhaseRootSignatureEvaluator oxygenSignatures(
        critical, kij, oxygenSpecies, 2.0e-12);
    const PRPhaseRootSignature oxygenSignature =
        oxygenSignatures.evaluate(oxygenRecovered);
    check(
        oxygenSignature.valid()
            && oxygenSignature.physicalRootCount == 1U
            && oxygenSignature.retainedOrdinal == 0U
            && !oxygenSignature.hasAdjacentRoot
            && std::isinf(oxygenSignature.adjacentRootGap)
            && oxygenSignature.adjacentRootGap > 0.0,
        "unique supercritical root has ordinal zero and explicit no-adjacent gap");

    PRPhaseRootSignatureEvaluator invalidSignatureEvaluator(
        critical, kij, oxygenSpecies, 0.0);
    check(
        !invalidSignatureEvaluator.configured()
            && invalidSignatureEvaluator.evaluate(oxygenRecovered).status
                == PRPhaseRootSignatureStatus::INVALID_CONFIGURATION,
        "root matching tolerance is explicit and bounded");
    check(
        propaneSignatures.evaluate(PRPhaseRhoEResult{}).status
            == PRPhaseRootSignatureStatus::INVALID_RECOVERY_RESULT,
        "signature evaluator rejects partial recovery payloads");

    const std::vector<std::string> species = {"O2", "C3H8"};
    PRTPFlash flash(critical, kij, species);
    const auto parentFlash = flash.flashMoleFractions(
        2.0e6, 300.0, {0.2, 0.8});
    check(
        parentFlash.valid()
            && parentFlash.region == EquilibriumRegion::TWO_PHASE
            && parentFlash.vaporState.valid()
            && !parentFlash.vaporMassFractions.empty(),
        "two-phase parent supplies an explicit vapor anchor");
    if (!parentFlash.valid()
        || parentFlash.region != EquilibriumRegion::TWO_PHASE
        || !parentFlash.vaporState.valid()
        || parentFlash.vaporMassFractions.empty())
    {
        return 1;
    }

    PRMixtureKernel kernel(critical, kij, species);
    PRPhaseCaloric caloric(critical, kij, species);
    PengRobinsonEOS eos(critical, kij);
    eos.setSpeciesNames(species);
    const auto parentState = caloric.evaluate(
        parentFlash.pressure,
        parentFlash.temperature,
        parentFlash.vaporMoleFractions,
        PhaseRoot::VAPOR);
    check(
        parentState.success(),
        "phase caloric authority materializes the vapor parent state");
    if (!parentState.success()) return 1;

    std::vector<double> displacedY = parentFlash.vaporMassFractions;
    displacedY[0] += 1.0e-4;
    displacedY[1] -= 1.0e-4;
    const auto displacedState = fixedDensityVaporState(
        kernel, caloric, eos,
        parentState.state->rho, 301.0, displacedY);
    check(
        displacedState.has_value(),
        "nearby fixed-density rho-e-Y request is independently constructible");
    if (!displacedState) return 1;

    const PRPhaseRhoERequest parentRequest(
        parentState.state->rho,
        parentState.state->e,
        parentFlash.vaporMassFractions,
        parentState.state->temperature);
    const PRPhaseRhoERequest displacedRequest(
        displacedState->rho,
        displacedState->e,
        displacedState->massFractions,
        displacedState->temperature);

    PRPhaseDisplacementAuthority authority(
        critical, kij, species, recoveryOptions(), normalBounds());
    check(authority.configured(),
        "explicit bounded displacement authority configures");
    const PRPhaseDisplacementResult first = authority.evaluate(
        parentRequest, displacedRequest,
        PRPhaseIdentity::VAPOR,
        PRPhaseRootIdentity::PR_VAPOR_ROOT);
    const PRPhaseDisplacementResult second = authority.evaluate(
        parentRequest, displacedRequest,
        PRPhaseIdentity::VAPOR,
        PRPhaseRootIdentity::PR_VAPOR_ROOT);
    check(
        first.success()
            && first.diagnostics.attemptedWaypoints == 9U
            && first.diagnostics.acceptedWaypoints == 9U,
        "depth-three request-anchored homotopy accepts parent plus eight waypoints");
    check(sameSuccessfulDisplacement(first, second),
        "repeated dyadic homotopy is bitwise deterministic");
    check(
        first.requested && first.requested->state
            && close(first.requested->state->rho,
                displacedRequest.density, 2.0e-9)
            && close(first.requested->state->e,
                displacedRequest.specificInternalEnergy,
                2.0e-10, 1.0e-5)
            && first.diagnostics.requestDensityRelativeResidual
                <= normalBounds().densityClosureRelativeTolerance
            && first.diagnostics
                .requestCompositionMaximumAbsoluteResidual
                <= normalBounds().compositionClosureAbsoluteTolerance,
        "publication closes against requested rho-e-Y, not the parent payload");
    check(
        first.rootSignature
            && first.rootSignature->physicalRootCount > 0U
            && first.rootSignature->retainedOrdinal
                == first.rootSignature->physicalRootCount - 1U
            && first.diagnostics
                .minimumObservedDimensionlessDpDrhoMargin
                >= normalBounds().minimumDimensionlessDpDrhoMargin,
        "accepted path retains the vapor ordinal and derivative margin");

    PRPhaseRhoEResult committed = propaneRecovered;
    PRPhaseDisplacementDiagnostics commitDiagnostics;
    const PRPhaseDisplacementStatus committedStatus = authority.tryCommit(
        parentRequest, displacedRequest,
        PRPhaseIdentity::VAPOR,
        PRPhaseRootIdentity::PR_VAPOR_ROOT,
        committed, commitDiagnostics);
    check(
        committedStatus == PRPhaseDisplacementStatus::SUCCESS
            && first.requested
            && samePublishedResult(committed, *first.requested),
        "successful transaction publishes only the complete request state");

    const PRPhaseRhoEResult beforeOutOfBounds = committed;
    const PRPhaseRhoERequest outOfBounds(
        displacedRequest.density,
        displacedRequest.specificInternalEnergy + 1.0e9,
        displacedRequest.massFractions,
        displacedRequest.temperatureSeed);
    const PRPhaseDisplacementStatus outOfBoundsStatus = authority.tryCommit(
        parentRequest, outOfBounds,
        PRPhaseIdentity::VAPOR,
        PRPhaseRootIdentity::PR_VAPOR_ROOT,
        committed, commitDiagnostics);
    check(
        outOfBoundsStatus
                == PRPhaseDisplacementStatus::DISPLACEMENT_OUT_OF_BOUNDS
            && samePublishedResult(committed, beforeOutOfBounds),
        "out-of-envelope request fails before recovery and rolls back publication");

    PRPhaseDisplacementAuthority strictContinuity(
        critical, kij, species, recoveryOptions(), strictStepBounds());
    const PRPhaseRhoEResult beforeContinuityFailure = committed;
    const PRPhaseDisplacementStatus continuityStatus =
        strictContinuity.tryCommit(
            parentRequest, displacedRequest,
            PRPhaseIdentity::VAPOR,
            PRPhaseRootIdentity::PR_VAPOR_ROOT,
            committed, commitDiagnostics);
    check(
        continuityStatus
                == PRPhaseDisplacementStatus::ROOT_CONTINUITY_FAILED
            && commitDiagnostics.failedWaypoint !=
                std::numeric_limits<std::size_t>::max()
            && commitDiagnostics.failedFraction > 0.0
            && samePublishedResult(committed, beforeContinuityFailure),
        "continuity failure reports its dyadic waypoint and rolls back");

    PRPhaseDisplacementBounds invalidBounds(
        0.75,
        0.05,
        1.0e-3,
        3,
        1.0e-5,
        1.0e-4,
        0.05,
        0.05,
        0.05,
        2.0e-9,
        1.0e-5,
        2.0e-10,
        2.0e-12);
    PRPhaseDisplacementAuthority unbounded(
        critical, kij, species, recoveryOptions(), invalidBounds);
    check(
        !unbounded.configured()
            && unbounded.configurationStatus()
                == PRPhaseDisplacementStatus::INVALID_CONFIGURATION,
        "hard ceilings reject an unbounded density displacement envelope");

    std::vector<double> invalidY = displacedRequest.massFractions;
    invalidY[0] = -1.0e-6;
    invalidY[1] = 1.0 + 1.0e-6;
    const PRPhaseRhoERequest invalidRequest(
        displacedRequest.density,
        displacedRequest.specificInternalEnergy,
        invalidY,
        displacedRequest.temperatureSeed);
    const PRPhaseDisplacementResult invalid = authority.evaluate(
        parentRequest, invalidRequest,
        PRPhaseIdentity::VAPOR,
        PRPhaseRootIdentity::PR_VAPOR_ROOT);
    check(
        invalid.status
                == PRPhaseDisplacementStatus::INVALID_DISPLACED_REQUEST
            && !invalid.requested && !invalid.rootSignature,
        "negative requested composition fails without a partial publication");

    const PRPhaseDisplacementResult wrongIdentity = authority.evaluate(
        parentRequest, displacedRequest,
        PRPhaseIdentity::SUPERCRITICAL,
        PRPhaseRootIdentity::PR_UNIQUE_SUPERCRITICAL_ROOT);
    check(
        wrongIdentity.success()
            && wrongIdentity.diagnostics.recoveryStatus
                == Foam::phaseEquilibrium::PRPhaseRhoEStatus::SUCCESS,
        "raw phase and root enums remain caller assertions without parent lineage");

    if (failures == 0)
    {
        std::cout
            << "PR_PHASE_DISPLACEMENT_AUTHORITY_PASS checks=" << checks
            << " dyadic_depth=3 root_signature=1 request_closure=1"
            << " rollback=1 parent_lineage=0 capability=0"
            << " production_ready=0\n";
        return 0;
    }
    std::cout << "PR_PHASE_DISPLACEMENT_AUTHORITY_FAIL checks=" << checks
              << " failures=" << failures << '\n';
    return 1;
}
