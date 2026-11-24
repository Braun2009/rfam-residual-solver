#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "EquationOfState/PengRobinsonEOS.H"
#include "MixingRules/BinaryInteractionDB.H"
#include "MixingRules/SpeciesCriticalDB.H"
#include "PhaseEquilibrium/PRMixtureKernel.H"
#include "PhaseEquilibrium/PRPhaseCaloric.H"
#include "PhaseEquilibrium/PRPhaseRhoERecovery.H"
#include "PhaseEquilibrium/PRTPFlash.H"
#include "PhaseEquilibrium/TabulatedIdealCaloric.H"
#include "SpeciesRegistry.H"

namespace
{

using Foam::BinaryInteractionDB;
using Foam::PengRobinsonEOS;
using Foam::SpeciesCriticalDB;
using Foam::SpeciesRegistry;
using Foam::phaseEquilibrium::EquilibriumRegion;
using Foam::phaseEquilibrium::PhaseEquilibriumStatus;
using Foam::phaseEquilibrium::PhaseRoot;
using Foam::phaseEquilibrium::PRMixtureKernel;
using Foam::phaseEquilibrium::PRPhaseCaloric;
using Foam::phaseEquilibrium::PRPhaseCaloricProperties;
using Foam::phaseEquilibrium::PRPhaseIdentity;
using Foam::phaseEquilibrium::PRPhaseRhoEOptions;
using Foam::phaseEquilibrium::PRPhaseRhoERecovery;
using Foam::phaseEquilibrium::PRPhaseRhoEResult;
using Foam::phaseEquilibrium::PRPhaseRhoEStatus;
using Foam::phaseEquilibrium::PRPhaseRootIdentity;
using Foam::phaseEquilibrium::PRTPFlash;
using Foam::phaseEquilibrium::TPFlashResult;

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

bool closeVector(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs,
    double absoluteTolerance)
{
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        if (!close(lhs[i], rhs[i], 0.0, absoluteTolerance)) return false;
    }
    return true;
}

bool exposesNoState(const PRPhaseRhoEResult& result)
{
    return !result.success() && !result.state.has_value()
        && !result.derivatives.has_value();
}

PRPhaseRhoEOptions vaporOptions()
{
    PRPhaseRhoEOptions options;
    options.minimumTemperature = 250.0;
    options.maximumTemperature = 500.0;
    options.temperatureScanIntervals = 128;
    options.maximumIterations = 80;
    return options;
}

PRPhaseRhoEOptions supercriticalOptions()
{
    PRPhaseRhoEOptions options = vaporOptions();
    options.minimumTemperature = 160.0;
    return options;
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
    if (!roots.valid() || roots.physicalRoots.empty()) return std::nullopt;

    const auto evaluated = caloric.evaluate(
        pressure, temperature, moleFractions, PhaseRoot::VAPOR);
    if (!evaluated.success()
        || !close(evaluated.state->rho, density, 2.0e-9)
        || !close(
            evaluated.state->compressibility,
            roots.physicalRoots.back(), 2.0e-12))
    {
        return std::nullopt;
    }
    return *evaluated.state;
}

void checkFailure(
    const PRPhaseRhoEResult& result,
    PRPhaseRhoEStatus expectedStatus,
    const std::string& label)
{
    check(
        result.status == expectedStatus && exposesNoState(result),
        label);
}

} // End anonymous namespace

int main()
{
    using Foam::phaseEquilibrium::tabulatedIdealCaloric;

    std::cout << "=== Offline PR phase-local rho-e-Y recovery ===\n";

    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB kij = BinaryInteractionDB::createDefault();
    const std::vector<std::string> species = {"O2", "C3H8"};

    PRTPFlash flash(critical, kij, species);
    const TPFlashResult parent = flash.flashMoleFractions(
        2.0e6, 300.0, {0.2, 0.8});
    check(
        parent.valid() && parent.region == EquilibriumRegion::TWO_PHASE
            && parent.vaporState.valid()
            && !parent.vaporMassFractions.empty(),
        "two-phase parent supplies an explicit vapor payload");
    if (!parent.valid() || parent.region != EquilibriumRegion::TWO_PHASE
        || !parent.vaporState.valid()
        || parent.vaporMassFractions.empty())
    {
        return 1;
    }

    PRMixtureKernel kernel(critical, kij, species);
    PRPhaseCaloric caloric(critical, kij, species);
    PengRobinsonEOS eos(critical, kij);
    eos.setSpeciesNames(species);
    PRPhaseRhoERecovery recovery(
        critical, kij, species, vaporOptions());
    check(
        kernel.configured() && caloric.configured()
            && recovery.configured()
            && recovery.configurationStatus() == PRPhaseRhoEStatus::SUCCESS,
        "sourced O2-C3H8 phase-local backend configures");

    const auto parentVapor = caloric.evaluate(
        parent.pressure, parent.temperature,
        parent.vaporMoleFractions, PhaseRoot::VAPOR);
    check(parentVapor.success(),
        "two-phase parent vapor has a phase-caloric target");
    if (!parentVapor.success()) return 1;

    const PRPhaseRhoEResult recovered = recovery.recover(
        parentVapor.state->rho,
        parentVapor.state->e,
        parentVapor.state->massFractions,
        PRPhaseIdentity::VAPOR,
        PRPhaseRootIdentity::PR_VAPOR_ROOT,
        310.0);
    check(recovered.success(),
        "biased seed recovers the committed two-phase vapor rho-e-Y");
    check(
        recovered.state
            && close(recovered.state->pressure, parent.pressure, 2.0e-9)
            && close(recovered.state->temperature,
                parent.temperature, 2.0e-9)
            && close(recovered.state->rho,
                parentVapor.state->rho, 2.0e-9)
            && close(recovered.state->e,
                parentVapor.state->e, 2.0e-9)
            && close(recovered.state->compressibility,
                parent.vaporState.compressibility, 2.0e-9)
            && closeVector(recovered.state->massFractions,
                parent.vaporMassFractions, 2.0e-12),
        "recovery replays parent p-T-rho-e-Z and ordered composition");
    check(
        recovered.phaseIdentity == PRPhaseIdentity::VAPOR
            && recovered.rootIdentity
                == PRPhaseRootIdentity::PR_VAPOR_ROOT,
        "result retains the caller-bound vapor and PR-root identities");

    const PRPhaseRhoEResult alternateSeed = recovery.recover(
        parentVapor.state->rho,
        parentVapor.state->e,
        parentVapor.state->massFractions,
        PRPhaseIdentity::VAPOR,
        PRPhaseRootIdentity::PR_VAPOR_ROOT,
        470.0);
    check(
        alternateSeed.success() && alternateSeed.state && recovered.state
            && close(alternateSeed.state->pressure,
                recovered.state->pressure, 2.0e-9)
            && close(alternateSeed.state->temperature,
                recovered.state->temperature, 2.0e-9)
            && close(alternateSeed.state->compressibility,
                recovered.state->compressibility, 2.0e-9),
        "an alternate in-range seed retains the same vapor root");

    PRPhaseRhoEOptions narrowOptions = vaporOptions();
    narrowOptions.minimumTemperature = 299.999999;
    narrowOptions.maximumTemperature = 300.000001;
    PRPhaseRhoERecovery narrowRecovery(
        critical, kij, species, narrowOptions);
    const PRPhaseRhoEResult narrowRecovered = narrowRecovery.recover(
        parentVapor.state->rho,
        parentVapor.state->e,
        parentVapor.state->massFractions,
        PRPhaseIdentity::VAPOR,
        PRPhaseRootIdentity::PR_VAPOR_ROOT,
        300.0);
    check(
        narrowRecovered.success() && narrowRecovered.state
            && close(narrowRecovered.state->temperature, 300.0, 1.0e-12)
            && narrowRecovered.diagnostics.validScanSamples == 129,
        "one contiguous tolerance band is one energy-root candidate");

    std::vector<double> displacedY = parentVapor.state->massFractions;
    displacedY[0] += 1.0e-7;
    displacedY[1] -= 1.0e-7;
    const auto displacedTarget = fixedDensityVaporState(
        kernel, caloric, eos, parentVapor.state->rho, 302.0, displacedY);
    check(displacedTarget.has_value(),
        "nearby fixed-density composition-energy target exists");
    if (displacedTarget)
    {
        const PRPhaseRhoEResult displaced = recovery.recover(
            parentVapor.state->rho,
            displacedTarget->e,
            displacedY,
            PRPhaseIdentity::VAPOR,
            PRPhaseRootIdentity::PR_VAPOR_ROOT,
            295.0);
        check(displaced.success(),
            "normalized displaced rho-e-Y target recovers");
        check(
            displaced.state
                && close(displaced.state->temperature, 302.0, 2.0e-9)
                && close(displaced.state->pressure,
                    displacedTarget->pressure, 2.0e-9)
                && close(displaced.state->rho,
                    parentVapor.state->rho, 2.0e-9)
                && close(displaced.state->e,
                    displacedTarget->e, 2.0e-9)
                && closeVector(displaced.state->massFractions,
                    displacedY, 2.0e-12),
            "displaced recovery closes p-T-rho-e-Y on the vapor root");
    }

    check(recovered.derivatives && recovered.derivatives->valid(),
        "successful vapor recovery publishes finite phase derivatives");
    if (recovered.state && recovered.derivatives)
    {
        const auto& state = *recovered.state;
        const auto& derivatives = *recovered.derivatives;
        const double expectedDedrho =
            (state.pressure - state.temperature*derivatives.dpdT_rho)
            /(state.rho*state.rho);
        check(derivatives.dpdrho_T > 0.0,
            "vapor dp/drho|T is positive");
        check(derivatives.dedT_rho > 0.0,
            "vapor de/dT|rho is positive");
        check(close(derivatives.dedT_rho, derivatives.cv, 2.0e-6, 1.0e-4),
            "vapor de/dT|rho equals explicit phase Cv");
        check(close(derivatives.dedrho_T, expectedDedrho, 2.0e-6, 1.0e-4),
            "vapor de/drho|T satisfies the pressure-temperature identity");
        check(
            recovered.diagnostics.densityRelativeResidual <= 2.0e-9
                && recovered.diagnostics.eosRootRelativeResidual <= 2.0e-12
                && recovered.diagnostics.cvIdentityRelativeResidual
                    <= 2.0e-6
                && recovered.diagnostics.energyDensityIdentityRelativeResidual
                    <= 2.0e-6,
            "vapor diagnostics expose bounded closure residuals");
    }

    std::vector<double> invalidY = parentVapor.state->massFractions;
    invalidY[0] = -1.0e-8;
    checkFailure(
        recovery.recover(
            parentVapor.state->rho, parentVapor.state->e, invalidY,
            PRPhaseIdentity::VAPOR,
            PRPhaseRootIdentity::PR_VAPOR_ROOT, 300.0),
        PRPhaseRhoEStatus::INVALID_COMPOSITION,
        "negative composition fails without a phase state");
    checkFailure(
        recovery.recover(
            0.0, parentVapor.state->e, parentVapor.state->massFractions,
            PRPhaseIdentity::VAPOR,
            PRPhaseRootIdentity::PR_VAPOR_ROOT, 300.0),
        PRPhaseRhoEStatus::INVALID_DENSITY,
        "nonpositive density fails without a phase state");
    checkFailure(
        recovery.recover(
            parentVapor.state->rho,
            std::numeric_limits<double>::quiet_NaN(),
            parentVapor.state->massFractions,
            PRPhaseIdentity::VAPOR,
            PRPhaseRootIdentity::PR_VAPOR_ROOT, 300.0),
        PRPhaseRhoEStatus::INVALID_ENERGY,
        "nonfinite energy fails without a phase state");
    checkFailure(
        recovery.recover(
            parentVapor.state->rho, parentVapor.state->e,
            parentVapor.state->massFractions,
            PRPhaseIdentity::VAPOR,
            PRPhaseRootIdentity::PR_VAPOR_ROOT, 0.0),
        PRPhaseRhoEStatus::INVALID_TEMPERATURE_SEED,
        "nonpositive temperature seed fails without a phase state");
    checkFailure(
        recovery.recover(
            parentVapor.state->rho, parentVapor.state->e,
            parentVapor.state->massFractions,
            PRPhaseIdentity::VAPOR,
            PRPhaseRootIdentity::PR_UNIQUE_SUPERCRITICAL_ROOT, 300.0),
        PRPhaseRhoEStatus::INVALID_PHASE_IDENTITY,
        "incompatible vapor/root identity fails without a phase state");
    checkFailure(
        recovery.recover(
            parentVapor.state->rho, 1.0e12,
            parentVapor.state->massFractions,
            PRPhaseIdentity::VAPOR,
            PRPhaseRootIdentity::PR_VAPOR_ROOT, 300.0),
        PRPhaseRhoEStatus::NOT_BRACKETED,
        "out-of-envelope energy fails without a phase state");

    const std::vector<std::string> propaneSpecies = {"C3H8"};
    PRMixtureKernel propaneKernel(critical, kij, propaneSpecies);
    PRPhaseCaloric propaneCaloric(critical, kij, propaneSpecies);
    PRPhaseRhoEOptions propaneOptions = vaporOptions();
    propaneOptions.minimumTemperature = 240.0;
    propaneOptions.maximumTemperature = 330.0;
    PRPhaseRhoERecovery propaneRecovery(
        critical, kij, propaneSpecies, propaneOptions);
    check(
        propaneKernel.configured() && propaneCaloric.configured()
            && propaneRecovery.configured(),
        "pure-propane phase-local backend configures");

    const auto propaneNasaTarget = propaneCaloric.evaluate(
        1.0e6, 300.0, {1.0}, PhaseRoot::VAPOR);
    const auto propaneNasaRoots = propaneKernel.roots(
        1.0e6, 300.0, {1.0});
    check(
        propaneNasaTarget.success() && propaneNasaRoots.valid()
            && propaneNasaRoots.physicalRoots.size() == 3U,
        "pure propane supplies a three-root 300 K NASA target");
    if (propaneNasaTarget.success() && propaneNasaRoots.valid()
        && !propaneNasaRoots.physicalRoots.empty())
    {
        const PRPhaseRhoEResult propaneNasa = propaneRecovery.recover(
            propaneNasaTarget.state->rho,
            propaneNasaTarget.state->e,
            propaneNasaTarget.state->massFractions,
            PRPhaseIdentity::VAPOR,
            PRPhaseRootIdentity::PR_VAPOR_ROOT,
            305.0);
        check(
            propaneNasa.success() && propaneNasa.state
                && close(
                    propaneNasa.state->compressibility,
                    propaneNasaRoots.physicalRoots.back(), 2.0e-12),
            "three-root recovery retains the largest PR vapor root");
    }

    const auto propaneTableTarget = propaneCaloric.evaluate(
        1.0e5, 250.0, {1.0}, PhaseRoot::VAPOR);
    PRPhaseRhoEResult propaneTable;
    if (propaneTableTarget.success())
    {
        propaneTable = propaneRecovery.recover(
            propaneTableTarget.state->rho,
            propaneTableTarget.state->e,
            propaneTableTarget.state->massFractions,
            PRPhaseIdentity::VAPOR,
            PRPhaseRootIdentity::PR_VAPOR_ROOT,
            260.0);
    }
    check(
        propaneTableTarget.success() && propaneTable.success()
            && propaneTable.state
            && close(propaneTable.state->temperature, 250.0, 2.0e-9),
        "pure-propane recovery supports the 250 K table branch");

    const auto propaneHandoffWindowTarget = propaneCaloric.evaluate(
        1.0e6, 296.5, {1.0}, PhaseRoot::VAPOR);
    PRPhaseRhoEResult propaneHandoffWindow;
    if (propaneHandoffWindowTarget.success())
    {
        propaneHandoffWindow = propaneRecovery.recover(
            propaneHandoffWindowTarget.state->rho,
            propaneHandoffWindowTarget.state->e,
            propaneHandoffWindowTarget.state->massFractions,
            PRPhaseIdentity::VAPOR,
            PRPhaseRootIdentity::PR_VAPOR_ROOT,
            296.5);
    }
    check(
        propaneHandoffWindowTarget.success()
            && propaneHandoffWindow.success()
            && propaneHandoffWindow.state
            && propaneHandoffWindow.derivatives
            && close(
                propaneHandoffWindow.state->temperature,
                296.5, 2.0e-9),
        "C3H8 recovery uses corrected Cp inside the C1 handoff window");

    const auto propaneBridgeTarget = propaneCaloric.evaluate(
        1.0e6, 299.0, {1.0}, PhaseRoot::VAPOR);
    PRPhaseRhoEResult propaneBridge;
    if (propaneBridgeTarget.success())
    {
        propaneBridge = propaneRecovery.recover(
            propaneBridgeTarget.state->rho,
            propaneBridgeTarget.state->e,
            propaneBridgeTarget.state->massFractions,
            PRPhaseIdentity::VAPOR,
            PRPhaseRootIdentity::PR_VAPOR_ROOT,
            305.0);
    }
    check(
        propaneBridgeTarget.success() && propaneBridge.success()
            && propaneBridge.state
            && close(propaneBridge.state->temperature, 299.0, 2.0e-9),
        "pure-propane recovery supports the 299 K NASA bridge");

    const auto propaneKinkTarget = propaneCaloric.evaluate(
        1.0e6, 298.15, {1.0}, PhaseRoot::VAPOR);
    PRPhaseRhoEResult propaneKink;
    if (propaneKinkTarget.success())
    {
        propaneKink = propaneRecovery.recover(
            propaneKinkTarget.state->rho,
            propaneKinkTarget.state->e,
            propaneKinkTarget.state->massFractions,
            PRPhaseIdentity::VAPOR,
            PRPhaseRootIdentity::PR_VAPOR_ROOT,
            298.15);
    }
    check(
        propaneKinkTarget.success() && propaneKink.success()
            && propaneKink.state && propaneKink.derivatives
            && close(propaneKink.state->temperature, 298.15, 2.0e-9),
        "C3H8 C1 table-to-bridge handoff supports rho-e recovery");

    const std::vector<std::string> oxygenSpecies = {"O2"};
    PRTPFlash oxygenFlash(critical, kij, oxygenSpecies);
    const TPFlashResult oxygenParent = oxygenFlash.flashMoleFractions(
        6.0e6, 180.0, {1.0});
    const TPFlashResult oxygenHandoffParent = oxygenFlash.flashMoleFractions(
        6.0e6, 200.0, {1.0});
    check(
        oxygenParent.valid()
            && oxygenParent.region == EquilibriumRegion::SINGLE_SUPERCRITICAL
            && oxygenParent.singlePhaseState.valid(),
        "parent flash explicitly identifies low-temperature supercritical O2");
    check(
        oxygenHandoffParent.valid()
            && oxygenHandoffParent.region
                == EquilibriumRegion::SINGLE_SUPERCRITICAL
            && oxygenHandoffParent.singlePhaseState.valid(),
        "parent flash identifies the 200 K O2 caloric handoff state");
    if (oxygenParent.valid()
        && oxygenParent.region == EquilibriumRegion::SINGLE_SUPERCRITICAL)
    {
        PRMixtureKernel oxygenKernel(critical, kij, oxygenSpecies);
        PRPhaseCaloric oxygenCaloric(critical, kij, oxygenSpecies);
        PengRobinsonEOS oxygenEos(critical, kij);
        oxygenEos.setSpeciesNames(oxygenSpecies);
        const auto oxygenTarget = oxygenCaloric.evaluate(
            oxygenParent.pressure, oxygenParent.temperature,
            {1.0}, PhaseRoot::VAPOR);
        check(oxygenTarget.success(),
            "supercritical O2 has a low-temperature phase-caloric target");

        PRPhaseRhoERecovery oxygenRecovery(
            critical, kij, oxygenSpecies, supercriticalOptions());
        check(oxygenRecovery.configured(),
            "pure-O2 supercritical recovery configures");
        const auto oxygenHandoffWindowTarget = oxygenCaloric.evaluate(
            6.0e6, 197.5, {1.0}, PhaseRoot::VAPOR);
        PRPhaseRhoEResult oxygenHandoffWindow;
        if (oxygenHandoffWindowTarget.success())
        {
            oxygenHandoffWindow = oxygenRecovery.recover(
                oxygenHandoffWindowTarget.state->rho,
                oxygenHandoffWindowTarget.state->e,
                oxygenHandoffWindowTarget.state->massFractions,
                PRPhaseIdentity::SUPERCRITICAL,
                PRPhaseRootIdentity::PR_UNIQUE_SUPERCRITICAL_ROOT,
                197.5);
        }
        check(
            oxygenHandoffWindowTarget.success()
                && oxygenHandoffWindow.success()
                && oxygenHandoffWindow.state
                && oxygenHandoffWindow.derivatives
                && close(
                    oxygenHandoffWindow.state->temperature,
                    197.5, 2.0e-9),
            "O2 recovery uses corrected Cp inside the C1 handoff window");
        const auto oxygenHandoffTarget = oxygenCaloric.evaluate(
            6.0e6, 200.0, {1.0}, PhaseRoot::VAPOR);
        PRPhaseRhoEResult oxygenHandoff;
        if (oxygenHandoffTarget.success())
        {
            oxygenHandoff = oxygenRecovery.recover(
                oxygenHandoffTarget.state->rho,
                oxygenHandoffTarget.state->e,
                oxygenHandoffTarget.state->massFractions,
                PRPhaseIdentity::SUPERCRITICAL,
                PRPhaseRootIdentity::PR_UNIQUE_SUPERCRITICAL_ROOT,
                200.0);
        }
        check(
            oxygenHandoffParent.valid()
                && oxygenHandoffTarget.success()
                && oxygenHandoff.success()
                && oxygenHandoff.state && oxygenHandoff.derivatives
                && close(
                    oxygenHandoff.state->temperature, 200.0, 2.0e-9),
            "O2 C1 table-to-NASA handoff supports rho-e recovery");
        if (oxygenTarget.success())
        {
            const PRPhaseRhoEResult oxygenRecovered = oxygenRecovery.recover(
                oxygenTarget.state->rho,
                oxygenTarget.state->e,
                oxygenTarget.state->massFractions,
                PRPhaseIdentity::SUPERCRITICAL,
                PRPhaseRootIdentity::PR_UNIQUE_SUPERCRITICAL_ROOT,
                240.0);
            check(oxygenRecovered.success(),
                "parent-bound unique supercritical rho-e-Y recovers");
            check(
                oxygenRecovered.state
                    && close(oxygenRecovered.state->pressure,
                        oxygenParent.pressure, 2.0e-9)
                    && close(oxygenRecovered.state->temperature,
                        oxygenParent.temperature, 2.0e-9)
                    && oxygenRecovered.phaseIdentity
                        == PRPhaseIdentity::SUPERCRITICAL
                    && oxygenRecovered.rootIdentity
                        == PRPhaseRootIdentity::PR_UNIQUE_SUPERCRITICAL_ROOT,
                "supercritical result preserves parent identity and closes p-T");

            const auto roots = oxygenKernel.roots(
                oxygenParent.pressure, oxygenParent.temperature, {1.0});
            const auto liquidAlias = oxygenCaloric.evaluate(
                oxygenParent.pressure, oxygenParent.temperature,
                {1.0}, PhaseRoot::LIQUID);
            check(
                roots.valid() && roots.physicalRoots.size() == 1U
                    && liquidAlias.success()
                    && close(liquidAlias.state->compressibility,
                        oxygenTarget.state->compressibility, 1.0e-13),
                "low-level liquid/vapor labels alias one supercritical root");
            check(
                oxygenRecovered.phaseIdentity
                    == PRPhaseIdentity::SUPERCRITICAL,
                "supercritical identity remains parent-owned, not root-inferred");

            const auto table = tabulatedIdealCaloric().evaluate("O2", 180.0);
            check(table.success(),
                "180 K O2 ideal Cp comes from the frozen low-temperature table");
            if (table.success() && oxygenRecovered.state
                && oxygenRecovered.derivatives)
            {
                const double molecularWeight = oxygenEos.MW({1.0});
                const double explicitCv = table.state->cp
                    - PRMixtureKernel::universalGasConstant/molecularWeight
                    + oxygenEos.cvDeparture(
                        oxygenRecovered.state->pressure,
                        oxygenRecovered.state->temperature,
                        oxygenRecovered.state->compressibility,
                        oxygenRecovered.state->massFractions);
                check(close(oxygenRecovered.derivatives->cv,
                        explicitCv, 2.0e-6, 1.0e-4)
                        && close(oxygenRecovered.derivatives->dedT_rho,
                            explicitCv, 2.0e-6, 1.0e-4),
                    "low-temperature de/dT|rho uses table Cp plus PR departure Cv");
            }

            checkFailure(
                oxygenRecovery.recover(
                    oxygenTarget.state->rho,
                    oxygenTarget.state->e,
                    oxygenTarget.state->massFractions,
                    PRPhaseIdentity::SUPERCRITICAL,
                    PRPhaseRootIdentity::PR_VAPOR_ROOT,
                    180.0),
                PRPhaseRhoEStatus::INVALID_PHASE_IDENTITY,
                "supercritical parent with vapor-root identity fails closed");
        }
    }

    const SpeciesRegistry registry;
    check(
        registry.size() == 9
            && !kij.has("O2", "H2O")
            && !kij.has("C3H8", "H2O"),
        "canonical nine-species fixture contains real missing kij pairs");
    PRPhaseRhoERecovery nineSpecies(
        critical, kij, registry.names(), vaporOptions());
    check(
        !nineSpecies.configured()
            && nineSpecies.configurationStatus()
                == PRPhaseRhoEStatus::INVALID_CONFIGURATION,
        "full nine-species recovery rejects real missing kij pairs");
    std::vector<double> pureOxygen(
        static_cast<std::size_t>(registry.size()), 0.0);
    pureOxygen[static_cast<std::size_t>(registry.index("O2"))] = 1.0;
    checkFailure(
        nineSpecies.recover(
            1.0, 0.0, pureOxygen,
            PRPhaseIdentity::VAPOR,
            PRPhaseRootIdentity::PR_VAPOR_ROOT, 300.0),
        PRPhaseRhoEStatus::INVALID_CONFIGURATION,
        "unconfigured nine-species object exposes no lower-dimensional state");

    PRPhaseRhoEOptions excessiveScan = vaporOptions();
    excessiveScan.temperatureScanIntervals =
        std::numeric_limits<int>::max();
    PRPhaseRhoEOptions excessiveIterations = vaporOptions();
    excessiveIterations.maximumIterations =
        std::numeric_limits<int>::max();
    const PRPhaseRhoERecovery scanRejected(
        critical, kij, species, excessiveScan);
    const PRPhaseRhoERecovery iterationsRejected(
        critical, kij, species, excessiveIterations);
    check(
        !scanRejected.configured() && !iterationsRejected.configured(),
        "unbounded scan and iteration requests fail at configuration");

    PRPhaseRhoEOptions excessiveEnergyTolerance = vaporOptions();
    excessiveEnergyTolerance.energyAbsoluteTolerance =
        std::numeric_limits<double>::max();
    PRPhaseRhoEOptions excessiveDerivativeTolerance = vaporOptions();
    excessiveDerivativeTolerance.derivativeIdentityAbsoluteTolerance =
        std::numeric_limits<double>::max();
    const PRPhaseRhoERecovery energyToleranceRejected(
        critical, kij, species, excessiveEnergyTolerance);
    const PRPhaseRhoERecovery derivativeToleranceRejected(
        critical, kij, species, excessiveDerivativeTolerance);
    check(
        !energyToleranceRejected.configured()
            && !derivativeToleranceRejected.configured(),
        "unbounded absolute tolerances fail at configuration");

    PRPhaseRhoEOptions excessiveEnergyRelative = vaporOptions();
    excessiveEnergyRelative.energyRelativeTolerance =
        std::numeric_limits<double>::max();
    PRPhaseRhoEOptions excessiveDensityRelative = vaporOptions();
    excessiveDensityRelative.densityRelativeTolerance =
        std::numeric_limits<double>::max();
    PRPhaseRhoEOptions excessiveRootRelative = vaporOptions();
    excessiveRootRelative.rootRelativeTolerance =
        std::numeric_limits<double>::max();
    PRPhaseRhoEOptions excessiveDerivativeRelative = vaporOptions();
    excessiveDerivativeRelative.derivativeIdentityRelativeTolerance =
        std::numeric_limits<double>::max();
    const PRPhaseRhoERecovery energyRelativeRejected(
        critical, kij, species, excessiveEnergyRelative);
    const PRPhaseRhoERecovery densityRelativeRejected(
        critical, kij, species, excessiveDensityRelative);
    const PRPhaseRhoERecovery rootRelativeRejected(
        critical, kij, species, excessiveRootRelative);
    const PRPhaseRhoERecovery derivativeRelativeRejected(
        critical, kij, species, excessiveDerivativeRelative);
    check(
        !energyRelativeRejected.configured()
            && !densityRelativeRejected.configured()
            && !rootRelativeRejected.configured()
            && !derivativeRelativeRejected.configured(),
        "unbounded relative tolerances fail at configuration");

    std::cout << std::setprecision(17)
              << "vapor_status="
              << Foam::phaseEquilibrium::phaseRhoEStatusName(recovered.status)
              << '\n'
              << "vapor_iterations=" << recovered.diagnostics.iterations
              << '\n'
              << "vapor_valid_scan_samples="
              << recovered.diagnostics.validScanSamples << '\n'
              << "vapor_invalid_scan_samples="
              << recovered.diagnostics.invalidScanSamples << '\n';

    if (failures != 0)
    {
        std::cerr << "PR_PHASE_RHO_E_RECOVERY_FAIL failures="
                  << failures << " checks=" << checks << '\n';
        return 1;
    }
    std::cout << "PR_PHASE_RHO_E_RECOVERY_PASS checks=" << checks << '\n'
              << "SCOPE offline-nonproduction EquilibriumState=0 Domain=0 "
                 "quality=0 acoustics=0 capability=0 production_ready=0\n";
    return 0;
}
