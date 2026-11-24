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
#include "NASA7Thermo.H"
#include "PhaseEquilibrium/PRMixtureKernel.H"
#include "PhaseEquilibrium/PRPhaseCaloric.H"
#include "PhaseEquilibrium/PRTPFlash.H"

namespace
{

using Foam::BinaryInteractionDB;
using Foam::NASA7Thermo;
using Foam::PengRobinsonEOS;
using Foam::SpeciesCriticalDB;
using Foam::phaseEquilibrium::PhaseEquilibriumStatus;
using Foam::phaseEquilibrium::PhaseRoot;
using Foam::phaseEquilibrium::PRMixtureKernel;
using Foam::phaseEquilibrium::PRPhaseCaloric;
using Foam::phaseEquilibrium::PRPhaseCaloricProperties;
using Foam::phaseEquilibrium::PRPhaseCaloricResult;
using Foam::phaseEquilibrium::PRTPFlash;
using Foam::phaseEquilibrium::TPFlashResult;

constexpr double universalGasConstant = 8.314462618;
int failures = 0;

void check(bool condition, const std::string& label)
{
    std::cout << (condition ? "  PASS " : "  FAIL ") << label << '\n';
    if (!condition) ++failures;
}

bool finitePositive(double value)
{
    return std::isfinite(value) && value > 0.0;
}

bool close(
    double lhs,
    double rhs,
    double relativeTolerance,
    double absoluteTolerance = 0.0)
{
    return std::abs(lhs - rhs) <= absoluteTolerance
        + relativeTolerance*std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

struct FixedDensityEvaluation
{
    double requestedDensity = std::numeric_limits<double>::quiet_NaN();
    PRPhaseCaloricProperties state;
};

std::optional<FixedDensityEvaluation> evaluateFixedDensity(
    const PRMixtureKernel& kernel,
    const PRPhaseCaloric& caloric,
    const PengRobinsonEOS& eos,
    double density,
    double temperature,
    const std::vector<double>& massFractions,
    PhaseRoot root)
{
    if (!finitePositive(density) || !finitePositive(temperature))
        return std::nullopt;

    std::vector<double> moleFractions;
    if (kernel.massToMoleFractions(massFractions, moleFractions)
        != PhaseEquilibriumStatus::SUCCESS)
    {
        return std::nullopt;
    }

    const double pressure = eos.pressure(density, temperature, massFractions);
    if (!finitePositive(pressure)) return std::nullopt;

    const PRPhaseCaloricResult evaluated = caloric.evaluate(
        pressure, temperature, moleFractions, root);
    if (!evaluated.success() || evaluated.state->root != root)
        return std::nullopt;

    // Computing p(rho,T,Y) does not prove that rho is the requested cubic
    // branch. Reject when the explicitly selected root returns another volume.
    if (!close(evaluated.state->rho, density, 2.0e-9))
        return std::nullopt;

    FixedDensityEvaluation result;
    result.requestedDensity = density;
    result.state = *evaluated.state;
    return result;
}

struct DerivativePayload
{
    double dpdrhoT = std::numeric_limits<double>::quiet_NaN();
    double dpdTRho = std::numeric_limits<double>::quiet_NaN();
    double dedrhoT = std::numeric_limits<double>::quiet_NaN();
    double dedTRho = std::numeric_limits<double>::quiet_NaN();
};

std::optional<DerivativePayload> centralDerivatives(
    const PRMixtureKernel& kernel,
    const PRPhaseCaloric& caloric,
    const PengRobinsonEOS& eos,
    const PRPhaseCaloricProperties& base)
{
    constexpr double relativeStep = 2.0e-7;
    const double densityStep = relativeStep*base.rho;
    const double temperatureStep = relativeStep*base.temperature;
    if (!finitePositive(densityStep) || !finitePositive(temperatureStep)
        || base.rho <= densityStep || base.temperature <= temperatureStep)
    {
        return std::nullopt;
    }

    const auto rhoMinus = evaluateFixedDensity(
        kernel, caloric, eos, base.rho - densityStep, base.temperature,
        base.massFractions, base.root);
    const auto rhoPlus = evaluateFixedDensity(
        kernel, caloric, eos, base.rho + densityStep, base.temperature,
        base.massFractions, base.root);
    const auto temperatureMinus = evaluateFixedDensity(
        kernel, caloric, eos, base.rho, base.temperature - temperatureStep,
        base.massFractions, base.root);
    const auto temperaturePlus = evaluateFixedDensity(
        kernel, caloric, eos, base.rho, base.temperature + temperatureStep,
        base.massFractions, base.root);
    if (!rhoMinus || !rhoPlus || !temperatureMinus || !temperaturePlus)
        return std::nullopt;

    DerivativePayload result;
    result.dpdrhoT =
        (rhoPlus->state.pressure - rhoMinus->state.pressure)/(2.0*densityStep);
    result.dedrhoT =
        (rhoPlus->state.e - rhoMinus->state.e)/(2.0*densityStep);
    result.dpdTRho =
        (temperaturePlus->state.pressure - temperatureMinus->state.pressure)
        /(2.0*temperatureStep);
    result.dedTRho =
        (temperaturePlus->state.e - temperatureMinus->state.e)
        /(2.0*temperatureStep);
    if (!std::isfinite(result.dpdrhoT) || !std::isfinite(result.dpdTRho)
        || !std::isfinite(result.dedrhoT) || !std::isfinite(result.dedTRho))
    {
        return std::nullopt;
    }
    return result;
}

} // End anonymous namespace

int main()
{
    using namespace Foam;
    using phaseEquilibrium::EquilibriumRegion;

    std::cout << "=== Offline phase-caloric derivative feasibility ===\n";

    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB interactions = BinaryInteractionDB::createDefault();
    const std::vector<std::string> species = {"O2", "C3H8"};
    PRMixtureKernel kernel(critical, interactions, species);
    PRPhaseCaloric caloric(critical, interactions, species);
    PengRobinsonEOS eos(critical, interactions);
    eos.setSpeciesNames(species);

    check(kernel.configured() && caloric.configured(),
        "frozen O2-C3H8 phase backend configures");

    PRTPFlash flash(critical, interactions, species);
    const TPFlashResult parent =
        flash.flashMoleFractions(2.0e6, 300.0, {0.2, 0.8});
    check(parent.valid() && parent.region == EquilibriumRegion::TWO_PHASE
            && parent.vaporState.valid() && !parent.vaporMassFractions.empty(),
        "two-phase parent supplies a committed vapor composition");
    if (!parent.valid() || parent.region != EquilibriumRegion::TWO_PHASE
        || parent.vaporMassFractions.empty())
    {
        return 1;
    }

    const PRPhaseCaloricResult vapor = caloric.evaluate(
        parent.pressure, parent.temperature, parent.vaporMoleFractions,
        PhaseRoot::VAPOR);
    check(vapor.success(), "committed vapor root has a caloric payload");
    if (!vapor.success()) return 1;

    const auto fixedBase = evaluateFixedDensity(
        kernel, caloric, eos, vapor.state->rho, vapor.state->temperature,
        vapor.state->massFractions, PhaseRoot::VAPOR);
    check(fixedBase.has_value(),
        "direct p(rho,T,Y) replay retains the committed vapor root");
    if (!fixedBase) return 1;

    const auto derivatives = centralDerivatives(
        kernel, caloric, eos, fixedBase->state);
    check(derivatives.has_value(),
        "four fixed-composition phase-local derivatives are finite");
    if (!derivatives) return 1;

    check(derivatives->dpdrhoT > 0.0,
        "vapor branch has positive isothermal pressure-density derivative");
    check(derivatives->dedTRho > 0.0,
        "vapor branch has positive fixed-density caloric derivative");

    const double thermodynamicDedrho =
        (fixedBase->state.pressure
            - fixedBase->state.temperature*derivatives->dpdTRho)
        /(fixedBase->state.rho*fixedBase->state.rho);
    check(close(derivatives->dedrhoT, thermodynamicDedrho, 2.0e-7, 1.0e-5),
        "de/drho|T satisfies the pressure-temperature identity");

    const NASA7Thermo& nasa = nasa7Thermo();
    const double idealCp = nasa.cp_mass_mix(
        species, fixedBase->state.massFractions, fixedBase->state.temperature);
    const double mixtureMolecularWeight =
        eos.MW(fixedBase->state.massFractions);
    const double idealCv = idealCp - universalGasConstant/mixtureMolecularWeight;
    const double departureCv = eos.cvDeparture(
        fixedBase->state.pressure,
        fixedBase->state.temperature,
        fixedBase->state.compressibility,
        fixedBase->state.massFractions);
    const double explicitCv = idealCv + departureCv;
    check(finitePositive(explicitCv),
        "NASA ideal Cv plus PR departure Cv is positive");
    check(close(derivatives->dedTRho, explicitCv, 2.0e-7, 1.0e-5),
        "finite-difference de/dT|rho matches explicit phase Cv");

    const std::vector<double> pureOxygenMole = {1.0, 0.0};
    const PRPhaseCaloricResult supercriticalVapor = caloric.evaluate(
        6.0e6, 200.0, pureOxygenMole, PhaseRoot::VAPOR);
    const PRPhaseCaloricResult supercriticalLiquid = caloric.evaluate(
        6.0e6, 200.0, pureOxygenMole, PhaseRoot::LIQUID);
    check(supercriticalVapor.success() && supercriticalLiquid.success(),
        "single-root supercritical O2 accepts both low-level root requests");
    check(
        supercriticalVapor.success() && supercriticalLiquid.success()
            && close(supercriticalVapor.state->compressibility,
                supercriticalLiquid.state->compressibility, 1.0e-13)
            && close(supercriticalVapor.state->rho,
                supercriticalLiquid.state->rho, 1.0e-13)
            && close(supercriticalVapor.state->e,
                supercriticalLiquid.state->e, 1.0e-13),
        "current LIQUID/VAPOR root enum aliases one supercritical root");

    const auto wrongBranch = evaluateFixedDensity(
        kernel, caloric, eos, 500.0, 300.0, {0.0, 1.0}, PhaseRoot::VAPOR);
    check(!wrongBranch.has_value(),
        "fixed-density helper rejects a state represented by another cubic root");

    std::cout << std::setprecision(17)
              << "vapor_rho=" << fixedBase->state.rho << '\n'
              << "vapor_T=" << fixedBase->state.temperature << '\n'
              << "dpdrho_T=" << derivatives->dpdrhoT << '\n'
              << "dpdT_rho=" << derivatives->dpdTRho << '\n'
              << "dedrho_T=" << derivatives->dedrhoT << '\n'
              << "dedT_rho=" << derivatives->dedTRho << '\n'
              << "explicit_Cv=" << explicitCv << '\n';

    if (failures != 0)
    {
        std::cerr << "PR_PHASE_CALORIC_DERIVATIVE_FEASIBILITY_FAIL failures="
                  << failures << '\n';
        return 1;
    }
    std::cout
        << "PR_PHASE_CALORIC_DERIVATIVE_FEASIBILITY_PASS checks=13\n"
        << "SCOPE offline-diagnostic-only runtime_api=0 "
           "supercritical_identity_owner=0 production_ready=0\n";
    return 0;
}
