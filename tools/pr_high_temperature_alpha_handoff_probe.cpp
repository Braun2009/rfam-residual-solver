#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "EquationOfState/PengRobinsonEOS.H"
#include "MixingRules/BinaryInteractionDB.H"
#include "MixingRules/SpeciesCriticalDB.H"
#include "PhaseEquilibrium/PRMixtureKernel.H"
#include "PhaseEquilibrium/PRPhaseCaloric.H"
#include "pr_high_temperature_residual_handoff_prototype.H"

namespace diagnostic = rfamHighTemperatureDiagnostic;
using namespace Foam;
using namespace Foam::phaseEquilibrium;

namespace
{

int checks = 0;

void require(bool condition, const std::string& message)
{
    ++checks;
    if (!condition)
    {
        std::cerr << "PR_HIGH_TEMPERATURE_ALPHA_HANDOFF_FAIL "
                  << message << '\n';
        std::exit(2);
    }
}

struct BackendState
{
    bool valid = false;
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double energy = std::numeric_limits<double>::quiet_NaN();
    double enthalpy = std::numeric_limits<double>::quiet_NaN();
    double cvDeparture = std::numeric_limits<double>::quiet_NaN();
    double dpdrho = std::numeric_limits<double>::quiet_NaN();
    std::size_t roots = 0;
};

BackendState evaluateBackend(
    const PRMixtureKernel& kernel,
    const PRPhaseCaloric& caloric,
    const PengRobinsonEOS& eos,
    double density,
    double temperature,
    const std::vector<double>& massFractions,
    const std::vector<double>& moleFractions)
{
    BackendState result;
    result.pressure = eos.pressure(density, temperature, massFractions);
    const PRRootSet roots = kernel.roots(
        result.pressure, temperature, moleFractions);
    const PRPhaseCaloricResult state = caloric.evaluate(
        result.pressure, temperature, moleFractions, PhaseRoot::VAPOR);
    if (!roots.valid() || !state.success()) return result;

    result.energy = state.state->e;
    result.enthalpy = state.state->h;
    result.cvDeparture = eos.cvDeparture(
        result.pressure, temperature,
        state.state->compressibility, massFractions);
    const double densityStep = 1.0e-6*density;
    const double pressureLow = eos.pressure(
        density - densityStep, temperature, massFractions);
    const double pressureHigh = eos.pressure(
        density + densityStep, temperature, massFractions);
    result.dpdrho = (pressureHigh - pressureLow)/(2.0*densityStep);
    result.roots = roots.physicalRoots.size();
    result.valid = std::isfinite(result.pressure)
        && std::isfinite(result.energy) && std::isfinite(result.enthalpy)
        && std::isfinite(result.cvDeparture) && std::isfinite(result.dpdrho);
    return result;
}

diagnostic::Component component(
    const SpeciesCriticalDB& database,
    const std::string& name,
    double massFraction)
{
    const CriticalProperties& properties = database.get(name);
    return {name, massFraction, properties.MW, properties.Tc,
            properties.Pc, properties.omega, true};
}

diagnostic::InteractionMap interactions(const BinaryInteractionDB& database)
{
    const std::vector<std::string> species = {
        "O2", "N2", "H2", "C3H8", "H2O", "CO", "CO2"};
    diagnostic::InteractionMap result;
    for (std::size_t i = 0; i < species.size(); ++i)
        for (std::size_t j = i + 1; j < species.size(); ++j)
            result[diagnostic::interactionKey(species[i], species[j])] =
                database.get(species[i], species[j]);
    return result;
}

struct MatrixCase
{
    std::string name;
    std::vector<diagnostic::Component> components;
};

std::vector<MatrixCase> matrixCases(const SpeciesCriticalDB& database)
{
    const double trace = 1.0e-8;
    const double scale = 1.0 - trace;
    return {
        {"pure-CO", {component(database, "CO", 1.0)}},
        {"pure-N2", {component(database, "N2", 1.0)}},
        {"CO-N2", {
            component(database, "CO", 0.5),
            component(database, "N2", 0.5)}},
        {"CO2-N2", {
            component(database, "CO2", 0.5),
            component(database, "N2", 0.5)}},
        {"H2O-N2", {
            component(database, "H2O", 0.5),
            component(database, "N2", 0.5)}},
        {"O2-N2", {
            component(database, "O2", 0.5),
            component(database, "N2", 0.5)}},
        {"reference-crossing", {
            component(database, "O2", 0.70),
            component(database, "N2", 0.15),
            component(database, "C3H8", 0.15)}},
        {"full-products", {
            component(database, "O2", 0.0943289863862),
            component(database, "H2", 0.1252197516978),
            component(database, "H2O", 0.254843792564),
            component(database, "CO", 0.299197804074),
            component(database, "CO2", 0.226409665278)}},
        {"full-products-trace-OH", {
            component(database, "O2", scale*0.0943289863862),
            component(database, "H2", scale*0.1252197516978),
            component(database, "H2O", scale*0.254843792564),
            component(database, "CO", scale*0.299197804074),
            component(database, "CO2", scale*0.226409665278),
            {"OH-ideal-trace", trace, 0.01700734,
             std::numeric_limits<double>::quiet_NaN(),
             std::numeric_limits<double>::quiet_NaN(),
             std::numeric_limits<double>::quiet_NaN(), false}}}
    };
}

void printMatrix(
    const std::vector<MatrixCase>& cases,
    const diagnostic::InteractionMap& binaryInteractions)
{
    for (const MatrixCase& testCase : cases)
    {
        for (const double density : {1.0, 10.0})
        {
            double minimumCv = std::numeric_limits<double>::infinity();
            double minimumCp = std::numeric_limits<double>::infinity();
            double minimumMechanical = std::numeric_limits<double>::infinity();
            double minimumSoundSpeed2 = std::numeric_limits<double>::infinity();
            double minimumConstrainedCv = std::numeric_limits<double>::infinity();
            double minimumConstrainedSound2 =
                std::numeric_limits<double>::infinity();
            for (int index = 0; index <= 620; ++index)
            {
                const double temperature = 900.0 + 5.0*index;
                const diagnostic::ThermodynamicState state =
                    diagnostic::evaluateHandoff(
                        testCase.components, binaryInteractions,
                        density, temperature);
                const diagnostic::ThermodynamicState constrained =
                    diagnostic::evaluateDirect(
                        testCase.components, binaryInteractions,
                        density, temperature,
                        diagnostic::RootPolicy::ConstrainedPositive);
                require(state.valid && constrained.valid,
                        testCase.name + " matrix validity");
                minimumCv = std::min(minimumCv, state.cv);
                minimumCp = std::min(minimumCp, state.cp);
                minimumMechanical = std::min(
                    minimumMechanical, state.dpdrhoIsothermal);
                minimumSoundSpeed2 = std::min(
                    minimumSoundSpeed2, state.soundSpeedSquared);
                minimumConstrainedCv = std::min(
                    minimumConstrainedCv, constrained.cv);
                minimumConstrainedSound2 = std::min(
                    minimumConstrainedSound2,
                    constrained.soundSpeedSquared);
            }
            require(minimumCv > 0.0 && minimumCp > 0.0
                    && minimumMechanical > 0.0
                    && minimumSoundSpeed2 > 0.0,
                    testCase.name + " handoff stability");
            require(minimumConstrainedCv > 0.0
                    && minimumConstrainedSound2 > 0.0,
                    testCase.name + " constrained stability");
            std::cout << "matrix case=" << testCase.name
                      << " rho_kg_m3=" << density
                      << " handoff_min_cv_J_molK=" << minimumCv
                      << " handoff_min_cp_J_molK=" << minimumCp
                      << " handoff_min_dpdrho=" << minimumMechanical
                      << " handoff_min_c2_m2_s2=" << minimumSoundSpeed2
                      << " constrained_min_cv_J_molK="
                      << minimumConstrainedCv
                      << " constrained_min_c2_m2_s2="
                      << minimumConstrainedSound2 << '\n';
        }
    }
}

} // namespace

int main()
{
    std::cout << std::setprecision(17);
    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB binary = BinaryInteractionDB::createDefault();
    const diagnostic::InteractionMap binaryInteractions = interactions(binary);

    const std::vector<std::string> species = {"O2", "N2", "C3H8"};
    const std::vector<double> massFractions = {0.70, 0.15, 0.15};
    const double density = 10.0;
    PRMixtureKernel kernel(critical, binary, species);
    PRPhaseCaloric caloric(critical, binary, species);
    PengRobinsonEOS eos(critical, binary);
    eos.setSpeciesNames(species);
    std::vector<double> moleFractions;
    require(kernel.configured() && caloric.configured()
            && kernel.massToMoleFractions(massFractions, moleFractions)
                == PhaseEquilibriumStatus::SUCCESS,
            "backend configuration");

    const diagnostic::Component nitrogen = component(critical, "N2", 1.0);
    const diagnostic::Component oxygen = component(critical, "O2", 1.0);
    const diagnostic::Component carbonMonoxide = component(critical, "CO", 1.0);
    const double crossing = diagnostic::alphaZeroTemperature(nitrogen);
    const double lowestZero = diagnostic::alphaZeroTemperature(carbonMonoxide);
    const double epsilon = 1.0e-8;
    const BackendState below = evaluateBackend(
        kernel, caloric, eos, density, crossing - epsilon,
        massFractions, moleFractions);
    const BackendState above = evaluateBackend(
        kernel, caloric, eos, density, crossing + epsilon,
        massFractions, moleFractions);
    require(below.valid && above.valid && below.roots == 1 && above.roots == 1,
            "backend one-sided states");
    const double backendJump = above.energy - below.energy;
    require(std::abs(backendJump - 171.55024922639132) < 5.0e-3,
            "backend 171.550249 J/kg jump");

    const diagnostic::RootDerivatives nitrogenAtZero =
        diagnostic::originalAlphaRoot(nitrogen, crossing);
    const diagnostic::RootDerivatives oxygenAtZero =
        diagnostic::originalAlphaRoot(oxygen, crossing);
    const double diagonalHalfWidth = std::sqrt(
        diagnostic::principalRootDerivativeThreshold
        /diagnostic::alphaScale(nitrogen))
        /std::abs(nitrogenAtZero.first);
    const double crossHalfWidth =
        diagnostic::principalRootDerivativeThreshold
        /(std::sqrt(diagnostic::alphaScale(nitrogen)
                   *diagnostic::alphaScale(oxygen))
          *std::abs(oxygenAtZero.value)*std::abs(nitrogenAtZero.first));
    const double crossingUlp = std::nextafter(
        crossing, std::numeric_limits<double>::infinity()) - crossing;
    require(diagonalHalfWidth > crossingUlp && crossHalfWidth < crossingUlp,
            "threshold representability classification");

    const std::vector<diagnostic::Component> reference = {
        component(critical, "O2", 0.70),
        component(critical, "N2", 0.15),
        component(critical, "C3H8", 0.15)};
    const double molecularWeight = diagnostic::mixtureMolecularWeight(reference);
    const auto energyJump = [&](diagnostic::RootPolicy policy)
    {
        const diagnostic::ThermodynamicState low = diagnostic::evaluateDirect(
            reference, binaryInteractions, density, crossing - epsilon, policy);
        const diagnostic::ThermodynamicState high = diagnostic::evaluateDirect(
            reference, binaryInteractions, density, crossing + epsilon, policy);
        require(low.valid && high.valid, "route one-sided states");
        return (high.internalEnergy - low.internalEnergy)/molecularWeight;
    };
    const double prototypePrincipalJump =
        energyJump(diagnostic::RootPolicy::PrincipalAbsolute);
    const double signedJump = energyJump(diagnostic::RootPolicy::SignedDiagnostic);
    const double constrainedJump =
        energyJump(diagnostic::RootPolicy::ConstrainedPositive);
    const diagnostic::ThermodynamicState handoffLow =
        diagnostic::evaluateHandoff(
            reference, binaryInteractions, density, crossing - epsilon);
    const diagnostic::ThermodynamicState handoffHigh =
        diagnostic::evaluateHandoff(
            reference, binaryInteractions, density, crossing + epsilon);
    require(handoffLow.valid && handoffHigh.valid, "handoff one-sided states");
    const double handoffJump =
        (handoffHigh.internalEnergy - handoffLow.internalEnergy)
        /molecularWeight;
    require(std::abs(prototypePrincipalJump - backendJump) < 5.0e-3,
            "prototype/backend principal route agreement");
    require(std::abs(signedJump) < 1.0e-2
            && std::abs(constrainedJump) < 1.0e-2
            && std::abs(handoffJump) < 1.0e-2,
            "continuous diagnostic routes");

    const double signedTemperature = 1380.0;
    const diagnostic::RootDerivatives signedCO =
        diagnostic::originalAlphaRoot(carbonMonoxide, signedTemperature);
    const diagnostic::RootDerivatives signedN2 =
        diagnostic::originalAlphaRoot(nitrogen, signedTemperature);
    const double signedCrossAttraction =
        std::sqrt(diagnostic::alphaScale(carbonMonoxide)
                  *diagnostic::alphaScale(nitrogen))
        *signedCO.value*signedN2.value
        *(1.0 - binary.get("CO", "N2"));
    require(signedCrossAttraction < 0.0,
            "signed route admits negative cross attraction");

    std::cout << "schema=rfam.pr-high-temperature-alpha-handoff-probe/v1\n"
              << "scope=offline-nonproduction-diagnostic-only\n"
              << "capability_change=0\n"
              << "principal_threshold="
              << diagnostic::principalRootDerivativeThreshold << '\n'
              << "handoff_lower_K="
              << diagnostic::handoffTemperatureLower << '\n'
              << "handoff_upper_K="
              << diagnostic::handoffTemperatureUpper << '\n'
              << "handoff_continuity=C3-residual-potential;C1-Cp-Cv-sound\n"
              << "high_temperature_residual_authority=zero-ideal-limit\n"
              << "lowest_hazardous_species=CO\n"
              << "lowest_hazardous_zero_K=" << lowestZero << '\n'
              << "handoff_zero_margin_K="
              << lowestZero - diagnostic::handoffTemperatureUpper << '\n'
              << "n2_zero_K=" << crossing << '\n'
              << "n2_zero_ulp_K=" << crossingUlp << '\n'
              << "threshold_diagonal_half_width_K="
              << diagonalHalfWidth << '\n'
              << "threshold_n2_o2_cross_half_width_K="
              << crossHalfWidth << '\n'
              << "backend_below_p_Pa=" << below.pressure << '\n'
              << "backend_above_p_Pa=" << above.pressure << '\n'
              << "backend_below_e_J_kg=" << below.energy << '\n'
              << "backend_above_e_J_kg=" << above.energy << '\n'
              << "backend_jump_J_kg=" << backendJump << '\n'
              << "prototype_principal_jump_J_kg="
              << prototypePrincipalJump << '\n'
              << "handoff_jump_J_kg=" << handoffJump << '\n'
              << "constrained_alpha_jump_J_kg=" << constrainedJump << '\n'
              << "signed_root_jump_J_kg=" << signedJump << '\n'
              << "signed_root_cross_attraction_1380K="
              << signedCrossAttraction << '\n'
              << "signed_root_authorization=diagnostic-only\n"
              << "constrained_alpha_authorization=secondary-unvalidated\n";

    printMatrix(matrixCases(critical), binaryInteractions);

    const double highTemperature = 4000.0;
    const double volume = molecularWeight/density;
    const diagnostic::PotentialDerivatives handoffPotential =
        diagnostic::handoffResidualPotential(
            reference, binaryInteractions, highTemperature, volume);
    const diagnostic::PotentialDerivatives constrainedPotential =
        diagnostic::directResidualPotential(
            reference, binaryInteractions, highTemperature, volume,
            diagnostic::RootPolicy::ConstrainedPositive);
    require(handoffPotential.valid && constrainedPotential.valid
            && handoffPotential.value == 0.0,
            "high-temperature route comparison");
    std::cout << "route_comparison T_K=" << highTemperature
              << " handoff_Ares_J_mol=" << handoffPotential.value
              << " constrained_Ares_J_mol=" << constrainedPotential.value
              << '\n'
              << "PR_HIGH_TEMPERATURE_ALPHA_HANDOFF_PASS checks="
              << checks << '\n';
    return 0;
}
