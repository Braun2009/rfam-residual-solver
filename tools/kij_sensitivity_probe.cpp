#include "PhaseEquilibrium/ConservedPhaseMixture.H"
#include "PhaseEquilibrium/PRPhaseCaloric.H"
#include "PhaseEquilibrium/PRTPFlash.H"
#include "MixingRules/BinaryInteractionDB.H"
#include "MixingRules/SpeciesCriticalDB.H"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace Foam;
using namespace Foam::phaseEquilibrium;

namespace
{

constexpr const char* schema = "rfam.kij-sensitivity-probe/v1";

struct ProbeState
{
    const char* id;
    double pressure;
    double temperature;
    double oxygenMassFraction;
};

constexpr std::array<double, 5> kijValues = {
    -0.05, 0.0, 0.055, 0.10, 0.15};

constexpr std::array<ProbeState, 10> states = {{
    {"o2_rich_low_2mpa", 2.0e6, 180.0, 0.90},
    {"global_low_2mpa", 2.0e6, 220.0, 2.8/3.85},
    {"fuel_rich_low_2mpa", 2.0e6, 260.0, 0.20},
    {"global_transition_4mpa", 4.0e6, 300.0, 2.8/3.85},
    {"global_inlet_mix_4mpa", 4.0e6, 450.0, 2.8/3.85},
    {"global_background_4mpa", 4.0e6, 700.0, 2.8/3.85},
    {"global_ignition_4mpa", 4.0e6, 1200.0, 2.8/3.85},
    {"global_pilot_4mpa", 4.0e6, 2400.0, 2.8/3.85},
    {"o2_supercritical_control_6mpa", 6.0e6, 200.0, 1.0},
    {"global_dense_6mpa", 6.0e6, 220.0, 2.8/3.85},
}};

const char* regionName(EquilibriumRegion region)
{
    switch (region)
    {
        case EquilibriumRegion::UNKNOWN: return "UNKNOWN";
        case EquilibriumRegion::SINGLE_LIQUID: return "SINGLE_LIQUID";
        case EquilibriumRegion::SINGLE_VAPOR: return "SINGLE_VAPOR";
        case EquilibriumRegion::SINGLE_SUPERCRITICAL:
            return "SINGLE_SUPERCRITICAL";
        case EquilibriumRegion::TWO_PHASE: return "TWO_PHASE";
        case EquilibriumRegion::PURE_COMPONENT_COEXISTENCE:
            return "PURE_COMPONENT_COEXISTENCE";
    }
    return "UNKNOWN";
}

ConservedPhaseProperties conservedPhase(
    const PRPhaseCaloricProperties& phase)
{
    ConservedPhaseProperties result;
    result.molecularWeight = phase.molecularWeight;
    result.rho = phase.rho;
    result.h = phase.h;
    result.e = phase.e;
    result.s = phase.s;
    result.massFractions = phase.massFractions;
    return result;
}

struct ProbeResult
{
    std::string status;
    EquilibriumRegion region = EquilibriumRegion::UNKNOWN;
    bool qualityDefined = false;
    double vaporMassFraction = std::numeric_limits<double>::quiet_NaN();
    double vaporVolumeFraction = std::numeric_limits<double>::quiet_NaN();
    double rho = std::numeric_limits<double>::quiet_NaN();
    double e = std::numeric_limits<double>::quiet_NaN();
    double liquidXOxygen = std::numeric_limits<double>::quiet_NaN();
    double vaporXOxygen = std::numeric_limits<double>::quiet_NaN();
    double vaporYOxygen = std::numeric_limits<double>::quiet_NaN();
    double vaporYPropane = std::numeric_limits<double>::quiet_NaN();
    double vaporRhoOxygen = std::numeric_limits<double>::quiet_NaN();
    double vaporRhoPropane = std::numeric_limits<double>::quiet_NaN();
    double singleRhoOxygen = std::numeric_limits<double>::quiet_NaN();
    double singleRhoPropane = std::numeric_limits<double>::quiet_NaN();
    double materialResidual = std::numeric_limits<double>::quiet_NaN();
    double energyResidual = std::numeric_limits<double>::quiet_NaN();
};

ProbeResult evaluate(const ProbeState& input, double kijValue)
{
    ProbeResult output;
    BinaryInteractionDB kij = BinaryInteractionDB::createDefault();
    kij.add("O2", "C3H8", kijValue);

    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const bool pureOxygen = input.oxygenMassFraction == 1.0;
    const std::vector<std::string> species = pureOxygen
        ? std::vector<std::string>{"O2"}
        : std::vector<std::string>{"O2", "C3H8"};
    const std::vector<double> massFractions = pureOxygen
        ? std::vector<double>{1.0}
        : std::vector<double>{
            input.oxygenMassFraction, 1.0 - input.oxygenMassFraction};

    PRTPFlash flash(critical, kij, species);
    const TPFlashResult split = flash.flashMassFractions(
        input.pressure, input.temperature, massFractions);
    output.status = statusName(split.status);
    output.region = split.region;
    if (!split.valid() || !split.phaseSplitReady()) return output;

    PRPhaseCaloric caloric(critical, kij, species);
    if (!caloric.configured())
    {
        output.status = phaseCaloricStatusName(
            caloric.configurationStatus());
        return output;
    }

    if (split.region == EquilibriumRegion::TWO_PHASE)
    {
        const PRPhaseCaloricResult liquid = caloric.evaluate(
            input.pressure,
            input.temperature,
            split.liquidMoleFractions,
            PhaseRoot::LIQUID);
        const PRPhaseCaloricResult vapor = caloric.evaluate(
            input.pressure,
            input.temperature,
            split.vaporMoleFractions,
            PhaseRoot::VAPOR);
        if (!liquid.success() || !vapor.success())
        {
            output.status = !liquid.success()
                ? phaseCaloricStatusName(liquid.status)
                : phaseCaloricStatusName(vapor.status);
            return output;
        }

        const ConservedMixtureResult mixture = assembleConservedPhaseMixture(
            input.pressure,
            split.vaporMoleFraction,
            conservedPhase(*liquid.state),
            conservedPhase(*vapor.state),
            split.feedMassFractions);
        if (!mixture.success())
        {
            output.status = "CONSERVED_MIXTURE_FAILURE";
            return output;
        }

        output.status = "SUCCESS";
        output.qualityDefined = true;
        output.vaporMassFraction = mixture.state->vaporMassFraction;
        output.vaporVolumeFraction = mixture.state->vaporVolumeFraction;
        output.rho = mixture.state->rho;
        output.e = mixture.state->e;
        output.liquidXOxygen = split.liquidMoleFractions[0];
        output.vaporXOxygen = split.vaporMoleFractions[0];
        output.vaporYOxygen = vapor.state->massFractions[0];
        output.vaporYPropane = pureOxygen
            ? 0.0 : vapor.state->massFractions[1];
        output.vaporRhoOxygen = vapor.state->rho*output.vaporYOxygen;
        output.vaporRhoPropane = vapor.state->rho*output.vaporYPropane;
        output.materialResidual = mixture.materialResidual;
        output.energyResidual = std::max(
            mixture.phaseEnergyIdentityResidual,
            mixture.mixtureEnergyIdentityResidual);
        return output;
    }

    const PhaseRoot singleRoot =
        split.region == EquilibriumRegion::SINGLE_LIQUID
        ? PhaseRoot::LIQUID : PhaseRoot::VAPOR;
    const PRPhaseCaloricResult single = caloric.evaluate(
        input.pressure,
        input.temperature,
        split.singlePhaseMoleFractions,
        singleRoot);
    if (!single.success())
    {
        output.status = phaseCaloricStatusName(single.status);
        return output;
    }

    output.status = "SUCCESS";
    output.rho = single.state->rho;
    output.e = single.state->e;
    output.materialResidual = 0.0;
    output.energyResidual = std::max(
        single.pressureRelativeResidual,
        single.energyIdentityRelativeResidual);
    if (split.region == EquilibriumRegion::SINGLE_LIQUID)
    {
        output.qualityDefined = true;
        output.vaporMassFraction = 0.0;
        output.vaporVolumeFraction = 0.0;
        output.liquidXOxygen = split.singlePhaseMoleFractions[0];
    }
    else if (split.region == EquilibriumRegion::SINGLE_VAPOR)
    {
        output.qualityDefined = true;
        output.vaporMassFraction = 1.0;
        output.vaporVolumeFraction = 1.0;
        output.vaporXOxygen = split.singlePhaseMoleFractions[0];
        output.vaporYOxygen = single.state->massFractions[0];
        output.vaporYPropane = pureOxygen
            ? 0.0 : single.state->massFractions[1];
        output.vaporRhoOxygen = single.state->rho*output.vaporYOxygen;
        output.vaporRhoPropane = single.state->rho*output.vaporYPropane;
    }
    else if (split.region == EquilibriumRegion::SINGLE_SUPERCRITICAL)
    {
        // Quality and vapor composition do not exist in a supercritical state.
        output.singleRhoOxygen =
            single.state->rho*single.state->massFractions[0];
        output.singleRhoPropane = pureOxygen
            ? 0.0 : single.state->rho*single.state->massFractions[1];
    }
    return output;
}

void writeValue(double value)
{
    if (std::isfinite(value)) std::cout << value;
    else std::cout << "nan";
}

void writeRow(
    const ProbeState& input,
    double kijValue,
    const ProbeResult& output)
{
    std::cout << input.id << ',' << kijValue << ','
              << input.pressure << ',' << input.temperature << ','
              << input.oxygenMassFraction << ','
              << 1.0 - input.oxygenMassFraction << ','
              << output.status << ',' << regionName(output.region) << ','
              << (output.qualityDefined ? 1 : 0) << ',';
    writeValue(output.vaporMassFraction);
    std::cout << ',';
    writeValue(output.vaporVolumeFraction);
    std::cout << ',';
    writeValue(output.rho);
    std::cout << ',';
    writeValue(output.e);
    std::cout << ',';
    writeValue(output.liquidXOxygen);
    std::cout << ',';
    writeValue(output.vaporXOxygen);
    std::cout << ',';
    writeValue(output.vaporYOxygen);
    std::cout << ',';
    writeValue(output.vaporYPropane);
    std::cout << ',';
    writeValue(output.vaporRhoOxygen);
    std::cout << ',';
    writeValue(output.vaporRhoPropane);
    std::cout << ',';
    writeValue(output.singleRhoOxygen);
    std::cout << ',';
    writeValue(output.singleRhoPropane);
    std::cout << ',';
    writeValue(output.materialResidual);
    std::cout << ',';
    writeValue(output.energyResidual);
    std::cout << '\n';
}

} // End anonymous namespace

int main()
{
    std::cout.imbue(std::locale::classic());
    std::cout << std::setprecision(17);
    std::cout << "# " << schema << '\n';
    std::cout
        << "case_id,kij,pressure_pa,temperature_k,y_o2,y_c3h8,status,"
        << "region,quality_defined,vapor_mass_fraction,"
        << "vapor_volume_fraction,rho_kg_m3,e_j_kg,liquid_x_o2,"
        << "vapor_x_o2,vapor_y_o2,vapor_y_c3h8,"
        << "vapor_rho_o2_kg_m3,vapor_rho_c3h8_kg_m3,"
        << "single_rho_o2_kg_m3,single_rho_c3h8_kg_m3,"
        << "material_residual,energy_identity_residual\n";

    bool allSuccessful = true;
    for (const ProbeState& state : states)
    {
        for (double kijValue : kijValues)
        {
            const ProbeResult result = evaluate(state, kijValue);
            writeRow(state, kijValue, result);
            allSuccessful = allSuccessful && result.status == "SUCCESS";
        }
    }
    std::cout << "# RFAM_KIJ_SENSITIVITY_"
              << (allSuccessful ? "PASS" : "FAIL")
              << " rows=" << states.size()*kijValues.size() << '\n';
    return allSuccessful ? 0 : 1;
}
