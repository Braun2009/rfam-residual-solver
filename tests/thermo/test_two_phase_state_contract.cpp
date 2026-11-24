#include "PhaseEquilibrium/ConservedPhaseMixture.H"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace
{

int failures = 0;

void check(bool condition, const char* label)
{
    std::cout << (condition ? "  PASS " : "  FAIL ") << label << '\n';
    if (!condition) ++failures;
}

bool close(double lhs, double rhs, double relative = 1.0e-12)
{
    const double scale = std::max({std::abs(lhs), std::abs(rhs), 1.0});
    return std::abs(lhs - rhs) <= relative * scale;
}

} // End anonymous namespace

int main()
{
    using namespace Foam::phaseEquilibrium;

    std::cout << "=== Two-phase conserved-state algebra contract ===\n";

    constexpr double pressure = 1.0e6;
    constexpr double betaMole = 0.25;

    ConservedPhaseProperties liquid;
    liquid.molecularWeight = 0.032;
    liquid.rho = 1000.0;
    liquid.h = 100000.0;
    liquid.e = liquid.h - pressure / liquid.rho;
    liquid.s = 1200.0;
    liquid.massFractions = {0.95, 0.05};

    ConservedPhaseProperties vapor;
    vapor.molecularWeight = 0.028;
    vapor.rho = 10.0;
    vapor.h = 500000.0;
    vapor.e = vapor.h - pressure / vapor.rho;
    vapor.s = 4200.0;
    vapor.massFractions = {0.20, 0.80};

    const double betaMass = betaMole * vapor.molecularWeight /
        ((1.0 - betaMole) * liquid.molecularWeight
            + betaMole * vapor.molecularWeight);
    const std::vector<double> overall = {
        (1.0 - betaMass) * liquid.massFractions[0]
            + betaMass * vapor.massFractions[0],
        (1.0 - betaMass) * liquid.massFractions[1]
            + betaMass * vapor.massFractions[1]
    };

    const auto result = assembleConservedPhaseMixture(
        pressure, betaMole, liquid, vapor, overall);
    check(result.success(), "consistent phase state is accepted");
    check(result.state.has_value(), "success carries a committed state");

    if (result.state)
    {
        const auto& state = *result.state;
        const double expectedVolume =
            (1.0 - betaMass) / liquid.rho + betaMass / vapor.rho;
        const double expectedAlpha =
            (betaMass / vapor.rho) / expectedVolume;
        const double expectedH =
            (1.0 - betaMass) * liquid.h + betaMass * vapor.h;
        const double expectedE =
            (1.0 - betaMass) * liquid.e + betaMass * vapor.e;
        const double expectedS =
            (1.0 - betaMass) * liquid.s + betaMass * vapor.s;

        check(close(state.vaporMassFraction, betaMass),
            "molar vapor fraction converts to mass basis");
        check(close(state.vaporVolumeFraction, expectedAlpha),
            "vapor volume fraction uses phase specific volumes");
        check(close(state.rho, 1.0 / expectedVolume),
            "mixture density is the reciprocal mass-weighted volume");
        check(close(state.h, expectedH),
            "mixture enthalpy follows the mass lever rule");
        check(close(state.e, expectedE),
            "mixture internal energy follows the mass lever rule");
        check(close(state.s, expectedS),
            "mixture entropy follows the mass lever rule when phases provide it");
        check(close(state.e, state.h - pressure / state.rho),
            "mixture e=h-p/rho identity closes");
        check(close(state.overallMassFractions[0], overall[0]) &&
              close(state.overallMassFractions[1], overall[1]),
            "overall species material balance closes");
    }

    {
        auto inconsistentOverall = overall;
        inconsistentOverall[0] += 1.0e-3;
        inconsistentOverall[1] -= 1.0e-3;
        const auto rejected = assembleConservedPhaseMixture(
            pressure, betaMole, liquid, vapor, inconsistentOverall);
        check(rejected.status ==
              ConservedMixtureStatus::MATERIAL_BALANCE_FAILED,
            "inconsistent feed composition is rejected");
        check(!rejected.state.has_value(),
            "material-balance failure exposes no state");
    }

    {
        auto inconsistentEnergy = vapor;
        inconsistentEnergy.e += 1000.0;
        const auto rejected = assembleConservedPhaseMixture(
            pressure, betaMole, liquid, inconsistentEnergy, overall);
        check(rejected.status ==
              ConservedMixtureStatus::PHASE_ENERGY_IDENTITY_FAILED,
            "inconsistent phase energy reference is rejected");
        check(!rejected.state.has_value(),
            "phase-energy failure exposes no state");
    }

    {
        const auto rejected = assembleConservedPhaseMixture(
            pressure, std::numeric_limits<double>::quiet_NaN(),
            liquid, vapor, overall);
        check(rejected.status == ConservedMixtureStatus::INVALID_INPUT,
            "non-finite phase fraction is rejected");
        check(!rejected.state.has_value(),
            "invalid input exposes no state");
    }

    {
        const auto rejected = assembleConservedPhaseMixture(
            pressure, 1.01, liquid, vapor, overall);
        check(rejected.status == ConservedMixtureStatus::INVALID_INPUT,
            "out-of-range phase fraction is rejected without clipping");
    }

    {
        auto infiniteEntropy = liquid;
        infiniteEntropy.s = std::numeric_limits<double>::infinity();
        const auto rejected = assembleConservedPhaseMixture(
            pressure, betaMole, infiniteEntropy, vapor, overall);
        check(rejected.status == ConservedMixtureStatus::INVALID_INPUT,
            "infinite phase entropy is rejected, not marked unavailable");
        check(!rejected.state.has_value(),
            "invalid entropy exposes no state");
    }

    {
        auto unavailableEntropy = vapor;
        unavailableEntropy.s = std::numeric_limits<double>::quiet_NaN();
        const auto accepted = assembleConservedPhaseMixture(
            pressure, betaMole, liquid, unavailableEntropy, overall);
        check(accepted.success() && std::isnan(accepted.state->s),
            "canonical NaN preserves unavailable mixture entropy");
    }

    {
        ConservedMixtureTolerances unsafe;
        unsafe.compositionAbsolute = 1.0;
        unsafe.energyRelative = 1.0;
        const auto rejected = assembleConservedPhaseMixture(
            pressure, betaMole, liquid, vapor, overall, unsafe);
        check(rejected.status == ConservedMixtureStatus::INVALID_INPUT,
            "contract-weakening tolerances are rejected");
        check(!rejected.state.has_value(),
            "invalid tolerances expose no state");
    }

    if (failures != 0)
    {
        std::cerr << "TWO_PHASE_STATE_CONTRACT_FAIL failures="
                  << failures << '\n';
        return 1;
    }

    std::cout << "TWO_PHASE_STATE_ALGEBRA_PASS\n"
              << "NOTE: this validates only conserved mixture algebra; "
              << "TP/PH/rho-E flash and production capability remain gated.\n";
    return 0;
}
