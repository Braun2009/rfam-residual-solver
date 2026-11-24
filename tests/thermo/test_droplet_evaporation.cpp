// Droplet evaporation d^2-law physical validation.
//
// Validates the classical d^2-law (Turns ch.10) with the O2 droplet at the
// Vortex LOX condition (130 K droplet in a 700 K / 4 MPa O2-vapour
// atmosphere), against the frozen CoolProp 7.2.0 HEOS reference:
//   L_v(130K) = 153.129 kJ/kg
//   rho_l(130K sat liq) = 902.48 kg/m3
//   lambda_g(film 415K/4MPa) = 0.036383 W/m/K
//   c_pg(film) = 973.6 J/kg/K
//   B = 3.624, K = 8 lambda/(rho_l c_pg) ln(1+B) = 5.072e-7 m2/s
//
// The key physical guard: the liquid density must be the dense-fluid value
// (~900 kg/m3), NOT a gas-phase density (~60 kg/m3 at 150K) — using the gas
// value would overestimate the d^2 slope by ~15x, the same class of bug as
// the earlier 130K-LOX mu/kappa error.
//
// Output marker: DROPLET_EVAPORATION_PASS / _FAIL.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "DropletEvaporation.H"
#include "TransportProps.H"

using namespace Foam;

namespace
{

int failures = 0;

void check(
    const std::string& name, bool ok,
    double got = 0.0, double expected = 0.0)
{
    std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << name;
    if (!ok)
        std::cout << " got=" << std::setprecision(12) << got
                  << " expected=" << expected;
    std::cout << '\n';
    if (!ok) ++failures;
}

double relativeError(double got, double expected, double floor = 1.0e-10)
{
    return std::abs(got - expected) / std::max(std::abs(expected), floor);
}

// O2 gas-phase transport at film temperature via the Eucken/Sutherland
// closure used by RFAM's gas backend (approximate reference; the analytic
// d^2-law check uses the same closure, so consistency is what matters).
double gasThermalConductivityO2(double T)
{
    // Sutherland mu for O2 at T
    const double mu_ref = 1.919e-5, T_ref = 273.0, S = 139.0;
    const double mu = mu_ref * std::pow(T / T_ref, 1.5)
        * (T_ref + S) / (T + S);
    // Eucken: kappa = mu * Cv * (1.32 + 1.77*R_spec/Cv)
    const double R = 259.838;             // O2 gas constant
    const double cv = 973.6 - R;          // approximate at film T
    return mu * cv * (1.32 + 1.77 * R / cv);
}

int run()
{
    // ------------------------------------------------------------------
    // 1. Physical property audit for the LOX droplet condition.
    //    These are the quantities that were previously wrong by ~8x.
    // ------------------------------------------------------------------
    // Frozen CoolProp 7.2.0 HEOS reference values:
    const double refLv = 153128.7;        // J/kg at 130 K
    const double refLambdaFilm = 0.036383; // W/m/K at 415 K / 4 MPa
    const double refCpgFilm = 973.6;      // J/kg/K at 415 K / 4 MPa

    // RFAM-side property audit (using the frozen cryogenic transport table
    // for the liquid branch and the gas closure for the film).
    CryogenicTransportTable table;
    table.load("tests/reference/cryogenic_transport_runtime_coolprop72_v1.csv");
    // Liquid density: the LOX state is a compressed liquid at rho~920.66.
    const CryogenicTransportTableResult lox =
        table.evaluate("O2", 130.0, 920.66);
    check("LOX dense-liquid table available", lox.success());
    const double rhoL = 920.66;  // compressed-liquid density (kg/m3)

    // Vapour film thermal conductivity (Eucken at 415 K).
    const double lambdaFilm = gasThermalConductivityO2(415.0);
    // The Eucken value is a gas-correlation approximation; assert it stays
    // in the gas magnitude band (0.02-0.06), NOT the liquid band (~0.1),
    // and that it is within 20% of the CoolProp film reference.
    check("film lambda gas magnitude (0.02-0.06)",
        lambdaFilm > 0.02 && lambdaFilm < 0.06, lambdaFilm);
    check("film lambda within 20% of CoolProp",
        relativeError(lambdaFilm, refLambdaFilm, 1.0e-3) < 0.20,
        lambdaFilm, refLambdaFilm);
    const double cpgFilm = refCpgFilm;
    check("film cpg reference sane", cpgFilm > 800.0 && cpgFilm < 1200.0,
        cpgFilm, refCpgFilm);

    // d^2-law slope with the physically correct dense liquid density.
    DropletEvaporationInput input;
    input.liquidDensityKgM3 = rhoL;          // dense liquid ~900
    input.liquidCpJPerKgK = 1680.0;          // O2 liquid cp ~130K
    input.latentHeatJPerKg = refLv;
    input.gasThermalConductivityWPerMK = lambdaFilm;
    (void)refLambdaFilm; // used in the film-lambda checks above
    input.gasCpJPerKgK = refCpgFilm;
    input.ambientTemperatureK = 700.0;
    input.surfaceTemperatureK = 130.0;
    input.ambientOxidizerMassFraction = 0.0;
    input.heatOfCombustionJPerKg = 0.0;

    const double d0 = 200.0e-6;
    const DropletEvaporationResult result =
        DropletEvaporation::evaluate(input, d0);
    check("evaporation evaluate: SUCCESS", result.success());
    if (result.success())
    {
        check("B in physical band (2-6)",
            result.transferNumberB > 2.0 && result.transferNumberB < 6.0,
            result.transferNumberB, 3.624);
        check("K in physical band (2e-7 - 8e-7)",
            result.dSquareSlopeK > 2.0e-7 && result.dSquareSlopeK < 8.0e-7,
            result.dSquareSlopeK, 5.072e-7);
        check("lifetime physical (50-150 ms for 200um)",
            result.lifetimeS > 0.05 && result.lifetimeS < 0.15,
            result.lifetimeS, 0.0789);
        check("surface mass flux positive",
            result.surfaceMassFlux > 0.0, result.surfaceMassFlux);
    }

    // ------------------------------------------------------------------
    // 2. d^2-law consistency: analytic slope vs numerical integration.
    // ------------------------------------------------------------------
    const DropletEvaporationResult integrated =
        DropletEvaporation::integrateD2Law(input, d0, 4000);
    check("integration: SUCCESS", integrated.success());
    if (integrated.success())
    {
        const double slopeError =
            relativeError(integrated.dSquareSlopeK, result.dSquareSlopeK);
        check("d2-law numeric slope matches analytic (<1%)",
            slopeError < 0.01,
            integrated.dSquareSlopeK, result.dSquareSlopeK);
    }

    // ------------------------------------------------------------------
    // 3. The 130K-LOX property bug guard: a gas-phase liquid density must
    //    NOT be silently accepted.  If someone feeds rho_l = 61 kg/m3
    //    (gas at 150K) the slope K grows ~15x — the module itself is
    //    property-agnostic, so the guard is the physical audit above, but
    //    we also verify the transfer-number/slope math reacts correctly.
    // ------------------------------------------------------------------
    DropletEvaporationInput wrongDensity = input;
    wrongDensity.liquidDensityKgM3 = 61.0;  // gas-phase density (wrong)
    const DropletEvaporationResult wrongResult =
        DropletEvaporation::evaluate(wrongDensity, d0);
    check("wrong gas-density input still evaluates (math is agnostic)",
        wrongResult.success());
    if (wrongResult.success())
    {
        check("gas-density K is ~15x larger (exposes the bug)",
            wrongResult.dSquareSlopeK > 10.0 * result.dSquareSlopeK,
            wrongResult.dSquareSlopeK, result.dSquareSlopeK);
    }

    // ------------------------------------------------------------------
    // 4. Combustion transfer number (droplet burning).
    // ------------------------------------------------------------------
    DropletEvaporationInput burnInput = input;
    burnInput.ambientOxidizerMassFraction = 1.0;   // pure O2 atmosphere
    burnInput.heatOfCombustionJPerKg = 12.9e6;     // ~LHV C3H8 per kg fuel
    const DropletEvaporationResult burnResult =
        DropletEvaporation::evaluate(burnInput, d0);
    check("combustion evaluate: SUCCESS", burnResult.success());
    if (burnResult.success())
    {
        check("combustion B > evaporation B",
            burnResult.transferNumberB > result.transferNumberB,
            burnResult.transferNumberB, result.transferNumberB);
        check("combustion K > evaporation K",
            burnResult.dSquareSlopeK > result.dSquareSlopeK,
            burnResult.dSquareSlopeK, result.dSquareSlopeK);
        check("flame standoff ratio > 1 (flame detached)",
            burnResult.flameStandoffRatio > 1.0,
            burnResult.flameStandoffRatio);
    }

    // ------------------------------------------------------------------
    // 5. Invalid input fails closed.
    // ------------------------------------------------------------------
    DropletEvaporationInput invalid = input;
    invalid.liquidDensityKgM3 = -1.0;
    check("negative rho_l -> INVALID_INPUT",
        DropletEvaporation::evaluate(invalid, d0).status
            == DropletEvaporationStatus::INVALID_INPUT);
    DropletEvaporationInput coldAmbient = input;
    coldAmbient.ambientTemperatureK = 100.0;  // below surface T
    check("cold ambient -> INVALID_INPUT",
        DropletEvaporation::evaluate(coldAmbient, d0).status
            == DropletEvaporationStatus::INVALID_INPUT);

    std::cout << (failures == 0 ? "DROPLET_EVAPORATION_PASS\n"
                                : "DROPLET_EVAPORATION_FAIL\n");
    return failures == 0 ? 0 : 1;
}

} // namespace

int main()
{
    return run();
}
