// Strict runtime dense-fluid transport table focused contract tests.
//
// Verifies:
//   - table load and schema/provenance binding;
//   - the Vortex LOX injector state (O2 130 K, rho ~ 920.66 kg/m3,
//     4 MPa compressed liquid) returns finite mu/kappa close to the frozen
//     CoolProp 7.2.0 HEOS reference;
//   - the propane injector state (C3H8 450 K, 4 MPa) returns finite values;
//   - a two-phase dome query fails closed with TWO_PHASE;
//   - out-of-table and unknown-species queries fail closed;
//   - never silently clamps or extrapolates.
//
// Output marker: CRYOGENIC_TRANSPORT_TABLE_PASS / _FAIL.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "CryogenicTransportTable.H"

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

const char* TABLE_PATH =
    "tests/reference/cryogenic_transport_runtime_coolprop72_v1.csv";

int run()
{
    CryogenicTransportTable table;
    if (!table.load(TABLE_PATH))
    {
        check("table load", false);
        std::cout << "CRYOGENIC_TRANSPORT_TABLE_FAIL\n";
        return 1;
    }
    check("table load", table.has("O2") && table.has("C3H8"));
    check("provenance recorded",
        table.provenanceSha256().find("loaded:") == 0);

    // 1. LOX injector state: O2 130 K, rho=920.66 kg/m3 (compressed liquid).
    //    CoolProp HEOS reference: mu=8.2105376e-5, kappa=9.6354290e-2.
    CryogenicTransportTableResult lox = table.evaluate("O2", 130.0, 920.66);
    check("LOX 130K/920.66: SUCCESS",
        lox.success(), static_cast<double>(lox.status), 0.0);
    if (lox.success())
    {
        const double muError = relativeError(lox.properties->mu, 8.2105376e-5);
        const double kappaError =
            relativeError(lox.properties->kappa, 9.6354290e-2);
        check("LOX mu within 5% of HEOS", muError < 0.05,
            lox.properties->mu, 8.2105376e-5);
        check("LOX kappa within 5% of HEOS", kappaError < 0.05,
            lox.properties->kappa, 9.6354290e-2);
        check("LOX mu physical (liquid magnitude)",
            lox.properties->mu > 5.0e-5 && lox.properties->mu < 2.0e-4,
            lox.properties->mu);
        check("LOX kappa physical (liquid magnitude)",
            lox.properties->kappa > 0.05 && lox.properties->kappa < 0.2,
            lox.properties->kappa);
    }

    // 2. Propane injector state: C3H8 450 K, 4 MPa (supercritical gas,
    //    rho ~ 57.5 kg/m3).  CoolProp: mu=1.3652e-5, kappa=4.3427e-2.
    CryogenicTransportTableResult propane =
        table.evaluate("C3H8", 450.0, 57.5);
    check("C3H8 450K/57.5: SUCCESS",
        propane.success(), static_cast<double>(propane.status), 0.0);
    if (propane.success())
    {
        const double muError =
            relativeError(propane.properties->mu, 1.3652e-5);
        const double kappaError =
            relativeError(propane.properties->kappa, 4.3427e-2);
        check("C3H8 mu within 5% of HEOS", muError < 0.05,
            propane.properties->mu, 1.3652e-5);
        check("C3H8 kappa within 5% of HEOS", kappaError < 0.05,
            propane.properties->kappa, 4.3427e-2);
    }

    // 3. Two-phase dome: O2 130 K, rho=500 kg/m3 is inside
    //    (rho_vap=68.4, rho_liq=902.5).
    CryogenicTransportTableResult twoPhase =
        table.evaluate("O2", 130.0, 500.0);
    check("O2 130K/500 two-phase: TWO_PHASE",
        twoPhase.status == CryogenicTransportTableStatus::TWO_PHASE,
        static_cast<double>(twoPhase.status), 0.0);

    // 4. Gas branch: O2 130 K, rho=50 kg/m3 is a stable gas.
    CryogenicTransportTableResult gas = table.evaluate("O2", 130.0, 50.0);
    check("O2 130K/50 gas: SUCCESS", gas.success());

    // 5. Out of table temperature.
    CryogenicTransportTableResult hot =
        table.evaluate("O2", 500.0, 920.66);
    check("O2 500K: OUT_OF_RANGE",
        hot.status == CryogenicTransportTableStatus::OUT_OF_RANGE,
        static_cast<double>(hot.status), 0.0);

    // 6. Out of table density.
    CryogenicTransportTableResult dense =
        table.evaluate("O2", 130.0, 5000.0);
    check("O2 rho=5000: OUT_OF_RANGE",
        dense.status == CryogenicTransportTableStatus::OUT_OF_RANGE,
        static_cast<double>(dense.status), 0.0);

    // 7. Unknown species.
    CryogenicTransportTableResult unknown =
        table.evaluate("N2", 300.0, 1.0);
    check("N2: UNKNOWN_SPECIES",
        unknown.status == CryogenicTransportTableStatus::UNKNOWN_SPECIES,
        static_cast<double>(unknown.status), 0.0);

    // 8. Invalid input.
    CryogenicTransportTableResult invalid =
        table.evaluate("O2", -1.0, 920.66);
    check("negative T: INVALID_INPUT",
        invalid.status == CryogenicTransportTableStatus::INVALID_INPUT,
        static_cast<double>(invalid.status), 0.0);

    // 9. Never silently clamps: a query just above the last node must be
    //    OUT_OF_RANGE, not the last node value.
    CryogenicTransportTableResult beyond =
        table.evaluate("O2", 130.0, 1300.0);
    check("O2 rho=1300 beyond max node: OUT_OF_RANGE",
        beyond.status == CryogenicTransportTableStatus::OUT_OF_RANGE,
        static_cast<double>(beyond.status), 0.0);

    std::cout << (failures == 0
        ? "CRYOGENIC_TRANSPORT_TABLE_PASS\n"
        : "CRYOGENIC_TRANSPORT_TABLE_FAIL\n");
    return failures == 0 ? 0 : 1;
}

} // namespace

int main()
{
    return run();
}
