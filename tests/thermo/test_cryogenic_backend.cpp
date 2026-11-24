// TransportProps dense-cryogenic backend selection focused tests.
//
// Verifies (docs/cryogenic_transport_backend.md section 5):
//   - gas-range states select LEGACY_GAS;
//   - pure/near-pure O2 or C3H8 dense states select CRYOGENIC_PURE_FLUID
//     and evaluate to the frozen table values (LOX injector state);
//   - dense multi-component mixtures select UNSUPPORTED (no silent fallback);
//   - two-phase dome queries fail closed;
//   - an unregistered table makes every dense query fail closed.
//
// Output marker: CRYOGENIC_BACKEND_PASS / _FAIL.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "CryogenicTransportTable.H"
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

const char* TABLE_PATH =
    "tests/reference/cryogenic_transport_runtime_coolprop72_v1.csv";

int run()
{
    // Reset registration state at start.
    TransportProps::registerCryogenicTable(nullptr);

    // Load and register the strict table.
    static CryogenicTransportTable table;
    if (!table.load(TABLE_PATH))
    {
        check("table load", false);
        std::cout << "CRYOGENIC_BACKEND_FAIL\n";
        return 1;
    }
    check("table load", true);
    TransportProps::registerCryogenicTable(&table);
    check("table registered", TransportProps::cryogenicTable() == &table);

    // 1. Gas-range state -> LEGACY_GAS (O2 700 K, rho=20).
    const std::vector<std::string> nine = {
        "O2","N2","H2","CH4","C2H6","C3H8","H2O","CO","CO2"};
    const std::vector<double> gasY = {
        0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.25};
    check("gas range -> LEGACY_GAS",
        TransportProps::selectBackend(nine, gasY, 700.0, 20.0)
            == TransportProps::Backend::LEGACY_GAS);

    // 2. Pure O2 dense compressed liquid (LOX injector state,
    //    130 K / 4 MPa / rho ~ 920.66) -> CRYOGENIC_PURE_FLUID.
    const std::vector<double> pureO2 = {1.0, 0,0,0,0,0,0,0,0};
    check("pure O2 dense -> CRYOGENIC_PURE_FLUID",
        TransportProps::selectBackend(nine, pureO2, 130.0, 920.66)
            == TransportProps::Backend::CRYOGENIC_PURE_FLUID);

    // 3. Pure C3H8 dense state -> CRYOGENIC_PURE_FLUID.
    const std::vector<double> pureC3H8 = {0,0,0,0,0,1.0,0,0,0};
    check("pure C3H8 dense -> CRYOGENIC_PURE_FLUID",
        TransportProps::selectBackend(nine, pureC3H8, 300.0, 500.0)
            == TransportProps::Backend::CRYOGENIC_PURE_FLUID);

    // 4. Dense multi-component mixture -> UNSUPPORTED (no silent fallback).
    const std::vector<double> denseMix = {
        0.4, 0.0, 0.0, 0.0, 0.0, 0.3, 0.3, 0.0, 0.0};
    check("dense mixture -> UNSUPPORTED",
        TransportProps::selectBackend(nine, denseMix, 130.0, 500.0)
            == TransportProps::Backend::UNSUPPORTED);

    // 5. LOX state evaluation matches the frozen table.
    const TransportProps::Result lox =
        TransportProps::evaluateCryogenicPureFluid(
            nine, pureO2, 4.0e6, 130.0, 920.66, 2023.27, 829.24);
    check("LOX cryogenic evaluate: SUCCESS", lox.success(),
        static_cast<double>(lox.status), 0.0);
    if (lox.success())
    {
        check("LOX mu liquid magnitude",
            lox.point.mu > 5.0e-5 && lox.point.mu < 2.0e-4, lox.point.mu);
        check("LOX kappa liquid magnitude",
            lox.point.kappa > 0.05 && lox.point.kappa < 0.2,
            lox.point.kappa);
    }

    // 6. Two-phase dome evaluation fails closed (O2 130 K, rho=500).
    const TransportProps::Result twoPhase =
        TransportProps::evaluateCryogenicPureFluid(
            nine, pureO2, 4.0e6, 130.0, 500.0, 1000.0, 800.0);
    check("two-phase -> UNSUPPORTED_PHASE",
        twoPhase.status == TransportProps::Status::UNSUPPORTED_PHASE,
        static_cast<double>(twoPhase.status), 0.0);

    // 7. Dense mixture evaluation fails closed.
    const TransportProps::Result mix =
        TransportProps::evaluateCryogenicPureFluid(
            nine, denseMix, 4.0e6, 130.0, 500.0, 1000.0, 800.0);
    check("dense mixture evaluate: UNSUPPORTED_PHASE",
        mix.status == TransportProps::Status::UNSUPPORTED_PHASE,
        static_cast<double>(mix.status), 0.0);

    // 8. Unregistered table -> dense query fails closed.
    TransportProps::registerCryogenicTable(nullptr);
    const TransportProps::Result unregistered =
        TransportProps::evaluateCryogenicPureFluid(
            nine, pureO2, 4.0e6, 130.0, 920.66, 2023.27, 829.24);
    check("unregistered table -> UNSUPPORTED_TEMPERATURE_RANGE",
        unregistered.status
            == TransportProps::Status::UNSUPPORTED_TEMPERATURE_RANGE,
        static_cast<double>(unregistered.status), 0.0);
    check("unregistered backend -> UNSUPPORTED",
        TransportProps::selectBackend(nine, pureO2, 130.0, 920.66)
            == TransportProps::Backend::UNSUPPORTED);
    // Restore for any later tests in the same process.
    TransportProps::registerCryogenicTable(&table);

    // 9. Near-pure (>99.9%) dense O2 is still accepted as pure-fluid.
    const std::vector<double> nearPureO2 = {
        0.9995, 0.0, 0.0, 0.0, 0.0, 0.0005, 0.0, 0.0, 0.0};
    check("near-pure O2 dense -> CRYOGENIC_PURE_FLUID",
        TransportProps::selectBackend(nine, nearPureO2, 130.0, 920.66)
            == TransportProps::Backend::CRYOGENIC_PURE_FLUID);

    std::cout << (failures == 0 ? "CRYOGENIC_BACKEND_PASS\n"
                                : "CRYOGENIC_BACKEND_FAIL\n");
    return failures == 0 ? 0 : 1;
}

} // namespace

int main()
{
    return run();
}
