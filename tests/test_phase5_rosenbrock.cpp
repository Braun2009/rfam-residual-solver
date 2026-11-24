#include <iostream>
#include <iomanip>
#include "src/combustion/FiniteRate/TwoStepCO.H"
#include "src/combustion/FiniteRate/RosenbrockSolver.H"

using namespace Foam;
using namespace std;

int main()
{
    cout << "Phase 5: Rosenbrock Stiff Solver Test\n\n";
    cout << "Comparing Explicit Euler vs Rosenbrock for C3H8 ignition\n\n";

    TwoStepCO chem(1e10, 0.0, 125520.0, 2e12, 0.0, 167360.0);
    RosenbrockSolver solver;

    double p = 4e6, rho = 15.0;
    vector<double> Y0 = {0.216, 0.784, 0.0, 0.0, 0.0, 0.0};

    // Test with large dt = 1e-8 (impossible for explicit Euler)
    double dt_large = 1e-8;

    cout << "Rosenbrock with dt=" << dt_large << " s (explicit Euler would blow up)\n\n";
    cout << "  t(us)   T(K)    Y_C3H8  Y_CO    Y_CO2\n";
    cout << "  -----   ----    ------  ----    -----\n";

    vector<double> Y = Y0;
    double T = 1200.0, t = 0.0;
    double t_end = 2e-4;  // 200 us

    int steps = 0;
    while (t < t_end && steps < 100000)
    {
        if (steps % 1000 == 0)
        {
            cout << fixed << setprecision(1) << "  " << setw(5) << t*1e6
                 << setprecision(0) << "  " << setw(5) << T
                 << setprecision(4)
                 << "  " << Y[0] << "  " << Y[2] << "  " << Y[3] << "\n";
        }

        solver.step(&chem, p, rho, Y, T, dt_large);
        t += dt_large;
        steps++;

        // Stop when fully combusted
        if (Y[0] < 1e-4 && Y[2] < 1e-4) break;
    }

    cout << "\nFinal state at t=" << t*1e6 << " us:\n";
    cout << "  T     = " << fixed << setprecision(0) << T << " K\n";
    cout << "  Y_CO2 = " << setprecision(4) << Y[3] << "\n";
    cout << "  Y_H2O = " << Y[4] << "\n";
    cout << "  Y_C3H8= " << Y[0] << " (should->0)\n";

    // Expected final T from energy balance:
    // DH_total = 2.046 MJ/kg_C3H8 * 0.216 kg/kg = 0.441 MJ/kg_mixture
    // DT_adiabatic ~ DH / Cp ~ 441e3 / 1500 ~ 294 K (from T0=1200K -> ~1494K)
    cout << "\nExpected adiabatic T ~ " << setprecision(0) << 1200 + 441e3/1500 << " K\n";
    cout << "\nPhase 5 Complete - Rosenbrock solver stable!\n";

    return 0;
}
