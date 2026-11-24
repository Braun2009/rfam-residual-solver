#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "src/combustion/FiniteRate/TwoStepCO.H"
#include "src/combustion/FiniteRate/CVODEWrapper.H"

using namespace Foam;
using namespace std;

// Species order {C3H8, O2, CO, CO2, H2O, N2}
static void elementMoles(const vector<double>& Y, const vector<double>& MW,
                         double& C, double& H, double& O)
{
    static const double nC[6] = {3, 0, 1, 1, 0, 0};
    static const double nH[6] = {8, 0, 0, 0, 2, 0};
    static const double nO[6] = {0, 2, 1, 2, 1, 0};
    C = H = O = 0.0;
    for (int i = 0; i < 6; ++i)
    {
        double n = Y[i] / MW[i];          // mol per kg mixture
        C += n * nC[i]; H += n * nH[i]; O += n * nO[i];
    }
}

int main()
{
    cout << "=== Phase 6: SUNDIALS CVODE Stiff Solver (0-D, 2-step C3H8) ===\n\n";

    // Same operating point as Phase 5 (test_phase5_rosenbrock)
    TwoStepCO chem(1e10, 0.0, 125520.0, 2e12, 0.0, 167360.0);
    const auto& MW = chem.molecularWeights();

    const double p = 4e6, rho = 15.0;
    vector<double> Y0 = {0.216, 0.784, 0.0, 0.0, 0.0, 0.0};   // near-stoichiometric
    const double T0 = 1200.0;
    const double t_end = 2e-4;        // 200 us

    CVODEWrapper cv(&chem, 1500.0);   // Cv=1500, matches RosenbrockSolver/ZeroD
    cv.setTolerances(1e-8, 1e-12);
    cv.init(0.0, p, rho, Y0, T0);

    double C0, H0, O0; elementMoles(Y0, MW, C0, H0, O0);

    cout << fixed;
    cout << "  t(us)    T(K)     Y_C3H8   Y_O2     Y_CO     Y_CO2    Y_H2O      sumY\n";
    cout << "  -----   ------    -------  -------  -------  -------  -------   --------\n";

    vector<double> Y = Y0; double T = T0;
    const int nout = 100;
    for (int k = 0; k <= nout; ++k)
    {
        double tout = t_end * double(k) / nout;
        if (k > 0 && !cv.advanceTo(tout, Y, T))
        {
            cout << "  CVODE FAILED at t=" << tout << "\n";
            return 2;
        }
        if (k % 10 == 0)
        {
            double s = 0; for (double y : Y) s += y;
            cout << "  " << setw(5) << setprecision(1) << tout * 1e6
                 << "   " << setw(6) << setprecision(0) << T
                 << "   " << setprecision(4)
                 << setw(8) << Y[0] << " " << setw(8) << Y[1] << " "
                 << setw(8) << Y[2] << " " << setw(8) << Y[3] << " "
                 << setw(8) << Y[4]
                 << "   " << setw(8) << setprecision(6) << s << "\n";
        }
    }

    cout << "\n  Final: T = " << setprecision(1) << T << " K\n";
    cout << "  CVODE stats: steps=" << cv.numSteps()
         << "  rhs=" << cv.numRhsEvals()
         << "  jac=" << cv.numJacEvals() << "\n\n";

    auto relerr = [](double a, double b){ return fabs(a - b) / (fabs(b) + 1e-30); };

    // ---- Test 1: element conservation (C/H/O) ----
    double Cf, Hf, Of; elementMoles(Y, MW, Cf, Hf, Of);
    double eC = relerr(Cf, C0), eH = relerr(Hf, H0), eO = relerr(Of, O0);
    bool t1 = (eC < 1e-4 && eH < 1e-4 && eO < 1e-4);
    cout << scientific << setprecision(2);
    cout << "[Test 1] element conservation  relerr  C=" << eC
         << "  H=" << eH << "  O=" << eO
         << "  => " << (t1 ? "PASS" : "FAIL") << "\n";

    // ---- Test 2: heat-release consistency  qdot == -sum(hf_i * omega_i/MW_i) ----
    // formation enthalpies [J/mol], order {C3H8,O2,CO,CO2,H2O,N2}
    const double hf[6] = {-103.85e3, 0.0, -110.53e3, -393.51e3, -241.82e3, 0.0};
    vector<double> Ymid = {0.10, 0.50, 0.05, 0.05, 0.10, 0.0};
    vector<double> omega; double qdot;
    chem.computeSourceTerms(2500.0, p, rho, Ymid, omega, qdot);
    double qdot_hf = 0.0;
    for (int i = 0; i < 6; ++i) qdot_hf += -hf[i] * (omega[i] / MW[i]);   // W/m^3
    double e2 = relerr(qdot_hf, qdot);
    bool t2 = (e2 < 1e-3);
    cout << "[Test 2] heat-release consistency  qdot(mech)=" << qdot
         << "  qdot(-sum hf*omega)=" << qdot_hf << "  relerr=" << e2
         << "  => " << (t2 ? "PASS" : "FAIL") << "\n";

    // ---- Test 3: realizability (non-negative + sumY=1 + T>0) ----
    double sumY = 0; bool nonneg = true;
    for (double y : Y) { sumY += y; if (y < -1e-8) nonneg = false; }
    bool t3 = nonneg && fabs(sumY - 1.0) < 1e-4 && T > 0.0;
    cout << "[Test 3] realizability  min(Y)>=0:" << (nonneg ? "yes" : "no")
         << "  |sumY-1|=" << fabs(sumY - 1.0)
         << "  T=" << fixed << setprecision(1) << T
         << "  => " << (t3 ? "PASS" : "FAIL") << "\n\n";

    bool all = t1 && t2 && t3;
    cout << (all ? "Phase 6 COMPLETE - all three conservation tests PASS"
                 : "Phase 6 - SOME TESTS FAILED") << "\n";
    return all ? 0 : 1;
}
