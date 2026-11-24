#include <iostream>
#include <iomanip>
#include <cmath>
#include "src/cht/SolidEnergyEqn.H"

using namespace Foam;
using namespace std;

int main()
{
    cout << "=== Stage 4 CHT: 1D Solid Heat Conduction Validation ===\n\n";

    // Steady composite: left convection (h, T_inf) | solid (L,k) | right fixed T
    const int    N   = 50;
    const double L   = 0.01;     // 10 mm wall
    const double k   = 20.0;     // W/m/K
    const double rho = 8000.0, cp = 400.0;
    const double h   = 1000.0, Tinf = 1000.0, TR = 300.0;

    SolidEnergyEqn solid(N, L, rho, cp, k, 300.0);
    auto left  = SolidBC::convection(h, Tinf);
    auto right = SolidBC::fixedT(TR);
    solid.solveSteady(left, right);

    // Analytic (series resistance): q = (Tinf-TR)/(1/h + L/k); T(x)=T0face - q x/k
    const double q_an   = (Tinf - TR) / (1.0/h + L/k);
    const double T0face = Tinf - q_an/h;

    cout << fixed << setprecision(2);
    cout << "Analytic: q = " << q_an << " W/m2,  inner-wall T = " << T0face << " K\n\n";
    cout << "  x(mm)   T_num(K)   T_analytic(K)    err(K)\n";
    cout << "  -----   --------   -------------    ------\n";
    double maxerr = 0.0;
    for (int i = 0; i < N; i += N/10)
    {
        double xi  = solid.x(i);
        double Tan = T0face - q_an*xi/k;
        double e   = fabs(solid.T(i) - Tan);
        maxerr = max(maxerr, e);
        cout << "  " << setw(5) << xi*1000.0
             << "   " << setw(8) << solid.T(i)
             << "   " << setw(11) << Tan
             << "    " << setprecision(3) << e << setprecision(2) << "\n";
    }
    double q_num = solid.leftFlux(left);
    double q_rel = fabs(q_num - q_an)/q_an;
    cout << "\nInterface flux: num=" << q_num << "  analytic=" << q_an
         << "  rel.err=" << scientific << setprecision(2) << q_rel << fixed << "\n";
    bool steady_ok = (maxerr < 2.0) && (q_rel < 0.02);
    cout << "[Steady CHT] " << (steady_ok ? "PASS" : "FAIL")
         << "  (max T err = " << setprecision(3) << maxerr << " K)\n\n";

    // Transient must relax to the steady solution (implicit BE, unconditionally stable)
    SolidEnergyEqn tr(N, L, rho, cp, k, 300.0);
    const double dt = 1e-3, tend = 50.0;   // alpha=k/(rho cp)=6.25e-6, L^2/alpha~16s < 50s
    int nstep = (int)(tend/dt);
    for (int s = 0; s < nstep; ++s) tr.step(dt, left, right);
    double maxdiff = 0.0;
    for (int i = 0; i < N; ++i) maxdiff = max(maxdiff, fabs(tr.T(i) - solid.T(i)));
    bool trans_ok = (maxdiff < 0.5);
    cout << "[Transient->Steady] max|T_trans - T_steady| = " << setprecision(4) << maxdiff
         << " K  => " << (trans_ok ? "PASS" : "FAIL") << "\n\n";

    bool all = steady_ok && trans_ok;
    cout << (all ? "Stage 4 SolidEnergyEqn VALIDATED" : "Stage 4 - CHECK NEEDED") << "\n";
    return all ? 0 : 1;
}
