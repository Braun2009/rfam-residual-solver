#include <iostream>
#include <iomanip>
#include <cmath>
#include "src/cht/CHTCoupler.H"

using namespace Foam;
using namespace std;

int main()
{
    cout << "=== Stage 4 CHT: Monolithic Conjugate Heat Transfer (2-layer wall) ===\n\n";

    // Composite wall: [convection h,Tinf] | A (fluid-side layer) | B (solid wall) | [fixed TR]
    const int    NA = 40, NB = 80;
    const double LA = 0.002, kA = 0.5,  rcA = 1.0e6;   // low-k "fluid-side" layer
    const double LB = 0.008, kB = 20.0, rcB = 3.2e6;   // higher-k solid wall
    const double h  = 1000.0, Tinf = 1000.0, TR = 300.0;

    CHTCoupler cht(NA, LA, kA, rcA, NB, LB, kB, rcB, 300.0);
    auto left  = SolidBC::convection(h, Tinf);
    auto right = SolidBC::fixedT(TR);
    cht.solveSteady(left, right);

    // Analytic series resistance
    const double R     = 1.0/h + LA/kA + LB/kB;
    const double q_an  = (Tinf - TR)/R;
    const double T0    = Tinf - q_an/h;        // inner surface (x=0+)
    const double T_int = T0   - q_an*LA/kA;    // A|B interface
    auto Tan = [&](double x){ return (x <= LA) ? (T0 - q_an*x/kA)
                                               : (T_int - q_an*(x-LA)/kB); };

    cout << fixed << setprecision(2);
    cout << "Analytic: q=" << q_an << " W/m2, inner-wall T=" << T0
         << " K, interface T=" << T_int << " K\n\n";
    cout << "  x(mm)   T_num(K)   T_analytic(K)   err(K)   [domain]\n";
    cout << "  -----   --------   -------------   ------   --------\n";
    double maxerr = 0.0;
    for (int i = 0; i < cht.M(); i += cht.M()/12)
    {
        double xi = cht.x(i), Ta = Tan(xi), e = fabs(cht.T(i)-Ta);
        maxerr = max(maxerr, e);
        cout << "  " << setw(5) << xi*1000.0 << "   " << setw(8) << cht.T(i)
             << "   " << setw(11) << Ta << "   " << setprecision(3) << e << setprecision(2)
             << "   " << (i < cht.NA() ? "A(fluid)" : "B(solid)") << "\n";
    }

    // Interface continuity & flux balance
    double q_num   = cht.interfaceFlux();
    double q_rel   = fabs(q_num - q_an)/q_an;
    // interface temperature reconstructed from A-last cell + its half-cell drop
    double Tint_num = cht.T(cht.NA()-1) - q_num*(0.5*(cht.LA()/cht.NA())/kA);
    double Tint_err = fabs(Tint_num - T_int);

    cout << "\nInterface flux:  num=" << q_num << "  analytic=" << q_an
         << "  rel.err=" << scientific << setprecision(2) << q_rel << fixed << "\n";
    cout << "Interface temp:  num=" << setprecision(2) << Tint_num
         << "  analytic=" << T_int << "  err=" << setprecision(3) << Tint_err << " K\n\n";

    bool ok = (maxerr < 2.0) && (q_rel < 0.02) && (Tint_err < 1.0);
    cout << "[Monolithic CHT] field maxerr=" << setprecision(3) << maxerr
         << " K  => " << (ok ? "PASS" : "FAIL") << "\n\n";

    cout << (ok ? "Stage 4 CHTCoupler VALIDATED (interface continuity + flux balance)"
                : "Stage 4 CHTCoupler - CHECK NEEDED") << "\n";
    return ok ? 0 : 1;
}
