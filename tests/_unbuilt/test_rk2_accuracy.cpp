/*---------------------------------------------------------------------------*\
    test_rk2_accuracy - TEMPORAL order-of-accuracy via Richardson self-convergence.

    The previous version compared RK2 vs Euler TOTAL error against an exact
    Gaussian. That was wrong: on a fixed grid the AUSM SPATIAL dissipation
    dominates the total error and swamps the temporal term, so RK2 and Euler
    looked identical. The temporal order is invisible in the total error.

    Correct method (isolates temporal error): run the SAME spatial grid at three
    timesteps dt, dt/2, dt/4 to the SAME end time with FIXED dt (no CFL). Because
    the spatial operator is byte-identical across the three runs, the spatial
    error cancels in the successive differences:

        u_dt   - u*  = C dt^p                 (u* = exact-in-time, discrete-in-space)
        u_dt   - u_dt2 = C dt^p (1 - 2^-p)
        u_dt2  - u_dt4 = C (dt/2)^p (1 - 2^-p)
        ratio  = |u_dt - u_dt2| / |u_dt2 - u_dt4| = 2^p

    So observed order p = log2(ratio):
        Euler  -> ratio ~ 2   (p ~ 1)
        RK2    -> ratio ~ 4   (p ~ 2)

    Quantitative gate: RK2 observed order > 1.6 (clearly super-linear, near 2),
    and RK2 order meaningfully higher than Euler order.
\*---------------------------------------------------------------------------*/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

#include "src/variables/ConservedVars.H"
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"
#include "src/framework/TimeIntegrator.H"
#include "src/framework/SolverBase.H"

using namespace Foam;
using namespace std;

static int g_pass=0, g_fail=0;
static void check(const string& n, bool ok, double got=0, double ref=0)
{
    cout << "  [" << (ok?"PASS":"FAIL") << "] " << n;
    if (!ok) cout << fixed << setprecision(4)
                  << "  (got=" << got << " ref=" << ref << ")";
    cout << "\n";
    ok ? ++g_pass : ++g_fail;
}

static double gaussian(double x, double x0, double sigma)
{
    double dx = x - x0;
    return std::exp(-dx*dx / (2*sigma*sigma));
}

// Build the SAME initial condition every time so the three dt runs are
// identical except for the timestep.
static void initField(Domain& dom)
{
    const int N = dom.N();
    const double p0=1e5, T0=300.0, u0=100.0, x0=0.5, sigma=0.05;
    for (int i = 0; i < N; ++i)
    {
        PrimitiveVars& w = dom.stagedW(i);
        w.T = T0; w.U[0] = u0; w.U[1]=w.U[2]=0.0;
        w.Y[0]=1.0; w.Y[1]=0.0;
        double x = dom.x(i);
        double rho_base = p0 / (287.0 * T0);
        double rho_pert = 0.5 * rho_base * gaussian(x, x0, sigma);
        w.p = (rho_base + rho_pert) * 287.0 * T0;   // hold T fixed
    }
    dom.primitiveToConserved();
}

// Run to tEnd with a FIXED timestep (nSteps of size dt), return final rho field.
static vector<double> runFixedDt(bool useRK2, int N, double L,
                                 double dt, int nSteps)
{
    Domain dom(N, L, 2);
    initField(dom);

    SolverBase solver(dom);
    if (useRK2)
        solver.setIntegrator(unique_ptr<TimeIntegrator>(new RK2Integrator()));
    solver.addOperator(unique_ptr<Operator>(
        new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));

    for (int s = 0; s < nSteps && !dom.tripped(); ++s)
        solver.step(dt);

    vector<double> rho(N);
    for (int i = 0; i < N; ++i) rho[i] = dom.U(i).rho;
    return rho;
}

static double l2diff(const vector<double>& a, const vector<double>& b)
{
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) { double d=a[i]-b[i]; s += d*d; }
    return std::sqrt(s / a.size());
}

// Observed temporal order from a 3-level dt refinement (same grid, same tEnd).
static double observedOrder(bool useRK2, int N, double L,
                            double dt0, int nSteps0, const char* label)
{
    // three runs: dt0 (nSteps0), dt0/2 (2*nSteps0), dt0/4 (4*nSteps0)
    vector<double> u1 = runFixedDt(useRK2, N, L, dt0,      nSteps0);
    vector<double> u2 = runFixedDt(useRK2, N, L, dt0/2.0,  nSteps0*2);
    vector<double> u4 = runFixedDt(useRK2, N, L, dt0/4.0,  nSteps0*4);

    double d12 = l2diff(u1, u2);
    double d24 = l2diff(u2, u4);
    double ratio = d12 / (d24 + 1e-30);
    double order = std::log2(ratio);

    cout << "  [" << label << "]  |u_dt - u_dt2| = " << scientific << setprecision(3) << d12
         << "   |u_dt2 - u_dt4| = " << d24 << "\n";
    cout << "            ratio = " << fixed << setprecision(3) << ratio
         << "   observed order p = " << order << "\n";
    return order;
}

int main()
{
    cout << "=== RK2 temporal order (Richardson self-convergence) ===\n\n";

    const int    N   = 200;
    const double L   = 1.0;
    // dx = 0.005; max wavespeed ~ u+a ~ 100+330 = 430; CFL~0.34 at dt0=4e-6.
    const double dt0 = 4.0e-6;
    const int    nSteps0 = 100;    // tEnd = 4e-4 s (well before pulse exits)

    cout << "  Grid N=" << N << "  dt0=" << scientific << setprecision(1) << dt0
         << "  tEnd=" << dt0*nSteps0 << "s (same grid, 3 dt levels)\n\n";

    double pEuler = observedOrder(false, N, L, dt0, nSteps0, "Euler");
    cout << "\n";
    double pRK2   = observedOrder(true,  N, L, dt0, nSteps0, "RK2  ");

    cout << "\n";
    cout << fixed << setprecision(3);
    cout << "  Euler observed order = " << pEuler << "  (expect ~1)\n";
    cout << "  RK2   observed order = " << pRK2   << "  (expect ~2)\n\n";

    check("Euler order in [0.7, 1.4] (first-order)",
          pEuler > 0.7 && pEuler < 1.4, pEuler, 1.0);
    check("RK2 order > 1.6 (clearly second-order)", pRK2 > 1.6, pRK2, 1.6);
    check("RK2 order exceeds Euler by > 0.5", (pRK2 - pEuler) > 0.5, pRK2-pEuler, 0.5);

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail==0 ? "RK2_ACCURACY_PASS" : "RK2_ACCURACY_FAIL") << "\n";
    return g_fail==0 ? 0 : 1;
}
