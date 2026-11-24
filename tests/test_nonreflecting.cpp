/*---------------------------------------------------------------------------*\
    test_nonreflecting - Stage 5+: LODI non-reflecting outlet vs zero-gradient.

    Setup: a quiescent 2D channel (uniform p0, T0, small mean outflow u0) with
    a Gaussian OVER-pressure pulse initialised near the middle. The pulse splits
    into a left- and right-running acoustic wave. The right-runner reaches the
    +x outlet and (partly) reflects. We measure the reflected wave that comes
    back into the interior AFTER the pulse has left, for two outlet BCs:

        1. zero-gradient  (W[Nx-1] = W[Nx-2])  — reflects strongly
        2. NonReflectingOutlet (LODI)          — should reflect far less

    METRIC: peak |p - p0| in the interior at a probe column, sampled in a time
    window AFTER the outgoing pulse has passed the outlet. A good non-reflecting
    BC gives a much smaller reflected peak than zero-gradient. The test asserts
    the LODI reflected amplitude is a fraction of the zero-gradient one.

    This isolates the acoustic reflection property (no chemistry, no CHT).
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
#include "src/framework/NonReflectingOutlet.H"

using namespace Foam;
using namespace std;

static int g_pass=0, g_fail=0;
static void check(const string& n, bool ok, double got=0, double ref=0)
{
    cout << "  [" << (ok?"PASS":"FAIL") << "] " << n;
    if (!ok) cout << "  (got=" << scientific << setprecision(4) << got
                  << " ref=" << ref << ")";
    cout << "\n";
    ok ? ++g_pass : ++g_fail;
}

static PrimitiveVars prim(double p, double T, double u)
{
    PrimitiveVars w(2);
    w.p=p; w.T=T; w.U[0]=u; w.U[1]=w.U[2]=0.0;
    w.Y[0]=1.0; w.Y[1]=0.0;   // pure O2
    return w;
}

// Run the pulse case. useLODI selects the outlet BC. Returns the peak reflected
// |p-p0| at a probe column (Nx*0.3) during the post-pulse window.
static double runPulse(bool useLODI)
{
    const int Nx=120, Ny=4;
    const double Lx=0.12, Ly=0.004;
    const double p0=1.0e6, T0=300.0, u0=20.0;   // subsonic mean outflow
    const double dp=2.0e4;                       // 2% pressure pulse

    Domain dom(Nx, Ny, Lx, Ly, 2);

    // Uniform base flow + Gaussian pressure bump at x_c
    const double x_c = 0.4*Lx, w_pulse = 0.03*Lx;
    for (int i=0;i<Nx;++i) for (int j=0;j<Ny;++j)
    {
        double x = dom.x(i,j);
        double bump = dp * std::exp(-(x-x_c)*(x-x_c)/(w_pulse*w_pulse));
        dom.stagedW(i,j) = prim(p0 + bump, T0, u0);
    }
    dom.primitiveToConserved();

    NonReflectingOutlet outlet(p0, Lx);   // p_inf = p0

    SolverBase solver(dom);
    solver.addOperator(unique_ptr<Operator>(
        new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));

    const int probe = int(0.3*Nx);   // interior probe, left of pulse origin
    const double CFL=0.4;

    // sound speed ~ 330 m/s at these conditions; transit Lx/a ~ 3.6e-4 s.
    // Run ~2.5 transit times so the reflected wave has time to return.
    const double tEnd = 9.0e-4;
    double t=0.0; int iter=0;
    double reflectedPeak = 0.0;
    bool windowOpen = false;

    while (t < tEnd && !dom.tripped())
    {
        // inlet (left): fix base flow
        for (int j=0;j<Ny;++j) dom.setCell(0,j, prim(p0,T0,u0));
        // top/bottom: zero-gradient (1-cell channel behaviour; Ny small)
        for (int i=0;i<Nx;++i){ dom.setCell(i,0,dom.W(i,1)); dom.setCell(i,Ny-1,dom.W(i,Ny-2)); }
        // outlet (right)
        if (useLODI) outlet.apply(dom, solver.computeTimestep(CFL));
        else for (int j=0;j<Ny;++j) dom.setCell(Nx-1,j, dom.W(Nx-2,j));

        double dt = solver.computeTimestep(CFL);
        solver.step(dt);
        t += dt; ++iter;

        // Open the measurement window after the outgoing pulse has reached the
        // outlet (t > ~1 transit from x_c to Lx). Then track the reflected peak
        // at the probe.
        if (t > 4.0e-4) windowOpen = true;
        if (windowOpen)
        {
            double pk = std::abs(dom.W(probe, Ny/2).p - p0);
            reflectedPeak = std::max(reflectedPeak, pk);
        }
    }
    return reflectedPeak;
}

int main()
{
    cout << "=== Non-reflecting outlet (LODI) vs zero-gradient ===\n\n";
    cout << "  Acoustic pulse (2% over-pressure) hits the +x outlet;\n";
    cout << "  measure reflected |p-p0| returning to an interior probe.\n\n";

    double refl_zg   = runPulse(false);   // zero-gradient
    double refl_lodi = runPulse(true);    // LODI non-reflecting

    cout << fixed << setprecision(1);
    cout << "  zero-gradient reflected peak = " << refl_zg   << " Pa\n";
    cout << "  LODI          reflected peak = " << refl_lodi << " Pa\n";
    double ratio = refl_lodi / (refl_zg + 1e-30);
    cout << setprecision(3);
    cout << "  ratio LODI/zero-grad = " << ratio << "  (smaller = less reflection)\n\n";

    check("zero-gradient reflects measurably", refl_zg > 100.0, refl_zg);
    check("LODI reflects less than zero-gradient", refl_lodi < refl_zg, refl_lodi, refl_zg);
    check("LODI reflection < 60% of zero-gradient", ratio < 0.6, ratio, 0.6);
    check("LODI outlet stayed finite/bounded", refl_lodi < 1e5, refl_lodi);

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail==0 ? "NONREFLECTING_PASS" : "NONREFLECTING_FAIL") << "\n";
    return g_fail==0 ? 0 : 1;
}
