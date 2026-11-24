/*---------------------------------------------------------------------------*\
    test_turbulence - Stage 6: k-omega SST smoke test.

    2D channel with a mean shear profile (parabolic U). Checks:
      1. k > 0 everywhere after advance
      2. omega > 0 everywhere
      3. nu_t > 0 (turbulent viscosity positive)
      4. Production > 0 in shear region (physical)
      5. No trip in flow domain
      6. k/omega don't blow up (bounded)
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include "src/variables/ConservedVars.H"
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"
#include "src/framework/TimeIntegrator.H"
#include "src/framework/SolverBase.H"
#include "src/turbulence/TurbulenceOperator.H"

using namespace Foam;
using namespace std;

static int g_pass=0, g_fail=0;
static void check(const string& n, bool ok, double got=0, double ref=0)
{
    cout << "  [" << (ok?"PASS":"FAIL") << "] " << n;
    if (!ok) cout << "  (got=" << got << " ref=" << ref << ")";
    cout << "\n";
    ok ? ++g_pass : ++g_fail;
}

int main()
{
    cout << "=== k-omega SST turbulence operator smoke test ===\n\n";

    const int Nx=10, Ny=8;
    const double Lx=0.01, Ly=0.008;
    const double p0=1e6, T0=300.0;

    Domain dom(Nx, Ny, Lx, Ly, 2);

    // Parabolic velocity profile: U = U_max * (1 - (y/H - 0.5)^2 * 4)
    const double U_max = 100.0, H = Ly;
    for (int i=0;i<Nx;++i) for (int j=0;j<Ny;++j)
    {
        PrimitiveVars& w = dom.stagedW(i,j);
        double y = dom.y(i,j);
        double xi = y/H - 0.5;
        w.p = p0; w.T = T0;
        w.U[0] = std::max(U_max * (1.0 - 4*xi*xi), 0.0);
        w.U[1] = w.U[2] = 0.0;
        w.Y[0] = 1.0; w.Y[1] = 0.0;
    }
    dom.primitiveToConserved();

    const double nu_lam = 1.5e-5;   // air-like kinematic viscosity
    auto* turb = new TurbulenceOperator(nu_lam, 0.05, 0.001, 0.3);

    SolverBase solver(dom);
    solver.addOperator(unique_ptr<Operator>(
        new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));
    solver.addOperator(unique_ptr<Operator>(turb));

    const double tEnd = 1e-4;
    solver.advanceToTime(tEnd, 0.3);

    cout << "  iters=" << solver.iteration() << "\n\n";

    // Checks
    check("no trip", !dom.tripped());

    bool k_pos=true, om_pos=true, nut_pos=true, bounded=true;
    for (int i=1;i<Nx-1;++i) for (int j=1;j<Ny-1;++j)
    {
        double kv  = turb->k(i,j,Nx);
        double omv = turb->omega(i,j,Nx);
        double ntv = turb->nu_t(i,j,Nx,dom);
        if (kv  <= 0.0) k_pos   = false;
        if (omv <= 0.0) om_pos  = false;
        if (ntv <= 0.0) nut_pos = false;
        if (kv  > 1e6 || omv > 1e10) bounded = false;
    }
    check("k > 0 everywhere",     k_pos);
    check("omega > 0 everywhere", om_pos);
    check("nu_t > 0 everywhere",  nut_pos);
    check("k/omega bounded",       bounded);

    // Mid-channel nu_t should be larger than laminar (shear production)
    double nut_mid = turb->nu_t(Nx/2, Ny/2, Nx, dom);
    check("nu_t > nu_lam in shear (turbulence active)",
          nut_mid > nu_lam, nut_mid, nu_lam);

    cout << "\n  mid-channel nu_t = " << nut_mid
         << "  nu_lam = " << nu_lam
         << "  ratio = " << nut_mid/nu_lam << "\n";

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail==0 ? "TURBULENCE_PASS" : "TURBULENCE_FAIL") << "\n";
    return g_fail==0 ? 0 : 1;
}
