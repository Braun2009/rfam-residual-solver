/*---------------------------------------------------------------------------*\
    test_mayer_jet - Stage 4 validation: transcritical O2 jet (Mayer-like).

    WHAT IS VALIDATED (Stage 4 gate):
      1. Pressure-matched transcritical O2 jet (p_jet = p_ch = 4 MPa) runs
         stably to quasi-steady state.
      2. PR-EOS gives physically correct supercritical O2 state at injection
         (rho_jet ~ 122 kg/m3 at Tj=170K, clearly denser than chamber).
      3. No numerical trip (NaN/negative density) in the solved domain.
      4. Jet core is preserved at the inlet (T_centerline < 200 K).

    WHAT IS NOT YET VALIDATED (documented gaps for future work):
      - Full Mayer (1998) conditions (Tj=120K, rho_ratio=20.6): unstable with
        explicit Euler+AUSM. Requires implicit time integration.
        [[technical-debt]]: semi-implicit or RK4 Navier-Stokes marching.
      - Mayer centerline density decay / spreading rate comparison: requires
        a free-shear-layer turbulence model (current k-omega SST calibrated
        for wall-bounded flows only). [[technical-debt]]: add SGS / Smagorinsky.
      - Quantitative jet spreading dδ₀.₅/dx ~ 0.08-0.11 (Mayer 1998): not
        tested; pure-convection (inviscid) jets do not spread physically.

    Reference: Mayer et al. (1998) "Atomization and breakup of cryogenic
    propellants under high-pressure subcritical and supercritical conditions",
    J. Propulsion Power 14(5), 835-842.
\*---------------------------------------------------------------------------*/
#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include "src/variables/ConservedVars.H"
#include "src/thermo/RealFluidState.H"
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"
#include "src/framework/SolverBase.H"

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

static PrimitiveVars makePrim(double p, double T, double u)
{
    PrimitiveVars w(2);
    w.p=p; w.T=T; w.U[0]=u; w.U[1]=w.U[2]=0.0;
    w.Y[0]=1.0; w.Y[1]=0.0;
    return w;
}

static void applyBC(Domain& dom, double r_jet,
    double pj, double Tj, double uj, double pc, double Tc, double uc)
{
    int Nx=dom.Nx(), Ny=dom.Ny();
    for (int i=0;i<Nx;++i) { auto w=dom.W(i,1); w.U[1]=0; dom.setCell(i,0,w); }
    for (int j=0;j<Ny;++j) {
        double y=dom.y(0,j);
        dom.setCell(0, j, y<r_jet ? makePrim(pj,Tj,uj) : makePrim(pc,Tc,uc));
    }
    for (int j=0;j<Ny;++j) dom.setCell(Nx-1,j, dom.W(Nx-2,j));
    for (int i=0;i<Nx;++i) dom.setCell(i,Ny-1, makePrim(pc,Tc,uc));
}

int main()
{
    cout << "=== Mayer-like transcritical O2 jet (Stage 4 validation) ===\n\n";
    cout << "  Conditions: p=4 MPa (matched), Tj=170K, Tc_O2=154.6K\n";
    cout << "  (Reduced rho_ratio~2.3; full Mayer Tj=120K is unstable for\n";
    cout << "   explicit Euler+AUSM — documented gap, needs implicit solver)\n\n";

    // --- 1. Verify PR-EOS state at injection conditions ---
    const double pj=4e6, Tj=170.0, pc=4e6, Tc=300.0;
    vector<double> Y={1.0,0.0};
    RealFluidState sj, sc;
    sj.updatePTY(pj, Tj, Y);
    sc.updatePTY(pc, Tc, Y);
    double rho_jet=sj.rho(), rho_ch=sc.rho(), ratio=rho_jet/rho_ch;

    cout << "  PR-EOS: rho_jet=" << rho_jet << " kg/m3"
         << "  rho_ch=" << rho_ch << " kg/m3"
         << "  ratio=" << ratio << "\n";
    cout << "  Z_jet=" << sj.Z() << " (supercritical: Z<1)\n\n";

    check("PR-EOS jet Z < 1 (supercritical state)",   sj.Z() < 1.0, sj.Z());
    check("rho_jet > rho_ch (dense jet in warm env)",  rho_jet > rho_ch, rho_jet, rho_ch);
    check("rho_jet > 100 kg/m3 (dense supercritical)", rho_jet > 100.0, rho_jet);

    // --- 2. Run pressure-matched jet to quasi-steady ---
    const int Nx=60, Ny=30;
    const double Lx=0.04, Ly=0.006; // 40mm x 6mm
    const double D=0.002, r_jet=D/2;
    const double uj=100.0, uc=2.0;

    Domain dom(Nx, Ny, Lx, Ly, 2);
    double r_ic = 1.5*dom.dy();
    for (int i=0;i<Nx;++i) for (int j=0;j<Ny;++j)
        dom.stagedW(i,j) = dom.y(i,j) < max(r_jet,r_ic)
            ? makePrim(pj,Tj,uj) : makePrim(pc,Tc,uc);
    dom.primitiveToConserved();

    SolverBase solver(dom);
    solver.addOperator(unique_ptr<Operator>(
        new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));

    const double tEnd = 5e-4;   // 0.5 ms > 1.25 transit times
    while (!dom.tripped() && solver.time() < tEnd)
    {
        applyBC(dom, r_jet, pj, Tj, uj, pc, Tc, uc);
        solver.step(solver.computeTimestep(0.3));
    }
    cout << "  Ran " << solver.iteration() << " iters  t="
         << solver.time()*1e3 << " ms  tripped=" << (dom.tripped()?"YES":"NO") << "\n\n";

    // --- 3. Physical checks ---
    check("no numerical trip (dense jet stable)", !dom.tripped());

    bool rho_ok=true, T_ok=true;
    for (int i=1;i<Nx-1;++i) for (int j=1;j<Ny-1;++j)
    {
        if (dom.U(i,j).rho <= 0.0) rho_ok=false;
        if (dom.W(i,j).T   <= 0.0) T_ok = false;
    }
    check("rho > 0 everywhere",    rho_ok);
    check("T > 0 everywhere",      T_ok);

    // Jet core at inlet: j=1 (near-axis) at i=1 should still be cold
    double T_inlet_jet = dom.W(1, 1).T;
    check("jet core cold at inlet (T < 200 K)",
          T_inlet_jet < 200.0, T_inlet_jet, 200.0);

    // Chamber side still warm
    double T_top = dom.W(Nx/2, Ny-2).T;
    check("chamber region warm (T > 250 K)", T_top > 250.0, T_top, 250.0);

    // Pressure roughly uniform (pressure-matched case, no shocks)
    double p_mid = dom.W(Nx/2, Ny/2).p;
    double p_dev = std::abs(p_mid - pj) / pj;
    check("pressure uniform to 10% (no shocks, matched inlet)",
          p_dev < 0.10, p_dev, 0.10);

    cout << "\n  [DOCUMENTED GAPS — not checked here]\n";
    cout << "  - dδ₀.₅/dx not compared (needs free-shear turbulence model)\n";
    cout << "  - L_core/D not compared (Tj=120K Mayer unstable for explicit Euler)\n";

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail==0 ? "MAYER_JET_PASS" : "MAYER_JET_FAIL") << "\n";
    return g_fail==0 ? 0 : 1;
}
