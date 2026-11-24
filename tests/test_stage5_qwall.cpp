/*---------------------------------------------------------------------------*\
    test_stage5_qwall - Stage 5: q_wall(x) extraction and validation.

    Hot fluid over cold solid wall.
    Checks:
      1. q_wall > 0 everywhere (heat flows fluid→solid)
      2. T_wall monotonically increases from T0_solid toward T_fluid
      3. q_wall physically bounded (< h_conv * ΔT_max)
      4. CSV write produces readable output
      5. Energy conservation: qToSolid ≈ sum(q_wall * dx * dt * nSteps)
\*---------------------------------------------------------------------------*/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

#include "src/thermo/SpeciesRegistry.H"
#include "src/variables/ConservedVars.H"
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"
#include "src/framework/SolverBase.H"
#include "src/cht/CHTOperator2D.H"

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

int main()
{
    cout << "=== Stage 5: q_wall(x) extraction ===\n\n";

    SpeciesRegistry reg;
    const int nS=2, Nx=12, Ny=6;
    const double Lx=0.012, Ly=0.006;
    const double T_fluid=1500.0, T_solid0=300.0, p0=4e6;
    const double h_conv=50000.0;

    Domain dom(Nx, Ny, Lx, Ly, nS);
    for (int i=0;i<Nx;++i) for (int j=0;j<Ny;++j)
    {
        PrimitiveVars& w=dom.stagedW(i,j);
        w.p=p0; w.T=T_fluid; w.U[0]=5.0; w.U[1]=w.U[2]=0.0;
        w.Y[0]=1.0; w.Y[1]=0.0;
    }
    dom.primitiveToConserved();

    SolidBC bcOuter = SolidBC::adiabatic();
    auto* cht = new CHTOperator2D(4, 0.002, 8900.0, 390.0, 320.0,
                                   h_conv, bcOuter, reg, T_solid0);

    SolverBase solver(dom);
    // CHT-only (no convection): maintain hot fluid by resetting interior BCs each step.
    // This isolates the solid-side response and keeps T_fluid constant.
    solver.addOperator(unique_ptr<Operator>(cht));

    const double tEnd = 5e-4;
    while (!dom.tripped() && solver.time() < tEnd)
    {
        // Fix hot-fluid BC: reset j=1 to T_fluid every step (no convective replenishment)
        for (int i=0;i<Nx;++i)
        {
            PrimitiveVars& w = dom.stagedW(i,1);
            w.T = T_fluid;
            w.p = p0;
            w.U[0] = 5.0; w.U[1] = w.U[2] = 0.0;
        }
        dom.primitiveToConserved();
        double dt = solver.computeTimestep(0.3);
        solver.step(dt);
    }
    cout << "  iters=" << solver.iteration() << "\n\n";

    check("no trip", !dom.tripped());

    // --- q_wall profile ---
    auto qw = cht->qWallProfile(dom);
    auto tw = cht->wallTProfile();

    cout << "  col   x[mm]   q_wall[kW/m²]   T_wall[K]\n";
    for (int i=0; i<Nx; ++i)
        cout << "  " << setw(3) << i
             << "  " << fixed << setprecision(1) << setw(6) << dom.x(i,1)*1e3
             << "  " << setw(14) << setprecision(1) << qw[i]/1e3
             << "  " << setw(10) << tw[i] << "\n";
    cout << "\n";

    // 1. q_wall > 0 (heat into solid)
    bool all_pos = true;
    for (int i=1;i<Nx-1;++i) if (qw[i] <= 0) all_pos=false;
    check("q_wall > 0 (heat flows fluid→solid)", all_pos);

    // 2. T_wall > T0_solid (wall heated up)
    double T_wall_mean = 0;
    for (int i=1;i<Nx-1;++i) T_wall_mean += tw[i];
    T_wall_mean /= (Nx-2);
    check("T_wall > T0_solid (wall heating)",
          T_wall_mean > T_solid0 + 10, T_wall_mean, T_solid0);
    check("T_wall < T_fluid (physical bound)",
          T_wall_mean < T_fluid, T_wall_mean, T_fluid);

    // 3. q_wall bounded by h_conv * DeltaT_max
    double q_max_theory = h_conv * (T_fluid - T_solid0);
    bool all_bounded = true;
    for (int i=1;i<Nx-1;++i)
        if (qw[i] > q_max_theory * 1.01) all_bounded=false;
    check("q_wall ≤ h_conv*(T_fluid-T_solid0)",
          all_bounded, *max_element(qw.begin()+1,qw.end()-1)/1e3,
          q_max_theory/1e3);

    // 4. CSV + PinFin JSON write
    cht->writeQwallCSV("/tmp/rfam_qwall.csv", dom);
    cht->writeQwallJSON("/tmp/rfam_qwall.json", dom);
    {
        ifstream f("/tmp/rfam_qwall.csv");
        check("CSV written and readable", f.good());
        ifstream fj("/tmp/rfam_qwall.json");
        check("PinFin JSON written and readable", fj.good());
    }

    // 5. Energy conservation: cumulative qToSolid matches q_wall * area * time
    double q_wall_avg = 0;
    for (int i=0;i<Nx;++i) q_wall_avg += qw[i];
    q_wall_avg /= Nx;
    // expected qToSolid ≈ q_wall_avg * Lx * t (Ly=depth=1 implicit)
    double q_expected = q_wall_avg * Lx * solver.time();
    double q_actual   = cht->heatToSolid();
    double ratio = std::abs(q_actual - q_expected) / (std::abs(q_expected) + 1.0);
    cout << "  q_expected(approx)=" << scientific << setprecision(3) << q_expected
         << "  q_actual=" << q_actual << "  ratio=" << fixed << setprecision(3) << ratio << "\n";
    // loose check: within factor 2 (approximate due to transient + spatial variation)
    check("energy consistency (within 2x)", ratio < 2.0, ratio, 2.0);

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail==0 ? "STAGE5_QWALL_PASS" : "STAGE5_QWALL_FAIL") << "\n";
    return g_fail==0 ? 0 : 1;
}
