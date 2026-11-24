/*---------------------------------------------------------------------------*\
    test_cht2d - Stage 4: 2D conjugate heat transfer smoke test.

    Hot fluid (1500K) flowing over a solid wall (initially 300K).
    Checks:
      1. Wall heats up (T_wall > T0_solid after time advance)
      2. Fluid-solid interface temperature is between fluid and initial solid T
      3. No trip in fluid domain
      4. Heat flows from hot fluid into cold solid (wall > 300K, < 1500K)
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <cmath>
#include <memory>
#include <vector>
#include <string>

#include "SpeciesRegistry.H"
#include "ConservedVars.H"
#include "Domain.H"
#include "FluxScheme.H"
#include "Operator.H"
#include "TimeIntegrator.H"
#include "SolverBase.H"
#include "cht/CHTOperator2D.H"

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
    cout << "=== 2D CHT: fluid-solid coupling smoke test ===\n\n";

    SpeciesRegistry reg;
    const int nS = 2;   // legacy 2-species (pure O2)
    const int Nx=8, Ny=6;
    const double Lx=0.01, Ly=0.006;
    const double T_fluid=1500.0, T_solid0=300.0;
    const double p0=4e6;

    Domain dom(Nx, Ny, Lx, Ly, nS);

    // Hot uniform fluid IC
    for (int i=0;i<Nx;++i) for (int j=0;j<Ny;++j)
    {
        PrimitiveVars& w=dom.stagedW(i,j);
        w.p=p0; w.T=T_fluid; w.U[0]=5.0; w.U[1]=w.U[2]=0.0;
        w.Y[0]=1.0; w.Y[1]=0.0;
    }
    dom.primitiveToConserved();

    // CHT operator: CuCrZr-like properties
    // rho=8900, cp=390, k=320 W/mK, h=50000 W/m2K (high-pressure convection)
    SolidBC bcOuter = SolidBC::fixedT(300.0);   // cold outer wall
    SolverBase solver(dom);
    solver.addOperator(unique_ptr<Operator>(
        new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));
    auto cht=std::make_unique<CHTOperator2D>(
        4,0.002,8900.0,390.0,320.0,50000.0,bcOuter,reg,T_solid0);
    const CHTOperator2D* chtView=cht.get();
    solver.addOperator(std::move(cht));

    const double tEnd=1e-4;
    solver.advanceToTime(tEnd, 0.3);

    cout << "  iters=" << solver.iteration() << "\n";

    const double T_wall_mean=chtView->meanWallT();
    const double heat_to_solid=chtView->heatToSolid();

    cout << "  mean solid interface T = " << T_wall_mean << "K\n";
    cout << "  fluid mid T = " << dom.W(Nx/2,Ny/2).T << "K\n\n";

    check("no trip", !dom.tripped());
    check("solid interface temperature is finite", std::isfinite(T_wall_mean),
          T_wall_mean);
    check("wall heated above initial solid T", T_wall_mean > T_solid0,
          T_wall_mean, T_solid0);
    check("wall remains below the hot fluid", T_wall_mean < T_fluid,
          T_wall_mean, T_fluid);
    check("heat flows from fluid into solid", heat_to_solid > 0.0,
          heat_to_solid);
    check("fluid interior T > 0", dom.W(Nx/2,Ny/2).T > 0.0,
          dom.W(Nx/2,Ny/2).T);

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail==0 ? "CHT2D_PASS" : "CHT2D_FAIL") << "\n";
    return g_fail==0 ? 0 : 1;
}
