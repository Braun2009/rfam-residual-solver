/*---------------------------------------------------------------------------*\
    test_cht2d_conservation - Stage 5: CHT2D interface ENERGY CONSERVATION.

    The Stage-4 CHT2D operator only relaxed the fluid halo temperature toward
    the wall; it did NOT remove the absorbed heat from the fluid energy, so the
    fluid-solid interface did not conserve energy. Stage 5 makes CHT2D a proper
    flux-conservative source operator:
        - heat q_int enters the solid at the interface (SolidEnergyEqn),
        - the SAME q_int is removed from the fluid interface cell's rhoE,
        - bookkeeping tracks cumulative heat in/out for a closed balance.

    Two conservation checks (both to tight tolerance):

      A. ADIABATIC OUTER WALL  -> all heat that leaves the fluid is stored in
         the solid:  solidStored == heatToSolid,  heatOutOuter == 0.

      B. CLOSED BALANCE (fixedT outer) -> solid stored-energy change equals
         (heat in at interface) - (heat out at outer face):
            solidStored == heatToSolid - heatOutOuter.

    These are the multi-region analogue of the 1D CHTCoupler machine-precision
    interface balance, now expressed as an explicit source on the fluid energy.
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
#include "src/framework/TimeIntegrator.H"
#include "src/framework/SolverBase.H"
#include "src/cht/CHTOperator2D.H"

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

// Run a hot-fluid / cold-solid case and return the CHT operator's balance.
// outerBC: adiabatic (case A) or fixedT (case B).
struct Balance { double stored, in, out, dxCol; };

static Balance runCase(const SolidBC& outerBC, double tEnd)
{
    SpeciesRegistry reg;
    const int nS = 2;                       // legacy 2-species (pure O2)
    const int Nx=8, Ny=6;
    const double Lx=0.01, Ly=0.006;
    const double T_fluid=1500.0, T_solid0=300.0, p0=4e6;

    Domain dom(Nx, Ny, Lx, Ly, nS);
    for (int i=0;i<Nx;++i) for (int j=0;j<Ny;++j)
    {
        PrimitiveVars& w=dom.stagedW(i,j);
        w.p=p0; w.T=T_fluid; w.U[0]=0.0; w.U[1]=w.U[2]=0.0;  // quiescent: isolate CHT
        w.Y[0]=1.0; w.Y[1]=0.0;
    }
    dom.primitiveToConserved();

    // CuCrZr-like: rho=8900, cp=390, k=320; strong interface HTC
    auto* cht = new CHTOperator2D(4, 0.002, 8900.0, 390.0, 320.0, 50000.0,
                                  outerBC, reg, T_solid0);

    SolverBase solver(dom);
    // NO convection operator: quiescent fluid, so the ONLY energy change in the
    // fluid interface cell is the CHT sink -> a clean conservation ledger.
    solver.addOperator(unique_ptr<Operator>(cht));
    solver.advanceToTime(tEnd, 0.3);

    return { cht->solidStoredDelta(Lx/Nx), cht->heatToSolid(),
             cht->heatOutOuter(), Lx/Nx };
}

int main()
{
    cout << "=== 2D CHT interface energy conservation (Stage 5) ===\n\n";

    // --- Case A: adiabatic outer wall -> stored == heat in, nothing out ---
    cout << "-- A. Adiabatic outer wall: solidStored == heatToSolid --\n";
    {
        Balance b = runCase(SolidBC::adiabatic(), 5e-4);
        double relErr = std::abs(b.stored - b.in) / (std::abs(b.in) + 1e-30);
        cout << "  heatToSolid = " << scientific << setprecision(4) << b.in
             << "  solidStored = " << b.stored
             << "  heatOutOuter = " << b.out << " J/m\n";
        check("A: heat entered solid (in > 0)", b.in > 0.0, b.in);
        check("A: outer wall adiabatic (out == 0)", std::abs(b.out) < 1e-9, b.out);
        check("A: stored == heat in (rel < 1e-6)", relErr < 1e-6, relErr, 0.0);
    }

    // --- Case B: fixedT outer wall -> closed balance stored == in - out ---
    cout << "\n-- B. FixedT outer wall: solidStored == heatToSolid - heatOutOuter --\n";
    {
        Balance b = runCase(SolidBC::fixedT(300.0), 2e-3);
        double residual = b.stored - (b.in - b.out);
        double relErr = std::abs(residual) / (std::abs(b.in) + 1e-30);
        cout << "  heatToSolid = " << scientific << setprecision(4) << b.in
             << "  heatOutOuter = " << b.out
             << "  solidStored = " << b.stored << " J/m\n";
        cout << "  residual = " << residual << " J/m  (rel " << relErr << ")\n";
        check("B: heat entered solid (in > 0)", b.in > 0.0, b.in);
        check("B: heat left outer wall (out > 0)", b.out > 0.0, b.out);
        check("B: closed balance stored==in-out (rel < 1e-6)", relErr < 1e-6, relErr, 0.0);
    }

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail==0 ? "CHT2D_CONSERVATION_PASS" : "CHT2D_CONSERVATION_FAIL") << "\n";
    return g_fail==0 ? 0 : 1;
}
