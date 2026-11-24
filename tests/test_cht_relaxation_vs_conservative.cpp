/*---------------------------------------------------------------------------*\
    test_cht_relaxation_vs_conservative - Stage 5 confirmation study.

    PURPOSE (user request): prove that the Stage-4 CHT2D coupling was
    NON-conservative because it only RELAXED the fluid halo temperature toward
    the wall, WITHOUT removing the absorbed heat from the fluid energy. The
    Stage-5 operator fixes this by subtracting the exact interface heat from the
    fluid rhoE (flux-conservative).

    Both methods run on the SAME quiescent hot-fluid / cold-solid setup with an
    ADIABATIC outer wall, so the closed system {fluid + solid} must conserve
    total energy: whatever the fluid loses, the solid must store. We measure

        imbalance = Δ(fluid energy) + Δ(solid stored energy)

    (per unit depth [J/m]). For a closed adiabatic system this must be ~0.

      - RELAXATION method  : fluid energy barely changes (only the thin halo T
        is overwritten and even that is not integrated into rhoE via commit),
        while the solid absorbs a large amount -> LARGE imbalance = spurious
        energy created out of nothing. THIS IS THE BUG.
      - CONSERVATIVE method : fluid loses exactly what the solid gains ->
        imbalance ~ 0 to solver accuracy.

    The relaxation method is reproduced inline here (it is no longer in the
    codebase); the conservative one is the real CHTOperator2D.
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

#include "src/thermo/SpeciesRegistry.H"
#include "src/thermo/RealFluidState.H"
#include "src/variables/ConservedVars.H"
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"
#include "src/framework/TimeIntegrator.H"
#include "src/framework/SolverBase.H"
#include "src/cht/SolidEnergyEqn.H"
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

// --- shared geometry / physics -------------------------------------------
static const int    Nx=8, Ny=6;
static const double Lx=0.01, Ly=0.006;
static const double T_fluid0=1500.0, T_solid0=300.0, p0=4e6;
static const int    Nsolid=4;
static const double Lsolid=0.002, rho_s=8900.0, cp_s=390.0, k_s=320.0, h_conv=50000.0;
static const double tEnd=5e-4;

// Total fluid energy per unit depth [J/m] = sum rhoE * dx*dy over the REAL
// fluid cells (j=1..Ny-1), EXCLUDING the wall halo row j=0. Rationale:
//   - the conservative CHT debits the interface row j=1 (a real fluid cell),
//   - the relaxation bug only overwrites the halo row j=0 (a ghost cell),
// so summing j>=1 makes the ledger see the conservative debit but NOT the
// relaxation's cosmetic halo change. With convection removed, nothing else
// crosses a cell face, so this is a genuinely closed fluid ledger.
static double fluidEnergy(const Domain& d)
{
    double E=0.0; const double cellV = d.dx()*d.dy();
    for (int i=0;i<d.Nx();++i) for (int j=1;j<d.Ny();++j)
        E += d.U(i,j).rhoE * cellV;
    return E;
}

static void initFluid(Domain& d)
{
    for (int i=0;i<Nx;++i) for (int j=0;j<Ny;++j)
    {
        PrimitiveVars& w=d.stagedW(i,j);
        w.p=p0; w.T=T_fluid0; w.U[0]=w.U[1]=w.U[2]=0.0;   // quiescent
        w.Y[0]=1.0; w.Y[1]=0.0;
    }
    d.primitiveToConserved();
}

// ============================================================
// METHOD 1: RELAXATION (the Stage-4 bug, reproduced inline)
//   Solid absorbs heat via SolidEnergyEqn; fluid halo T is set to the wall T
//   but the fluid rhoE is NOT debited. Closed-system energy is NOT conserved.
// ============================================================
static double runRelaxation()
{
    SpeciesRegistry reg; (void)reg;
    Domain dom(Nx,Ny,Lx,Ly,2);
    initFluid(dom);

    vector<vector<double>> solidT(Nx, vector<double>(Nsolid, T_solid0));
    SolidBC bcOuter = SolidBC::adiabatic();

    double Efluid0 = fluidEnergy(dom);

    // NO convection: isolate the coupling. With a quiescent field the only
    // energy change comes from the CHT coupling itself, so the ledger is clean.
    double t=0.0; double solidStored=0.0;
    const double dt = 1e-6;
    while (t<tEnd)
    {
        // relaxation coupling: solid gets fluid-T convection BC; fluid halo T
        // is overwritten to wall T, but NO rhoE debit on the fluid.
        for (int i=0;i<Nx;++i)
        {
            double Tf = dom.W(i,1).T;
            SolidBC bcTop = SolidBC::convection(h_conv, Tf);
            SolidEnergyEqn solid(Nsolid,Lsolid,rho_s,cp_s,k_s,T_solid0);
            solid.setField(solidT[i]);
            solid.step(dt, bcTop, bcOuter);
            for (int k=0;k<Nsolid;++k) solidT[i][k]=solid.T(k);
            // halo relaxation only (the bug):
            PrimitiveVars wbc = dom.W(i,0); wbc.T = solidT[i][0]; wbc.p = dom.W(i,1).p;
            dom.setCell(i,0,wbc);
        }
        t += dt;
    }

    double Efluid1 = fluidEnergy(dom);
    // solid stored energy per unit depth [J/m]
    const double cellThk = Lsolid/Nsolid, dxCol = Lx/Nx;
    for (int i=0;i<Nx;++i) for (int k=0;k<Nsolid;++k)
        solidStored += rho_s*cp_s*(solidT[i][k]-T_solid0)*cellThk*dxCol;

    double dFluid = Efluid1 - Efluid0;
    double imbalance = dFluid + solidStored;   // closed adiabatic -> should be 0
    cout << "  [relaxation] dEfluid=" << scientific << setprecision(3) << dFluid
         << "  solidStored=" << solidStored
         << "  imbalance=" << imbalance << " J/m\n";
    return imbalance;
}

// ============================================================
// METHOD 2: CONSERVATIVE (real CHTOperator2D)
// ============================================================
static double runConservative()
{
    SpeciesRegistry reg;
    Domain dom(Nx,Ny,Lx,Ly,2);
    initFluid(dom);
    SolidBC bcOuter = SolidBC::adiabatic();

    double Efluid0 = fluidEnergy(dom);

    auto* cht = new CHTOperator2D(Nsolid,Lsolid,rho_s,cp_s,k_s,h_conv,
                                  bcOuter,reg,T_solid0);
    SolverBase solver(dom);
    // NO convection: isolate the CHT coupling ledger (mirror of runRelaxation).
    solver.addOperator(unique_ptr<Operator>(cht));

    double t=0.0;
    while (t<tEnd && !dom.tripped())
    {
        double dt = solver.computeTimestep(0.3);
        solver.step(dt);
        t += dt;
    }

    double Efluid1 = fluidEnergy(dom);
    double dFluid = Efluid1 - Efluid0;
    double solidStored = cht->solidStoredDelta(Lx/Nx);
    double imbalance = dFluid + solidStored;   // should be ~0
    cout << "  [conservative] dEfluid=" << scientific << setprecision(3) << dFluid
         << "  solidStored=" << solidStored
         << "  imbalance=" << imbalance << " J/m\n";
    return imbalance;
}

int main()
{
    cout << "=== CHT coupling: relaxation (Stage-4 bug) vs conservative (Stage-5) ===\n";
    cout << "    Closed adiabatic system {fluid+solid}: total energy must be conserved.\n";
    cout << "    imbalance = d(fluid E) + d(solid stored). |imbalance| ~ 0 iff conservative.\n\n";

    double imb_relax = std::abs(runRelaxation());
    double imb_cons  = std::abs(runConservative());

    cout << "\n";
    cout << fixed << setprecision(3);
    double ratio = imb_relax / (imb_cons + 1e-30);
    cout << "  relaxation |imbalance| = " << scientific << imb_relax << " J/m\n";
    cout << "  conservative |imbalance| = " << imb_cons << " J/m\n";
    cout << "  ratio = " << scientific << ratio << " (relaxation is this many x worse)\n\n";

    // The confirmation: relaxation leaks a LOT of energy; conservative does not.
    check("relaxation method leaks energy (imbalance large)", imb_relax > 1.0, imb_relax);
    check("conservative method conserves (imbalance small)", imb_cons < 1e-3*imb_relax + 1e-6,
          imb_cons, 1e-3*imb_relax);
    check("conservative is >=100x better than relaxation", ratio > 100.0, ratio, 100.0);

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail==0 ? "CHT_RELAX_VS_CONS_PASS" : "CHT_RELAX_VS_CONS_FAIL") << "\n";
    return g_fail==0 ? 0 : 1;
}
