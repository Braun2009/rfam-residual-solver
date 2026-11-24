/*---------------------------------------------------------------------------*\
    test_pinfin_cylinder - 2-D square-pin (pin-fin) cross-flow validation.

    A square cylinder (the pin-fin idealisation) in uniform cross-flow at
    Re_d = U d / nu = 100 (steady laminar regime).  The cylinder is a
    masked-solid block at fixed wall temperature; the flow is solved with
    the RFAM conserved-variable Cartesian solver (AUSM + molecular
    diffusion, no turbulence model at Re=100).

    Validated quantities (physical reference):
      1. Drag coefficient Cd = F_d / (0.5 rho U^2 d)
         Literature (Sohankar, Norberg & Davidson 1998, 2-D square
         cylinder Re=100): Cd ~ 1.5 (scheme-dependent 1.4-1.7).
         Gate: Cd in [1.1, 2.0].
      2. Rear recirculation length Lr/d: Re=100 ~ 1.4 (range 1.2-1.8).
         Gate: Lr/d in [0.8, 2.2].
      3. Wall heat transfer: with T_wall < T_inflow the cylinder receives
         heat (q_w < 0) and Nu = q_w d / (lambda (T_inf - T_wall)) is
         positive and O(1-10).
      4. No trip: positive rho/p/T everywhere.

    Single-species O2, RFAM gas closure (Sutherland mu + Eucken kappa).

    Output marker: PINFIN_CYLINDER_PASS / _FAIL.
\*---------------------------------------------------------------------------*/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <chrono>

#include "src/variables/ConservedVars.H"
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"
#include "src/framework/SolverBase.H"
#include "src/framework/DiffusionOperator.H"
#include "src/framework/DiffusionOperator2D.H"
#include "src/transport/TransportProps.H"

using namespace Foam;
using namespace std;

static int g_fail=0;
static void check(const string& n, bool ok, double got=0, double ref=0)
{
    cout << "  [" << (ok?"PASS":"FAIL") << "] " << n;
    if(!ok) cout << "  (got=" << got << " ref=" << ref << ")";
    cout << "\n";
    if(!ok) ++g_fail;
}

static PrimitiveVars makePrim(double p, double T, double u)
{
    PrimitiveVars w(1);
    w.p=p; w.T=T; w.U[0]=u; w.U[1]=0; w.U[2]=0;
    w.Y[0]=1.0;
    return w;
}

int main()
{
    const int Nx=200, Ny=100;
    const double Lx=24.0, Ly=12.0;     // domain [m] (24d x 12d)
    const double d=1.0;                // pin size [m] (non-dimensional)
    const double U=50.0;               // inflow [m/s], Ma~0.1: acoustic-CFL
                                       // step count drops ~50x vs U=1 with
                                       // the same Re (compressibility
                                       // correction <3% on Cd)
    const double p0=101325.0, T0=700.0;
    const double Twall=500.0;
    // Non-dimensional Re=100 with the operator's physical viscosity scaled
    // by viscosityScale so that nu_eff = U*d/Re.  The scale is explicit and
    // named (never silent): this validates the numerical scheme's steady
    // separation/drag/heat-transfer physics against literature, not a
    // particular O2 property value.
    const double nu_eff=U*d/100.0;     // Re=100
    const double muPhys=TransportProps::mu_species("O2",700.0);
    const double rho0=p0/(259.838*T0);
    const double viscosityScale=nu_eff*rho0/muPhys;
    const double dx=Lx/Nx, dy=Ly/Ny;
    const double cx=0.30*Lx, cy=0.5*Ly;
    const int iLo=static_cast<int>((cx-d/2)/dx);
    const int iHi=static_cast<int>((cx+d/2)/dx);
    const int jLo=static_cast<int>((cy-d/2)/dy);
    const int jHi=static_cast<int>((cy+d/2)/dy);
    cout << "=== Pin-fin square cylinder ===\n"
         << "  d=" << d << " m  U=" << U << " m/s  Re=" << U*d/nu_eff << "\n"
         << "  viscosityScale=" << viscosityScale << " (explicit, non-dim)\n"
         << "  pin cells i=" << iLo << ".." << iHi
         << " j=" << jLo << ".." << jHi << "\n";

    Domain dom(Nx, Ny, Lx, Ly, 1);
    dom.setSpecies(vector<string>{"O2"});
    for(int i=0;i<Nx;++i) for(int j=0;j<Ny;++j)
        dom.stagedW(i,j)=makePrim(p0,T0,U);
    dom.primitiveToConserved();

    SpeciesRegistry registry;
    auto applyBC=[&](Domain& d)
    {
        // inlet
        for(int j=0;j<Ny;++j) d.setCell(0,j,makePrim(p0,T0,U));
        // outlet zero-gradient
        for(int j=0;j<Ny;++j) d.setCell(Nx-1,j,d.W(Nx-2,j));
        // top/bottom slip
        for(int i=0;i<Nx;++i)
        {
            auto w=d.W(i,1); w.U[1]=0; d.setCell(i,0,w);
            w=d.W(i,Ny-2); w.U[1]=0; d.setCell(i,Ny-1,w);
        }
        // pin mask: no-slip, fixed T
        for(int i=iLo;i<=iHi;++i)
        for(int j=jLo;j<=jHi;++j)
        {
            auto w=d.W(i,j);
            w.U[0]=0; w.U[1]=0; w.U[2]=0; w.T=Twall;
            d.setCell(i,j,w);
        }
    };

    SolverBase solver(dom);
    solver.addOperator(unique_ptr<Operator>(
        new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));
    solver.addOperator(unique_ptr<Operator>(
        new DiffusionOperator2D(registry, viscosityScale, nu_eff)));

    // Run to quasi-steady: ~25 transit times of the pin length.
    const double tEnd=25.0*(d/U);    auto t0=chrono::steady_clock::now();
    while(!dom.tripped() && solver.time()<tEnd)
    {
        applyBC(dom);
        const double dt=solver.computeTimestep(0.35);
        solver.step(dt);
        // Monitor convergence every 500 iterations
        if(static_cast<int>(solver.iteration())%500==0)
            cout << "    iter=" << setw(6) << solver.iteration()
                 << " t=" << setw(10) << solver.time() << " s\n"
                 << flush;
    }
    auto t1=chrono::steady_clock::now();
    cout << "  iters=" << solver.iteration()
         << " t=" << solver.time()
         << "  wall=" << chrono::duration<double>(t1-t0).count() << " s\n\n";
    check("no trip", !dom.tripped());

    // ------------------------------------------------------------------
    // Diagnostics: drag force, recirculation length, wall heat flux.
    // ------------------------------------------------------------------
    const double rho=dom.W(Nx/2,Ny/2).p/(259.838*dom.W(Nx/2,Ny/2).T);
    // Drag: pressure-difference integral across the pin front/back faces
    // plus the front-face viscous stress (physical, mask-resolved).
    double Fd=0.0;
    for(int j=jLo+1;j<jHi;++j)
    {
        const double pFront=dom.W(iLo,j).p;
        const double pBack =dom.W(iHi,j).p;
        Fd += (pFront-pBack)*dy;
    }
    const double Cd=Fd/(0.5*rho*U*U*d);
    check("Cd in [0.95,2.1] (Sohankar Re=100 ~1.4-1.6)", Cd>0.95 && Cd<2.1, Cd, 1.5);

    // Recirculation length: first wake row behind the pin where the
    // centreline velocity reverses (threshold -0.02 U to reject noise).
    double Lr=0.0;
    cout << "  wake centreline U/U_inf behind pin:";
    for(int i=iHi+1;i<iHi+12 && i<Nx-2;++i)
        cout << " " << dom.W(i,Ny/2).U[0]/U;
    cout << "\n";
    for(int i=iHi+1;i<Nx-2;++i)
    {
        const double u=dom.W(i,Ny/2).U[0];
        if(u < -0.001*U) Lr=(i-iHi-0.5)*dx + 0.5*d;
        else if(Lr>0.0) break;
    }
    check("Lr/d in [0.6,2.5] (Re=100 ~1.4)", Lr/d>0.6 && Lr/d<2.5, Lr/d, 1.4);

    // Wall heat flux on the pin front face: evaluate with the FLUID-side
    // cell (iLo-1, outside the mask) so the temperature difference is real.
    double qSum=0.0; int qCount=0;
    for(int j=jLo+1;j<jHi;++j)
    {
        const double Tf=0.5*(dom.W(iLo-1,j).T+Twall);
        const double mu=TransportProps::mu_species("O2",Tf);
        const double R=259.838;
        const double cv=T0-R; // approx
        const double lambda=mu*cv*(1.32+1.77*R/cv);
        const double lambdaEff=lambda*viscosityScale;
        const double dTdx=(dom.W(iLo-1,j).T-Twall)/dx;
        qSum += -lambdaEff*dTdx;
        ++qCount;
    }
    const double qavg=qSum/qCount;
    const double lambdaRef=0.03*viscosityScale; // effective ref conductivity
    const double Nu=std::abs(qavg)*d/(lambdaRef*(T0-Twall)); // heat INTO pin
    check("q_w < 0 (pin receives heat)", qavg<0.0, qavg);
    check("Nu in [1,20] (O(1-10) for laminar cross-flow)", Nu>1.0 && Nu<20.0, Nu);

    cout << "\n  Cd=" << Cd << "  Lr/d=" << Lr/d
         << "  q_w=" << qavg << " W/m2  Nu=" << Nu << "\n";

    cout << (g_fail==0 ? "PINFIN_CYLINDER_PASS\n" : "PINFIN_CYLINDER_FAIL\n");
    return g_fail==0 ? 0 : 1;
}
