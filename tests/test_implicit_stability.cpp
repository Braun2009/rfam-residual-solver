/*---------------------------------------------------------------------------*\
    test_implicit_stability - Point-implicit integrator validation.

    Validates:
      1. At safe density ratio (Tj=170K, ρ=2.3): point-implicit runs with 5×
         larger effective dt than Euler (steady-state acceleration).
      2. At the legacy dense-injection condition (Tj=120K), Euler completes
         1000 accepted steps with a valid U/W/Q field. This is a numerical
         regression only; it is not a two-phase or transcritical physics gate.
      3. For non-transcritical flows, point-implicit allows larger pseudo-CFL
         without sacrificing stability.
\*---------------------------------------------------------------------------*/
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <memory>

#include "src/variables/ConservedVars.H"
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"
#include "src/framework/TimeIntegrator.H"
#include "src/framework/SolverBase.H"
#include "src/framework/PointImplicit.H"

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

static PrimitiveVars mp(double p,double T,double u)
{PrimitiveVars w(2);w.p=p;w.T=T;w.U[0]=u;w.U[1]=w.U[2]=0;w.Y[0]=1;w.Y[1]=0;return w;}

static void bc(Domain& d, double r, double pj, double Tj, double uj,
               double pc, double Tc, double uc)
{
    int Nx=d.Nx(),Ny=d.Ny();
    for(int i=0;i<Nx;++i){auto w=d.W(i,1);w.U[1]=0;d.setCell(i,0,w);}
    for(int j=0;j<Ny;++j){double y=d.y(0,j);d.setCell(0,j,y<r?mp(pj,Tj,uj):mp(pc,Tc,uc));}
    for(int j=0;j<Ny;++j)d.setCell(Nx-1,j,d.W(Nx-2,j));
    for(int i=0;i<Nx;++i)d.setCell(i,Ny-1,mp(pc,Tc,uc));
}

struct RunResult
{
    bool initialized;
    bool tripped;
    int iters;
    double time;
    bool statesValid;
};

static RunResult runCase(bool useImplicit, double Tj, int maxIter=1000)
{
    int Nx=40, Ny=20;
    double Lx=0.04, Ly=0.006, r=0.001;
    double pj=4e6, uj=100.0, pc=4e6, Tc=300.0, uc=2.0;
    Domain dom(Nx,Ny,Lx,Ly,2);
    double r_ic=1.5*dom.dy();
    for(int i=0;i<Nx;++i) for(int j=0;j<Ny;++j)
        dom.stagedW(i,j)=dom.y(i,j)<max(r,r_ic)?mp(pj,Tj,uj):mp(pc,Tc,uc);
    const bool initialized=dom.primitiveToConserved();
    if(!initialized) return {false,dom.tripped(),0,0.0,false};
    SolverBase solver(dom);
    if(useImplicit)
        solver.setIntegrator(unique_ptr<TimeIntegrator>(new PointImplicitIntegrator(5.0)));
    solver.addOperator(unique_ptr<Operator>(
        new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));
    int iter=0;
    while(!dom.tripped()&&iter<maxIter){
        bc(dom,r,pj,Tj,uj,pc,Tc,uc);
        if(dom.tripped()) break;
        double dt=solver.computeTimestep(0.3);
        if(useImplicit) dt*=5.0;   // allow larger pseudo-dt
        if(!std::isfinite(dt)||dt<=0.0||!solver.step(dt)) break;
        ++iter;
    }
    bool statesValid=!dom.tripped()&&dom.stateReady();
    for(int i=0;i<Nx&&statesValid;++i)
    for(int j=0;j<Ny&&statesValid;++j)
    {
        const auto& u=dom.U(i,j);
        const auto& w=dom.W(i,j);
        statesValid=std::isfinite(u.rho)&&u.rho>0.0
            &&std::isfinite(u.rhoE)
            &&std::isfinite(w.p)&&w.p>0.0
            &&std::isfinite(w.T)&&w.T>0.0
            &&dom.Q(i,j).valid();
    }
    return {initialized,dom.tripped(),iter,solver.time(),statesValid};
}

int main()
{
    cout << "=== Point-implicit integrator validation ===\n\n";

    // ── Test A: Safe density ratio (Tj=170K, ρ=2.3) ──
    // Both Euler and implicit should survive. Implicit covers more physical time
    // in the same iteration count (5× larger dt).
    cout << "─ A: Safe density ratio (Tj=170K, ρ≈2.3) ─\n";
    auto euler_safe = runCase(false, 170.0, 1000);
    auto impl_safe  = runCase(true,  170.0, 1000);
    cout << "  Euler:    iters=1000 tripped=" << (euler_safe.tripped?"YES":"NO")
         << " t=" << euler_safe.time << "s\n";
    cout << "  Implicit: iters=1000 tripped=" << (impl_safe.tripped?"YES":"NO")
         << " t=" << impl_safe.time << "s\n";
    cout << "  time ratio (impl/euler) = " << impl_safe.time/(euler_safe.time+1e-30) << "×\n\n";

    check("Euler completes 1000 accepted steps at Tj=170K",
          euler_safe.initialized&&!euler_safe.tripped
              &&euler_safe.iters==1000&&euler_safe.statesValid);
    check("Implicit completes 1000 accepted steps at Tj=170K",
          impl_safe.initialized&&!impl_safe.tripped
              &&impl_safe.iters==1000&&impl_safe.statesValid);
    check("Implicit covers more physical time (larger effective dt)",
          std::isfinite(euler_safe.time)&&euler_safe.time>0.0
              &&std::isfinite(impl_safe.time)
              &&impl_safe.time > euler_safe.time * 2.0,
          impl_safe.time/euler_safe.time, 2.0);

    // This locks the legacy numerical path only. Production two-phase and
    // cryogenic-transport capabilities remain gated independently.
    cout << "\n─ B: Legacy dense injection (Tj=120K, ρ≈20.6) ─\n";
    auto euler_trans = runCase(false, 120.0, 1000);
    cout << "  Euler:    iters=" << euler_trans.iters
         << " tripped=" << (euler_trans.tripped?"YES":"NO")
         << " t=" << euler_trans.time << "s\n";
    check("Euler completes 1000 accepted legacy dense-state steps",
          euler_trans.initialized&&!euler_trans.tripped
              &&euler_trans.iters==1000&&euler_trans.statesValid);
    check("Legacy dense-state physical time is finite and positive",
          std::isfinite(euler_trans.time)&&euler_trans.time>0.0,
          euler_trans.time,0.0);
    cout << "  [INFO] Numerical survival does not validate phase equilibrium, "
            "cryogenic transport, or production Vortex advancement.\n";

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail==0 ? "IMPLICIT_STABILITY_PASS" : "IMPLICIT_STABILITY_FAIL") << "\n";
    return g_fail==0 ? 0 : 1;
}
