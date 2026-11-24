/*---------------------------------------------------------------------------*\
    test_2d_reacting - Stage 4: 2D reacting flow with chemistry + diffusion.

    Uniform premixed C3H8/air box (1200K, 4MPa) on a coarse 10x5 2D grid.
    Operators stacked: Convection + ChemistrySource + Diffusion.
    Verifies that the 3 operators compose correctly in 2D (no framework change
    needed — each operator already dispatches on d.Ny()>1 or is cell-local).

    Checks:
      1. No trip (NaN/negative rho)
      2. rho > 0 everywhere
      3. a > 0 everywhere
      4. Fuel C3H8 consumed (species coupling works in 2D)
      5. Carbon element conserved (< 1e-4)
      6. Chemical heat added (rhoE increased)
      7. Every committed state remains inside the transport caloric range
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include "src/thermo/SpeciesRegistry.H"
#include "src/thermo/NASA7Thermo.H"
#include "src/transport/TransportProps.H"
#include "src/variables/ConservedVars.H"
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"
#include "src/framework/DiffusionOperator.H"
#include "src/framework/TimeIntegrator.H"
#include "src/framework/SolverBase.H"
#include "src/combustion/ChemistrySourceOperator.H"
#include "src/combustion/FiniteRate/TwoStepCO.H"

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
    cout << "=== 2D reacting flow: Convection + Chemistry + Diffusion ===\n\n";

    SpeciesRegistry reg;
    const int nS = reg.size();   // 9 species (global registry)
    const int Nx = 10, Ny = 5;
    const double Lx = 0.01, Ly = 0.005;  // 1 cm x 0.5 cm
    const int nSteps = 32;
    const double CFL = 0.3;

    Domain dom(Nx, Ny, Lx, Ly, nS);
    dom.setSpecies(reg.names());

    // Uniform stoichiometric propane-air IC. N2 dilution keeps the reacting
    // state inside the active-species NASA-7 range required by gas transport.
    const double p0=4e6, T0=1200.0;
    vector<double> Y0(nS, 0.0);
    Y0[reg.index("C3H8")] = 0.0603436613169685;
    Y0[reg.index("O2")]   = 0.218947454389961;
    Y0[reg.index("N2")]   = 0.72070888429307;

    for (int i=0;i<Nx;++i) for (int j=0;j<Ny;++j)
    {
        PrimitiveVars& w = dom.stagedW(i,j);
        w.p=p0; w.T=T0;
        w.U[0]=w.U[1]=w.U[2]=0.0;
        w.Y=Y0;
    }
    dom.primitiveToConserved();

    // Record initial element and energy
    auto boxAtoms = [&](double& C, double& H, double& O)
    {
        C=H=O=0.0;
        for (int i=1;i<Nx-1;++i) for (int j=1;j<Ny-1;++j)
        {
            double rr=dom.U(i,j).rho;
            for (int g=0;g<nS;++g)
            {
                double mr=rr*dom.W(i,j).Y[g]/reg.molecularWeight(g);
                C+=mr*reg.nC(g); H+=mr*reg.nH(g); O+=mr*reg.nO(g);
            }
        }
    };
    double C0,H0,O0; boxAtoms(C0,H0,O0);
    double rhoE0=0.0;
    for (int i=1;i<Nx-1;++i) for (int j=1;j<Ny-1;++j)
        rhoE0+=dom.U(i,j).rhoE;

    // Stack operators: Convection (assign interior) -> Chemistry (+= rhoY/rhoE) -> Diffusion (+= grad terms)
    TwoStepCO chem;
    SolverBase solver(dom);
    solver.addOperator(unique_ptr<Operator>(
        new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));
    solver.addOperator(unique_ptr<Operator>(
        new ChemistrySourceOperator(&chem, reg)));
    solver.addOperator(unique_ptr<Operator>(
        new DiffusionOperator(reg)));

    // A bounded number of coupled steps is sufficient for this composition
    // smoke test. Marching toward a distant physical end time stalls on the
    // explicit chemistry limiter as reversible R2 approaches equilibrium.
    for (int step=0;step<nSteps && !dom.tripped();++step)
    {
        const double dt=solver.computeTimestep(CFL);
        if (!isfinite(dt) || dt<=0.0 || !solver.step(dt)) break;
    }

    // Final measurements
    double C1,H1,O1; boxAtoms(C1,H1,O1);
    double rhoE1=0.0;
    for (int i=1;i<Nx-1;++i) for (int j=1;j<Ny-1;++j)
        rhoE1+=dom.U(i,j).rhoE;

    double thr=C0+H0+O0+1e-30;
    const PrimitiveVars& wmid=dom.W(Nx/2, Ny/2);

    cout << "  iters=" << solver.iteration() << "\n";
    cout << "  mid: Y_C3H8=" << wmid.Y[reg.index("C3H8")]
         << " T=" << wmid.T << "K\n\n";

    check("no trip", !dom.tripped());

    bool rho_ok=true, a_ok=true;
    for (int i=1;i<Nx-1;++i) for (int j=1;j<Ny-1;++j)
    {
        if (dom.U(i,j).rho<=0.0) rho_ok=false;
        RealFluidState s; s.updatePTY(dom.W(i,j).p, dom.W(i,j).T, dom.W(i,j).Y);
        if (s.a()<=0.0) a_ok=false;
    }
    check("rho > 0 all interior cells", rho_ok);
    check("a > 0 all interior cells", a_ok);

    bool transport_ok=true, temperature_range_ok=true;
    const vector<string>& names=dom.speciesNames();
    for (int i=0;i<dom.N();++i)
    {
        const PrimitiveVars& w=dom.W(i);
        const TransportProps::Result transport=
            TransportProps::evaluateCommittedLegacyGas(
                names,w,dom.Q(i),dom.thermoModel(),reg);
        if (!transport.success()) transport_ok=false;

        for (int g=0;g<nS;++g)
        {
            if (w.Y[g]<=1e-14) continue;
            double minimum=0.0, maximum=0.0;
            if (!nasa7Thermo().temperatureRange(
                    names[g],minimum,maximum)
                || w.T<minimum || w.T>maximum)
                temperature_range_ok=false;
        }
    }
    check("committed direct-Q transport succeeds in all cells", transport_ok);
    check("all active species remain in NASA-7 temperature range",
          temperature_range_ok);
    check("C3H8 consumed (Y < 0.01)", wmid.Y[reg.index("C3H8")] < 0.01,
          wmid.Y[reg.index("C3H8")]);
    check("element C conserved (< 1e-4 rel)", fabs(C1-C0)/thr < 1e-4);
    check("element H conserved (< 1e-4 rel)", fabs(H1-H0)/thr < 1e-4);
    check("element O conserved (< 1e-4 rel)", fabs(O1-O0)/thr < 1e-4);
    check("chemical heat added (rhoE up)", rhoE1 > rhoE0);

    cout << "\n--------------------------------------\n";
    cout << "  RESULT: " << g_pass << " passed, " << g_fail << " failed\n";
    cout << (g_fail==0 ? "REACTING2D_PASS" : "REACTING2D_FAIL") << "\n";
    return g_fail==0 ? 0 : 1;
}
