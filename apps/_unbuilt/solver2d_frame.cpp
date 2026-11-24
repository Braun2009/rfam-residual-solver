/*---------------------------------------------------------------------------*\
    solver2d_frame - Stage 4: 2D transcritical O2 injection via framework.

    Reproduces solver2d_stable IC using Domain(2D)+ConvectionOperator(2D).
    IC: cold O2 jet (100K, 6MPa, 10m/s) into warm coflow (300K, 0.5MPa, 1m/s).
    BCs applied via Domain::setCell before each beginStep.
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
#include "src/framework/VTKWriter.H"

using namespace Foam;
using namespace std;

// Build a PrimitiveVars from jet or coflow params (2-species: O2/C3H8)
static PrimitiveVars makePrim(double p, double T, double ux, double uy)
{
    PrimitiveVars w(2);
    w.p = p; w.T = T;
    w.U[0] = ux; w.U[1] = uy; w.U[2] = 0.0;
    w.Y[0] = 1.0; w.Y[1] = 0.0;   // pure O2
    return w;
}

// Apply 2D BCs. Symmetry FIRST, then inlet so inlet wins on j=0.
static void applyBC(Domain& d,
                    double r_jet, double p_jet, double T_jet, double u_jet,
                    double p_co,  double T_co,  double u_co)
{
    const int Nx = d.Nx(), Ny = d.Ny();

    // 1. symmetry (bottom): copy j=1, zero V
    for (int i = 0; i < Nx; ++i)
    {
        PrimitiveVars w = d.W(i, 1);
        w.U[1] = 0.0;
        d.setCell(i, 0, w);
    }
    // 2. inlet (left) — runs after symmetry, overrides j=0 at x=0
    for (int j = 0; j < Ny; ++j)
    {
        double y = d.y(0, j);
        PrimitiveVars w = (y < r_jet)
            ? makePrim(p_jet, T_jet, u_jet, 0.0)
            : makePrim(p_co,  T_co,  u_co,  0.0);
        d.setCell(0, j, w);
    }
    // 3. outlet (right): zero-gradient
    for (int j = 0; j < Ny; ++j)
        d.setCell(Nx-1, j, d.W(Nx-2, j));
    // 4. top (coflow)
    for (int i = 0; i < Nx; ++i)
        d.setCell(i, Ny-1, makePrim(p_co, T_co, u_co, 0.0));
}

int main(int argc, char* argv[])
{
    const int    Nx   = (argc > 1) ? atoi(argv[1]) : 20;
    const int    Ny   = (argc > 2) ? atoi(argv[2]) : 10;
    const double Lx   = 0.020, Ly = 0.010;
    const double tEnd = 1.0e-5;
    const double CFL  = 0.3;

    Domain dom(Nx, Ny, Lx, Ly, 2);

    const double r_jet = 0.0005;
    const double p_jet=6e6, T_jet=100.0, u_jet=10.0;
    const double p_co =0.5e6, T_co=300.0, u_co=1.0;

    // IC — put jet core in j=0 explicitly (r_jet < dy, so y(i,0)=0.5*dy > r_jet
    // for a 10-cell domain; enlarge jet to fill j=0 row to ensure a visible core)
    const double r_jet_ic = 1.5 * dom.dy();   // fill first 1-2 rows as jet
    for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
    {
        double y = dom.y(i, j);
        dom.stagedW(i,j) = (y < r_jet_ic)
            ? makePrim(p_jet, T_jet, u_jet, 0.0)
            : makePrim(p_co,  T_co,  u_co,  0.0);
    }
    dom.primitiveToConserved();

    SolverBase solver(dom);
    solver.addOperator(unique_ptr<Operator>(
        new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));

    cout << "=== solver2d_frame: 2D transcritical O2 injection ===\n";
    cout << "    " << Nx << "x" << Ny << " cells  "
         << Lx*100 << "x" << Ly*100 << " cm  t_end=" << tEnd << "s\n\n";

    double t = 0.0; int iter = 0;
    while (t < tEnd && !dom.tripped())
    {
    applyBC(dom, r_jet_ic, p_jet, T_jet, u_jet, p_co, T_co, u_co);
        double dt = solver.computeTimestep(CFL);
        if (t + dt > tEnd) dt = tEnd - t;
        if (!solver.step(dt)) break;
        t += dt; ++iter;
        if (iter % 20 == 0)
        {
            double rmin=1e30, rmax=0, Tmin=1e30, Tmax=0;
            for (int i=1;i<Nx-1;++i) for (int j=1;j<Ny-1;++j)
            {
                rmin=min(rmin, dom.U(i,j).rho); rmax=max(rmax, dom.U(i,j).rho);
                Tmin=min(Tmin, dom.W(i,j).T);   Tmax=max(Tmax, dom.W(i,j).T);
            }
            cout << "  iter=" << setw(4) << iter
                 << " t=" << scientific << setprecision(2) << t
                 << fixed << " rho=[" << setprecision(1)
                 << rmin << "," << rmax << "] T=[" << Tmin << "," << Tmax << "]\n";
        }
    }

    // Final summary
    double rmin=1e30, rmax=0, Tmin=1e30, Tmax=0, amin=1e30;
    for (int i=1;i<Nx-1;++i) for (int j=1;j<Ny-1;++j)
    {
        rmin=min(rmin,dom.U(i,j).rho); rmax=max(rmax,dom.U(i,j).rho);
        Tmin=min(Tmin,dom.W(i,j).T);   Tmax=max(Tmax,dom.W(i,j).T);
        RealFluidState s; s.updatePTY(dom.W(i,j).p, dom.W(i,j).T, dom.W(i,j).Y);
        amin=min(amin,s.a());
    }
    cout << (dom.tripped() ? "TRIP\n" : "Done\n");
    cout << "  " << iter << " iters  t=" << t << "s\n";
    cout << "  rho=[" << rmin << "," << rmax << "]  T=["
         << Tmin << "," << Tmax << "]  a_min=" << amin << "\n";

    bool ok = !dom.tripped() && rmin>0 && Tmin>0 && amin>0 && Tmax<1e6;

    // Verify inlet BC is maintaining the jet at the inlet face (i=0)
    double T_inlet_jet  = dom.W(0, 0).T;    // i=0, j=0: jet core (should be cold)
    double T_inlet_top  = dom.W(0, Ny-2).T; // i=0, j=Ny-2: coflow (should be warm)
    cout << "  inlet i=0 j=0:    T=" << T_inlet_jet  << "K (jet, expect ~" << T_jet << ")\n";
    cout << "  inlet i=0 j=Ny-2: T=" << T_inlet_top  << "K (coflow, expect ~" << T_co << ")\n";
    // Jet at inlet should be significantly colder than coflow
    bool jet_present = (T_inlet_jet < 200.0) && (T_inlet_top > 250.0);
    ok = ok && jet_present;

    cout << (ok ? "SOLVER2D_FRAME_PASS" : "SOLVER2D_FRAME_FAIL") << "\n";

    // Write final VTK output (for ParaView / post-processing)
    {
        VTKWriter vtk("solver2d_final.vtk", dom);
        vtk.writeScalar("density",     [&](int i,int j){ return dom.U(i,j).rho; });
        vtk.writeScalar("pressure",    [&](int i,int j){ return dom.W(i,j).p;   });
        vtk.writeScalar("temperature", [&](int i,int j){ return dom.W(i,j).T;   });
        vtk.writeVelocity();
        cout << "  VTK written: solver2d_final.vtk\n";
    }

    return ok ? 0 : 1;
}
