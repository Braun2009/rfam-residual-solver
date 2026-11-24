/*---------------------------------------------------------------------------*\
    solver_mayer_o2 - Mayer-like transcritical O2 jet (pressure-matched).

    Approximates Mayer (1998) coaxial jet experiment in O2 (vs LN2 original):
        p_jet = p_ch  = 4.0 MPa   (pressure-matched — no shock cells)
        T_jet         = 130 K     (Vortex LOX inlet design point: 130 K /
                                   6.0 MPa pump outlet, Pc = 4.0 MPa;
                                   SUBCRITICAL compressed liquid —
                                   Tc=154.6K, Pc=5.04MPa, so the jet is a
                                   dense subcritical liquid, NOT supercritical)
        T_ch          = 300 K     (warm O2 environment)
        D_jet         = 2 mm      (r_jet = 1 mm)
        u_jet         = 5 m/s
        Domain        = 40mm × 10mm (axisymmetric half-plane, y=0 is axis)

    PR-EOS conditions (from check):
        rho_jet ~ 1091 kg/m3    (dense supercritical liquid-like)
        rho_ch  ~   53 kg/m3    (ideal-gas-like)
        density ratio ~ 20.6

    Post-processing output (for Mayer comparison):
        rho_centerline(x)   — axial density decay
        half-width y_05(x)  — jet spreading rate
        potential core L_c  — where rho_cl / rho_jet < 0.95

    Reference: Mayer et al. (1998) "Atomization and breakup of cryogenic
    propellants under high-pressure subcritical and supercritical conditions."
    J. Propulsion Power 14(5).
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

#include "src/variables/ConservedVars.H"
#include "src/thermo/RealFluidState.H"
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"
#include "src/framework/TimeIntegrator.H"
#include "src/framework/SolverBase.H"
#include "src/framework/NonReflectingOutlet.H"
#include "src/framework/VTKWriter.H"
#include "src/framework/DiffusionOperator2D.H"
#include "src/thermo/GERGResidualTable.H"
#include "src/thermo/GERGSaturationTable.H"

using namespace Foam;
using namespace std;

static GERGResidualTable g_mayerGerg;
static GERGSaturationTable g_mayerSat;

static void registerGERG(RealFluidState& st)
{
    st.setGERGTable(&g_mayerGerg);
    st.setSaturationTable(&g_mayerSat);
}

static PrimitiveVars makePrim2sp(double p, double T, double ux)
{
    PrimitiveVars w(2);
    w.p = p; w.T = T;
    w.U[0] = ux; w.U[1] = w.U[2] = 0.0;
    w.Y[0] = 1.0; w.Y[1] = 0.0;  // pure O2
    return w;
}

static void applyBC(Domain& d,
    double r_jet, double p_jet, double T_jet, double u_jet,
    double p_co,  double T_co,  double u_co)
{
    const int Nx = d.Nx(), Ny = d.Ny();

    // 1. axis symmetry first
    for (int i = 0; i < Nx; ++i)
    {
        PrimitiveVars w = d.W(i, 1);
        w.U[1] = 0.0;
        d.setCell(i, 0, w);
    }
    // 2. inlet (left) — jet / coflow split
    for (int j = 0; j < Ny; ++j)
    {
        double y = d.y(0, j);
        PrimitiveVars w = (y < r_jet)
            ? makePrim2sp(p_jet, T_jet, u_jet)
            : makePrim2sp(p_co,  T_co,  u_co);
        d.setCell(0, j, w);
    }
    // 3. top: far-field coflow
    for (int i = 0; i < Nx; ++i)
        d.setCell(i, Ny-1, makePrim2sp(p_co, T_co, u_co));
}

// ---- Post-processing helpers ----

// Axial density profile at j=0..1 (axis / near-axis)
static vector<pair<double,double>> centerlineDensity(const Domain& d, int Nx)
{
    vector<pair<double,double>> out;
    for (int i = 1; i < Nx-1; ++i)
    {
        double x = d.x(i, 1);      // near-axis row j=1
        double rho = d.U(i, 1).rho;
        out.push_back({x, rho});
    }
    return out;
}

// Jet half-width: radial distance where rho drops to 50% of centerline value
static double halfWidth(const Domain& d, int i, int Ny, double rho_cl)
{
    const double half_target = 0.5 * rho_cl;
    for (int j = 1; j < Ny-1; ++j)
    {
        if (d.U(i, j).rho < half_target)
        {
            // linear interpolation
            double rho_jm = d.U(i, j-1).rho;
            double rho_j  = d.U(i, j).rho;
            double frac = (rho_jm - half_target) / (rho_jm - rho_j + 1e-30);
            double y_jm = d.y(i, j-1);
            double y_j  = d.y(i, j);
            return y_jm + frac * (y_j - y_jm);
        }
    }
    return d.Ly();  // not found — return domain edge
}

int main(int argc, char* argv[])
{
    const int    Nx   = (argc > 1) ? atoi(argv[1]) : 60;
    const int    Ny   = (argc > 2) ? atoi(argv[2]) : 30;
    // Narrow domain: Ly=3mm=1.5D, so 10 cells per radius at Ny=30
    const double Lx   = 0.030, Ly = 0.003;
    const double D    = 0.002;    // jet diameter [m]
    const double r_jet= 0.5*D;

    // Mayer-like conditions (pressure matched)
    // u_jet = 100 m/s (Ma~0.18, subsonic): flow-through time Lx/u = 0.4ms
    // Run to t=1ms (~2.5 transit times) to reach quasi-steady state.
    const double p_jet = 4.0e6, T_jet = 130.0, u_jet = 100.0;
    const double p_co  = 4.0e6, T_co  = 300.0, u_co  = 2.0;

    // Load the GERG-2008 residual table (dense-cryogenic backend).
    {
        const string gergPath = string(RFAM_SOURCE_DIR)
            + "/tests/reference/gerg_residual_table_v3.csv";
        if (!g_mayerGerg.load(gergPath))
        {
            cerr << "FATAL: GERG residual table missing: " << gergPath << "\n";
            return 2;
        }
        cout << "  GERG-2008 residual table loaded (v3)\n";
        const string satPath = string(RFAM_SOURCE_DIR)
            + "/tests/reference/gerg_saturation_table_v4.csv";
        if (!g_mayerSat.load(satPath))
        {
            cerr << "FATAL: saturation table missing: " << satPath << "\n";
            return 2;
        }
        cout << "  GERG saturation table loaded (170 rows)\n";
    }

    // Estimate rho_jet from GERG (PR shown for contrast)
    {
        RealFluidState sj; vector<double> Y={1,0};
        registerGERG(sj);
        sj.updatePTY(p_jet, T_jet, Y);
        cout << "=== Mayer-like O2 jet (GERG-2008) ===\n";
        cout << "  rho_jet  = " << fixed << setprecision(1) << sj.rho() << " kg/m3\n";
        RealFluidState sc; registerGERG(sc); sc.updatePTY(p_co, T_co, Y);
        cout << "  rho_ch   = " << sc.rho() << " kg/m3\n";
        cout << "  rho_ratio= " << sj.rho()/sc.rho() << "\n";
        cout << "  a_jet    = " << sj.a() << " m/s   u_jet=" << u_jet << " m/s  Ma_jet=" << u_jet/sj.a() << "\n";
        cout << "  Domain   = " << Lx*1000 << " x " << Ly*1000 << " mm  (" << Lx/D << " x " << Ly/D << " D)\n";
        cout << "  Grid     = " << Nx << " x " << Ny << "\n\n";
    }

    Domain dom(Nx, Ny, Lx, Ly, 2);
    dom.setSpecies({"O2", "C3H8"});
    dom.thermoModel().setGERGTable(&g_mayerGerg);
    dom.thermoModel().setSaturationTable(&g_mayerSat);
    const double r_ic = 1.5 * dom.dy();  // ensure j=0 gets jet IC

    for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
    {
        double y = dom.y(i, j);
        dom.stagedW(i,j) = (y < max(r_jet, r_ic))
            ? makePrim2sp(p_jet, T_jet, u_jet)
            : makePrim2sp(p_co,  T_co,  u_co);
    }
    dom.primitiveToConserved();

    NonReflectingOutlet nro(p_co, Lx);

    SolverBase solver(dom);
    solver.addOperator(unique_ptr<Operator>(
        new ConvectionOperator(unique_ptr<FluxScheme>(new AUSMFluxAdapter()))));
    // Molecular + SGS-like mixing: DiffusionOperator2D with the legacy
    // gas transport plus an explicit effective kinematic viscosity
    // (turbulent mixing placeholder, DIAGNOSTIC-ONLY, named not silent).
    {
        const double nuEff = (argc > 4) ? atof(argv[4]) : 1.0e-5;
        SpeciesRegistry reg;
        solver.addOperator(unique_ptr<Operator>(
            new DiffusionOperator2D(reg, 1.0, nuEff, 0.3)));
        cout << "  Diffusion: legacy gas transport + nu_eff="
             << scientific << setprecision(3) << nuEff
             << " m2/s (SGS placeholder)\n";
    }

    const double CFL = 0.3;
    const double tEnd = (argc > 3) ? atof(argv[3]) : 1.0e-3;

    cout << "Running to t=" << scientific << setprecision(3)
         << tEnd << " s ...\n";
    double t = 0.0; int iter = 0;
    while (t < tEnd && !dom.tripped())
    {
        applyBC(dom, r_jet, p_jet, T_jet, u_jet, p_co, T_co, u_co);
        nro.apply(dom, solver.computeTimestep(CFL));
        double dt = solver.computeTimestep(CFL);
        if (t + dt > tEnd) dt = tEnd - t;
        if (!solver.step(dt)) break;
        t += dt; ++iter;
        if (iter % 200 == 0)
            cout << "  iter=" << setw(5) << iter
                 << " t=" << scientific << setprecision(2) << t << "s\n";
    }
    cout << "Done: " << iter << " iters  t=" << t << "\n\n";

    if (dom.tripped())
    {
        cout << "TRIP - case diverged. code="
             << Domain::failureCodeName(dom.failure().code)
             << " cell=" << dom.failure().cell
             << " detail=" << dom.failure().detail
             << " backend=" << dom.failure().thermoDiagnostics.backendStatus
             << " eqstatus=" << static_cast<int>(
                 dom.failure().equilibriumStatus) << "\n";
        {
            // diagnostic: the conserved state of the failed cell
            const int ci = dom.failure().cell;
            const auto& C = dom.Unew(ci);
            const double rho = C.rho;
            const double rhoE = C.rhoE;
            const double e = rhoE / rho;
            const auto& UY = C.rhoY;
            double yO2 = -1.0;
            if (UY.size() >= 1) yO2 = UY[0] / rho;
            // Nx/Ny are the main() grid dimensions (in scope)
            cout << "  TRIP state: rho=" << rho << " e=" << e
                 << " Y[O2]=" << yO2
                 << " cell=" << ci << " (i=" << ci % Nx
                 << " j=" << (ci / Nx) % Ny << ")"
                 << " x=" << dom.x(ci % Nx) << " y=" << dom.y(0, (ci / Nx) % Ny)
                 << "\n";
        }
        return 1;
    }

    // ---- Post-processing: extract Mayer comparison profiles ----
    // (2026-08-04 audit: rho_jet must use the SAME GERG backend as the
    // solve; previously it fell back to PR and normalised against a
    // different density.)
    RealFluidState sj; { vector<double> Y={1,0}; registerGERG(sj); sj.updatePTY(p_jet,T_jet,Y); }
    const double rho_jet = sj.rho();

    // Write VTK
    {
        VTKWriter vtk("mayer_o2_final.vtk", dom);
        vtk.writeScalar("density",     [&](int i,int j){ return dom.U(i,j).rho; });
        vtk.writeScalar("temperature", [&](int i,int j){ return dom.W(i,j).T;   });
        vtk.writeScalar("pressure",    [&](int i,int j){ return dom.W(i,j).p;   });
        vtk.writeVelocity();
        cout << "VTK written: mayer_o2_final.vtk\n";
    }

    // Write centerline density profile
    cout << "\n=== Centerline density (axial profile) ===\n";
    cout << "  x[mm]   x/D    rho[kg/m3]   rho/rho_jet\n";

    double L_core = -1.0;
    ofstream fout("mayer_centerline.dat");
    fout << "# x[m]  x/D  rho[kg/m3]  rho/rho_jet\n";

    for (int i = 1; i < Nx-1; ++i)
    {
        double x = dom.x(i, 1);
        double rho = dom.U(i, 1).rho;
        double xD = x / D;
        double rnorm = rho / rho_jet;

        // Potential core: first x where normalized density < 0.95
        if (L_core < 0 && rnorm < 0.95) L_core = x;

        fout << x << " " << xD << " " << rho << " " << rnorm << "\n";
        if (i % 5 == 0)
            cout << "  " << fixed << setprecision(1) << x*1e3
                 << "   " << setprecision(2) << xD
                 << "   " << setprecision(1) << rho
                 << "       " << setprecision(3) << rnorm << "\n";
    }
    cout << "\n  Potential core length: "
         << (L_core > 0 ? to_string(L_core/D) + " D" : ">domain") << "\n";

    // Write spreading rate
    cout << "\n=== Jet spreading (half-width profile) ===\n";
    cout << "  x[mm]   x/D    y_05[mm]   y_05/D\n";
    ofstream fspr("mayer_spreading.dat");
    fspr << "# x[m]  x/D  y_half[m]  y_half/D\n";

    for (int i = 2; i < Nx-2; i += 3)
    {
        double x = dom.x(i, 1);
        double rho_cl = dom.U(i, 1).rho;
        double y05 = halfWidth(dom, i, Ny, rho_cl);
        double xD = x / D, y05D = y05 / D;

        fspr << x << " " << xD << " " << y05 << " " << y05D << "\n";
        cout << "  " << fixed << setprecision(1) << x*1e3
             << "   " << setprecision(2) << xD
             << "   " << setprecision(2) << y05*1e3
             << "     " << setprecision(3) << y05D << "\n";
    }

    RealFluidState sc;
    {
        vector<double> Y={1,0};
        registerGERG(sc);
        sc.updatePTY(p_co,T_co,Y);
    }
    cout << "\nMayer (1998) reference values (ρ_jet/ρ_ch ~ 4-10):\n";
    cout << "  L_core/D   ~ 8-12\n";
    cout << "  dδ₀.₅/dx  ~ 0.08-0.11\n";
    cout << "  our ρ_jet/ρ_ch = " << fixed << setprecision(1)
         << rho_jet/sc.rho() << "\n";

    cout << "\nSOLVER_MAYER_DONE\n";
    return 0;
}
