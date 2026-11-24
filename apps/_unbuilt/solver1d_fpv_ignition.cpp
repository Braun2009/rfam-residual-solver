#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "../variables/ConservedVars.H"
#include "../numerics/flux/AUSMPlusUpFlux.H"
#include "../combustion/FPV/FPVTable.H"
#include "../thermo/RealFluidState.H"

using namespace Foam;
using namespace std;

class Solver1DFPVIgnition
{
private:
    int N_;
    double L_, dx_;
    vector<ConservedVars> U_;
    vector<PrimitiveVars> W_;
    AUSMPlusUpFlux fluxScheme_;
    FPVTable fpvTable_;
    double t_;
    int iter_;
    bool tripped_;
    double D_Z_, D_C_;

public:
    Solver1DFPVIgnition(int N, double L, const string& fpvTableFile)
        : N_(N), L_(L), t_(0.0), iter_(0), tripped_(false), D_Z_(1e-5), D_C_(1e-5)
    {
        dx_ = L_ / N_;
        U_.resize(N_);
        W_.resize(N_);
        if (!fpvTable_.loadFromCSV(fpvTableFile)) { cerr << "ERROR: Failed to load FPV table" << endl; exit(1); }
    }

    void initialize(double p0, double T0, double Z0, double C0)
    {
        cout << "Initializing uniform mixture: Z=" << Z0 << ", C=" << C0 << endl;
        for (int i = 0; i < N_; i++)
        {
            W_[i].p = p0; W_[i].T = T0; W_[i].Z = Z0; W_[i].C = C0;
            W_[i].U[0] = 0.0; W_[i].U[1] = 0.0; W_[i].U[2] = 0.0;
            W_[i].Y[0] = 0.0; W_[i].Y[1] = 1.0;
            StateConversion::primitiveToConserved(W_[i], U_[i]);
        }
    }

    void setIgnitionPatch(double x_center, double radius, double C_ignite = 1.0)
    {
        cout << "Ignition at x=" << x_center << ", r=" << radius << endl;
        int count = 0;
        for (int i = 0; i < N_; i++)
        {
            double x = (i + 0.5) * dx_;
            if (abs(x - x_center) < radius)
            {
                W_[i].C = C_ignite;
                StateConversion::primitiveToConserved(W_[i], U_[i]);
                count++;
            }
        }
        cout << "  Ignited " << count << " cells" << endl;
    }

    double computeTimestep(double CFL)
    {
        double dt_min = 1e10;
        for (int i = 0; i < N_; i++)
        {
            RealFluidState state;
            state.updatePTY(W_[i].p, W_[i].T, W_[i].Y);
            double a = state.a();
            double u = abs(W_[i].U[0]);
            double dt_conv = CFL * dx_ / (u + a + 1e-10);
            double dt_diff = 0.5 * dx_ * dx_ / max(D_Z_, D_C_);
            dt_min = min(dt_min, min(dt_conv, dt_diff));
        }
        return dt_min;
    }

    void step(double dt)
    {
        auto U_new = U_;
        vector<double> normal(3, 0.0);
        normal[0] = 1.0;

        for (int i = 1; i < N_-1; i++)
        {
            ConservedVars F_left, F_right;
            fluxScheme_.computeFlux(W_[i-1], W_[i], normal, F_left);
            fluxScheme_.computeFlux(W_[i], W_[i+1], normal, F_right);

            U_new[i].rho = U_[i].rho - dt/dx_ * (F_right.rho - F_left.rho);
            U_new[i].rhoE = U_[i].rhoE - dt/dx_ * (F_right.rhoE - F_left.rhoE);
            for (size_t j = 0; j < 3; j++) U_new[i].rhoU[j] = U_[i].rhoU[j] - dt/dx_ * (F_right.rhoU[j] - F_left.rhoU[j]);
            for (size_t j = 0; j < U_[i].rhoY.size(); j++) U_new[i].rhoY[j] = U_[i].rhoY[j] - dt/dx_ * (F_right.rhoY[j] - F_left.rhoY[j]);
            U_new[i].rhoZ = U_[i].rhoZ - dt/dx_ * (F_right.rhoZ - F_left.rhoZ);
            U_new[i].rhoC = U_[i].rhoC - dt/dx_ * (F_right.rhoC - F_left.rhoC);
        }

        for (int i = 1; i < N_-1; i++)
        {
            double Z_l = U_[i-1].rhoZ / (U_[i-1].rho + 1e-10);
            double Z_c = U_[i].rhoZ / (U_[i].rho + 1e-10);
            double Z_r = U_[i+1].rhoZ / (U_[i+1].rho + 1e-10);
            double dZ_dx2 = (Z_r - 2.0*Z_c + Z_l) / (dx_*dx_);
            U_new[i].rhoZ += dt * U_[i].rho * D_Z_ * dZ_dx2;

            double C_l = U_[i-1].rhoC / (U_[i-1].rho + 1e-10);
            double C_c = U_[i].rhoC / (U_[i].rho + 1e-10);
            double C_r = U_[i+1].rhoC / (U_[i+1].rho + 1e-10);
            double dC_dx2 = (C_r - 2.0*C_c + C_l) / (dx_*dx_);

            RealFluidState state;
            state.updatePTY(W_[i].p, W_[i].T, W_[i].Y);
            double h = state.h();
            double omega_C = fpvTable_.getOmegaDot(W_[i].Z, W_[i].C, h);

            U_new[i].rhoC += dt * (U_[i].rho * D_C_ * dC_dx2 + omega_C);
            U_new[i].rhoC = max(0.0, min(U_new[i].rho * 1.0, U_new[i].rhoC));
        }

        for (int i = 0; i < N_; i++)
        {
            if (!isfinite(U_new[i].rho) || U_new[i].rho <= 0.0) { tripped_ = true; return; }
        }

        U_ = U_new;
        for (int i = 0; i < N_; i++) StateConversion::conservedToPrimitive(U_[i], W_[i]);
        t_ += dt; iter_++;
    }

    void solve(double t_end, double CFL, int freq = 100)
    {
        cout << "Starting t_end=" << t_end << ", CFL=" << CFL << endl;
        while (t_ < t_end && iter_ < 100000 && !tripped_)
        {
            step(computeTimestep(CFL));
            if (iter_ % freq == 0)
                cout << "iter=" << iter_ << ", t=" << t_ << ", C_max=" << getMaxC() << ", C_avg=" << getAvgC() << endl;
        }
        cout << (tripped_ ? "TRIPPED" : "COMPLETE") << ": " << iter_ << " iterations" << endl;
    }

    double getMaxC() const { double m = 0.0; for (int i = 0; i < N_; i++) m = max(m, W_[i].C); return m; }
    double getAvgC() const { double s = 0.0; for (int i = 0; i < N_; i++) s += W_[i].C; return s / N_; }
};

int main()
{
    cout << "1D FPV Ignition Test\n" << endl;
    Solver1DFPVIgnition solver(200, 0.1, "fpv_table_proof.csv");
    solver.initialize(1.0e6, 300.0, 0.273, 0.0);
    solver.setIgnitionPatch(0.05, 0.005, 1.0);
    solver.solve(0.001, 0.3, 100);
    return 0;
}
