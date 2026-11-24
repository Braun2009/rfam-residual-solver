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

class Solver1DFPV
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
    
    // Diffusion coefficients
    double D_Z_;   // Mixture fraction diffusivity
    double D_C_;   // Progress variable diffusivity

public:
    Solver1DFPV(int N, double L, const string& fpvTableFile)
        : N_(N), L_(L), t_(0.0), iter_(0), tripped_(false),
          D_Z_(1e-5), D_C_(1e-5)
    {
        dx_ = L_ / N_;
        U_.resize(N_);
        W_.resize(N_);
        
        // Load FPV table
        if (!fpvTable_.loadFromCSV(fpvTableFile))
        {
            cerr << "ERROR: Failed to load FPV table" << endl;
            exit(1);
        }
    }

    void initialize(double pL, double TL, double ZL, double CL,
                    double pR, double TR, double ZR, double CR, double x0)
    {
        cout << "Initializing 1D FPV domain..." << endl;
        cout << "  N = " << N_ << ", L = " << L_ << " m" << endl;

        for (int i = 0; i < N_; i++)
        {
            double x = (i + 0.5) * dx_;
            
            if (x < x0)
            {
                W_[i].p = pL;
                W_[i].T = TL;
                W_[i].Z = ZL;
                W_[i].C = CL;
            }
            else
            {
                W_[i].p = pR;
                W_[i].T = TR;
                W_[i].Z = ZR;
                W_[i].C = CR;
            }
            
            W_[i].U[0] = 0.0;
            W_[i].U[1] = 0.0;
            W_[i].U[2] = 0.0;
            W_[i].Y[0] = 0.0;
            W_[i].Y[1] = 1.0;  // N2
            
            StateConversion::primitiveToConserved(W_[i], U_[i]);
        }
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
            double dt_conv = CFL * dx_ / (u + a);
            
            // Diffusion CFL
            double dt_diff = 0.5 * dx_ * dx_ / max(D_Z_, D_C_);
            
            double dt = min(dt_conv, dt_diff);
            dt_min = min(dt_min, dt);
        }
        
        return dt_min;
    }

    void step(double dt)
    {
        auto U_new = U_;
        vector<double> normal(3, 0.0);
        normal[0] = 1.0;
        
        // Convective fluxes
        for (int i = 1; i < N_-1; i++)
        {
            ConservedVars F_left, F_right;
            fluxScheme_.computeFlux(W_[i-1], W_[i], normal, F_left);
            fluxScheme_.computeFlux(W_[i], W_[i+1], normal, F_right);
            
            // Conservative update
            U_new[i].rho = U_[i].rho - dt/dx_ * (F_right.rho - F_left.rho);
            U_new[i].rhoE = U_[i].rhoE - dt/dx_ * (F_right.rhoE - F_left.rhoE);
            
            for (size_t j = 0; j < 3; j++)
            {
                U_new[i].rhoU[j] = U_[i].rhoU[j] - dt/dx_ * (F_right.rhoU[j] - F_left.rhoU[j]);
            }
            
            for (size_t j = 0; j < U_[i].rhoY.size(); j++)
            {
                U_new[i].rhoY[j] = U_[i].rhoY[j] - dt/dx_ * (F_right.rhoY[j] - F_left.rhoY[j]);
            }
            
            // FPV convection
            U_new[i].rhoZ = U_[i].rhoZ - dt/dx_ * (F_right.rhoZ - F_left.rhoZ);
            U_new[i].rhoC = U_[i].rhoC - dt/dx_ * (F_right.rhoC - F_left.rhoC);
        }
        
        // Diffusion + Source terms for Z and C
        for (int i = 1; i < N_-1; i++)
        {
            // Z diffusion (passive scalar, no source)
            double Z_left = U_[i-1].rhoZ / U_[i-1].rho;
            double Z_center = U_[i].rhoZ / U_[i].rho;
            double Z_right = U_[i+1].rhoZ / U_[i+1].rho;
            
            double dZ_dx2 = (Z_right - 2.0*Z_center + Z_left) / (dx_*dx_);
            U_new[i].rhoZ += dt * U_[i].rho * D_Z_ * dZ_dx2;
            
            // C diffusion + source
            double C_left = U_[i-1].rhoC / U_[i-1].rho;
            double C_center = U_[i].rhoC / U_[i].rho;
            double C_right = U_[i+1].rhoC / U_[i+1].rho;
            
            double dC_dx2 = (C_right - 2.0*C_center + C_left) / (dx_*dx_);
            
            // Compute enthalpy for FPV lookup
            RealFluidState state;
            state.updatePTY(W_[i].p, W_[i].T, W_[i].Y);
            double h = state.h();
            
            // FPV source term
            double omega_C = fpvTable_.getOmegaDot(W_[i].Z, W_[i].C, h);
            
            U_new[i].rhoC += dt * (U_[i].rho * D_C_ * dC_dx2 + omega_C);
        }
        
        // Check for NaN/negative
        for (int i = 0; i < N_; i++)
        {
            bool bad = !std::isfinite(U_new[i].rho) || U_new[i].rho <= 0.0
                    || !std::isfinite(U_new[i].rhoE);
            if (bad)
            {
                cout << "  [TRIP] cell " << i << " rho=" << U_new[i].rho << endl;
                tripped_ = true;
                return;
            }
        }
        
        U_ = U_new;
        
        // Update primitive variables
        for (int i = 0; i < N_; i++)
        {
            StateConversion::conservedToPrimitive(U_[i], W_[i]);
        }
        
        t_ += dt;
        iter_++;
    }

    void solve(double t_end, double CFL)
    {
        cout << "\nStarting FPV time integration..." << endl;
        cout << "  t_end = " << t_end << " s, CFL = " << CFL << endl;
        
        while (t_ < t_end && iter_ < 50000 && !tripped_)
        {
            double dt = computeTimestep(CFL);
            step(dt);
            
            if (iter_ % 1000 == 0)
            {
                cout << "  iter=" << iter_ << ", t=" << t_ 
                     << ", C_max=" << getMaxC() << endl;
            }
        }
        
        if (tripped_)
        {
            cout << "\n❌ Simulation tripped!" << endl;
        }
        else
        {
            cout << "\n✅ Complete: " << iter_ << " iterations" << endl;
        }
    }
    
    double getMaxC() const
    {
        double maxC = 0.0;
        for (int i = 0; i < N_; i++)
        {
            maxC = max(maxC, W_[i].C);
        }
        return maxC;
    }
};

int main()
{
    cout << "1D FPV Combustion Solver\n" << endl;
    
    // Test: mixture fraction wave with ignition
    Solver1DFPV solver(200, 1.0, "fpv_table_proof.csv");
    
    // Initialize: left (fuel-rich), right (oxidizer)
    solver.initialize(
        6.0e6, 300.0, 0.95, 0.0,  // Left: Z=0.95 (fuel), C=0 (unburnt)
        6.0e6, 300.0, 0.05, 0.0,  // Right: Z=0.05 (oxidizer), C=0
        0.5                         // Interface at x=0.5
    );
    
    solver.solve(0.001, 0.5);
    
    return 0;
}
