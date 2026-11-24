#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "../variables/ConservedVars.H"
#include "../numerics/flux/AUSMPlusUpFlux.H"
#include "../numerics/preconditioning/ExperimentalCflAcousticScaler.H"
#include "../thermo/RealFluidState.H"

using namespace Foam;
using namespace std;

class Solver1DCflScalerDiagnostic
{
private:
    int N_;
    double L_, dx_;
    vector<ConservedVars> U_;
    vector<PrimitiveVars> W_;
    AUSMPlusUpFlux fluxScheme_;
    ExperimentalCflAcousticScaler cflScaler_;
    double t_;
    int iter_;

public:
    Solver1DCflScalerDiagnostic(int N, double L, bool useCflScaler = true)
        : N_(N), L_(L), t_(0.0), iter_(0)
    {
        dx_ = L_ / N_;
        U_.resize(N_);
        W_.resize(N_);
        cflScaler_.setEnabled(useCflScaler);
    }

    void initialize(double pL, double TL, double pR, double TR, double x0)
    {
        cout << "Initializing 1D CFL-scaler diagnostic..." << endl;
        cout << "  N = " << N_ << ", L = " << L_ << " m" << endl;
        cout << "  Experimental CFL acoustic scaler: "
             << (cflScaler_.isEnabled() ? "ON" : "OFF") << endl;
        cout << "  Residual/Jacobian preconditioning: NOT IMPLEMENTED" << endl;

        for (int i = 0; i < N_; i++)
        {
            double x = (i + 0.5) * dx_;
            
            if (x < x0)
            {
                W_[i].p = pL;
                W_[i].T = TL;
            }
            else
            {
                W_[i].p = pR;
                W_[i].T = TR;
            }
            
            W_[i].U[0] = 0.0;
            W_[i].U[1] = 0.0;
            W_[i].U[2] = 0.0;
            W_[i].Y[0] = 0.0;
            W_[i].Y[1] = 1.0;  // Pure C3H8
            
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
            
            // Update the explicit-timestep heuristic from local Mach.
            double M = state.getMach(W_[i].U);
            cflScaler_.updateForMach(M);
            
            double aCfl = cflScaler_.cflAcousticSpeed(state, W_[i].U);
            double u = abs(W_[i].U[0]);
            
            double dt = CFL * dx_ / (u + aCfl);
            dt_min = min(dt_min, dt);
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
            
            for (size_t j = 0; j < 3; j++)
            {
                U_new[i].rhoU[j] = U_[i].rhoU[j] - dt/dx_ * (F_right.rhoU[j] - F_left.rhoU[j]);
            }
            
            for (size_t j = 0; j < U_[i].rhoY.size(); j++)
            {
                U_new[i].rhoY[j] = U_[i].rhoY[j] - dt/dx_ * (F_right.rhoY[j] - F_left.rhoY[j]);
            }
        }
        
        U_ = U_new;
        
        for (int i = 0; i < N_; i++)
        {
            StateConversion::conservedToPrimitive(U_[i], W_[i]);
        }
        
        t_ += dt;
        iter_++;
    }

    void solve(double t_end, double CFL)
    {
        cout << "\nStarting time integration..." << endl;
        cout << "  t_end = " << t_end << " s, CFL = " << CFL << endl;
        
        while (t_ < t_end && iter_ < 10000)
        {
            double dt = computeTimestep(CFL);
            step(dt);
            
            if (iter_ % 1000 == 0)
            {
                cout << "  iter=" << iter_ << ", t=" << t_ << " s, dt=" << dt << " s" << endl;
            }
        }
        
        cout << "\nComplete: " << iter_ << " iterations" << endl;
    }
};

int main()
{
    cout << "1D Experimental CFL Acoustic-Scaler Diagnostic\n" << endl;
    cout << "This is not Turkel/Weiss all-Mach preconditioning.\n" << endl;
    
    // Test case: low-speed pressure wave
    Solver1DCflScalerDiagnostic solver(200, 1.0, true);
    solver.initialize(6.0e6, 300.0, 5.9e6, 300.0, 0.5);
    cout << "CFL_SCALER_NO_ADVANCE dt_estimate="
         << solver.computeTimestep(0.5) << endl;
    cout << "Advancing with this enlarged heuristic timestep is not "
            "authorized without a consistent residual/Jacobian transform."
         << endl;
    
    return 0;
}
