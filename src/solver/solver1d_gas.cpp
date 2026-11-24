/*---------------------------------------------------------------------------*\
  realFluidAllMach - Simple 1D Euler Solver

  Purpose: Verify RealFluidState integration with conservative solver
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "ConservedVars.H"
#include "AUSMPlusUpFlux.H"
#include "RealFluidState.H"

using namespace Foam;
using namespace std;

class Solver1D
{
private:
    int N_;                              // Number of cells
    double L_;                           // Domain length
    double dx_;                          // Cell size
    vector<ConservedVars> U_;            // Conserved variables
    vector<PrimitiveVars> W_;            // Primitive variables
    AUSMPlusUpFlux fluxScheme_;
    bool tripped_ = false;

public:
    Solver1D(int N, double L) : N_(N), L_(L)
    {
        dx_ = L_ / N_;
        U_.resize(N_);
        W_.resize(N_);
    }

    void initialize(double pL, double TL, double pR, double TR, double x0)
    {
        cout << "Initializing 1D domain..." << endl;
        cout << "  N = " << N_ << " cells" << endl;
        cout << "  L = " << L_ << " m" << endl;
        cout << "  dx = " << dx_ << " m" << endl;

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

        cout << "  Left state:  p=" << pL/1e6 << " MPa, T=" << TL << " K" << endl;
        cout << "  Right state: p=" << pR/1e6 << " MPa, T=" << TR << " K" << endl;
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
            double dt = CFL * dx_ / (u + a);

            dt_min = min(dt_min, dt);
        }

        return dt_min;
    }

    void step(double dt)
    {
        vector<ConservedVars> U_new = U_;
        vector<double> normal(3, 0.0);
        normal[0] = 1.0;

        // Compute fluxes and update
        for (int i = 1; i < N_-1; i++)
        {
            // Left interface flux
            ConservedVars F_left;
            fluxScheme_.computeFlux(W_[i-1], W_[i], normal, F_left);

            // Right interface flux
            ConservedVars F_right;
            fluxScheme_.computeFlux(W_[i], W_[i+1], normal, F_right);

            // Update (Forward Euler)
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

        // --- tripwire: report first non-finite / non-positive cell ---
        for (int i = 0; i < N_; i++)
        {
            bool bad = !std::isfinite(U_new[i].rho) || U_new[i].rho <= 0.0
                    || !std::isfinite(U_new[i].rhoE);
            if (bad)
            {
                std::cout << "  [TRIP] cell " << i
                          << " rho=" << U_new[i].rho
                          << " rhoUx=" << U_new[i].rhoU[0]
                          << " rhoE=" << U_new[i].rhoE << std::endl;
                std::cout << "         neighbors rho: "
                          << (i>0?U_[i-1].rho:0) << " | " << U_[i].rho
                          << " | " << (i<N_-1?U_[i+1].rho:0) << std::endl;
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
    }

    void solve(double t_end, double CFL)
    {
        cout << endl;
        cout << "Starting time integration..." << endl;
        cout << "  t_end = " << t_end << " s" << endl;
        cout << "  CFL = " << CFL << endl;

        double t = 0.0;
        int iter = 0;

        while (t < t_end)
        {
            double dt = computeTimestep(CFL);
            if (t + dt > t_end) dt = t_end - t;

            step(dt);
            if (tripped_) { std::cout << "  STOPPED at iter " << iter << ", t=" << t << std::endl; break; }

            t += dt;
            iter++;

            if (iter % 100 == 0)
            {
                cout << "  iter=" << iter << ", t=" << t << " s, dt=" << dt << " s" << endl;
            }
        }

        cout << "Completed: " << iter << " iterations" << endl;
    }

    void output(const string& filename)
    {
        ofstream file(filename);
        file << "# x rho p T u Cv" << endl;

        for (int i = 0; i < N_; i++)
        {
            double x = (i + 0.5) * dx_;

            RealFluidState state;
            state.updatePTY(W_[i].p, W_[i].T, W_[i].Y);

            file << x << " "
                 << state.rho() << " "
                 << W_[i].p << " "
                 << W_[i].T << " "
                 << W_[i].U[0] << " "
                 << state.Cv() << endl;
        }

        file.close();
        cout << "Output written to " << filename << endl;
    }
};

int main()
{
    cout << "========================================" << endl;
    cout << "  realFluidAllMach 1D Solver Test" << endl;
    cout << "========================================" << endl;
    cout << endl;

    // Test case: Pressure wave with real fluid
    // Left: high pressure, Right: low pressure

    Solver1D solver(200, 1.0);

    double pL = 2.0e5;
    double TL = 450.0;
    double pR = 2.0e4;
    double TR = 450.0;

    solver.initialize(pL, TL, pR, TR, 0.5);
    solver.solve(0.0015, 0.4);
    solver.output("solution_gas.dat");

    cout << endl;
    cout << "========================================" << endl;
    cout << "  Test completed successfully" << endl;
    cout << "========================================" << endl;

    return 0;
}
