#ifndef RosenbrockSolver_H
#define RosenbrockSolver_H

/*---------------------------------------------------------------------------*\
    RosenbrockSolver - L-stable stiff ODE solver

    Rosenbrock-Euler (1st order, L-stable):
        (I - γ*dt*J) * k = f(y)
        y_{n+1} = y + dt * k

    where γ = 1 (simplest L-stable variant)

    Handles chemical source terms without blow-up at large dt.
    Numerical Jacobian via finite differences.
\*---------------------------------------------------------------------------*/

#include "ChemistryModel.H"
#include <vector>
#include <cmath>

namespace Foam
{

class RosenbrockSolver
{
private:
    static constexpr double GAMMA = 1.0;   // L-stable
    static constexpr double EPS_FD = 1e-6; // FD step size

    // Gaussian elimination with partial pivoting
    bool solveLinear(
        std::vector<std::vector<double>>& A,
        std::vector<double>& b,
        std::vector<double>& x
    ) const
    {
        int n = b.size();
        for (int col = 0; col < n; col++)
        {
            // Pivot
            int maxRow = col;
            for (int row = col+1; row < n; row++)
                if (std::abs(A[row][col]) > std::abs(A[maxRow][col]))
                    maxRow = row;
            std::swap(A[col], A[maxRow]);
            std::swap(b[col], b[maxRow]);

            if (std::abs(A[col][col]) < 1e-30) return false;

            for (int row = col+1; row < n; row++)
            {
                double f = A[row][col] / A[col][col];
                for (int j = col; j < n; j++) A[row][j] -= f * A[col][j];
                b[row] -= f * b[col];
            }
        }
        // Back substitution
        x.resize(n);
        for (int i = n-1; i >= 0; i--)
        {
            x[i] = b[i];
            for (int j = i+1; j < n; j++) x[i] -= A[i][j] * x[j];
            x[i] /= A[i][i];
        }
        return true;
    }

public:
    // One Rosenbrock step for {Y_0,...,Y_{N-1}, T}
    bool step(
        ChemistryModel* chem,
        double p, double rho,
        std::vector<double>& Y,
        double& T,
        double dt
    )
    {
        int nS = chem->nSpecies();
        int n  = nS + 1;  // state = [Y..., T]

        // Current state
        std::vector<double> y(n);
        for (int i = 0; i < nS; i++) y[i] = Y[i];
        y[nS] = T;

        // f(y): source terms
        auto f_func = [&](const std::vector<double>& yy) -> std::vector<double>
        {
            std::vector<double> Y_tmp(yy.begin(), yy.begin()+nS);
            double T_tmp = yy[nS];
            // Clamp
            for (auto& yi : Y_tmp) yi = std::max(0.0, yi);
            T_tmp = std::max(100.0, T_tmp);

            std::vector<double> omega_i;
            double qdot;
            chem->computeSourceTerms(T_tmp, p, rho, Y_tmp, omega_i, qdot);

            std::vector<double> dy(n);
            for (int i = 0; i < nS; i++)
                dy[i] = omega_i[i] / rho;
            dy[nS] = qdot / (rho * 1500.0);  // dT/dt = qdot/(rho*Cv)
            return dy;
        };

        std::vector<double> f0 = f_func(y);

        // Numerical Jacobian J_{ij} = df_i/dy_j
        std::vector<std::vector<double>> J(n, std::vector<double>(n, 0.0));
        for (int j = 0; j < n; j++)
        {
            auto y_perturb = y;
            double eps = EPS_FD * (std::abs(y[j]) + 1e-8);
            y_perturb[j] += eps;
            auto f_perturb = f_func(y_perturb);
            for (int i = 0; i < n; i++)
                J[i][j] = (f_perturb[i] - f0[i]) / eps;
        }

        // Build (I - γ*dt*J)
        std::vector<std::vector<double>> M(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                M[i][j] = -GAMMA * dt * J[i][j];
            M[i][i] += 1.0;
        }

        // Solve M * k = f0
        std::vector<double> rhs = f0;
        std::vector<double> k;
        if (!solveLinear(M, rhs, k)) return false;

        // Update
        for (int i = 0; i < nS; i++)
        {
            Y[i] += dt * k[i];
            if (Y[i] < 0.0) Y[i] = 0.0;
        }
        T += dt * k[nS];

        return true;
    }
};

} // End namespace Foam

#endif
