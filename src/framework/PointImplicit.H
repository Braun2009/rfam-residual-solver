#ifndef PointImplicit_H
#define PointImplicit_H

/*---------------------------------------------------------------------------*\
    PointImplicitIntegrator - local-Jacobian implicit time advancement.

    Standard in production CFD (STAR-CCM+, Fluent, SU2). Each cell's update is
    damped by the local spectral radius of the flux Jacobian:

        ΔU = dt · R(U) / (1 + dt · |∂F/∂U|_diag)

    where |∂F/∂U|_diag ≈ (|u|+a)/dx + (|v|+a)/dy (spectral radius per cell).

    This is equivalent to a backward-Euler linearisation with a diagonal
    approximation of the implicit matrix (no global matrix inversion needed).

    WHY THIS WORKS for transcritical jets:
    - Explicit Euler crashes at density ratio >5 because the CFL-limited dt in
      the light gas allows the dense liquid jet to over-update in one step.
    - Point-implicit damps the dense-jet cells MORE (they have smaller spectral
      radius → larger denominator), naturally stabilising the interface.
    - The dt input becomes a CFL-like PSEUDO time step (can be much larger than
      the explicit stability limit) because each cell's actual update is damped.

    Workflow:
      1. beginStep()
      2. Apply all operators (explicit residuals into Unew)
      3. For each cell: Unew = U + (Unew - U) / (1 + dt · λ_i)
         where λ_i is the local spectral radius [1/s]
      4. commit()

    Pure std C++, header-only. No linear algebra library needed.
\*---------------------------------------------------------------------------*/

#include "TimeIntegrator.H"
#include "AcousticStateAdmission.H"
#include "Domain.H"
#include "Operator.H"
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>

namespace Foam
{

class PointImplicitIntegrator : public TimeIntegrator
{
private:
    double CFL_implicit_;   // effective CFL (can be >> 1 for steady-state)

public:
    explicit PointImplicitIntegrator(double CFL_impl = 5.0) : CFL_implicit_(CFL_impl) {}

    bool step(
        Domain& d,
        const std::vector<std::unique_ptr<Operator>>& ops,
        double dt_global   // used as reference; actual per-cell dt is implicit
    ) const override
    {
        for (int cell=0; cell<d.N(); ++cell)
            if (!solverConsumerAdmitsAcousticState(d.Q(cell)))
            {
                d.recordFailure(Domain::FailureCode::STATE_RECOVERY_FAILED,
                    cell,name(),
                    "solver acoustic state is not admitted");
                return false;
            }

        detail::StepTransaction transaction(d,ops);
        if (!transaction.active()) return false;

        for (const auto& op : ops)
            if (op->ownsPersistentState())
            {
                d.recordFailure(Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                    name(),"point-implicit damping is not defined for operator-owned state: "
                        + op->name());
                return false;
            }

        const int Nx = d.Nx(), Ny = d.Ny(), Nz = d.Nz();
        const double dx = d.dx();
        const double dy = (Ny > 1) ? d.dy() : dx;
        const double dz = (Nz > 1) ? d.dz() : dx;

        // Step 1-2: explicit residual computation (same as Euler)
        d.beginStep();   // Unew <- U
        if (!transaction.apply(dt_global)) return false;
        // Now Unew[k] = U[k] + dt_global * R[k]  (from the operators)
        // The residual increment is: dU[k] = Unew[k] - U[k]

        // Step 3: point-implicit damping
        for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
        for (int k = 0; k < Nz; ++k)
        {
            // Get local state for spectral radius estimate
            const PrimitiveVars& w = d.W(i,j,k);
            double ux = std::abs(w.U[0]);
            double uy = std::abs(w.U[1]);
            double uz = std::abs(w.U[2]);

            const thermo::EquilibriumState& q = d.Q(i,j,k);
            const double a = q.acoustic.soundSpeed();

            // Local spectral radius [1/s]
            double lambda = (ux + a) / dx;
            if (Ny > 1) lambda += (uy + a) / dy;
            if (Nz > 1) lambda += (uz + a) / dz;

            // Implicit damping factor
            double imp_factor = 1.0 / (1.0 + dt_global * lambda);

            // Apply damping: Unew = U + imp_factor * (Unew - U)
            ConservedVars& un = d.Unew(i,j,k);
            const ConservedVars& u = d.U(i,j,k);

            un.rho  = u.rho  + imp_factor * (un.rho  - u.rho);
            un.rhoE = u.rhoE + imp_factor * (un.rhoE - u.rhoE);
            for (size_t k = 0; k < un.rhoU.size(); ++k)
                un.rhoU[k] = u.rhoU[k] + imp_factor * (un.rhoU[k] - u.rhoU[k]);
            for (size_t k = 0; k < un.rhoY.size(); ++k)
                un.rhoY[k] = u.rhoY[k] + imp_factor * (un.rhoY[k] - u.rhoY[k]);
            un.rhoZ = u.rhoZ + imp_factor * (un.rhoZ - u.rhoZ);
            un.rhoC = u.rhoC + imp_factor * (un.rhoC - u.rhoC);
        }

        // Step 4: commit with NaN check + primitive refresh
        if (!d.commit()) return false;
        transaction.accept();
        return true;
    }

    std::string name() const override { return "PointImplicit (diagonal Jacobian)"; }
};

} // End namespace Foam
#endif
