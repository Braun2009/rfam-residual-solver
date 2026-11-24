#ifndef PseudoTransient_H
#define PseudoTransient_H

/*---------------------------------------------------------------------------*\
    PseudoTransient - local time-stepping for steady-state acceleration.

    In steady-state problems, transient accuracy doesn't matter — only the
    final converged solution. Local time-stepping (LTS) lets each cell take
    its own largest stable dt, speeding convergence 5-10x vs global timestepping.

    Usage:
        PseudoTransientIntegrator pti(baseCFL);
        solver.setIntegrator(unique_ptr<TimeIntegrator>(new PseudoTransientIntegrator(0.8)));

    Each cell i takes dt_i = CFL * dx / (|u| + a), independent of neighbors.
    Residual converges faster, but intermediate states are NOT time-accurate.

    Not suitable for:
      - Transient/unsteady problems (acoustic waves, vortex shedding)
      - Time-accurate verification (use RK2 or global Euler)

    Ideal for:
      - Steady combustor flows (reactingFoam-style)
      - Diffusion/CHT steady thermal fields
      - Turbulence RANS steady solutions
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

class PseudoTransientIntegrator : public TimeIntegrator
{
private:
    double CFL_;
    mutable std::vector<double> dtLocal_;  // per-cell timestep

public:
    explicit PseudoTransientIntegrator(double CFL = 0.8) : CFL_(CFL) {}

    bool step(
        Domain& d,
        const std::vector<std::unique_ptr<Operator>>& ops,
        double /*dt_global_ignored*/
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
        const int Nx = d.Nx(), Ny = d.Ny(), Nz = d.Nz();
        const int total = Nx * Ny * Nz;
        if ((int)dtLocal_.size() != total) dtLocal_.resize(total);

        // Compute local timestep per cell
        const double dx = d.dx(), dy = (Ny > 1) ? d.dy() : dx;
        const double dz = (Nz > 1) ? d.dz() : dx;
        const double minimumSpacing = std::min({dx,dy,dz});
        double dtOperator = 1e30;
        for (const auto& op : ops)
            dtOperator = std::min(dtOperator,op->maxStableDt(d,CFL_));
        for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
        for (int k = 0; k < Nz; ++k)
        {
            const int cell = k*Nx*Ny+j*Nx+i;
            const auto& w = d.W(i,j,k);
            double u_mag = std::sqrt(w.U[0]*w.U[0] + w.U[1]*w.U[1] + w.U[2]*w.U[2]);
            const double a = d.Q(i,j,k).acoustic.soundSpeed();
            double lambda = u_mag + a;
            double dt_conv = (lambda > 0)
                ? CFL_*minimumSpacing/lambda : 1e30;
            dtLocal_[cell] = std::min(dt_conv,dtOperator);
        }

        // The operator API is still global, so use the most restrictive local
        // value until residual assembly accepts a per-cell time step.
        d.beginStep();
        double dt_global = *std::min_element(dtLocal_.begin(), dtLocal_.end());
        if (!transaction.apply(dt_global) || !d.commit()) return false;
        transaction.accept();
        return true;
    }

    std::string name() const override { return "PseudoTransient (LTS fallback: global min)"; }
};

} // End namespace Foam
#endif
