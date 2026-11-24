#ifndef DiffusionOperator_H
#define DiffusionOperator_H

/*---------------------------------------------------------------------------*\
    DiffusionOperator - species diffusion + thermal conduction (Stage 3).

    Implements the diffusive fluxes as a framework Operator (additive update):

      rhoY_i contribution:  d/dx (rho D_i dY_i/dx)  [Fick species diffusion]
      rhoE   contribution:  d/dx (kappa dT/dx)       [Fourier heat conduction]

    SCHEME: 2nd-order central differences on interior cells (i=1..N-2).
    Boundary cells (i=0, i=N-1) are left at their current Unew values (zero-
    gradient / adiabatic BC is the default; a BC abstraction layer is future
    work). This is consistent with ConvectionOperator's interior-only scope.

    STABILITY: explicit diffusion stability limit
      dt_diff <= CFL_diff * dx^2 / (2 * max(D, alpha_T))
    where alpha_T = kappa / (rho * Cv) is thermal diffusivity and Cv comes
    from the same real-fluid state surface used by the energy equation.
    CFL_diff = 0.3
    is conservative for explicit first-order time integration.

    PROPERTIES: TransportProps (Wilke mu, Eucken kappa, Sc=0.7 diffusivity).
    See [[technical-debt]]: per-species D_i and Soret/Dufour cross effects
    are deferred.

    Pure std C++, no SUNDIALS. Depends on: Operator.H, TransportProps.H,
    RealFluidState.H, SpeciesRegistry.H.
\*---------------------------------------------------------------------------*/

#include "Operator.H"
#include "TransportProps.H"
#include "SpeciesRegistry.H"
#include <memory>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace Foam
{

class DiffusionOperator : public Operator
{
private:
    const SpeciesRegistry& reg_;
    // diffusive CFL safety factor (explicit stability)
    double cflDiff_;

    std::vector<TransportProps::Point> transportField(
        const Domain& d) const
    {
        const std::vector<std::string>& names=d.speciesNames();
        std::vector<TransportProps::Point> points(
            static_cast<size_t>(d.N()));
        for(int i=0;i<d.N();++i)
        {
            const TransportProps::Result result=
                TransportProps::evaluateCommittedLegacyGas(
                    names,d.W(i),d.Q(i),d.thermoModel(),reg_);
            if(!result.success())
                throw std::runtime_error(
                    "DiffusionOperator: transport cell "
                    +std::to_string(i)+" status="
                    +TransportProps::statusName(result.status));
            points[static_cast<size_t>(i)]=result.point;
        }
        return points;
    }

public:
    explicit DiffusionOperator(const SpeciesRegistry& reg, double cflDiff = 0.3)
        : reg_(reg), cflDiff_(cflDiff) {}

    // ------------------------------------------------------------------
    void apply(Domain& d, double dt) const override
    {
        const int N   = d.N();
        const int nS  = d.nSpecies();
        const double dx = d.dx();

        const std::vector<TransportProps::Point> points=transportField(d);

        // Interior cells only (same scope as ConvectionOperator)
        for (int i = 1; i < N-1; ++i)
        {
            const PrimitiveVars& wL = d.W(i-1);
            const PrimitiveVars& w  = d.W(i);
            const PrimitiveVars& wR = d.W(i+1);

            const double rho=d.Q(i).rho;
            const TransportProps::Point& transport=
                points[static_cast<size_t>(i)];
            const double kap = transport.kappa;
            const double D = transport.diffusivity;

            // ---- species diffusion: d/dx(rho D dY/dx) ----
            // face fluxes (1st-order upwind in diffusivity, 2nd-order in Y)
            ConservedVars& un = d.Unew(i);

            for (int s = 0; s < nS; ++s)
            {
                // d/dx(rho D dY/dx) ~ (rho D / dx^2) * (Y_{i+1} - 2Y_i + Y_{i-1})
                double lapY = (wR.Y[s] - 2.0*w.Y[s] + wL.Y[s]) / (dx*dx);
                un.rhoY[s] += rho * D * lapY * dt;
            }

            // ---- thermal conduction: d/dx(kappa dT/dx) ----
            double lapT = (wR.T - 2.0*w.T + wL.T) / (dx*dx);
            un.rhoE += kap * lapT * dt;
        }
    }

    // ------------------------------------------------------------------
    double maxStableDt(const Domain& d, double /*CFL*/) const override
    {
        const std::vector<TransportProps::Point> points=transportField(d);

        const double dx   = d.dx();
        const double dx2  = dx * dx;
        double dtMin = 1e30;
        for (int i = 0; i < d.N(); ++i)
        {
            const double rho=d.Q(i).rho;
            if (rho < 1e-10) continue;

            const TransportProps::Point& transport=
                points[static_cast<size_t>(i)];
            const double D = transport.diffusivity;
            const double alpha_T = transport.kappa
                /(rho*transport.Cv);

            double diff_max = std::max(D, alpha_T);
            if (diff_max > 1e-30)
            {
                double dt_i = cflDiff_ * dx2 / (2.0 * diff_max);
                if (dt_i < dtMin) dtMin = dt_i;
            }
        }
        return dtMin;
    }

    std::string name() const override
    { return "Diffusion (species Fick + thermal Fourier)"; }
};

} // End namespace Foam

#endif
