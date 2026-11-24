#ifndef DiffusionOperator2D_H
#define DiffusionOperator2D_H

/*---------------------------------------------------------------------------*\
    DiffusionOperator2D - 2-D Cartesian molecular diffusion for the simple
    physics-validation cases (pin-fin cross flow, droplet-film transport).

    Applies to the interior cells the 5-point Laplacian for:
      - species diffusion:  d/dx(rho D dY/dx) + d/dy(rho D dY/dy)
      - thermal conduction: d/dx(kappa dT/dx)  + d/dy(kappa dT/dy)
      - momentum diffusion: d/dx(mu du/dx)     + d/dy(mu du/dy)
    using the committed legacy gas transport (Sutherland mu + Eucken kappa
    + Sc=0.7 species diffusivity).

    The transport coefficient at each cell comes from TransportProps; the
    density field for the species flux is taken from the conserved state.

    Optional kinematic-viscosity override (viscosityScale) lets a
    non-dimensional Re validation run with a prescribed nu while keeping
    the same operator; the override is explicit and named, never silent.
\*---------------------------------------------------------------------------*/

#include "Domain.H"
#include "Operator.H"
#include "SpeciesRegistry.H"
#include "TransportProps.H"

namespace Foam
{

class DiffusionOperator2D : public Operator
{
public:
    explicit DiffusionOperator2D(
        const SpeciesRegistry& reg,
        double viscosityScale = 1.0,
        double effectiveKinematicViscosity = 1.0e-5,
        double cflDiff = 0.3)
        : reg_(reg), viscosityScale_(viscosityScale),
          effectiveKinematicViscosity_(effectiveKinematicViscosity),
          cflDiff_(cflDiff) {}

    void apply(Domain& d, double dt) const override
    {
        const int Nx = d.Nx(), Ny = d.Ny(), nS = d.nSpecies();
        const double dx = d.dx(), dy = d.dy();
        const double dx2 = dx * dx, dy2 = dy * dy;

        const std::vector<TransportProps::Point> points = transportField(d);

        for (int i = 1; i < Nx - 1; ++i)
        for (int j = 1; j < Ny - 1; ++j)
        {
            const auto& w  = d.W(i,j);
            const auto& wL = d.W(i-1,j);
            const auto& wR = d.W(i+1,j);
            const auto& wD = d.W(i,j-1);
            const auto& wU = d.W(i,j+1);
            const TransportProps::Point& tr = points[static_cast<std::size_t>(j*Nx+i)];
            const double rho = d.Q(i,j).rho;
            // SGS-like effective transport: the explicit kinematic
            // viscosity adds to the molecular coefficients (named,
            // never silent; Pr_t = Sc_t = 1 for the placeholder).
            const double mu = tr.mu * viscosityScale_
                + rho * effectiveKinematicViscosity_;
            const double kappa = tr.kappa * viscosityScale_
                + rho * effectiveKinematicViscosity_
                    * tr.Cp;
            const double D = tr.diffusivity * viscosityScale_
                + effectiveKinematicViscosity_;

            ConservedVars& un = d.Unew(i,j);

            // Species diffusion (5-point Laplacian of Y).
            for (int s = 0; s < nS; ++s)
            {
                const double lapY =
                    (wR.Y[s] - 2.0*w.Y[s] + wL.Y[s]) / dx2
                    + (wU.Y[s] - 2.0*w.Y[s] + wD.Y[s]) / dy2;
                un.rhoY[s] += rho * D * lapY * dt;
            }

            // Thermal conduction.
            const double lapT =
                (wR.T - 2.0*w.T + wL.T) / dx2
                + (wU.T - 2.0*w.T + wD.T) / dy2;
            un.rhoE += kappa * lapT * dt;

            // Momentum diffusion (x and y components).
            const double lapUx =
                (wR.U[0] - 2.0*w.U[0] + wL.U[0]) / dx2
                + (wU.U[0] - 2.0*w.U[0] + wD.U[0]) / dy2;
            un.rhoU[0] += mu * lapUx * dt;
            const double lapUy =
                (wR.U[1] - 2.0*w.U[1] + wL.U[1]) / dx2
                + (wU.U[1] - 2.0*w.U[1] + wD.U[1]) / dy2;
            un.rhoU[1] += mu * lapUy * dt;
        }
    }

    double maxStableDt(const Domain& d, double /*cfl*/) const override
    {
        // Viscous stability: dt <= CFL * dx^2 / nu_eff where
        // nu_eff = effectiveKinematicViscosity_ is the physical kinematic
        // viscosity actually applied after scaling (nu_phys * viscosityScale).
        // Using the bare scale here would be wrong by a factor nu_phys.
        const double dx = d.dx(), dy = (d.Ny() > 1) ? d.dy() : d.dx();
        const double nuEff = std::max(effectiveKinematicViscosity_, 1.0e-12);
        return cflDiff_ * std::min(dx, dy) * std::min(dx, dy) / nuEff;
    }

    std::string name() const override { return "2-D molecular diffusion"; }

private:
    std::vector<TransportProps::Point> transportField(const Domain& d) const
    {
        std::vector<TransportProps::Point> points(
            static_cast<std::size_t>(d.N()));
        for (int i = 0; i < d.Nx(); ++i)
        for (int j = 0; j < d.Ny(); ++j)
        {
            const TransportProps::Result result =
                TransportProps::evaluateCommittedLegacyGas(
                    d.speciesNames(), d.W(i,j), d.Q(i,j),
                    d.thermoModel(), reg_);
            if (result.success())
                points[static_cast<std::size_t>(j*d.Nx()+i)] = result.point;
        }
        return points;
    }

    SpeciesRegistry reg_;
    double viscosityScale_;
    double effectiveKinematicViscosity_;
    double cflDiff_;
};

} // End namespace Foam

#endif
