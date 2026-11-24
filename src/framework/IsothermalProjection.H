#ifndef IsothermalProjection_H
#define IsothermalProjection_H

#include "Operator.H"
#include "RealFluidState.H"
#include <algorithm>
#include <cmath>
#include <vector>

namespace Foam
{

// Driver-layer projection for non-reacting isothermal baselines. It operates
// after convection: mass, momentum and species remain untouched, while rhoE is
// reset to the EOS-consistent internal energy at the configured temperature.
class IsothermalProjectionOperator : public Operator
{
    double temperature_;

public:
    explicit IsothermalProjectionOperator(double temperature)
        : temperature_(temperature) {}

    void apply(Domain& d, double /*dt*/) const override
    {
        const std::vector<std::string>& names = d.speciesNames();
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < d.N(); ++i)
        {
            ConservedVars& u = d.Unew(i);
            if (u.rho <= 1e-12) continue;

            std::vector<double> Y(u.rhoY.size(), 0.0);
            double sumY = 0.0;
            for (size_t s = 0; s < Y.size(); ++s)
            {
                Y[s] = std::max(0.0, u.rhoY[s] / u.rho);
                sumY += Y[s];
            }
            if (sumY > 0.0)
                for (double& y : Y) y /= sumY;

            RealFluidState state;
            if (!names.empty()) state.setSpecies(names);
            if (state.updateRhoTY(u.rho, temperature_, Y) != RealFluidState::Status::SUCCESS)
                continue;

            double kinetic = 0.0;
            for (int m = 0; m < 3; ++m)
            {
                const double velocity = u.rhoU[m] / u.rho;
                kinetic += 0.5 * velocity * velocity;
            }
            u.rhoE = u.rho * (state.e() + kinetic);
        }
    }

    std::string name() const override { return "Isothermal projection"; }
};

} // End namespace Foam

#endif
