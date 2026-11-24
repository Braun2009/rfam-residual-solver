/*---------------------------------------------------------------------------*\
    solver_frame_precond - legacy filename for a no-advance CFL-scaler audit.

    The experimental scaler has no residual/Jacobian transformation. This
    diagnostic therefore must not march with the enlarged heuristic timestep.
    It proves only that the CFL estimate changes while the same-dt AUSM
    conservative update remains identical.
\*---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>

#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"

using namespace Foam;

namespace
{
void initialize(Domain& domain)
{
    for (int i = 0; i < domain.N(); ++i)
    {
        PrimitiveVars& w = domain.stagedW(i);
        w.p = i < domain.N()/2 ? 6.0e6 : 5.9e6;
        w.T = 300.0;
        w.U = {2.0, 0.0, 0.0};
        w.Y[0] = 1.0;
        w.Y[1] = 0.0;
    }
}

double difference(const ConservedVars& a, const ConservedVars& b)
{
    double result = std::max(std::abs(a.rho-b.rho), std::abs(a.rhoE-b.rhoE));
    for (std::size_t i = 0; i < a.rhoU.size(); ++i)
        result = std::max(result, std::abs(a.rhoU[i]-b.rhoU[i]));
    for (std::size_t i = 0; i < a.rhoY.size(); ++i)
        result = std::max(result, std::abs(a.rhoY[i]-b.rhoY[i]));
    return result;
}
}

int main()
{
    std::cout << "=== Experimental CFL acoustic-scaler no-advance audit ===\n";
    std::cout << "CFL_SCALER_ONLY: residual/Jacobian preconditioning is absent\n";

    Domain physicalDomain(200, 1.0, 2);
    Domain scaledDomain(200, 1.0, 2);
    initialize(physicalDomain);
    initialize(scaledDomain);
    if (!physicalDomain.primitiveToConserved()
        || !scaledDomain.primitiveToConserved())
    {
        std::cerr << "FRAME_CFL_SCALER_DIAGNOSTIC_FAIL initialization\n";
        return 1;
    }

    ConvectionOperator physicalOperator(
        std::unique_ptr<FluxScheme>(new AUSMFluxAdapter()), false);
    ConvectionOperator scaledCflOperator(
        std::unique_ptr<FluxScheme>(new AUSMFluxAdapter()), true);
    const double physicalDt = physicalOperator.maxStableDt(physicalDomain, 0.5);
    const double scaledDt = scaledCflOperator.maxStableDt(scaledDomain, 0.5);

    // Deliberately use the physical acoustic timestep for both updates. The
    // enlarged heuristic timestep is not stability-authorized.
    const double commonDt = 0.2 * physicalDt;
    physicalDomain.beginStep();
    scaledDomain.beginStep();
    physicalOperator.apply(physicalDomain, commonDt);
    scaledCflOperator.apply(scaledDomain, commonDt);

    double maxDifference = 0.0;
    for (int i = 0; i < physicalDomain.N(); ++i)
        maxDifference = std::max(maxDifference,
            difference(physicalDomain.Unew(i), scaledDomain.Unew(i)));

    const bool passed = std::isfinite(physicalDt) && physicalDt > 0.0
        && std::isfinite(scaledDt) && scaledDt > physicalDt
        && maxDifference == 0.0;
    std::cout << "physical_dt=" << physicalDt
              << " scaled_dt=" << scaledDt
              << " dt_ratio=" << scaledDt/physicalDt
              << " same_dt_update_difference=" << maxDifference << '\n';
    std::cout << (passed
        ? "FRAME_CFL_SCALER_DIAGNOSTIC_PASS"
        : "FRAME_CFL_SCALER_DIAGNOSTIC_FAIL") << '\n';
    return passed ? 0 : 1;
}
