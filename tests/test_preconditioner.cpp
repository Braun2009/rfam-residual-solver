#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "src/numerics/preconditioning/ExperimentalCflAcousticScaler.H"
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"

using namespace Foam;

namespace
{
int failures = 0;

void check(const std::string& label, bool condition)
{
    std::cout << "  [" << (condition ? "PASS" : "FAIL") << "] "
              << label << '\n';
    if (!condition) ++failures;
}

void initializeLowMachDomain(Domain& domain)
{
    for (int i = 0; i < domain.N(); ++i)
    {
        PrimitiveVars& w = domain.stagedW(i);
        w.p = (i < domain.N()/2) ? 5.0e5 : 4.95e5;
        w.T = 300.0;
        w.U = {2.0, 0.0, 0.0};
        w.Y[0] = 1.0;
        w.Y[1] = 0.0;
    }
    check("domain initialization", domain.primitiveToConserved());
}

double conservedDifference(const ConservedVars& a, const ConservedVars& b)
{
    double result = std::max(std::abs(a.rho-b.rho), std::abs(a.rhoE-b.rhoE));
    for (std::size_t i = 0; i < a.rhoU.size(); ++i)
        result = std::max(result, std::abs(a.rhoU[i]-b.rhoU[i]));
    for (std::size_t i = 0; i < a.rhoY.size(); ++i)
        result = std::max(result, std::abs(a.rhoY[i]-b.rhoY[i]));
    result = std::max(result, std::abs(a.rhoZ-b.rhoZ));
    result = std::max(result, std::abs(a.rhoC-b.rhoC));
    return result;
}
}

int main()
{
    std::cout << "=== Experimental CFL acoustic scaler contract ===\n";

    static_assert(!ExperimentalCflAcousticScaler::modifiesFluxResidual(),
        "A CFL-only scaler must not claim to modify the residual");
    static_assert(!ExperimentalCflAcousticScaler::modifiesTimeDerivative(),
        "A CFL-only scaler must not claim a time-derivative transform");
    static_assert(!ExperimentalCflAcousticScaler::modifiesJacobian(),
        "A CFL-only scaler must not claim a Jacobian transform");

    check("capability is labeled CFL-only",
        std::string(ExperimentalCflAcousticScaler::capabilityLabel())
            == "experimental-cfl-acoustic-scaler-only");

    RealFluidState state;
    const std::vector<double> composition{1.0, 0.0};
    check("thermodynamic reference state",
        state.updatePTY(5.0e5, 300.0, composition)
            == RealFluidState::Status::SUCCESS);
    const double physicalSoundSpeed = state.a();

    ExperimentalCflAcousticScaler scaler;
    const std::vector<double> velocity{2.0, 0.0, 0.0};
    scaler.updateForMach(state.getMach(velocity));
    check("disabled by default",
        scaler.cflAcousticSpeed(state, velocity) == physicalSoundSpeed);

    scaler.setEnabled(true);
    scaler.updateForMach(state.getMach(velocity));
    const double cflSpeed = scaler.cflAcousticSpeed(state, velocity);
    check("low-Mach CFL speed is reduced", cflSpeed < physicalSoundSpeed);
    check("physical EOS sound speed is unchanged",
        state.a() == physicalSoundSpeed);
    check("scale is finite and bounded",
        std::isfinite(scaler.acousticScaleSquared())
        && scaler.acousticScaleSquared() > 0.0
        && scaler.acousticScaleSquared() <= 1.0);

    bool rejectedNonFiniteMach = false;
    try
    {
        scaler.updateForMach(std::numeric_limits<double>::quiet_NaN());
    }
    catch (const std::invalid_argument&)
    {
        rejectedNonFiniteMach = true;
    }
    check("non-finite Mach fails closed", rejectedNonFiniteMach);

    bool rejectedNonFiniteVelocity = false;
    try
    {
        scaler.updateForMach(0.01);
        scaler.cflAcousticSpeed(state,
            {std::numeric_limits<double>::infinity(), 0.0, 0.0});
    }
    catch (const std::invalid_argument&)
    {
        rejectedNonFiniteVelocity = true;
    }
    check("non-finite velocity fails closed", rejectedNonFiniteVelocity);

    Domain physicalDomain(8, 1.0, 2);
    Domain scaledDomain(8, 1.0, 2);
    initializeLowMachDomain(physicalDomain);
    initializeLowMachDomain(scaledDomain);

    ConvectionOperator physicalOperator(
        std::unique_ptr<FluxScheme>(new AUSMFluxAdapter()), false);
    ConvectionOperator scaledCflOperator(
        std::unique_ptr<FluxScheme>(new AUSMFluxAdapter()), true);

    const double physicalDt = physicalOperator.maxStableDt(physicalDomain, 0.4);
    const double scaledDt = scaledCflOperator.maxStableDt(scaledDomain, 0.4);
    check("experimental switch changes only the CFL estimate",
        std::isfinite(physicalDt) && std::isfinite(scaledDt)
        && scaledDt > physicalDt);
    check("operator name discloses experimental CFL-only scope",
        scaledCflOperator.name()
            == "Convection (AUSM; experimental CFL scaler only)");

    const double commonDt = 0.2 * std::min(physicalDt, scaledDt);
    physicalDomain.beginStep();
    scaledDomain.beginStep();
    physicalOperator.apply(physicalDomain, commonDt);
    scaledCflOperator.apply(scaledDomain, commonDt);

    double residualDifference = 0.0;
    for (int i = 0; i < physicalDomain.N(); ++i)
        residualDifference = std::max(residualDifference,
            conservedDifference(
                physicalDomain.Unew(i), scaledDomain.Unew(i)));
    check("same dt gives identical AUSM conservative update",
        residualDifference == 0.0);

    std::cout << std::setprecision(12)
              << "  physical_a=" << physicalSoundSpeed
              << " cfl_a=" << cflSpeed
              << " dt_ratio=" << scaledDt/physicalDt << '\n';
    std::cout << (failures == 0
        ? "EXPERIMENTAL_CFL_SCALER_CONTRACT_PASS"
        : "EXPERIMENTAL_CFL_SCALER_CONTRACT_FAIL") << '\n';
    return failures == 0 ? 0 : 1;
}
