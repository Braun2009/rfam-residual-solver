#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "DiffusionOperator.H"
#include "RealFluidState.H"
#include "SpeciesRegistry.H"
#include "TransportProps.H"

using namespace Foam;

namespace
{

int failures = 0;

void check(const std::string& name, bool ok, double value = 0.0)
{
    std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << name;
    if (!ok) std::cout << " (value=" << value << ")";
    std::cout << '\n';
    if (!ok) ++failures;
}

} // namespace

int main()
{
    std::cout << "=== Diffusion real-fluid Cv stability closure ===\n";
    SpeciesRegistry registry;
    const std::vector<std::string> names = registry.names();
    Domain domain(3, 0.03, registry.size());
    domain.setSpecies(names);

    std::vector<double> oxygen(names.size(), 0.0);
    oxygen[registry.index("O2")] = 1.0;
    for (int i = 0; i < domain.N(); ++i)
    {
        PrimitiveVars& state = domain.stagedW(i);
        state.p = 4.0e6;
        state.T = 300.0;
        state.U[0] = state.U[1] = state.U[2] = 0.0;
        state.Y = oxygen;
    }
    domain.primitiveToConserved();
    check("initial dense O2 field converts", !domain.tripped());

    RealFluidState thermo;
    thermo.setSpecies(names);
    check("reference high-pressure O2 gas state succeeds",
        thermo.updatePTY(4.0e6, 300.0, oxygen)
            == RealFluidState::Status::SUCCESS);

    TransportProps::Point transport;
    const bool transportOk = TransportProps::evaluateSinglePhaseGas(
        names, oxygen, 4.0e6, 300.0, thermo.rho(), registry, transport);
    check("single transport point succeeds", transportOk);
    if (!transportOk) return 1;
    check("validated gas closure accepts reference state",
        TransportProps::gasClosureSupported(
            4.0e6, 300.0, thermo.rho()));

    RealFluidState cryogenic;
    cryogenic.setSpecies(names);
    check("cryogenic O2 thermo state succeeds for capability audit",
        cryogenic.updatePTY(4.0e6, 130.0, oxygen)
            == RealFluidState::Status::SUCCESS);
    check("gas transport rejects dense cryogenic O2",
        !TransportProps::gasClosureSupported(
            4.0e6, 130.0, cryogenic.rho()));
    const double speciesDiffusivity = transport.diffusivity;
    const double thermalDiffusivity = transport.kappa
        /(thermo.rho()*transport.Cv);
    const double controllingDiffusivity = std::max(
        speciesDiffusivity, thermalDiffusivity);
    const double expected = 0.3*domain.dx()*domain.dx()
        /(2.0*controllingDiffusivity);

    DiffusionOperator diffusion(registry);
    const double actual = diffusion.maxStableDt(domain, 1.0);
    const double relativeError = std::abs(actual-expected)
        /std::max(std::abs(expected), 1.0e-30);
    check("diffusion dt uses EOS Cv", relativeError < 1.0e-12,
        relativeError);

    std::cout << (failures == 0
        ? "DIFFUSION_REAL_CV_PASS\n"
        : "DIFFUSION_REAL_CV_FAIL\n");
    return failures == 0 ? 0 : 1;
}
