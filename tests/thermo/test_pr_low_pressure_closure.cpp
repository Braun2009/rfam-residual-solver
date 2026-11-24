#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "BinaryInteractionDB.H"
#include "PengRobinsonEOS.H"
#include "SpeciesCriticalDB.H"

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
    std::cout << "=== PR low-pressure state-surface closure ===\n";
    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB interactions = BinaryInteractionDB::createDefault();
    PengRobinsonEOS eos(critical, interactions);
    const std::vector<double> oxygen{1.0, 0.0};
    const double criticalPressure = eos.pseudoCriticalPressure(oxygen);
    // Stay on the single gas branch. At 90 K the wider scan crosses the real
    // oxygen saturation boundary, which must not be confused with the removed
    // artificial Pr=0.01 ideal-gas switch.
    const double temperature = 130.0;

    double previousZ = 0.0;
    bool havePrevious = false;
    for (double reducedPressure : {0.002, 0.005, 0.0099, 0.0101, 0.02})
    {
        const double pressure = reducedPressure*criticalPressure;
        const auto root = eos.solveZ(pressure, temperature, oxygen);
        check("finite physical root at Pr=" + std::to_string(reducedPressure),
            root.valid && std::isfinite(root.Z) && root.Z > 0.0, root.Z);
        if (!root.valid) continue;

        const double density = eos.density(pressure, temperature, oxygen);
        const double recoveredPressure = eos.pressure(
            density, temperature, oxygen);
        const double relativeClosure = std::abs(recoveredPressure-pressure)
            /std::max(pressure, 1.0);
        check("p-rho-T closure at Pr=" + std::to_string(reducedPressure),
            std::isfinite(relativeClosure) && relativeClosure < 1.0e-9,
            relativeClosure);

        if (havePrevious)
            check("root branch remains continuous",
                std::abs(root.Z-previousZ) < 0.1,
                std::abs(root.Z-previousZ));
        previousZ = root.Z;
        havePrevious = true;
    }

    std::cout << (failures == 0
        ? "PR_LOW_PRESSURE_CLOSURE_PASS\n"
        : "PR_LOW_PRESSURE_CLOSURE_FAIL\n");
    return failures == 0 ? 0 : 1;
}
