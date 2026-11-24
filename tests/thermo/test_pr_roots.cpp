#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "BinaryInteractionDB.H"
#include "PengRobinsonEOS.H"
#include "SpeciesCriticalDB.H"

using namespace Foam;

namespace
{

int failures = 0;

void check(
    const std::string& name, bool ok,
    double got = 0.0, double expected = 0.0)
{
    std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << name;
    if (!ok)
        std::cout << " got=" << std::setprecision(12) << got
                  << " expected=" << expected;
    std::cout << '\n';
    if (!ok) ++failures;
}

bool sameRoot(double lhs, double rhs)
{
    return std::abs(lhs - rhs)
        <= 1.0e-11*std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

void checkBoundaryContinuity(
    const std::string& label,
    const PengRobinsonEOS& eos,
    const std::vector<double>& Y,
    double pBelow, double TBelow,
    double pAbove, double TAbove)
{
    const auto below = eos.solveZ(pBelow, TBelow, Y);
    const auto above = eos.solveZ(pAbove, TAbove, Y);
    check(label + " states valid", below.valid && above.valid);
    if (below.valid && above.valid)
    {
        check(label + " has no artificial Z jump",
            std::abs(above.Z - below.Z) < 1.0e-4,
            above.Z, below.Z);
    }
}

} // namespace

int main()
{
    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB kij = BinaryInteractionDB::createDefault();
    PengRobinsonEOS eos(critical, kij);
    eos.setSpeciesNames({"C3H8"});

    const double p = 1.0e6;
    const double T = 300.0;
    const std::vector<double> Y = {1.0};

    std::cout << "=== Peng-Robinson root and phase-hint contract ===\n";

    const auto gas = eos.solveZ(
        p, T, Y, PengRobinsonEOS::PhaseHint::GAS);
    const auto liquid = eos.solveZ(
        p, T, Y, PengRobinsonEOS::PhaseHint::LIQUID);
    const auto automatic = eos.solveZ(
        p, T, Y, PengRobinsonEOS::PhaseHint::AUTO);

    check("three real roots found", gas.numRoots == 3, gas.numRoots, 3.0);
    check("gas root valid", gas.valid);
    check("liquid root valid", liquid.valid);
    check("automatic root valid", automatic.valid);
    check("GAS selects larger physical root",
        gas.valid && liquid.valid && gas.Z > liquid.Z, gas.Z, liquid.Z);

    if (gas.valid && liquid.valid && automatic.valid)
    {
        const auto& stable = gas.G < liquid.G ? gas : liquid;
        check("AUTO selects minimum-Gibbs physical root",
            sameRoot(automatic.Z, stable.Z), automatic.Z, stable.Z);
    }

    const double nan = std::numeric_limits<double>::quiet_NaN();
    check("NaN pressure returns invalid root", !eos.solveZ(nan, T, Y).valid);
    check("NaN temperature returns invalid root", !eos.solveZ(p, nan, Y).valid);
    check("negative pressure returns invalid root", !eos.solveZ(-p, T, Y).valid);
    check("bad composition size returns invalid root",
        !eos.solveZ(p, T, {0.5, 0.5}).valid);
    check("non-normalized composition returns invalid root",
        !eos.solveZ(p, T, {0.5}).valid);
    check("density fails closed on invalid state",
        std::isnan(eos.density(nan, T, Y)));

    const CriticalProperties& propane = critical.get("C3H8");
    constexpr double epsilon = 1.0e-7;
    checkBoundaryContinuity(
        "former lower-Tr boundary", eos, Y,
        1.2*propane.Pc, 0.95*propane.Tc*(1.0 - epsilon),
        1.2*propane.Pc, 0.95*propane.Tc*(1.0 + epsilon));
    checkBoundaryContinuity(
        "former upper-Tr boundary", eos, Y,
        1.2*propane.Pc, 1.05*propane.Tc*(1.0 - epsilon),
        1.2*propane.Pc, 1.05*propane.Tc*(1.0 + epsilon));
    checkBoundaryContinuity(
        "former lower-Pr boundary", eos, Y,
        0.9*propane.Pc*(1.0 - epsilon), propane.Tc,
        0.9*propane.Pc*(1.0 + epsilon), propane.Tc);
    checkBoundaryContinuity(
        "former upper-Pr boundary", eos, Y,
        1.5*propane.Pc*(1.0 - epsilon), propane.Tc,
        1.5*propane.Pc*(1.0 + epsilon), propane.Tc);

    std::cout << (failures == 0 ? "PR_ROOTS_PASS\n" : "PR_ROOTS_FAIL\n");
    return failures == 0 ? 0 : 1;
}
