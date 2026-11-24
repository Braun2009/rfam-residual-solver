/* Cell-level test of the 9<->59 species bridge and the San Diego sector
   chemistry operator: expansion/projection round-trip, and a hot-cell
   autoignition push (fuel depletion, product formation, exothermic rhoE
   increase) while cold cells stay untouched. */

#include "Domain.H"
#include "SandiegoSectorChemistry.H"
#include "SandiegoNasa7Data.H"

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using namespace Foam;

namespace
{
constexpr double Runiv = 8.314462618;
int failures = 0;

void check(bool condition, const std::string& label)
{
    std::printf("  [%s] %s\n", condition ? "PASS" : "FAIL", label.c_str());
    if (!condition) ++failures;
}

PrimitiveVars primitive(
    double p, double T, const std::vector<double>& Y)
{
    PrimitiveVars w(static_cast<int>(Y.size()));
    w.p = p; w.T = T;
    w.U = {0.0, 0.0, 0.0};
    w.Y = Y;
    w.Z = 0.0; w.C = 0.0;
    return w;
}
} // namespace

int main()
{
    setvbuf(stdout, NULL, _IONBF, 0);
    std::printf("=== San Diego sector chemistry bridge ===\n");

    SpeciesRegistry registry;
    const std::vector<std::string> names = registry.names();
    check(names.size() == 9, "9-species registry");

    std::vector<double> Y9(9, 0.0);
    Y9[registry.index("CO2")] = 0.30;
    Y9[registry.index("H2O")] = 0.40;
    Y9[registry.index("N2")]  = 0.25;
    Y9[registry.index("O2")]  = 0.05;
    SandiegoBridge::SpeciesAdapter adapter(registry);
    std::vector<double> Y59, Y9back;
    adapter.expand(Y9, Y59);
    adapter.project(Y59, Y9back);
    check(std::abs(Y9back[registry.index("CO2")] - 0.30) < 1e-12
            && std::abs(Y9back[registry.index("H2O")] - 0.40) < 1e-12,
        "9->59->9 round trip preserves transported species");

    // Domain: 5x5x5, 9 species; cold background (700 K, N2), one hot cell
    // (2000 K, stoichiometric C3H8/O2).
    const int iC3 = registry.index("C3H8");
    const int iO2 = registry.index("O2");
    auto legacy = thermo::makeLegacySinglePhaseThermoModel(names);
    Domain dom(5, 5, 5, 0.05, 0.05, 0.05, 9, legacy);
    dom.setSpecies(names);
    std::vector<double> background(9, 0.0);
    background[registry.index("N2")] = 1.0;
    std::vector<double> hot(9, 0.0);
    hot[iC3] = 0.2222;
    hot[iO2] = 0.7778;
    for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
    for (int k = 0; k < 5; ++k)
    {
        const bool isHot = (i == 2 && j == 2 && k == 2);
        dom.stagedW(i, j, k) = primitive(
            4.0e6, isHot ? 2000.0 : 700.0, isHot ? hot : background);
    }
    check(dom.primitiveToConserved(), "domain initializes");
    dom.beginStep();

    SandiegoSectorChemistryOperator chemistry(registry, 700.0);
    bool threw = false;
    std::string error;
    try { chemistry.apply(dom, 2.0e-6); }
    catch (const std::runtime_error& e) { threw = true; error = e.what(); }
    check(!threw, "operator apply succeeds" + (threw ? ": " + error : ""));

    const auto& hotU = dom.Unew(2, 2, 2);
    const auto& coldU = dom.Unew(0, 0, 0);
    const double dFuel = hotU.rhoY[iC3] - dom.U(2, 2, 2).rhoY[iC3];
    const double dRhoE = hotU.rhoE - dom.U(2, 2, 2).rhoE;
    check(dFuel < -1e-6, "hot cell fuel mass decreases (autoignition progress)");
    check(dRhoE > 0.0, "hot cell rhoE increases (exothermic)");
    check(std::abs(coldU.rhoY[registry.index("N2")] - dom.U(0, 0, 0).rhoY[registry.index("N2")]) < 1e-12
            && std::abs(coldU.rhoE - dom.U(0, 0, 0).rhoE) < 1e-9,
        "cold cell stays untouched");
    check(chemistry.lastActive() >= 1 && chemistry.lastFailures() == 0,
        "operator reports active cells and zero failures");

    std::printf("%s\n", failures == 0
        ? "SANDIEGO_SECTOR_CHEM_BRIDGE_PASS"
        : "SANDIEGO_SECTOR_CHEM_BRIDGE_FAIL");
    return failures == 0 ? 0 : 1;
}
