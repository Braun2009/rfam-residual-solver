/* Structural test for the generated San Diego 59-species NASA-7 data table.
   Verifies counts, index lookup, thermodynamic self-consistency (dh/dT == cp
   via finite differences), and monotonicity of sensible enthalpy across the
   valid range. The coefficient fidelity to the Cantera source is locked by
   tools/audit_sandiego_nasa7_header.py (offline, WSL). */

#include "SandiegoNasa7Data.H"
#include <cmath>
#include <cstdio>
#include <string>

using namespace Foam;

namespace
{
int failures = 0;

void check(bool condition, const std::string& label)
{
    std::printf("  [%s] %s\n", condition ? "PASS" : "FAIL", label.c_str());
    if (!condition) ++failures;
}
} // namespace

int main()
{
    std::printf("=== Sandiego 59-species NASA-7 data ===\n");

    check(SandiegoNasa7::nSpecies == 59,
        "exactly 59 species are registered");
    check(SandiegoNasa7::index("N2") >= 0,
        "N2 is present");
    check(SandiegoNasa7::index("C3H8") >= 0,
        "C3H8 is present");
    check(SandiegoNasa7::index("OH") >= 0,
        "OH is present");
    check(SandiegoNasa7::index("NOT_A_SPECIES") < 0,
        "unknown species resolves to -1");

    // Composition spot checks.
    const int iN2 = SandiegoNasa7::index("N2");
    const int iC3 = SandiegoNasa7::index("C3H8");
    const int iOH = SandiegoNasa7::index("OH");
    check(SandiegoNasa7::species[iN2].nN == 2 && SandiegoNasa7::species[iN2].nC == 0,
        "N2 element composition");
    check(SandiegoNasa7::species[iC3].nC == 3 && SandiegoNasa7::species[iC3].nH == 8,
        "C3H8 element composition");
    check(SandiegoNasa7::species[iOH].nO == 1 && SandiegoNasa7::species[iOH].nH == 1,
        "OH element composition");

    // cp > 0 everywhere and sensible-enthalpy monotonic across each species'
    // valid range: dh/dT == cp > 0 (finite-difference self-consistency).
    double maxDhdTError = 0.0;
    for (int i = 0; i < SandiegoNasa7::nSpecies; ++i)
    {
        const SandiegoNasa7::SandiegoSpeciesData& s = SandiegoNasa7::species[i];
        for (double T = s.Tlow + 1.0; T < s.Thigh - 1.0; T += 50.0)
        {
            if (std::abs(T - s.Tcommon) < 1.0) continue;
            const double cp = SandiegoNasa7::cpMole(i, T);
            if (!(cp > 0.0)) { check(false, std::string("cp>0 for ") + s.name); goto done; }
            const double dT = 0.01;
            const double dhdT = (SandiegoNasa7::hMole(i, T + dT)
                - SandiegoNasa7::hMole(i, T - dT)) / (2.0 * dT);
            maxDhdTError = std::max(maxDhdTError,
                std::abs(dhdT - cp) / std::max(std::abs(cp), 1.0));
        }
    }
done:
    check(maxDhdTError < 1.0e-6,
        "finite-difference dh/dT matches cp to 1e-6 (monotone, self-consistent)");

    std::printf("%s\n", failures == 0
        ? "SANDIEGO_NASA7_DATA_PASS"
        : "SANDIEGO_NASA7_DATA_FAIL");
    return failures == 0 ? 0 : 1;
}
