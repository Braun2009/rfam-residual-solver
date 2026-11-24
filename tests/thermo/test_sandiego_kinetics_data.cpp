/* Structural + rate cross-check test for the generated San Diego kinetics
   table. For every reaction the high-rate Arrhenius parameters are evaluated
   at the three bundle reference-state temperatures and compared (after the
   kmol<->mol basis conversion) against the Cantera-embedded high rate
   constants. Falloff reactions additionally cross-check the low (k0) rate.
   Three-body/falloff effective-rate evaluation with concentrations is the
   next integration step and is not exercised here. */

#include "SandiegoKineticsData.H"
#include "SandiegoNasa7Data.H"
#include <cmath>
#include <cstdio>
#include <string>

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

double kArrhenius(const SandiegoKinetics::ArrheniusParams& p, double T)
{
    return p.A * std::pow(T, p.n) * std::exp(-p.Ea / (Runiv * T));
}
} // namespace

int main()
{
    std::printf("=== Sandiego 271-reaction kinetics data ===\n");

    constexpr int referenceTemperatures[3] = {1800, 2500, 3800};
    check(SandiegoKinetics::nReactions == 271,
        "exactly 271 reactions are registered");

    int nArr = 0, nTb = 0, nFt = 0, nDup = 0;
    for (int i = 0; i < SandiegoKinetics::nReactions; ++i)
    {
        const auto& r = SandiegoKinetics::reactions[i];
        if (r.type == 0) ++nArr;
        else if (r.type == 1) ++nTb;
        else ++nFt;
        if (r.duplicate) ++nDup;
    }
    check(nArr == 233 && nTb == 15 && nFt == 23,
        "reaction-type census matches the Cantera audit (233/15/23)");
    check(nDup == 6, "duplicate reaction census matches the Cantera audit (6)");

    double maxHighError = 0.0, maxLowError = 0.0;
    for (int i = 0; i < SandiegoKinetics::nReactions; ++i)
    {
        const auto& r = SandiegoKinetics::reactions[i];
        for (int ref = 0; ref < 3; ++ref)
        {
            const double T = referenceTemperatures[ref];
            // A_mol = A_kmol*1000^(1-order), so the kmol-basis reported rate
            // is recovered by multiplying by 1000^(order-1).
            const double kMol = kArrhenius(r.high, T);
            const double kKmol = kMol * std::pow(1000.0, r.high.order - 1.0);
            const double error =
                std::abs(kKmol - r.highKRef[ref]) / std::max(std::abs(r.highKRef[ref]), 1.0);
            maxHighError = std::max(maxHighError, error);
            if (r.hasLow)
            {
                const double k0Mol = kArrhenius(r.low, T);
                const double k0Kmol = k0Mol * std::pow(1000.0, r.low.order - 1.0);
                const double lowError =
                    std::abs(k0Kmol - r.lowKRef[ref]) / std::max(std::abs(r.lowKRef[ref]), 1.0);
                maxLowError = std::max(maxLowError, lowError);
            }
        }
        // Element balance across reactants/products must be exact using the
        // San Diego element compositions (C/H/O/N).
        long double c = 0, h = 0, o = 0, n = 0;
        for (int k = 0; k < r.nReactants; ++k)
        {
            const auto& sp = SandiegoNasa7::species[r.reactants[k].species];
            c -= r.reactants[k].coeff * sp.nC;
            h -= r.reactants[k].coeff * sp.nH;
            o -= r.reactants[k].coeff * sp.nO;
            n -= r.reactants[k].coeff * sp.nN;
        }
        for (int k = 0; k < r.nProducts; ++k)
        {
            const auto& sp = SandiegoNasa7::species[r.products[k].species];
            c += r.products[k].coeff * sp.nC;
            h += r.products[k].coeff * sp.nH;
            o += r.products[k].coeff * sp.nO;
            n += r.products[k].coeff * sp.nN;
        }
        if (std::abs(c) > 1e-10 || std::abs(h) > 1e-10
            || std::abs(o) > 1e-10 || std::abs(n) > 1e-10)
        {
            check(false, std::string("element balance for ") + r.equation);
        }
    }
    check(maxHighError < 1.0e-9,
        "all 271 high-rate constants match Cantera to 1e-9 "
        "(max_rel=" + std::to_string(maxHighError) + ")");
    check(maxLowError < 1.0e-9,
        "all 23 falloff low-rate constants match Cantera to 1e-9 "
        "(max_rel=" + std::to_string(maxLowError) + ")");

    std::printf("%s\n", failures == 0
        ? "SANDIEGO_KINETICS_DATA_PASS"
        : "SANDIEGO_KINETICS_DATA_FAIL");
    return failures == 0 ? 0 : 1;
}
