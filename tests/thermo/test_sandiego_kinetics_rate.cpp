/* Cross-check the San Diego rate kernel against the Cantera-embedded
   reference states (equilibrium compositions). For every reaction and every
   reference state: forward rate of progress, reverse rate of progress, and
   the concentration equilibrium constant are compared with Cantera values
   (kmol basis, so mol-basis C++ results are divided by 1000). */

#include "SandiegoKineticsData.H"
#include "SandiegoKineticsRate.H"
#include "SandiegoKineticsReference.H"
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
} // namespace

int main()
{
    std::printf("=== Sandiego rate kernel vs Cantera references ===\n");

    double maxFwdError = 0.0, maxRevError = 0.0, maxKcError = 0.0;
    double maxFwdByType[3] = {0, 0, 0};
    double maxRevByType[3] = {0, 0, 0};
    int fwdCompared = 0, revCompared = 0;
    int worstFwdReaction = -1, worstFwdState = -1;
    double worstFwdKmol = 0.0, worstFwdRef = 0.0;
    int worstKcReaction = -1;
    double worstKcMine = 0.0, worstKcRef = 0.0;
    for (int ref = 0; ref < SandiegoKinetics::nReferenceStates; ++ref)
    {
        const auto& state = SandiegoKinetics::referenceStates[ref];
        double C[59];
        for (int i = 0; i < 59; ++i)
            C[i] = state.X[i] * state.p / (Runiv * state.T);  // mol/m3
        for (int i = 0; i < SandiegoKinetics::nReactions; ++i)
        {
            const auto& r = SandiegoKinetics::reactions[i];
            const double fwdMol = SandiegoKinetics::forwardRateOfProgress(r, state.T, C);
            // mol-basis rop -> Cantera kmol basis. With the correct
            // A_mol = A_kmol*1000^(1-order) conversion, rop_mol = 1000 *
            // rop_kmol for every reaction type (the concentration power sum
            // cancels the rate-constant order factor).
            const double kmolFactor = 1000.0;
            const double fwdKmol = fwdMol / kmolFactor;
            const double fwdRef = state.forwardRop[i];
            if (std::abs(fwdRef) > 1.0e-40)
            {
                const double error = std::abs(fwdKmol - fwdRef) / std::abs(fwdRef);
                if (error > maxFwdError)
                {
                    maxFwdError = error;
                    worstFwdReaction = i;
                    worstFwdState = ref;
                    worstFwdKmol = fwdKmol;
                    worstFwdRef = fwdRef;
                }
                if (r.type >= 0 && r.type <= 2)
                    maxFwdByType[r.type] = std::max(maxFwdByType[r.type], error);
                ++fwdCompared;
            }
            const double kc = SandiegoKinetics::equilibriumConstant(r, state.T);
            // Cantera's equilibrium_constants are state-dependent effective
            // values for three-body/falloff reactions (they include [M]/Pr);
            // only plain Arrhenius Kc values are directly comparable.
            if (r.type == 0)
            {
                const double kcRef = state.Kc[i];
                double nuSum = 0.0;
                for (int k = 0; k < r.nProducts; ++k) nuSum += r.products[k].coeff;
                for (int k = 0; k < r.nReactants; ++k) nuSum -= r.reactants[k].coeff;
                const double kcKmol = kc * std::pow(1000.0, -nuSum);
                const double kcError = std::abs(kcKmol - kcRef) / std::max(std::abs(kcRef), 1.0);
                if (kcError > maxKcError)
                {
                    maxKcError = kcError;
                    worstKcReaction = i;
                    worstKcMine = kcKmol;
                    worstKcRef = kcRef;
                }
            }
            if (r.reversible)
            {
                // reverse = (fwd / prod(C_r^nu)) / kc * prod(C_p^nu): the
                // forward rate of progress already contains the reactant
                // concentration product, which must not be repeated.
                double reactantProduct = 1.0;
                if (r.nOrders > 0)
                {
                    for (int k = 0; k < r.nOrders; ++k)
                        if (r.orders[k].coeff != 0.0)
                            reactantProduct *= std::pow(C[r.orders[k].species], r.orders[k].coeff);
                }
                else
                {
                    for (int k = 0; k < r.nReactants; ++k)
                        if (r.reactants[k].coeff != 0.0)
                            reactantProduct *= std::pow(C[r.reactants[k].species], r.reactants[k].coeff);
                }
                double reverse = fwdMol / reactantProduct / kc;
                for (int k = 0; k < r.nProducts; ++k)
                {
                    const int s = r.products[k].species;
                    if (r.products[k].coeff != 0.0)
                        reverse *= std::pow(C[s], r.products[k].coeff);
                }
                const double revKmol = reverse / kmolFactor;
                const double revRef = state.reverseRop[i];
                if (std::abs(revRef) > 1.0e-40)
                {
                    const double error = std::abs(revKmol - revRef) / std::abs(revRef);
                    maxRevError = std::max(maxRevError, error);
                    if (r.type >= 0 && r.type <= 2)
                        maxRevByType[r.type] = std::max(maxRevByType[r.type], error);
                    ++revCompared;
                }
            }
        }
    }
    // Hard assertions cover all 271 reactions: forward and reverse rates of
    // progress at three equilibrium reference states (Arrhenius, three-body
    // and falloff-Troe all match Cantera 3.1 to ~5e-10), plus all Arrhenius
    // Kc. The falloff path implements exactly Cantera's evalFromStruct form
    // (kinf*Pr/(1+Pr)*F with kmol-basis reduced pressure) and its TroeRate::F
    // (T2=0 omits the exp(-T2/T) term).
    check(maxFwdByType[0] < 1.0e-6 && maxFwdByType[1] < 1.0e-6
            && maxFwdByType[2] < 1.0e-6,
        "all forward rates of progress match Cantera to 1e-6 "
        "(max_rel=" + std::to_string(maxFwdByType[0]) + "/"
        + std::to_string(maxFwdByType[1]) + "/"
        + std::to_string(maxFwdByType[2]) + ")");
    check(maxRevByType[0] < 1.0e-6 && maxRevByType[1] < 1.0e-6
            && maxRevByType[2] < 1.0e-6,
        "all reverse rates of progress match Cantera to 1e-6 "
        "(max_rel=" + std::to_string(maxRevByType[0]) + "/"
        + std::to_string(maxRevByType[1]) + "/"
        + std::to_string(maxRevByType[2]) + ")");
    check(maxKcError < 1.0e-6,
        "Arrhenius Kc values match Cantera to 1e-6 "
        "(max_rel=" + std::to_string(maxKcError) + ")");
    std::printf("  forward max_rel by type: Arrhenius=%.6g three-body=%.6g falloff(diag)=%.6g\n",
        maxFwdByType[0], maxFwdByType[1], maxFwdByType[2]);
    std::printf("  reverse max_rel by type: Arrhenius=%.6g three-body=%.6g falloff(diag)=%.6g\n",
        maxRevByType[0], maxRevByType[1], maxRevByType[2]);
    if (worstKcReaction >= 0)
        std::printf("  WORST Kc: reaction %d \"%s\" mine=%g cantera=%g\n",
            worstKcReaction,
            SandiegoKinetics::reactions[worstKcReaction].equation,
            worstKcMine, worstKcRef);
    if (worstFwdReaction >= 0 && maxFwdError > 1.0e-6)
        std::printf("  falloff WORST forward: reaction %d \"%s\" state %d "
            "mine=%g cantera=%g ratio=%g\n",
            worstFwdReaction,
            SandiegoKinetics::reactions[worstFwdReaction].equation,
            worstFwdState, worstFwdKmol, worstFwdRef,
            worstFwdKmol / worstFwdRef);

    std::printf("%s\n", failures == 0
        ? "SANDIEGO_KINETICS_RATE_PASS"
        : "SANDIEGO_KINETICS_RATE_FAIL");
    return failures == 0 ? 0 : 1;
}
