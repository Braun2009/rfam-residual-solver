#include "PengRobinsonEOS.H"
#include "CubicRoots.H"
#include "PengRobinsonAlpha.H"
#include <cmath>
#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>

namespace Foam
{

static constexpr double R_UNIVERSAL = 8.314462618;

PengRobinsonEOS::PengRobinsonEOS(
    const SpeciesCriticalDB& criticalData,
    const BinaryInteractionDB& kij
)
:
    criticalData_(criticalData),
    kij_(kij),
    mixingRule_(criticalData, kij)
{
    // Set default species names
    std::vector<std::string> names;
    names.push_back("O2");
    names.push_back("C3H8");
    mixingRule_.setSpeciesNames(names);
}

PengRobinsonEOS::RootResult PengRobinsonEOS::solveZ(
    double p, double T, const std::vector<double>& Y, PhaseHint hint
) const
{
    RootResult result;

    if (!std::isfinite(p) || !std::isfinite(T) || p <= 0.0 || T <= 0.0 ||
        !mixingRule_.isValidComposition(Y))
        return result;

    double a, b;
    if (mixingRule_.multiSpeciesEnabled())
    {
        // Stage 3 step 3: true quadratic VdW mixing (a_mix = sum_i sum_j x_i x_j sqrt(ai aj)(1-kij))
        a  = mixingRule_.a_quad(Y, T);
        b  = mixingRule_.b_lin(Y);
    }
    else
    {
        // Legacy: pseudo-critical a/b from Tc_mix/Pc_mix/omega_mix (bit-identical MVP path)
        double Tc    = mixingRule_.Tc(Y);
        double Pc    = mixingRule_.Pc(Y);
        double omega = mixingRule_.omega(Y);
        a = evaluatePRAlpha(Tc, Pc, omega, T).a;
        b = 0.07780 * R_UNIVERSAL*Tc / Pc;
    }

    if (!std::isfinite(a) || !std::isfinite(b) || b <= 0.0)
        return result;

    double A = a * p / (R_UNIVERSAL*R_UNIVERSAL*T*T);
    double B = b * p / (R_UNIVERSAL * T);

    double a2 = B - 1.0;
    double a1 = A - 2.0*B - 3.0*B*B;
    double a0 = -A*B + B*B + B*B*B;

    double r1, r2, r3;
    int numRoots = CubicRoots::solve(a2, a1, a0, r1, r2, r3);

    result.numRoots = numRoots;

    if (numRoots != 1 && numRoots != 3) return result;

    if (numRoots == 1)
    {
        result.Z = r1;
        result.G = computeGibbs(r1, A, B, T);
        result.valid = std::isfinite(r1) && std::isfinite(result.G) &&
            (r1 > B * 1.01);
    }
    else
    {
        result.Z = selectRootGibbs(r1, r2, r3, numRoots, A, B, T, hint);
        result.G = computeGibbs(result.Z, A, B, T);
        result.valid = std::isfinite(result.Z) && std::isfinite(result.G) &&
            (result.Z > B * 1.01);
    }

    return result;
}

double PengRobinsonEOS::computeGibbs(double Z, double A, double B, double T) const
{
    const double invalid = std::numeric_limits<double>::infinity();
    if (!std::isfinite(Z) || !std::isfinite(A) || !std::isfinite(B) ||
        !std::isfinite(T) || B <= 0.0 || Z <= B * 1.01)
        return invalid;

    double RT = R_UNIVERSAL * T;
    double sqrt2 = std::sqrt(2.0);
    double term = (Z + (1.0 + sqrt2)*B) / (Z - (sqrt2 - 1.0)*B);

    if (!std::isfinite(term) || term <= 0.0) return invalid;

    double G = RT * (Z - 1.0 - std::log(Z - B) - A/(2.0*sqrt2*B)*std::log(term));
    return std::isfinite(G) ? G : invalid;
}

double PengRobinsonEOS::selectRootGibbs(
    double r1, double r2, double r3, int numRoots,
    double A, double B, double T, PhaseHint hint
) const
{
    const std::array<double, 3> roots = {r1, r2, r3};
    double selected = std::numeric_limits<double>::quiet_NaN();
    double best = hint == PhaseHint::GAS
        ? -std::numeric_limits<double>::infinity()
        : std::numeric_limits<double>::infinity();

    for (int i = 0; i < numRoots && i < 3; ++i)
    {
        const double root = roots[static_cast<size_t>(i)];
        if (!std::isfinite(root) || root <= B * 1.01) continue;

        if (hint == PhaseHint::GAS)
        {
            if (root > best) { best = root; selected = root; }
        }
        else if (hint == PhaseHint::LIQUID)
        {
            if (root < best) { best = root; selected = root; }
        }
        else
        {
            const double gibbs = computeGibbs(root, A, B, T);
            if (gibbs < best) { best = gibbs; selected = root; }
        }
    }

    return selected;
}

double PengRobinsonEOS::density(double p, double T, const std::vector<double>& Y) const
{
    auto result = solveZ(p, T, Y);
    if (!result.valid || !std::isfinite(result.Z) || result.Z <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    double MW = mixingRule_.MW(Y);
    double rho = p * MW / (result.Z * R_UNIVERSAL * T);
    return std::isfinite(rho) && rho > 0.0
        ? rho : std::numeric_limits<double>::quiet_NaN();
}

double PengRobinsonEOS::pressure(double rho, double T, const std::vector<double>& Y) const
{
    if (!std::isfinite(rho) || !std::isfinite(T) ||
        rho <= 0.0 || T <= 0.0 || !mixingRule_.isValidComposition(Y))
        return std::numeric_limits<double>::quiet_NaN();

    double MW = mixingRule_.MW(Y);
    double a, b;

    if (mixingRule_.multiSpeciesEnabled())
    {
        a = mixingRule_.a_quad(Y, T);
        b = mixingRule_.b_lin(Y);
    }
    else
    {
        double Tc    = mixingRule_.Tc(Y);
        double Pc    = mixingRule_.Pc(Y);
        double omega = mixingRule_.omega(Y);
        a = evaluatePRAlpha(Tc, Pc, omega, T).a;
        b = 0.07780 * R_UNIVERSAL*Tc / Pc;
    }

    double v = MW / rho;
    if (!std::isfinite(MW) || !std::isfinite(a) || !std::isfinite(b) ||
        !std::isfinite(v) || b <= 0.0 || v <= b)
        return std::numeric_limits<double>::quiet_NaN();

    double p = R_UNIVERSAL * T / (v - b) - a / (v*v + 2.0*b*v - b*b);

    return std::isfinite(p) ? p : std::numeric_limits<double>::quiet_NaN();
}

double PengRobinsonEOS::enthalpyDeparture(
    double p, double T, double Z, const std::vector<double>& Y
) const
{
    double MW = mixingRule_.MW(Y);
    double a, b, dadT;

    if (mixingRule_.multiSpeciesEnabled())
    {
        a    = mixingRule_.a_quad(Y, T);
        b    = mixingRule_.b_lin(Y);
        dadT = mixingRule_.dadT_quad(Y, T);
    }
    else
    {
        double Tc    = mixingRule_.Tc(Y);
        double Pc    = mixingRule_.Pc(Y);
        double omega = mixingRule_.omega(Y);
        const PRAlphaDerivatives alpha = evaluatePRAlpha(Tc, Pc, omega, T);
        b    = 0.07780 * R_UNIVERSAL*Tc / Pc;
        a    = alpha.a;
        dadT = alpha.dadT;
    }

    if (!std::isfinite(MW) || !std::isfinite(a) || !std::isfinite(b) ||
        !std::isfinite(dadT) || !std::isfinite(p) || !std::isfinite(T) ||
        !std::isfinite(Z) || MW <= 0.0 || b <= 0.0 || p <= 0.0 || T <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    double B = b * p / (R_UNIVERSAL * T);

    double sqrt2 = std::sqrt(2.0);
    double term = (Z + (1.0 + sqrt2)*B) / (Z - (sqrt2 - 1.0)*B);
    if (!std::isfinite(term) || term <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    // h_dep = R*T*(Z-1) + (T*da/dT - a)/(2*sqrt(2)*b) * ln(term)
    double h_dep = R_UNIVERSAL * T * (Z - 1.0)
                 + (T * dadT - a) / (2.0 * sqrt2 * b) * std::log(term);

    return h_dep / MW;
}

double PengRobinsonEOS::cvDeparture(
    double p, double T, double Z, const std::vector<double>& Y
) const
{
    double MW = mixingRule_.MW(Y);
    double b, d2adT2;

    if (mixingRule_.multiSpeciesEnabled())
    {
        b       = mixingRule_.b_lin(Y);
        d2adT2  = mixingRule_.d2adT2_quad(Y, T);
    }
    else
    {
        double Tc    = mixingRule_.Tc(Y);
        double Pc    = mixingRule_.Pc(Y);
        double omega = mixingRule_.omega(Y);
        const PRAlphaDerivatives alpha = evaluatePRAlpha(Tc, Pc, omega, T);
        b       = 0.07780 * R_UNIVERSAL*Tc / Pc;
        d2adT2  = alpha.d2adT2;
    }

    if (!std::isfinite(MW) || !std::isfinite(b) ||
        !std::isfinite(d2adT2) || !std::isfinite(p) || !std::isfinite(T) ||
        !std::isfinite(Z) || MW <= 0.0 || b <= 0.0 || p <= 0.0 || T <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    double B = b * p / (R_UNIVERSAL * T);

    double sqrt2 = std::sqrt(2.0);
    double term = (Z + (1.0 + sqrt2)*B) / (Z - (sqrt2 - 1.0)*B);
    if (!std::isfinite(term) || term <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    // cv_dep = T*d2a/dT2/(2*sqrt(2)*b) * ln(term)
    double cv_dep = T * d2adT2 / (2.0 * sqrt2 * b) * std::log(term);
    return cv_dep / MW;
}

double PengRobinsonEOS::entropyDeparture(
    double p, double T, double Z, const std::vector<double>& Y
) const
{
    const double MW = mixingRule_.MW(Y);
    double b = std::numeric_limits<double>::quiet_NaN();
    double dadT = std::numeric_limits<double>::quiet_NaN();

    if (mixingRule_.multiSpeciesEnabled())
    {
        b = mixingRule_.b_lin(Y);
        dadT = mixingRule_.dadT_quad(Y, T);
    }
    else
    {
        const double Tc = mixingRule_.Tc(Y);
        const double Pc = mixingRule_.Pc(Y);
        const double omega = mixingRule_.omega(Y);
        const PRAlphaDerivatives alpha = evaluatePRAlpha(Tc, Pc, omega, T);
        b = 0.07780*R_UNIVERSAL*Tc/Pc;
        dadT = alpha.dadT;
    }

    if (!std::isfinite(MW) || !std::isfinite(b) ||
        !std::isfinite(dadT) || !std::isfinite(p) || !std::isfinite(T) ||
        !std::isfinite(Z) || MW <= 0.0 || b <= 0.0 || p <= 0.0 || T <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    const double B = b*p/(R_UNIVERSAL*T);
    const double sqrt2 = std::sqrt(2.0);
    const double term = (Z + (1.0 + sqrt2)*B)
        /(Z - (sqrt2 - 1.0)*B);
    if (!std::isfinite(B) || B <= 0.0 || !std::isfinite(term) ||
        term <= 0.0 || Z <= B)
        return std::numeric_limits<double>::quiet_NaN();

    // s_dep = R ln(Z-B) + da/dT/(2 sqrt(2) b) ln(term).
    // This convention is paired with enthalpyDeparture and computeGibbs:
    // h_dep - g_dep = T*s_dep.
    const double entropy = R_UNIVERSAL*std::log(Z - B)
        + dadT/(2.0*sqrt2*b)*std::log(term);
    return std::isfinite(entropy)
        ? entropy/MW : std::numeric_limits<double>::quiet_NaN();
}

bool PengRobinsonEOS::isNearCritical(double p, double T, const std::vector<double>& Y) const
{
    double Tc = mixingRule_.Tc(Y);
    double Pc = mixingRule_.Pc(Y);
    double Tr = T / Tc;
    double Pr = p / Pc;
    return (Tr > Tr_lower && Tr < Tr_upper && Pr > Pr_lower && Pr < Pr_upper);
}

double PengRobinsonEOS::pseudoCriticalTemperature(
    const std::vector<double>& Y) const
{
    return mixingRule_.Tc(Y);
}

double PengRobinsonEOS::pseudoCriticalPressure(
    const std::vector<double>& Y) const
{
    return mixingRule_.Pc(Y);
}

double PengRobinsonEOS::MW(const std::vector<double>& Y) const
{
    return mixingRule_.MW(Y);
}

} // End namespace Foam
