/* San Diego rate-evaluation kernel (effective forward rate constants, Troe
   falloff, reverse rates via NASA-7 equilibrium constants).

   Concentrations are mol/m3; rate constants are on the mol,m,s basis as
   stored in SandiegoKineticsData.H. Forward rate of progress is
   k_eff * prod_i C_i^nu_i where nu is the explicit order set (falling back
   to the reactant stoichiometry), exactly as Cantera defines it. */
#ifndef SandiegoKineticsRate_H
#define SandiegoKineticsRate_H

#include "SandiegoKineticsData.H"
#include "SandiegoNasa7Data.H"

#include <cmath>
#include <stdexcept>

namespace Foam
{
namespace SandiegoKinetics
{

constexpr double Runiv = 8.314462618;       // J/(mol K)
constexpr double StandardPressure = 101325.0;  // Pa, Cantera OneAtm

inline double arrheniusRate(const ArrheniusParams& p, double T)
{
    return p.A * std::pow(T, p.n) * std::exp(-p.Ea / (Runiv * T));
}

// Cantera 3.1 TroeRate::F (standard Troe form):
//   logFcent = log10((1-a)*exp(-T/T3) + a*exp(-T/T1) [+ exp(-T2/T) if T2!=0])
//   cc = -0.4 - 0.67*logFcent;  nn = 0.75 - 1.27*logFcent
//   f1 = (log10(pr)+cc) / (nn - 0.14*(log10(pr)+cc))
//   F = 10^(logFcent / (1 + f1*f1))
inline double troeFactor(const SandiegoReactionData& r, double T, double pr)
{
    double fcent = (1.0 - r.troe[0]) * std::exp(-T / r.troe[1])
        + r.troe[0] * std::exp(-T / r.troe[2]);
    if (r.troe[3] != 0.0)
        fcent += std::exp(-r.troe[3] / T);
    const double logFcent = std::log10(std::max(fcent, 1.0e-300));
    const double cc = -0.4 - 0.67 * logFcent;
    const double nn = 0.75 - 1.27 * logFcent;
    const double lpr = std::log10(std::max(pr, 1.0e-300));
    const double f1 = (lpr + cc) / (nn - 0.14 * (lpr + cc));
    const double logF = logFcent / (1.0 + f1 * f1);
    return std::pow(10.0, logF);
}

// Effective forward rate constant for one reaction at T [K] and third-body
// concentration m [mol/m3] (ignored by plain Arrhenius).
inline double effectiveForwardRate(const SandiegoReactionData& r, double T, double m)
{
    if (r.type == 0)
        return arrheniusRate(r.high, T);
    if (r.type == 1)
        return arrheniusRate(r.high, T) * m;
    // falloff-Troe: kf = kinf * Pr/(1+Pr) * F = k0*[M]/(1+Pr)*F, exactly the
    // Cantera 3.1 FalloffRate::evalFromStruct form. Pr must be formed from
    // the kmol-basis reported constants (k0 and kinf carry different unit
    // orders), so the stored mol-basis values are converted back first.
    const double k0 = arrheniusRate(r.low, T);
    const double kinf = arrheniusRate(r.high, T);
    if (!(k0 > 0.0) || !(kinf > 0.0) || !(m >= 0.0))
        return 0.0;
    // A_mol = A_kmol*1000^(1-order), so the kmol-basis reported rate is
    // recovered by multiplying by 1000^(order-1).
    const double k0Kmol = k0 * std::pow(1000.0, r.low.order - 1.0);
    const double kinfKmol = kinf * std::pow(1000.0, r.high.order - 1.0);
    const double mKmol = m / 1000.0;
    const double pr = mKmol * k0Kmol / kinfKmol;
    const double F = troeFactor(r, T, pr);
    return kinf * pr / (1.0 + pr) * F;
}

// Third-body concentration [mol/m3] for reaction r given molar concentrations
// C[i] [mol/m3]; unlisted species have efficiency 1.
inline double thirdBodyConcentration(
    const SandiegoReactionData& r, const double* C)
{
    double m = 0.0;
    for (int i = 0; i < 59; ++i) m += C[i];
    for (int k = 0; k < r.nEfficiencies; ++k)
    {
        const int s = r.efficiencies[k].species;
        m += (r.efficiencies[k].coeff - 1.0) * C[s];
    }
    return m;
}

// Forward rate of progress [mol/m3/s]: k_eff * prod C^nu where nu is the
// explicit order set if present, otherwise the reactant stoichiometry.
inline double forwardRateOfProgress(
    const SandiegoReactionData& r, double T, const double* C)
{
    double result = effectiveForwardRate(r, T, thirdBodyConcentration(r, C));
    if (r.nOrders > 0)
    {
        for (int k = 0; k < r.nOrders; ++k)
        {
            const int s = r.orders[k].species;
            const double nu = r.orders[k].coeff;
            if (nu != 0.0)
                result *= std::pow(C[s], nu);
        }
    }
    else
    {
        for (int k = 0; k < r.nReactants; ++k)
        {
            const int s = r.reactants[k].species;
            const double nu = r.reactants[k].coeff;
            if (nu != 0.0)
                result *= std::pow(C[s], nu);
        }
    }
    return result;
}

// Standard molar Gibbs energy difference for a reaction at T [K]:
//   sum_i nu_i * (h_i - T s_i)   [J/mol]
inline double gibbsDelta(const SandiegoReactionData& r, double T)
{
    double result = 0.0;
    for (int k = 0; k < r.nProducts; ++k)
    {
        const int s = r.products[k].species;
        result += r.products[k].coeff
            * (SandiegoNasa7::hMole(s, T) - T * SandiegoNasa7::sMole(s, T));
    }
    for (int k = 0; k < r.nReactants; ++k)
    {
        const int s = r.reactants[k].species;
        result -= r.reactants[k].coeff
            * (SandiegoNasa7::hMole(s, T) - T * SandiegoNasa7::sMole(s, T));
    }
    return result;
}

// Concentration equilibrium constant Kc [mol,m,s basis] at T [K].
inline double equilibriumConstant(const SandiegoReactionData& r, double T)
{
    double nuSum = 0.0;
    for (int k = 0; k < r.nProducts; ++k)
        nuSum += r.products[k].coeff;
    for (int k = 0; k < r.nReactants; ++k)
        nuSum -= r.reactants[k].coeff;
    const double kp = std::exp(-gibbsDelta(r, T) / (Runiv * T));
    return kp * std::pow(StandardPressure / (Runiv * T), nuSum);
}

// Net rate of progress [mol/m3/s] with reverse rate kf/Kc.  forward already
// contains the reactant concentration product, so the reverse rate
// kr*prod(C_p^nu) is obtained by dividing it back out before applying the
// product side.
inline double netRateOfProgress(
    const SandiegoReactionData& r, double T, const double* C)
{
    const double forward = forwardRateOfProgress(r, T, C);
    if (!r.reversible) return forward;
    const double kc = equilibriumConstant(r, T);
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
    if (reactantProduct <= 0.0) return forward;   // no reactants -> no reverse
    double reverse = forward / reactantProduct / kc;
    for (int k = 0; k < r.nProducts; ++k)
    {
        const int s = r.products[k].species;
        const double nu = r.products[k].coeff;
        if (nu != 0.0)
            reverse *= std::pow(C[s], nu);
    }
    return forward - reverse;
}

} // namespace SandiegoKinetics
} // namespace Foam

#endif
