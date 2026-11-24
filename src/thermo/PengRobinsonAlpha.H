#ifndef PengRobinsonAlpha_H
#define PengRobinsonAlpha_H

#include <cmath>

namespace Foam
{

struct PRAlphaDerivatives
{
    double a;
    double dadT;
    double d2adT2;
};

inline PRAlphaDerivatives evaluatePRAlpha(
    double Tc, double Pc, double omega, double T)
{
    constexpr double R_UNIVERSAL = 8.314462618;
    const double kappa =
        0.37464 + 1.54226*omega - 0.26992*omega*omega;
    const double sqrtTr = std::sqrt(T / Tc);
    const double f = 1.0 + kappa*(1.0 - sqrtTr);
    const double a0 =
        0.45724 * R_UNIVERSAL*R_UNIVERSAL * Tc*Tc / Pc;
    const double sqrtTTc = std::sqrt(T * Tc);

    PRAlphaDerivatives result;
    result.a = a0 * f * f;
    result.dadT = -a0 * kappa * f / sqrtTTc;
    result.d2adT2 =
        a0 * kappa * (1.0 + kappa) / (2.0 * T * sqrtTTc);
    return result;
}

} // End namespace Foam

#endif
