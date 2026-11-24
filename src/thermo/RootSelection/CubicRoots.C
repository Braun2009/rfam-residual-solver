/*---------------------------------------------------------------------------*\
  realFluidAllMach - Cubic Equation Solver Implementation
\*---------------------------------------------------------------------------*/

#include "CubicRoots.H"
#include <algorithm>

namespace Foam
{

int CubicRoots::solve(
    double a2,
    double a1,
    double a0,
    double& r1,
    double& r2,
    double& r3
)
{
    // Solve Z³ + a2*Z² + a1*Z + a0 = 0
    // Using Cardano's formula

    // Transform to depressed cubic: t³ + p*t + q = 0
    // where Z = t - a2/3
    double p = a1 - a2*a2/3.0;
    double q = a0 + 2.0*a2*a2*a2/27.0 - a1*a2/3.0;

    // Discriminant
    double disc = q*q/4.0 + p*p*p/27.0;

    if (disc > SMALL)
    {
        // One real root
        double u = std::cbrt(-q/2.0 + std::sqrt(disc));
        double v = std::cbrt(-q/2.0 - std::sqrt(disc));
        double t = u + v;
        r1 = t - a2/3.0;
        r2 = r1;
        r3 = r1;
        return 1;
    }
    else if (std::abs(disc) <= SMALL)
    {
        // Three real roots, at least two equal
        double t1, t2;
        if (std::abs(q) < SMALL)
        {
            t1 = 0.0;
            t2 = 0.0;
        }
        else
        {
            t1 = 3.0*q/p;
            t2 = -1.5*q/p;
        }

        r1 = t1 - a2/3.0;
        r2 = t2 - a2/3.0;
        r3 = t2 - a2/3.0;
        return 3;
    }
    else
    {
        // Three distinct real roots
        double r = std::sqrt(-p*p*p/27.0);
        double theta = std::acos(-q/(2.0*r));

        double t1 = 2.0 * std::cbrt(r) * std::cos(theta/3.0);
        double t2 = 2.0 * std::cbrt(r) * std::cos((theta + 2.0*PI)/3.0);
        double t3 = 2.0 * std::cbrt(r) * std::cos((theta + 4.0*PI)/3.0);

        r1 = t1 - a2/3.0;
        r2 = t2 - a2/3.0;
        r3 = t3 - a2/3.0;

        // Sort roots
        if (r1 > r2) std::swap(r1, r2);
        if (r2 > r3) std::swap(r2, r3);
        if (r1 > r2) std::swap(r1, r2);

        return 3;
    }
}

} // End namespace Foam
