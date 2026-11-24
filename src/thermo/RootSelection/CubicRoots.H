/*---------------------------------------------------------------------------*\
  realFluidAllMach - Cubic Equation Solver
\*---------------------------------------------------------------------------*/

#ifndef CubicRoots_H
#define CubicRoots_H

#include <cmath>

namespace Foam
{

class CubicRoots
{
public:
    // Solve Z³ + a2*Z² + a1*Z + a0 = 0
    // Returns number of real roots (1 or 3)
    static int solve(
        double a2,
        double a1,
        double a0,
        double& r1,
        double& r2,
        double& r3
    );

private:
    static constexpr double SMALL = 1e-14;
    static constexpr double PI = 3.14159265358979323846;
};

} // End namespace Foam

#endif
