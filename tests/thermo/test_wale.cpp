// WALE subgrid-scale model focused contract tests.
//
// Verifies the physical properties of the WALE closure (Nicoud & Ducros
// 1999):
//   - a two-dimensional pure shear flow produces EXACTLY zero SGS viscosity
//     (the model is built on the traceless symmetric part of the square of
//     the velocity gradient tensor, which vanishes for 2-D shear);
//   - a genuine three-dimensional straining/shearing field produces a
//     positive SGS viscosity that scales linearly with the shear magnitude
//     and quadratically with the filter width;
//   - pure dilatation produces exactly zero SGS viscosity;
//   - degenerate/non-finite inputs fail closed to zero;
//   - the SGS thermal conductivity and diffusivity helpers are finite and
//     positive for valid states.
//
// Output marker: WALE_CONTRACT_PASS / WALE_CONTRACT_FAIL.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "WALE.H"

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

double relativeError(double got, double expected, double floor = 1.0e-10)
{
    return std::abs(got - expected)
        / std::max(std::abs(expected), floor);
}

void zeroGradient(double g[3][3])
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            g[i][j] = 0.0;
}

// Two-dimensional pure shear: u = S*y, v = 0, w = 0 (only du/dy = S).
// The square of this gradient tensor is identically zero, hence nu_sgs = 0.
void twoDShearGradient(double S, double g[3][3])
{
    zeroGradient(g);
    g[0][1] = S;
}

// Three-dimensional shearing field (dominant du/dy plus a weaker secondary
// shear dv/dx): the gradient square is nonzero and nu_sgs > 0.
void threeDShearGradient(double S, double g[3][3])
{
    zeroGradient(g);
    g[0][1] = S;
    g[1][0] = 0.5 * S;
    g[0][2] = 0.25 * S;
}

// Pure dilatation: u = a*x, v = a*y, w = a*z.
void dilatationGradient(double a, double g[3][3])
{
    zeroGradient(g);
    g[0][0] = a;
    g[1][1] = a;
    g[2][2] = a;
}

int run()
{
    double g[3][3];
    const double delta = 1.0e-4;

    // 1. 2-D pure shear: exactly zero SGS viscosity (WALE signature).
    twoDShearGradient(1000.0, g);
    check("2-D shear: nu_sgs == 0", WALE::nuSgs(g, delta) == 0.0,
          WALE::nuSgs(g, delta), 0.0);

    // 2. 3-D shearing field: positive SGS viscosity.
    threeDShearGradient(1000.0, g);
    const double nu3d = WALE::nuSgs(g, delta);
    check("3-D shear: nu_sgs > 0", nu3d > 0.0, nu3d);

    // 3. Linear scaling with shear magnitude S (fixed tensor structure).
    double g2[3][3];
    threeDShearGradient(2000.0, g2);
    const double nu3d2 = WALE::nuSgs(g2, delta);
    const double ratio = nu3d2 / nu3d;
    check("3-D shear: linear scaling with |S|", relativeError(ratio, 2.0) < 1.0e-6,
          ratio, 2.0);

    // 4. Quadratic scaling with filter width delta.
    const double nu3dDelta2 = WALE::nuSgs(g, 2.0 * delta);
    const double deltaRatio = nu3dDelta2 / nu3d;
    check("3-D shear: delta^2 scaling", relativeError(deltaRatio, 4.0) < 1.0e-6,
          deltaRatio, 4.0);

    // 5. Pure dilatation: exactly zero.
    dilatationGradient(300.0, g);
    check("pure dilatation: nu_sgs == 0", WALE::nuSgs(g, delta) == 0.0,
          WALE::nuSgs(g, delta), 0.0);

    // 6. Rotational 2-D field (u = -omega y, v = omega x): the WALE kernel
    //    is built on the traceless symmetric square of the gradient, which
    //    does not vanish for this field; assert finite non-negative output.
    zeroGradient(g);
    g[0][1] = -500.0; g[1][0] = 500.0;
    const double nuRot = WALE::nuSgs(g, delta);
    check("rigid rotation: finite non-negative", std::isfinite(nuRot) && nuRot >= 0.0, nuRot);

    // 7. Non-finite input fails closed to zero.
    threeDShearGradient(1000.0, g);
    g[1][2] = std::numeric_limits<double>::quiet_NaN();
    check("NaN gradient: nu_sgs == 0", WALE::nuSgs(g, delta) == 0.0);
    threeDShearGradient(1000.0, g);
    check("zero delta: nu_sgs == 0", WALE::nuSgs(g, 0.0) == 0.0);
    check("negative delta: nu_sgs == 0", WALE::nuSgs(g, -1.0e-4) == 0.0);

    // 8. Filter width from cell volume.
    check("filterWidth: cbrt(V)", relativeError(WALE::filterWidth(1.0e-9), 1.0e-3) < 1.0e-12);

    // 9. SGS thermal conductivity and diffusivity helpers.
    threeDShearGradient(1000.0, g);
    const double nuS = WALE::nuSgs(g, delta);
    const double kappa = WALE::kappaSgs(1.0, 1000.0, nuS);
    check("kappaSgs: finite positive", kappa > 0.0 && std::isfinite(kappa), kappa);
    const double dSgs = WALE::diffusivitySgs(nuS);
    check("diffusivitySgs: equals nu_sgs/Sc_t",
          relativeError(dSgs, nuS) < 1.0e-12, dSgs, nuS);
    check("kappaSgs: invalid input -> 0", WALE::kappaSgs(0.0, 1000.0, nuS) == 0.0);
    check("diffusivitySgs: invalid input -> 0", WALE::diffusivitySgs(-1.0) == 0.0);

    // 10. Monotone in filter width.
    threeDShearGradient(1000.0, g);
    const double nuSmall = WALE::nuSgs(g, 1.0e-5);
    const double nuLarge = WALE::nuSgs(g, 1.0e-3);
    check("monotone in delta", nuLarge > nuSmall, nuLarge, nuSmall);

    std::cout << (failures == 0 ? "WALE_CONTRACT_PASS\n"
                                : "WALE_CONTRACT_FAIL\n");
    return failures == 0 ? 0 : 1;
}

} // namespace

int main()
{
    return run();
}
