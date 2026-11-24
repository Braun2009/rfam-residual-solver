#ifndef WALE_H
#define WALE_H

/*---------------------------------------------------------------------------*\
    WALE - Wall-Adapting Local Eddy-viscosity subgrid-scale model
    (Nicoud & Ducros, "Subgrid-scale stress modelling based on the square of
    the velocity gradient tensor", Flow, Turbulence and Combustion 62:183-200,
    1999).

    nu_sgs = (Cw * delta)^2
             * (S^d_ij S^d_ij)^(3/2)
             / ( (S_ij S_ij)^(5/2) + (S^d_ij S^d_ij)^(5/4) )

    where g_ij = du_i/dx_j is the resolved velocity gradient tensor,
          S_ij = 0.5 (g_ij + g_ji) is the strain-rate tensor,
          G2_ij = g_ik g_kj is the square of the velocity gradient tensor,
          S^d_ij = G2_ij - (1/3) delta_ij G2_mm is the traceless symmetric
                   part of G2.

    Properties (used for validation):
      - pure shear:        nu_sgs > 0 and scales with (Cw*delta)^2 * |S|
      - rigid rotation:    nu_sgs = 0  (S^d = 0)
      - pure dilatation:   nu_sgs = 0  (S^d = 0)
      - near wall:         nu_sgs ~ O(y^3) (automatic wall behavior, no
                           damping function required)

    Cw = 0.325 is the standard value calibrated for decaying isotropic
    turbulence (Nicoud & Ducros 1999).

    Pure std C++, no dependencies. All quantities in SI.
\*---------------------------------------------------------------------------*/

#include <cmath>
#include <cstddef>

namespace Foam
{

struct WALE
{
    // Standard WALE constant (Nicoud & Ducros 1999).
    static constexpr double Cw = 0.325;

    // Filter width.  Default construction from cell volume:
    //   delta = V^(1/3)
    static double filterWidth(double cellVolume)
    {
        return std::cbrt(cellVolume);
    }

    // Subgrid-scale eddy viscosity [m2/s].
    //   g   : resolved velocity gradient tensor g[i][j] = du_i/dx_j [1/s]
    //   delta : filter width [m]
    // Returns 0 for pure rotation/dilatation and for any non-finite input.
    static double nuSgs(const double g[3][3], double delta)
    {
        if (!finiteTensor(g) || !(delta > 0.0)
            || !std::isfinite(delta))
            return 0.0;

        // Strain-rate tensor S_ij = 0.5 (g_ij + g_ji)
        double S[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                S[i][j] = 0.5 * (g[i][j] + g[j][i]);

        // Square of the velocity gradient tensor G2_ij = g_ik g_kj
        double G2[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
            {
                double sum = 0.0;
                for (int k = 0; k < 3; ++k)
                    sum += g[i][k] * g[k][j];
                G2[i][j] = sum;
            }

        const double traceG2 = G2[0][0] + G2[1][1] + G2[2][2];

        // Traceless symmetric part of G2: S^d_ij
        double Sd[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                Sd[i][j] = G2[i][j]
                    - (i == j ? (1.0 / 3.0) * traceG2 : 0.0);

        // S_ij S_ij
        double SS = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                SS += S[i][j] * S[i][j];

        // S^d_ij S^d_ij
        double SdSd = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                SdSd += Sd[i][j] * Sd[i][j];

        // nu_sgs = (Cw delta)^2 * SdSd^(3/2) / (SS^(5/2) + SdSd^(5/4))
        // Numerical guards: if either invariant is degenerate the model
        // returns zero (rigid rotation, pure dilatation, solid body).
        if (SdSd <= 0.0 || !std::isfinite(SdSd))
            return 0.0;

        const double numerator = SdSd * std::sqrt(SdSd);           // ^(3/2)
        const double strainTerm = SS * SS * std::sqrt(SS);         // ^(5/2)
        const double deviatorTerm = std::pow(SdSd, 1.25);          // ^(5/4)
        const double denominator = strainTerm + deviatorTerm;
        if (!(denominator > 0.0) || !std::isfinite(denominator))
            return 0.0;

        const double coefficient = Cw * delta;
        return coefficient * coefficient * numerator / denominator;
    }

    // SGS thermal conductivity via a constant turbulent Prandtl number:
    //   kappa_sgs = rho * cp * nu_sgs / Pr_t
    // Pr_t = 0.9 matches the legacy RANS thermal closure.
    static double kappaSgs(
        double rho, double cp, double nuSgs, double turbulentPrandtl = 0.9)
    {
        if (!(rho > 0.0) || !(cp > 0.0) || !(nuSgs > 0.0)
            || !(turbulentPrandtl > 0.0))
            return 0.0;
        return rho * cp * nuSgs / turbulentPrandtl;
    }

    // SGS species diffusivity via constant turbulent Schmidt number:
    //   D_sgs = nu_sgs / Sc_t
    // Sc_t = 1.0 is the standard first-cut for passive scalars.
    static double diffusivitySgs(double nuSgs, double turbulentSchmidt = 1.0)
    {
        if (!(nuSgs > 0.0) || !(turbulentSchmidt > 0.0))
            return 0.0;
        return nuSgs / turbulentSchmidt;
    }

private:
    static bool finiteTensor(const double g[3][3])
    {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (!std::isfinite(g[i][j]))
                    return false;
        return true;
    }
};

} // End namespace Foam

#endif
