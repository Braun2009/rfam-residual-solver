#ifndef MUSCLReconstruction_H
#define MUSCLReconstruction_H

/*---------------------------------------------------------------------------*\
    MUSCLReconstruction - second-order spatial reconstruction with slope limiter.

    Converts cell-averaged primitive variables (W_L, W_C, W_R) into left/right
    face states (W_L+, W_R−) for the Riemann solver, achieving second-order
    spatial accuracy while suppressing spurious oscillations near discontinuities.

    Limiter choices:
      - minmod:     most dissipative, very robust (default)
      - vanLeer:    moderate dissipation, smooth limiter
      - vanAlbada:  similar to vanLeer, differentiable
      - none:       unlimited (oscillatory at shocks, for testing)

    Usage in ConvectionOperator:
      MUSCLReconstruction recon;
      auto [WL_plus, WR_minus] = recon.reconstruct(W[i-1], W[i], W[i+1], W[i+2]);
      flux(WL_plus, WR_minus, normal, F);   // instead of flux(W[i], W[i+1], ...)

    This is a PRIMITIVE-VARIABLE reconstruction (p, T, u, Y). More robust for
    real-gas flows than conservative-variable reconstruction (no need to invert
    rhoE → T at the reconstructed state).
\*---------------------------------------------------------------------------*/

#include <cmath>
#include <algorithm>

namespace Foam
{

enum class Limiter { MINMOD, VANLEER, VANALBADA, NONE };

class MUSCLReconstruction
{
private:
    Limiter limiter_;

    // slope limiters: φ(r) where r = (U_C - U_L)/(U_R - U_C) [consecutive gradient ratio]
    static double phi(double r, Limiter lim)
    {
        if (r <= 0.0) return 0.0;   // opposite slopes → zero gradient (TVD)
        switch (lim)
        {
            case Limiter::MINMOD:
                return std::min(r, 1.0);
            case Limiter::VANLEER:
                return (r + std::abs(r)) / (1.0 + std::abs(r));
            case Limiter::VANALBADA:
                return (r*r + r) / (r*r + 1.0);
            case Limiter::NONE:
                return 1.0;   // unlimited: central slope (2nd order but oscillatory)
        }
        return std::min(r, 1.0);
    }

    // Reconstruct a single scalar: given (left, centre, right), return
    // the extrapolated value at the RIGHT face of the centre cell.
    static double reconRight(double L, double C, double R, Limiter lim)
    {
        double dR = R - C;
        double dL = C - L;
        if (std::abs(dR) < 1e-30) return C;   // flat → no extrapolation
        double r = dL / dR;
        return C + 0.5 * phi(r, lim) * dR;
    }

    // Reconstruct a single scalar: given (left, centre, right), return
    // the extrapolated value at the LEFT face of the centre cell.
    static double reconLeft(double L, double C, double R, Limiter lim)
    {
        double dL = C - L;
        double dR = R - C;
        if (std::abs(dL) < 1e-30) return C;
        double r = dR / dL;
        return C - 0.5 * phi(r, lim) * dL;
    }

public:
    explicit MUSCLReconstruction(Limiter lim = Limiter::MINMOD) : limiter_(lim) {}

    void setLimiter(Limiter lim) { limiter_ = lim; }
    Limiter limiter() const { return limiter_; }

    // Reconstruct face states between cells i and i+1.
    // Input: W at cells i-1, i, i+1, i+2 (four-point stencil).
    // Output: W_left (right face of cell i), W_right (left face of cell i+1).
    //
    // For 1D: call for each face pair.
    // For 2D: call per-direction (x or y) with the appropriate 4-cell stencil.
    struct FacePair
    {
        double p_L, T_L, u_L, v_L;
        double p_R, T_R, u_R, v_R;
        // species: handled separately in the operator (loop over Y_k)
    };

    // Scalar reconstruction helper (public for per-species use)
    double leftFace(double Lval, double Cval, double Rval) const
    { return reconRight(Lval, Cval, Rval, limiter_); }   // C extrapolated to right face

    double rightFace(double Lval, double Cval, double Rval) const
    { return reconLeft(Lval, Cval, Rval, limiter_); }    // C extrapolated to left face
};

} // End namespace Foam
#endif
