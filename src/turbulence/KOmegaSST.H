#ifndef KOmegaSST_H
#define KOmegaSST_H

/*---------------------------------------------------------------------------*\
    KOmegaSST - Menter (1994) k-omega SST two-equation RANS model.

    Blends k-omega (near-wall) with k-epsilon (free stream) via the SST
    blending function F1. Produces physically robust eddy viscosity for
    rocket combustor CFD.

    Transport equations (per cell, operator-split with convection):
        dk/dt    = P_k - beta* k omega + d/dx[(nu + sigma_k nu_t) dk/dx]
        domega/dt = alpha/nu_t P_k - beta omega^2 + CDkw + diff_omega

    Eddy viscosity (Bradshaw limiter in SST):
        nu_t = a1 k / max(a1 omega, S F2)
    where S = |2 S_ij| is the strain-rate magnitude.

    Production (Boussinesq):
        P_k = nu_t * S^2  (clipped at 20 beta* k omega to prevent excess build-up)

    Pure std C++, no dependencies. All quantities in SI.
\*---------------------------------------------------------------------------*/

// ---------------------------------------------------------------------------
// Integration notes for the RFAM real-fluid path (read before wiring this in).
//
// 1. This file is the CLOSURE ONLY: constants, blending functions, and the
//    source terms. It holds no state and no mesh. The two equations are
//    advanced by an Operator (TurbulenceOperator for the legacy Cartesian
//    path, SectorTurbulenceOperator for the sector path). Nothing here
//    constrains the conserved variables, which is the point -- k and omega are
//    carried OUTSIDE ConservedVars so that adding turbulence cannot perturb
//    the byte-for-byte MVP/precond path.
//
// 2. The production files (k, omega) are participant state of the same
//    accepted/rejected step transaction as U/W, via
//    ownsPersistentState()/beginStepTransaction()/rollbackStepTransaction()
//    on the Operator interface. Do NOT advance k/omega outside that protocol:
//    an accepted flow state paired with a rolled-back turbulence state is an
//    inconsistent field that no downstream audit can detect.
//
// 3. The real-fluid coupling is one-way in this implementation. The closure
//    reads local (nu, S) from the committed Q payload and returns nu_t; it does
//    not receive rho or the thermodynamics. That is adequate for a Boussinesq
//    eddy viscosity but it is NOT a variable-density or Favre-averaged
//    formulation. In the injector, where rho spans ~3x across the LOX/fuel
//    interface, the turbulent stress contribution computed this way is a
//    modelling choice, not a derived term. Record it as such in any report.
//
// 4. wallDistance() currently returns a geometric estimate (first-cell
//    half-height). F1 and F2 both have 500*nu/(y^2*omega) terms, so a wrong y
//    changes the near-wall blend directly. On the tetra mesh with 0.22 um wall
//    cells this is the single most suspect input to the model and must be a
//    real wall-distance field before any wall-heat-flux number is quoted.
//
// 5. Cross-diffusion is the one term that CANNOT be computed locally, because
//    it needs grad(k).grad(omega). This closure therefore takes
//    cdKwPositive as an argument instead of estimating it. Without it, F1 and
//    the omega source disagree; see crossDiffusion() for how the clip is
//    applied.
// ---------------------------------------------------------------------------

#include <cmath>
#include <algorithm>

namespace Foam
{

struct KOmegaSST
{
    // ---- SST model constants (Menter 1994) ----
    static constexpr double sigma_k1  = 0.85;
    static constexpr double sigma_k2  = 1.00;
    static constexpr double sigma_w1  = 0.50;
    static constexpr double sigma_w2  = 0.856;
    static constexpr double beta1     = 0.075;
    static constexpr double beta2     = 0.0828;
    static constexpr double beta_star = 0.09;
    static constexpr double alpha1    = 5.0/9.0;
    static constexpr double alpha2    = 0.44;
    static constexpr double a1        = 0.31;
    static constexpr double c1        = 10.0;   // F3 limiter (optional)

    static constexpr double k_min   = 1e-12;   // floor to prevent NaN
    static constexpr double om_min  = 1e-10;

    // Blended constants (F1=1 -> k-omega set 1, F1=0 -> k-epsilon set 2)
    static double sigma_k(double F1) { return F1*sigma_k1 + (1-F1)*sigma_k2; }
    static double sigma_w(double F1) { return F1*sigma_w1 + (1-F1)*sigma_w2; }
    static double beta   (double F1) { return F1*beta1    + (1-F1)*beta2;    }
    static double alpha  (double F1) { return F1*alpha1   + (1-F1)*alpha2;   }

    // ---- SST blending function F1 ----------------------------------------
    // F1 = 1 near wall (k-omega), F1 = 0 far field (k-epsilon).
    // y = wall distance, nu = kinematic viscosity.
    // cdKwPositive is max(2*sigma_w2/omega * grad(k).grad(omega), eps),
    // with units 1/s2.  It must be supplied by the spatial operator because
    // no local scalar-only approximation can represent cross diffusion.
    static double F1(double k, double omega, double y, double nu,
                     double cdKwPositive)
    {
        k = std::max(k, k_min);
        omega = std::max(omega, om_min);
        y = std::max(y, 1e-12);
        nu = std::max(nu, 0.0);
        const double CDkw = std::max(cdKwPositive, 1e-20);
        double arg1 = std::min(
            std::max(std::sqrt(k)/(beta_star*omega*y),
                     500.0*nu/(y*y*omega)),
            4.0*sigma_w2*k/(CDkw*y*y));
        return std::tanh(std::pow(arg1, 4));
    }

    // ---- F2 (Bradshaw limiter) -------------------------------------------
    static double F2(double k, double omega, double y, double nu)
    {
        k = std::max(k, k_min);
        omega = std::max(omega, om_min);
        y = std::max(y, 1e-12);
        nu = std::max(nu, 0.0);
        double arg2 = std::max(2.0*std::sqrt(k)/(beta_star*omega*y),
                               500.0*nu/(y*y*omega));
        return std::tanh(arg2 * arg2);
    }

    static double crossDiffusion(double omega,
                                 const double gradK[3],
                                 const double gradOmega[3])
    {
        omega=std::max(omega,om_min);
        const double dot=gradK[0]*gradOmega[0]
                        +gradK[1]*gradOmega[1]
                        +gradK[2]*gradOmega[2];
        // Menter's CD_{k omega} is clipped to a non-negative value.  The
        // same clipped quantity is used in F1 and in the omega source; using
        // the signed gradient product in the source would create an
        // unphysical anti-diffusive sink whenever grad(k).grad(omega)<0.
        return std::max(2.0*sigma_w2*dot/omega, 0.0);
    }

    // ---- Eddy viscosity [m^2/s] ------------------------------------------
    // S = |strain-rate| = sqrt(2 S_ij S_ij). In 1D/2D this is |du/dy|.
    static double nu_t(double k, double omega, double S, double f2val)
    {
        k     = std::max(k,     k_min);
        omega = std::max(omega, om_min);
        return a1 * k / std::max(a1 * omega, S * f2val);
    }

    // ---- Production of k [J/m^3/s / rho] = nu_t * S^2 -------------------
    // Clipped at 20 beta* k omega to prevent excessive build-up in stagnation.
    static double Pk(double nut, double S, double k, double omega)
    {
        double prod = nut * S * S;
        return std::min(prod, 20.0 * beta_star * k * std::max(omega, om_min));
    }

    // ---- Source for dk/dt [m^2/s^3]: P_k - beta* k omega ----------------
    static double src_k(double pk, double k, double omega)
    {
        return pk - beta_star * std::max(k, k_min) * std::max(omega, om_min);
    }

    // ---- Source for domega/dt [1/s^2]: alpha/nu_t P_k - beta omega^2 + CDkw
    // CDkw = cross-diffusion term: 2 sigma_w2 (1-F1) / omega * dk/dx * domega/dx
    // Simplified here (CDkw passed in by caller from cell gradient product).
    static double src_omega(double pk, double nut, double k, double omega,
                            double F1val, double CDkw)
    {
        (void)k;
        omega = std::max(omega, om_min);
        double a  = alpha(F1val);
        double b  = beta (F1val);
        double al = (nut > 1e-30) ? a * pk / nut : 0.0;
        double ds = b * omega * omega;
        double cd = (1.0 - F1val) * CDkw;
        return al - ds + cd;
    }
};

} // End namespace Foam
#endif
