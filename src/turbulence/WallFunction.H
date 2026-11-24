#ifndef WallFunction_H
#define WallFunction_H

/*---------------------------------------------------------------------------*\
    WallFunction - Spalding unified wall function for k-omega SST.

    Maps the first off-wall cell's velocity to wall shear stress τ_w, without
    requiring y+ < 1 resolution. Valid from viscous sublayer (y+ < 5) through
    the log layer (y+ > 30). This is the standard approach for coarse-grid RANS.

    Spalding's law (implicit, solved by Newton):
      y+ = u+ + exp(-κB) * [exp(κu+) - 1 - κu+ - (κu+)²/2 - (κu+)³/6]

    Given (y, U_tangential, ν), iteratively solve for u_τ:
      u+ = U/u_τ,  y+ = y·u_τ/ν  →  Spalding(u+) = y+

    Usage in TurbulenceOperator:
      double tau_w = WallFunction::wallShear(y_1, U_1, nu_lam, rho);
      // replaces the explicit du/dy gradient at the wall

    Also provides:
      omega_wall = WallFunction::omegaWall(y_1, u_tau, nu_lam);
      // Menter (2003) wall BC for omega: omega_w = 6ν/(β₁·y₁²)

    Pure std C++, header-only.
\*---------------------------------------------------------------------------*/

#include <cmath>
#include <algorithm>

namespace Foam
{

class WallFunction
{
private:
    static constexpr double kappa_ = 0.41;    // von Kármán constant
    static constexpr double B_     = 5.0;     // smooth-wall intercept
    static constexpr double beta1_ = 0.075;   // k-omega SST β₁

    // Spalding's law: y+(u+)
    static double spalding_yplus(double uplus)
    {
        double ku = kappa_ * uplus;
        return uplus + std::exp(-kappa_*B_)
             * (std::exp(ku) - 1.0 - ku - 0.5*ku*ku - ku*ku*ku/6.0);
    }

    // Derivative d(y+)/d(u+) for Newton iteration
    static double spalding_dyp_dup(double uplus)
    {
        double ku = kappa_ * uplus;
        return 1.0 + std::exp(-kappa_*B_) * kappa_
             * (std::exp(ku) - 1.0 - ku - 0.5*ku*ku);
    }

public:
    // ── Core: compute friction velocity u_τ from first-cell (y, U, ν) ──
    // Newton iteration on Spalding's law.
    static double frictionVelocity(double y, double U, double nu)
    {
        if (U < 1e-10 || y < 1e-12) return std::sqrt(nu * U / (y + 1e-30));

        // Initial guess: use log-law u_τ ≈ U / (1/κ·ln(y·U/ν) + B)
        // (crude but converges in 5-10 Newton steps)
        double u_tau = std::max(U * 0.05, std::sqrt(nu * U / y));  // fallback

        for (int iter = 0; iter < 20; ++iter)
        {
            double uplus = U / u_tau;
            double yplus_target = y * u_tau / nu;
            double yplus_spalding = spalding_yplus(uplus);
            double residual = yplus_spalding - yplus_target;

            // d(residual)/d(u_tau):
            //   d(yplus_spalding)/d(u_tau) = d(yp)/d(up) · d(up)/d(u_tau)
            //     = dyp_dup · (-U/u_tau²)
            //   d(yplus_target)/d(u_tau) = y/ν
            double dup_dutau = -U / (u_tau * u_tau);
            double dres = spalding_dyp_dup(uplus) * dup_dutau - y / nu;

            if (std::abs(dres) < 1e-30) break;
            double delta = -residual / dres;
            u_tau += delta;
            u_tau = std::max(u_tau, 1e-6);   // prevent negative
            if (std::abs(delta) < 1e-8 * u_tau) break;
        }
        return u_tau;
    }

    // ── Wall shear stress τ_w = ρ · u_τ² ──
    static double wallShear(double y, double U, double nu, double rho)
    {
        double u_tau = frictionVelocity(y, U, nu);
        return rho * u_tau * u_tau;
    }

    // ── Wall y+ for diagnostics ──
    static double yPlus(double y, double U, double nu)
    {
        double u_tau = frictionVelocity(y, U, nu);
        return y * u_tau / nu;
    }

    // ── Wall u+ for diagnostics ──
    static double uPlus(double y, double U, double nu)
    {
        double u_tau = frictionVelocity(y, U, nu);
        return U / u_tau;
    }

    // ── Menter (2003) wall BC for omega: ω_w = 6ν/(β₁·y²) ──
    // This is the asymptotic near-wall limit of the exact ω equation solution.
    static double omegaWall(double y, double nu)
    {
        return 6.0 * nu / (beta1_ * y * y);
    }

    // Kader's all-y+ temperature wall law. The viscous limit is T+=Pr*y+;
    // the outer limit approaches the logarithmic thermal layer correlation.
    static double temperaturePlus(double yplus,double Pr)
    {
        if(yplus<=0.0||Pr<=0.0) return 0.0;
        const double beta=std::pow(3.85*std::cbrt(Pr)-1.3,2)
                         +2.12*std::log(Pr);
        const double prY=Pr*yplus;
        const double gamma=0.01*std::pow(prY,4)
            /(1.0+5.0*std::pow(Pr,3)*yplus);
        const double viscousWeight=std::exp(-gamma);
        const double logWeight=std::exp(-1.0/std::max(gamma,1e-30));
        const double logLaw=2.12*std::log(std::max(yplus,1e-30))+beta;
        return std::max(Pr*yplus*viscousWeight+logLaw*logWeight,1e-30);
    }

    // Positive heat flux is from gas to wall. This exactly recovers
    // k*(Tgas-Twall)/y as y+ -> 0.
    static double thermalHeatFlux(double y,double tangentialSpeed,double nu,
                                  double rho,double cp,double conductivity,
                                  double Tgas,double Twall)
    {
        if(y<=0.0||nu<=0.0||rho<=0.0||cp<=0.0||conductivity<=0.0)
            return 0.0;
        if(tangentialSpeed<=1e-10)
            return conductivity*(Tgas-Twall)/y;
        const double uTau=frictionVelocity(y,std::max(tangentialSpeed,0.0),nu);
        const double yplus=y*uTau/nu;
        const double Pr=rho*nu*cp/conductivity;
        const double tplus=temperaturePlus(yplus,Pr);
        return rho*cp*uTau*(Tgas-Twall)/tplus;
    }
};

} // End namespace Foam
#endif
