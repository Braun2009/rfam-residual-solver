#ifndef TurbulenceOperator_H
#define TurbulenceOperator_H

/*---------------------------------------------------------------------------*\
    TurbulenceOperator - k-omega SST RANS as a framework Operator (Stage 6).

    Transport k and omega as two extra scalar fields stored in the Domain's
    turbulence state (held in this operator, NOT in ConservedVars — avoids
    touching the existing MVP/precond path).

    Per step (additive-update contract):
      1. Read k, omega from internal state (previous step).
      2. Compute nu_t, S, F1, F2, Pk, source terms.
      3. Explicit Euler advance k and omega (implicit option: future work).
      4. Clamp k >= k_min, omega >= omega_min.
      5. Diffusion: central differences for (nu+sigma_k nu_t) dk/dx etc.
      6. Add turbulent stress contribution -d(nu_t dU/dx)/dx to rhoU (Unew).
      NOTE: does NOT touch rhoE (turbulent heating = small, deferred to Stage7).

    Stability: dt_turb = min_cell(safety * k / max(Pk, eps)) — mirrors the
    chemistry limiter pattern; prevents blow-up of k during ignition.

    Wall distance y_wall is currently set to dy/2 (first-cell half-height) —
    sufficient for log-law resolved grids. Full wall-distance field is a
    follow-up ([[technical-debt]]).

    Pure std C++, header-only. Depends on KOmegaSST.H, Domain.H.
\*---------------------------------------------------------------------------*/

#include "Operator.H"
#include "Domain.H"
#include "KOmegaSST.H"
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <stdexcept>

namespace Foam
{

// Legacy Cartesian diagnostic operator. Its constructor deliberately accepts
// one constant molecular nu; production real-fluid sector calculations use
// SectorTurbulenceOperator with local TransportProps instead.
class TurbulenceOperator : public Operator
{
private:
    // Turbulent fields (not in ConservedVars — opt-in, zero impact on MVP)
    mutable std::vector<double> k_;
    mutable std::vector<double> omega_;
    mutable bool initialized_;
    mutable std::vector<double> kStepStart_, omegaStepStart_;
    mutable bool initializedStepStart_ = false;
    mutable bool stepTransactionActive_ = false;

    double nu_lam_;     // laminar kinematic viscosity [m^2/s]
    double k_inlet_;    // inlet turbulence intensity k [m^2/s^2]
    double om_inlet_;   // inlet specific dissipation omega [1/s]
    double safety_;     // explicit CFL factor for turb source limiter

    int idx(int i, int j, int Nx) const { return j * Nx + i; }

    void lazyInit(const Domain& d) const
    {
        if (initialized_) return;
        int total = d.Nx() * d.Ny();
        k_.assign(total, k_inlet_);
        omega_.assign(total, om_inlet_);
        initialized_ = true;
    }

    // strain-rate magnitude |du/dy| (2D, Cartesian; placeholder for full S_ij)
    //
    // TODO(turbulence): replace with S = sqrt(2 S_ij S_ij) from the velocity
    // gradient tensor. The consequence of the 1-D form is not a small error:
    // on a shear layer with both dudy and dvdx, S is underestimated, so Pk is
    // underestimated, so k decays too fast and the modelled eddy viscosity
    // fades in exactly the region turbulence most needs to be produced.
    // SectorTurbulenceOperator already carries the full tensor; the fix is to
    // share that implementation rather than add a second one here.
    static double strainRate(const Domain& d, int i, int j)
    {
        int Nx = d.Nx(), Ny = d.Ny();
        double dy = d.dy() > 0 ? d.dy() : d.dx();
        double dudy = 0.0;
        if (Ny > 1)
        {
            int jlo = std::max(j-1, 0), jhi = std::min(j+1, Ny-1);
            double dj = (jhi - jlo) * dy;
            if (dj > 0)
                dudy = (d.W(i,jhi).U[0] - d.W(i,jlo).U[0]) / dj;
        }
        return std::abs(dudy);
    }

public:
    TurbulenceOperator(double nu_lam,
                       double turbIntensity = 0.05,   // 5% freestream Tu
                       double lengthScale   = 0.001,  // 1 mm integral length
                       double safety        = 0.3)
        : initialized_(false), nu_lam_(nu_lam), safety_(safety)
    {
        // k = 1.5*(U*Tu)^2, omega = k^0.5 / (C_mu^0.25 * L)
        double U_ref = 100.0;   // typical order-of-magnitude [m/s]
        k_inlet_  = 1.5 * std::pow(U_ref * turbIntensity, 2);
        om_inlet_ = std::sqrt(k_inlet_) / (std::pow(KOmegaSST::beta_star, 0.25) * lengthScale);
    }

    bool ownsPersistentState() const override { return true; }

    void beginStepTransaction(const Domain& d) const override
    {
        if (stepTransactionActive_)
            throw std::logic_error("nested turbulence step transaction");
        kStepStart_ = k_;
        omegaStepStart_ = omega_;
        initializedStepStart_ = initialized_;
        stepTransactionActive_ = true;
        try
        {
            lazyInit(d);
        }
        catch (...)
        {
            rollbackStepTransaction();
            throw;
        }
    }

    void acceptStepTransaction() const noexcept override
    { stepTransactionActive_ = false; }

    void rollbackStepTransaction() const noexcept override
    {
        if (!stepTransactionActive_) return;
        k_.swap(kStepStart_);
        omega_.swap(omegaStepStart_);
        initialized_ = initializedStepStart_;
        stepTransactionActive_ = false;
    }

    void averageWithStepStartState() const override
    {
        if (!stepTransactionActive_ || k_.size() != omega_.size())
            throw std::logic_error("turbulence RK2 state transaction is invalid");

        if (initializedStepStart_ &&
            (k_.size() != kStepStart_.size()
             || omega_.size() != omegaStepStart_.size()))
            throw std::logic_error("turbulence RK2 state size changed");

        for (std::size_t i = 0; i < k_.size(); ++i)
        {
            const double startK = initializedStepStart_ ? kStepStart_[i] : k_inlet_;
            const double startOmega = initializedStepStart_
                ? omegaStepStart_[i] : om_inlet_;
            k_[i] = 0.5*(startK + k_[i]);
            omega_[i] = 0.5*(startOmega + omega_[i]);
        }
    }

    // Accessors for post-processing / testing.
    // Returns 0 if not yet initialized (apply() not yet called).
    double k(int i, int j, int Nx) const
    {
        if (!initialized_ || (int)k_.size() <= idx(i,j,Nx)) return k_inlet_;
        return k_[idx(i,j,Nx)];
    }
    double omega(int i, int j, int Nx) const
    {
        if (!initialized_ || (int)omega_.size() <= idx(i,j,Nx)) return om_inlet_;
        return omega_[idx(i,j,Nx)];
    }
    double nu_t(int i, int j, int Nx, const Domain& d) const
    {
        double kv  = k(i,j,Nx);
        double omv = omega(i,j,Nx);
        int id = idx(i,j,Nx);
        (void)id;
        double dy = d.dy() > 0 ? d.dy() : d.dx();
        double y  = std::max((j + 0.5) * dy, 1e-6);
        double S  = strainRate(d, i, j);
        double f2 = KOmegaSST::F2(kv, omv, y, nu_lam_);
        return KOmegaSST::nu_t(kv, omv, S, f2);
    }

    void apply(Domain& d, double dt) const override
    {
        if (d.Ny() < 2) return;   // 1D domain: turbulence disabled
        lazyInit(d);

        const int Nx = d.Nx(), Ny = d.Ny();
        const double dx = d.dx(), dy = d.dy();

        std::vector<double> k_new   = k_;
        std::vector<double> om_new  = omega_;

        for (int i = 1; i < Nx-1; ++i)
        for (int j = 1; j < Ny-1; ++j)
        {
            int id   = idx(i,   j,   Nx);
            int idxm = idx(i-1, j,   Nx);
            int idxp = idx(i+1, j,   Nx);
            int idym = idx(i,   j-1, Nx);
            int idyp = idx(i,   j+1, Nx);

            double kc  = k_[id],    oc  = omega_[id];
            double kxm = k_[idxm],  kxp = k_[idxp];
            double kym = k_[idym],  kyp = k_[idyp];
            double oxm = omega_[idxm], oxp = omega_[idxp];
            double oym = omega_[idym], oyp = omega_[idyp];

            double y  = std::max((j + 0.5) * dy, 1e-6);
            double S  = strainRate(d, i, j);
            double dkdx  = (k_[idxp] - k_[idxm]) / (2*dx);
            double dodx  = (omega_[idxp] - omega_[idxm]) / (2*dx);
            double dkdy  = (k_[idyp] - k_[idym]) / (2*dy);
            double dody  = (omega_[idyp] - omega_[idym]) / (2*dy);
            const double gradK[3]={dkdx,dkdy,0.0};
            const double gradOmega[3]={dodx,dody,0.0};
            const double CDkw=KOmegaSST::crossDiffusion(oc,gradK,gradOmega);
            double f1 = KOmegaSST::F1(kc, oc, y, nu_lam_,
                                      std::max(CDkw,1e-20));
            double f2 = KOmegaSST::F2(kc, oc, y, nu_lam_);
            double nut = KOmegaSST::nu_t(kc, oc, S, f2);
            double pk  = KOmegaSST::Pk(nut, S, kc, oc);

            // Source terms
            double sk  = KOmegaSST::src_k(pk, kc, oc);
            double som = KOmegaSST::src_omega(pk, nut, kc, oc, f1, CDkw);

            // Diffusion coefficients
            double sigk  = KOmegaSST::sigma_k(f1);
            double sigom = KOmegaSST::sigma_w(f1);
            double Dk  = nu_lam_ + sigk  * nut;
            double Dom = nu_lam_ + sigom * nut;

            // Laplacians (central differences)
            double lapk  = (kxp  - 2*kc + kxm)/(dx*dx) + (kyp  - 2*kc + kym)/(dy*dy);
            double lapom = (oxp  - 2*oc + oxm)/(dx*dx) + (oyp  - 2*oc + oym)/(dy*dy);

            k_new [id] = kc + dt * (sk  + Dk  * lapk);
            om_new[id] = oc + dt * (som + Dom * lapom);

            // Clamp
            k_new [id] = std::max(k_new [id], KOmegaSST::k_min);
            om_new[id] = std::max(om_new[id], KOmegaSST::om_min);

            // Turbulent stress contribution to rhoU: -d(nu_t du/dx)/dx (x-momentum)
            // For now only x-direction; full tensor is follow-up.
            if (d.Unew(i,j).rhoU.size() > 0)
            {
                double rho = d.U(i,j).rho;
                double dudx = (d.W(i+1,j).U[0] - d.W(i-1,j).U[0]) / (2*dx);
                double dudy = (d.W(i,j+1).U[0] - d.W(i,j-1).U[0]) / (2*dy);
                // tau_t = rho*nu_t * (du/dx + du/dy) (simplified diagonal + shear)
                double tau_xx = -2.0 * rho * nut * dudx;
                double tau_xy = -rho * nut * dudy;
                d.Unew(i,j).rhoU[0] += dt * (-tau_xx / dx - tau_xy / dy);
            }
        }

        k_     = k_new;
        omega_ = om_new;
    }

    double maxStableDt(const Domain& d, double CFL) const override
    {
        if (d.Ny() < 2) return 1e30;
        const int Nx = d.Nx(), Ny = d.Ny();
        const double dy2 = d.dy() * d.dy();
        double dtmin = 1e30;
        for (int i=1; i<Nx-1; ++i)
        for (int j=1; j<Ny-1; ++j)
        {
            int id = idx(i,j,Nx);
            const double kval = initialized_ ? k_[id] : k_inlet_;
            const double omegaValue = initialized_ ? omega_[id] : om_inlet_;
            double S   = strainRate(d, i, j);
            double f2  = KOmegaSST::F2(kval, omegaValue,
                std::max((j+0.5)*d.dy(),1e-6), nu_lam_);
            double nut = KOmegaSST::nu_t(kval, omegaValue, S, f2);
            double Deff = nu_lam_ + nut;
            if (Deff > 1e-30) dtmin = std::min(dtmin, safety_ * dy2 / (2*Deff));
        }
        return dtmin;
    }

    std::string name() const override
    { return "Legacy diagnostic k-omega SST (constant molecular nu)"; }
};

} // End namespace Foam
#endif
