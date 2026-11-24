#ifndef HybridHelmholtzThermo_H
#define HybridHelmholtzThermo_H

/*---------------------------------------------------------------------------*\
    HybridHelmholtzThermo - unified NASA7-ideal + faded-residual Helmholtz
    potential (GERG adoption contract, 2026-08-02).

    a(T,rho,Y) = a_IG_NASA(T,rho,Y) + f(T) * a_res(T,rho,Y)

    - a_IG_NASA : ideal-gas Helmholtz from NASA7 sensible caloric data
                  (single authoritative ideal baseline + formation-energy
                  reference for RFAM);
    - a_res     : residual Helmholtz (PR placeholder; GERG-2008 for the
                  stable set in production);
    - f(T)      : C2 quintic smoothstep residual fade (1 in validated
                  low-T/dense/VLE region, 0 in high-T ideal-gas region,
                  f'=f''=0 at both ends).

    EVERY derived property (p, e, h, s, Cp, Cv, sound speed, ...) comes
    from differentiating this one potential (Mayer/Maxwell hold by
    construction).  Property blending is forbidden.

    Single-phase isochoric evaluate() fails CLOSED inside the PR two-phase
    dome (the isochoric cubic has only a spurious vapour-like root there).
    The isobaric path (isobaricRoots/bubblePressure/prLnFugacity) is the
    correct PR phase-boundary interface and is validated against CoolProp PR.

    Units: mass basis (J/kg).
\*---------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "GERGResidualTable.H"
#include "MixingRules/SpeciesCriticalDB.H"
#include "NASA7Thermo.H"
#include "SpeciesRegistry.H"

namespace Foam
{

struct ResidualFade
{
    double T0 = 600.0;
    double T1 = 900.0;

    double fade(double T) const
    {
        if (T <= T0) return 1.0;
        if (T >= T1) return 0.0;
        const double xi = (T - T0) / (T1 - T0);
        return 1.0 - (6.0*xi*xi*xi*xi*xi - 15.0*xi*xi*xi*xi
                      + 10.0*xi*xi*xi);
    }
    double derivative(double T) const
    {
        if (T <= T0 || T >= T1) return 0.0;
        const double xi = (T - T0) / (T1 - T0);
        const double dxi = 1.0 / (T1 - T0);
        return -(30.0*xi*xi*xi*xi - 60.0*xi*xi*xi + 30.0*xi*xi) * dxi;
    }
    double secondDerivative(double T) const
    {
        if (T <= T0 || T >= T1) return 0.0;
        const double xi = (T - T0) / (T1 - T0);
        const double dxi = 1.0 / (T1 - T0);
        return -(120.0*xi*xi*xi - 180.0*xi*xi + 60.0*xi) * dxi * dxi;
    }
};

class HybridHelmholtzThermo
{
public:
    enum class Status
    {
        SUCCESS = 0,
        INVALID_INPUT = 1,
        NONPHYSICAL_STATE = 2,
        NOT_CONVERGED = 3
    };

    struct IdealGasState
    {
        double e = std::numeric_limits<double>::quiet_NaN();
        double h = std::numeric_limits<double>::quiet_NaN();
        double s = std::numeric_limits<double>::quiet_NaN();
        double cv = std::numeric_limits<double>::quiet_NaN();
        double cp = std::numeric_limits<double>::quiet_NaN();
    };

    struct ResidualState
    {
        double a_res = std::numeric_limits<double>::quiet_NaN();
        double da_res_dT = std::numeric_limits<double>::quiet_NaN();
        double d2a_res_dT2 = std::numeric_limits<double>::quiet_NaN();
        double da_res_drho = std::numeric_limits<double>::quiet_NaN();
        double d2a_res_drho2 = std::numeric_limits<double>::quiet_NaN();
        double d2a_res_dTdrho = std::numeric_limits<double>::quiet_NaN();
    };

    struct State
    {
        Status status = Status::INVALID_INPUT;
        double T = std::numeric_limits<double>::quiet_NaN();
        double rho = std::numeric_limits<double>::quiet_NaN();
        double p = std::numeric_limits<double>::quiet_NaN();
        double e = std::numeric_limits<double>::quiet_NaN();
        double h = std::numeric_limits<double>::quiet_NaN();
        double s = std::numeric_limits<double>::quiet_NaN();
        double cv = std::numeric_limits<double>::quiet_NaN();
        double cp = std::numeric_limits<double>::quiet_NaN();
        double soundSpeed = std::numeric_limits<double>::quiet_NaN();
        double dpdrho_T = std::numeric_limits<double>::quiet_NaN();
        double dpdT_rho = std::numeric_limits<double>::quiet_NaN();
        double fade = std::numeric_limits<double>::quiet_NaN();
        double residualMagnitude = std::numeric_limits<double>::quiet_NaN();

        bool success() const { return status == Status::SUCCESS; }
    };

    struct IsobaricRoots
    {
        double rhoLiquidKgM3 = std::numeric_limits<double>::quiet_NaN();
        double rhoVaporKgM3 = std::numeric_limits<double>::quiet_NaN();
        double ZLiquid = std::numeric_limits<double>::quiet_NaN();
        double ZVapor = std::numeric_limits<double>::quiet_NaN();
        bool twoPhase = false;
        bool valid = false;
    };

    HybridHelmholtzThermo(const SpeciesRegistry& registry,
                          double fadeT0 = 600.0,
                          double fadeT1 = 900.0)
        : registry_(registry),
          criticalDB_(SpeciesCriticalDB::createDefault()),
          fadeT0_(fadeT0), fadeT1_(fadeT1) {}

    // Register an external GERG-2008 residual table.  When registered, the
    // residual path uses the table (defined for every rho, including the
    // dense liquid branch); otherwise the PR placeholder is used.
    void registerResidualTable(const GERGResidualTable* table)
    { residualTable_ = table; }

    ResidualFade fade() const
    {
        ResidualFade f;
        f.T0 = fadeT0_;
        f.T1 = fadeT1_;
        return f;
    }

    // ------------------------------------------------------------------ //
    // Ideal-gas Helmholtz from NASA7 (mass basis).  s includes the mixing
    // entropy and the density term -Rmix ln(rho/rho_ref) with rho_ref=1,
    // which keeps e = eIG(T) (rho-independent) and cv = cvIG.
    // ------------------------------------------------------------------ //
    IdealGasState idealGasState(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double T, double rho) const
    {
        IdealGasState out;
        if (names.size() != Y.size()) return out;
        double e = 0.0, s = 0.0, cv = 0.0, cp = 0.0;
        double Rmix = 0.0;
        for (std::size_t i = 0; i < names.size(); ++i)
        {
            if (Y[i] <= 0.0) continue;
            const double Mi = registry_.molecularWeight(
                registry_.index(names[i]));
            const double Ri = R_UNIVERSAL / Mi;
            const double hSens = nasa7Thermo().hsens_mole(names[i], T) / Mi;
            const double eSens = hSens - Ri * T;
            const double sSens = nasa7Thermo().s_mole(names[i], T) / Mi;
            const double cpi = nasa7Thermo().cp_mole(names[i], T) / Mi;
            e += Y[i] * eSens;
            s += Y[i] * (sSens - Ri * std::log(Y[i]));
            cv += Y[i] * (cpi - Ri);
            cp += Y[i] * cpi;
            Rmix += Y[i] * Ri;
        }
        s -= Rmix * std::log(rho / 1.0);
        out.e = e;
        out.cv = cv;
        out.cp = cp;
        out.s = s;
        out.h = e + Rmix * T;
        return out;
    }

    // ------------------------------------------------------------------ //
    // PR residual Helmholtz (placeholder for GERG residual), mass basis.
    // ------------------------------------------------------------------ //
    ResidualState prResidual(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double T, double rho) const
    {
        ResidualState out;
        if (names.size() != Y.size() || T <= 0.0 || rho <= 0.0) return out;
        // External GERG residual table (pure/near-pure, or O2/C3H8 binary
        // mixture at the nearest fixed-composition table level): defined for
        // every rho (including the dense liquid branch).
        if (residualTable_)
        {
            const std::string pure = activePureSpecies(names, Y);
            if (!pure.empty() && residualTable_->has(pure))
            {
                const GERGResidualTableResult res =
                    residualTable_->evaluate(pure, T, rho);
                if (res.success() && res.properties)
                {
                    out.a_res = res.properties->aRes;
                    out.da_res_dT = res.properties->daResDT;
                    out.d2a_res_dT2 = res.properties->d2aResDT2;
                    out.da_res_drho = res.properties->daResDrho;
                    out.d2a_res_drho2 = res.properties->d2aResDrho2;
                    out.d2a_res_dTdrho = res.properties->d2aResDTdrho;
                    return out;
                }
            }
            else
            {
                // O2/C3H8 binary mixture at the nearest table composition
                // level (GERG reducing-function mix, no departure terms in
                // the torch-flash parameter set).  Only the two stable
                // species participate; radicals are excluded by contract.
                const std::string mix = binaryMixKey(names, Y);
                if (!mix.empty() && residualTable_->has(mix))
                {
                    const GERGResidualTableResult res =
                        residualTable_->evaluate(mix, T, rho);
                    if (res.success() && res.properties)
                    {
                        out.a_res = res.properties->aRes;
                        out.da_res_dT = res.properties->daResDT;
                        out.d2a_res_dT2 = res.properties->d2aResDT2;
                        out.da_res_drho = res.properties->daResDrho;
                        out.d2a_res_drho2 = res.properties->d2aResDrho2;
                        out.d2a_res_dTdrho = res.properties->d2aResDTdrho;
                        return out;
                    }
                }
            }
            // Mixed composition with no table coverage: fail closed (do not
            // fall back silently to the PR placeholder).
            return out;
        }
        const double MW = mixtureMw(names, Y);
        const double rhoMolar = rho / MW;
        const double a = prAMix(names, Y, T);
        const double b = prBMix(names, Y);
        const double B = b * pIdeal(T, rhoMolar) / (R_UNIVERSAL*T);
        // Single-phase isochoric PR root.  In the dense two-phase dome the
        // isochoric cubic has only a spurious vapour-like root; the PR
        // residual is not single-phase-defined there.  Fail closed; the
        // two-phase flash (or the GERG residual, defined for every rho)
        // replaces this.
        const double Z = prRoot(names, Y, T, rho);
        if (!std::isfinite(Z) || Z <= 0.0) return out;
        if (rho > 300.0 && Z > 0.5) return out;  // two-phase dome
        const double sqrt2 = std::sqrt(2.0);
        const double term = (Z + (1.0+sqrt2)*B) / (Z + (1.0-sqrt2)*B);
        if (!std::isfinite(term) || term <= 0.0) return out;

        const double aResMolar = R_UNIVERSAL*T*(-std::log(Z-B))
            + a/(2.0*sqrt2*b)*std::log(term);
        out.a_res = aResMolar / MW;

        // Numeric derivatives on the selected (isochoric) root.
        {
            const double dT = 1.0e-4 * T;
            const double aTp = prResidualValue(names, Y, T+dT, rho, MW);
            const double aTm = prResidualValue(names, Y, T-dT, rho, MW);
            if (!std::isfinite(aTp) || !std::isfinite(aTm)) return out;
            out.da_res_dT = (aTp - aTm)/(2.0*dT)/MW;
            out.d2a_res_dT2 = (aTp - 2.0*aResMolar + aTm)/(dT*dT)/MW;
        }
        {
            const double drho = 1.0e-6 * rho;
            const double aRp = prResidualValue(names, Y, T, rho+drho, MW);
            const double aRm = prResidualValue(names, Y, T, rho-drho, MW);
            if (!std::isfinite(aRp) || !std::isfinite(aRm)) return out;
            out.da_res_drho = (aRp - aRm)/(2.0*drho)/MW;
            out.d2a_res_drho2 = (aRp - 2.0*aResMolar + aRm)/(drho*drho)/MW;
        }
        {
            const double dT = 1.0e-4 * T, drho = 1.0e-6 * rho;
            out.d2a_res_dTdrho =
                (prResidualValue(names, Y, T+dT, rho+drho, MW)
                 - prResidualValue(names, Y, T-dT, rho+drho, MW)
                 - prResidualValue(names, Y, T+dT, rho-drho, MW)
                 + prResidualValue(names, Y, T-dT, rho-drho, MW))
                /(4.0*dT*drho)/MW;
        }
        return out;
    }

    // ------------------------------------------------------------------ //
    // Unified state from the single potential.
    // ------------------------------------------------------------------ //
    State evaluate(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double T, double rho) const
    {
        State out;
        out.T = T;
        out.rho = rho;
        if (names.size() != Y.size() || !std::isfinite(T) || !std::isfinite(rho)
            || T <= 0.0 || rho <= 0.0)
            return out;

        const IdealGasState ig = idealGasState(names, Y, T, rho);
        if (!std::isfinite(ig.e) || !std::isfinite(ig.cv)) return out;

        const ResidualFade f = fade();
        const double fd = f.fade(T);
        const double fp = f.derivative(T);
        const double fpp = f.secondDerivative(T);
        out.fade = fd;

        const ResidualState rs = prResidual(names, Y, T, rho);
        if (!std::isfinite(rs.a_res)) return out;
        out.residualMagnitude = std::abs(rs.a_res);

        const double e = ig.e + fd*rs.a_res
            - T*(fp*rs.a_res + fd*rs.da_res_dT);
        const double s = ig.s - fp*rs.a_res - fd*rs.da_res_dT;
        const double Rmix = gasConstant(names, Y);
        const double p = rho*Rmix*T + rho*rho*fd*rs.da_res_drho;
        const double cv = ig.cv
            - T*(fpp*rs.a_res + 2.0*fp*rs.da_res_dT
                 + fd*rs.d2a_res_dT2);
        const double dpdrho = Rmix*T + fd*(2.0*rho*rs.da_res_drho
            + rho*rho*rs.d2a_res_drho2);
        const double dpdT = rho*Rmix + rho*rho*(fp*rs.da_res_drho
            + fd*rs.d2a_res_dTdrho);
        double cp = cv;
        if (dpdrho > 0.0)
            cp = cv + T*dpdT*dpdT/(rho*rho*dpdrho);
        double a2 = 0.0;
        if (cv > 0.0)
            a2 = dpdrho + T/(rho*rho*cv)*dpdT*dpdT;

        out.p = p;
        out.e = e;
        out.h = e + p/rho;
        out.s = s;
        out.cv = cv;
        out.cp = cp;
        out.dpdrho_T = dpdrho;
        out.dpdT_rho = dpdT;
        out.soundSpeed = (a2 > 0.0) ? std::sqrt(a2)
            : std::numeric_limits<double>::quiet_NaN();
        if (std::isfinite(p) && p > 0.0
            && std::isfinite(e) && std::isfinite(s)
            && std::isfinite(cv) && cv > 0.0
            && std::isfinite(cp) && cp > cv
            && dpdrho > 0.0 && a2 > 0.0)
            out.status = Status::SUCCESS;
        else
            out.status = Status::NONPHYSICAL_STATE;
        return out;
    }

    // Recovery cost counters (GERG main path): adaptive bracket-end scans
    // and bisection iterations of recoverRhoEY.
    struct RecoveryCounters
    {
        long long calls = 0;
        long long bisectionIters = 0;
        long long adaptiveScans = 0;
    };
    mutable RecoveryCounters recoveryCounters_;
    const RecoveryCounters& recoveryCounters() const { return recoveryCounters_; }
    void resetRecoveryCounters() const { recoveryCounters_ = {}; }

    // ------------------------------------------------------------------ //
    // rho-e-z recovery: F(T) = e(T,rho,Y) - e_target = 0 (bracketed).
    // ------------------------------------------------------------------ //
    State recoverRhoEY(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double rho, double eTarget,
        double tMin = 60.0, double tMax = 4000.0) const
    {
        State out;
        out.rho = rho;
        ++recoveryCounters_.calls;
        // Auto-find physical bracket ends: in the dense liquid branch the
        // low-temperature side lies inside the two-phase dome (no single-phase
        // state), so tMin/tMax are pushed inward until both ends succeed.
        const int samples = 64;
        double lo = tMin, hi = tMax;
        State sLo = evaluate(names, Y, tMin, rho);
        State sHi = evaluate(names, Y, tMax, rho);
        if (!sLo.success())
        {
            bool found = false;
            for (int i = 1; i <= samples; ++i)
            {
                ++recoveryCounters_.adaptiveScans;
                const double tq = tMin
                    + (tMax - tMin) * static_cast<double>(i) / samples;
                const State sq = evaluate(names, Y, tq, rho);
                if (sq.success()) { lo = tq; sLo = sq; found = true; break; }
            }
            if (!found) return out;
        }
        if (!sHi.success())
        {
            bool found = false;
            for (int i = 1; i <= samples; ++i)
            {
                ++recoveryCounters_.adaptiveScans;
                const double tq = tMax
                    - (tMax - tMin) * static_cast<double>(i) / samples;
                const State sq = evaluate(names, Y, tq, rho);
                if (sq.success()) { hi = tq; sHi = sq; found = true; break; }
            }
            if (!found) return out;
        }
        if (lo >= hi) return out;
        if (eTarget < sLo.e || eTarget > sHi.e)
        {
            out.status = Status::NONPHYSICAL_STATE;
            return out;
        }
        for (int iter = 0; iter < 100; ++iter)
        {
            ++recoveryCounters_.bisectionIters;
            const double Tmid = 0.5*(lo+hi);
            const State s = evaluate(names, Y, Tmid, rho);
            if (!s.success()) { hi = Tmid; continue; }
            const double F = s.e - eTarget;
            if (std::abs(F) < 1.0e-8*std::abs(eTarget) + 1.0e-10)
                return s;
            if (F > 0.0) hi = Tmid; else lo = Tmid;
        }
        out.status = Status::NOT_CONVERGED;
        return out;
    }

    // ------------------------------------------------------------------ //
    // Isobaric roots at (T, p, Y).
    // ------------------------------------------------------------------ //
    IsobaricRoots isobaricRoots(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double T, double p) const
    {
        IsobaricRoots out;
        if (names.size() != Y.size() || T <= 0.0 || p <= 0.0) return out;
        const double a = prAMix(names, Y, T);
        const double b = prBMix(names, Y);
        const double MW = mixtureMw(names, Y);
        const double Ap = a*p/(R_UNIVERSAL*R_UNIVERSAL*T*T);
        const double Bp = b*p/(R_UNIVERSAL*T);
        const double c2 = -(1.0-Bp);
        const double c1v = Ap - 3.0*Bp*Bp - 2.0*Bp;
        const double c0 = -(Ap*Bp - Bp*Bp - Bp*Bp*Bp);
        const double p3 = (3.0*c1v - c2*c2)/3.0;
        const double q3 = (2.0*c2*c2*c2 - 9.0*c2*c1v + 27.0*c0)/27.0;
        const double disc = q3*q3/4.0 + p3*p3*p3/27.0;
        std::vector<double> roots;
        if (disc >= 0.0)
        {
            const double u = std::cbrt(-q3/2.0 + std::sqrt(disc));
            const double v = std::cbrt(-q3/2.0 - std::sqrt(disc));
            roots.push_back(u + v - c2/3.0);
        }
        else
        {
            const double r = std::sqrt(-p3*p3*p3/27.0);
            const double theta = std::acos(-q3/(2.0*r));
            for (int k = 0; k < 3; ++k)
                roots.push_back(2.0*std::sqrt(-p3/3.0)
                    *std::cos(theta/3.0 + 2.0*M_PI*k/3.0) - c2/3.0);
        }
        std::vector<double> admissible;
        for (double r : roots)
            if (r > Bp*1.0001 && std::isfinite(r))
                admissible.push_back(r);
        if (admissible.empty()) return out;
        std::sort(admissible.begin(), admissible.end());
        out.valid = true;
        if (admissible.size() >= 3)
        {
            out.twoPhase = true;
            out.ZLiquid = admissible.front();
            out.ZVapor = admissible.back();
            out.rhoLiquidKgM3 = p*MW/(out.ZLiquid*R_UNIVERSAL*T);
            out.rhoVaporKgM3 = p*MW/(out.ZVapor*R_UNIVERSAL*T);
        }
        else
        {
            out.ZLiquid = out.ZVapor = admissible.front();
            out.rhoLiquidKgM3 = out.rhoVaporKgM3 =
                p*MW/(admissible.front()*R_UNIVERSAL*T);
        }
        return out;
    }

    // Standard PR fugacity-coefficient logarithm at an isobaric root Z.
    double prLnFugacity(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double T, double p, double Z) const
    {
        if (!std::isfinite(Z) || Z <= 0.0)
            return std::numeric_limits<double>::quiet_NaN();
        const double a = prAMix(names, Y, T);
        const double b = prBMix(names, Y);
        const double A = a*p/(R_UNIVERSAL*R_UNIVERSAL*T*T);
        const double B = b*p/(R_UNIVERSAL*T);
        if (Z - B <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        const double sqrt2 = std::sqrt(2.0);
        const double term = (Z + (1.0+sqrt2)*B)/(Z + (1.0-sqrt2)*B);
        if (term <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        return Z - 1.0 - std::log(Z - B)
            - (A/(2.0*sqrt2*B))*std::log(term);
    }

    // Single-component bubble pressure via ln(phi_L) = ln(phi_V).
    double bubblePressure(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double T) const
    {
        double pLo = 1.0e1, pHi = 1.0e8;
        for (int iter = 0; iter < 80; ++iter)
        {
            const double pMid = 0.5*(pLo + pHi);
            const IsobaricRoots roots = isobaricRoots(names, Y, T, pMid);
            if (!roots.valid || !roots.twoPhase)
            {
                if (roots.valid && roots.rhoLiquidKgM3 < 200.0)
                    pLo = pMid;
                else
                    pHi = pMid;
                continue;
            }
            const double lfL = prLnFugacity(names, Y, T, pMid, roots.ZLiquid);
            const double lfV = prLnFugacity(names, Y, T, pMid, roots.ZVapor);
            if (!std::isfinite(lfL) || !std::isfinite(lfV)) { pHi = pMid; continue; }
            if (lfL < lfV) pHi = pMid; else pLo = pMid;
        }
        return 0.5*(pLo + pHi);
    }

private:
    static constexpr double R_UNIVERSAL = 8.314462618;

    // Single active species above a mass-fraction threshold.
    static std::string activePureSpecies(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double threshold = 0.999)
    {
        if (names.size() != Y.size()) return std::string();
        std::string pure;
        for (std::size_t i = 0; i < names.size(); ++i)
        {
            if (!std::isfinite(Y[i]) || Y[i] < 0.0) return std::string();
            if (Y[i] > threshold)
            {
                if (!pure.empty()) return std::string();
                pure = names[i];
            }
        }
        return pure;
    }

    // O2/C3H8 binary mixture key at the nearest fixed-composition table
    // level.  Returns empty unless O2+C3H8 dominate the stable set (>99.9%
    // of mass) with both present; radicals never enter this route.
    static std::string binaryMixKey(
        const std::vector<std::string>& names,
        const std::vector<double>& Y)
    {
        if (names.size() != Y.size()) return std::string();
        double nO2 = 0.0, nC3 = 0.0, otherMass = 0.0;
        for (std::size_t i = 0; i < names.size(); ++i)
        {
            const std::string& n = names[i];
            const double y = Y[i];
            if (!std::isfinite(y) || y < 0.0) return std::string();
            if (n == "O2") nO2 += y / 0.0319988;
            else if (n == "C3H8") nC3 += y / 0.04409562;
            else otherMass += y;
        }
        if (nO2 <= 0.0 || nC3 <= 0.0) return std::string();
        if (otherMass > 1.0e-3) return std::string();
        const double xO2 = nO2 / (nO2 + nC3);
        // Table composition levels that survive the low-temperature GERG
        // health audit (2026-08-03): x0.05/x0.25 fail in the deep-cold
        // propane-rich region (both torch-flash and CoolProp root solvers
        // diverge there; the 20 propane-rich NIST points stay uncovered
        // until teqp C++ provides a healthier GERG root solver).
        static const double levels[] = {0.50, 0.75, 0.95};
        static const char* levelKeys[] =
            {"O2C3H8_x0.5", "O2C3H8_x0.75", "O2C3H8_x0.95"};
        int best = 0;
        double bestDist = std::abs(xO2 - levels[0]);
        for (int i = 1; i < 3; ++i)
        {
            const double dist = std::abs(xO2 - levels[i]);
            if (dist < bestDist) { bestDist = dist; best = i; }
        }
        return levelKeys[best];
    }

    double gasConstant(
        const std::vector<std::string>& names,
        const std::vector<double>& Y) const
    {
        double Rm = 0.0;
        for (std::size_t i = 0; i < names.size(); ++i)
            if (Y[i] > 0.0)
                Rm += Y[i] * R_UNIVERSAL
                    / registry_.molecularWeight(registry_.index(names[i]));
        return Rm;
    }

    double mixtureMw(
        const std::vector<std::string>& names,
        const std::vector<double>& Y) const
    {
        double inv = 0.0;
        for (std::size_t i = 0; i < names.size(); ++i)
            if (Y[i] > 0.0)
                inv += Y[i] / registry_.molecularWeight(registry_.index(names[i]));
        return 1.0 / inv;
    }

    double pIdeal(double T, double rhoMolar) const
    { return rhoMolar * R_UNIVERSAL * T; }

    double prAMix(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double T) const
    {
        double am = 0.0;
        const int n = static_cast<int>(names.size());
        std::vector<double> x(n), ai(n);
        for (int i = 0; i < n; ++i)
        {
            if (Y[i] <= 0.0) { x[i] = 0.0; ai[i] = 0.0; continue; }
            x[i] = (Y[i]/registry_.molecularWeight(registry_.index(names[i])))
                * mixtureMw(names, Y);
            const auto c = criticalDB_.get(names[i]);
            const double kappa = 0.37464 + 1.54226*c.omega - 0.26992*c.omega*c.omega;
            const double alpha = std::pow(
                1.0 + kappa*(1.0 - std::sqrt(T/c.Tc)), 2.0);
            ai[i] = 0.45724 * R_UNIVERSAL*R_UNIVERSAL*c.Tc*c.Tc/c.Pc * alpha;
        }
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (x[i] > 0.0 && x[j] > 0.0)
                    am += x[i]*x[j]*std::sqrt(ai[i]*ai[j])
                        * (1.0 - prKij(names[i], names[j]));
        return am;
    }

    double prDAdT(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double T) const
    {
        const int n = static_cast<int>(names.size());
        double am = 0.0;
        std::vector<double> x(n), ai(n), dai(n);
        for (int i = 0; i < n; ++i)
        {
            if (Y[i] <= 0.0) { x[i]=0; ai[i]=0; dai[i]=0; continue; }
            x[i] = (Y[i]/registry_.molecularWeight(registry_.index(names[i])))
                * mixtureMw(names, Y);
            const auto c = criticalDB_.get(names[i]);
            const double kappa = 0.37464 + 1.54226*c.omega - 0.26992*c.omega*c.omega;
            const double a0 = 0.45724*R_UNIVERSAL*R_UNIVERSAL*c.Tc*c.Tc/c.Pc;
            const double f = 1.0 + kappa*(1.0 - std::sqrt(T/c.Tc));
            ai[i] = a0*f*f;
            dai[i] = -a0*kappa*f/std::sqrt(T*c.Tc);
        }
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
            {
                if (x[i] <= 0.0 || x[j] <= 0.0) continue;
                const double sq = std::sqrt(ai[i]*ai[j]);
                if (sq <= 0.0) continue;
                const double dsq = 0.5*(dai[i]*ai[j] + ai[i]*dai[j])/sq;
                am += x[i]*x[j]*(1.0 - prKij(names[i],names[j]))*dsq;
            }
        return am;
    }

    double prD2AdT2(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double T) const
    {
        const int n = static_cast<int>(names.size());
        double am = 0.0;
        std::vector<double> x(n), ai(n), dai(n), d2ai(n);
        for (int i = 0; i < n; ++i)
        {
            if (Y[i] <= 0.0) { x[i]=0; ai[i]=0; dai[i]=0; d2ai[i]=0; continue; }
            x[i] = (Y[i]/registry_.molecularWeight(registry_.index(names[i])))
                * mixtureMw(names, Y);
            const auto c = criticalDB_.get(names[i]);
            const double kappa = 0.37464 + 1.54226*c.omega - 0.26992*c.omega*c.omega;
            const double a0 = 0.45724*R_UNIVERSAL*R_UNIVERSAL*c.Tc*c.Tc/c.Pc;
            const double f = 1.0 + kappa*(1.0 - std::sqrt(T/c.Tc));
            ai[i] = a0*f*f;
            dai[i] = -a0*kappa*f/std::sqrt(T*c.Tc);
            d2ai[i] = a0*kappa*(1.0+kappa)/(2.0*T*std::sqrt(T*c.Tc));
        }
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
            {
                if (x[i] <= 0.0 || x[j] <= 0.0) continue;
                const double sq = std::sqrt(ai[i]*ai[j]);
                if (sq <= 0.0) continue;
                const double dsq = 0.5*(dai[i]*ai[j]+ai[i]*dai[j])/sq;
                const double d2sq = 0.5*(d2ai[i]*ai[j] + 2.0*dai[i]*dai[j]
                    + ai[i]*d2ai[j])/sq - dsq*dsq/sq;
                am += x[i]*x[j]*(1.0 - prKij(names[i],names[j]))*d2sq;
            }
        return am;
    }

    double prBMix(
        const std::vector<std::string>& names,
        const std::vector<double>& Y) const
    {
        double bm = 0.0;
        const int n = static_cast<int>(names.size());
        for (int i = 0; i < n; ++i)
        {
            if (Y[i] <= 0.0) continue;
            const double x = (Y[i]/registry_.molecularWeight(registry_.index(names[i])))
                * mixtureMw(names, Y);
            const auto c = criticalDB_.get(names[i]);
            bm += x * 0.07780 * R_UNIVERSAL * c.Tc / c.Pc;
        }
        return bm;
    }

    double prKij(const std::string& a, const std::string& b) const
    {
        if ((a=="O2" && b=="C3H8") || (a=="C3H8" && b=="O2")) return 0.055;
        return 0.0;
    }

    // Largest admissible isochoric PR root.
    double prRoot(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double T, double rho) const
    {
        const double MW = mixtureMw(names, Y);
        const double rhoMolar = rho / MW;
        const double a = prAMix(names, Y, T);
        const double b = prBMix(names, Y);
        const double A = a*pIdeal(T,rhoMolar)/(R_UNIVERSAL*R_UNIVERSAL*T*T);
        const double B = b*pIdeal(T,rhoMolar)/(R_UNIVERSAL*T);
        const double c2 = -(1.0-B);
        const double c1 = A - 3.0*B*B - 2.0*B;
        const double c0 = -(A*B - B*B - B*B*B);
        const double p3 = (3.0*c1 - c2*c2)/3.0;
        const double q3 = (2.0*c2*c2*c2 - 9.0*c2*c1 + 27.0*c0)/27.0;
        const double disc = q3*q3/4.0 + p3*p3*p3/27.0;
        std::vector<double> roots;
        if (disc >= 0.0)
        {
            const double u = std::cbrt(-q3/2.0 + std::sqrt(disc));
            const double v = std::cbrt(-q3/2.0 - std::sqrt(disc));
            roots.push_back(u + v - c2/3.0);
        }
        else
        {
            const double r = std::sqrt(-p3*p3*p3/27.0);
            const double theta = std::acos(-q3/(2.0*r));
            for (int k = 0; k < 3; ++k)
                roots.push_back(2.0*std::sqrt(-p3/3.0)
                    *std::cos(theta/3.0 + 2.0*M_PI*k/3.0) - c2/3.0);
        }
        double Z = -1.0;
        for (double r : roots)
            if (r > B*1.01 && r > Z) Z = r;
        return Z;
    }

    // Residual Helmholtz value (molar) at (T,rho) for numeric derivatives.
    double prResidualValue(
        const std::vector<std::string>& names,
        const std::vector<double>& Y,
        double T, double rho, double MW) const
    {
        const double rhoMolar = rho / MW;
        const double a = prAMix(names, Y, T);
        const double b = prBMix(names, Y);
        const double B = b*pIdeal(T,rhoMolar)/(R_UNIVERSAL*T);
        const double Z = prRoot(names, Y, T, rho);
        if (!std::isfinite(Z) || Z <= 0.0 || Z-B <= 0.0)
            return std::numeric_limits<double>::quiet_NaN();
        const double sqrt2 = std::sqrt(2.0);
        const double term = (Z+(1.0+sqrt2)*B)/(Z+(1.0-sqrt2)*B);
        if (term <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        return R_UNIVERSAL*T*(-std::log(Z-B))
            + a/(2.0*sqrt2*b)*std::log(term);
    }

    const SpeciesRegistry& registry_;
    SpeciesCriticalDB criticalDB_;
    const GERGResidualTable* residualTable_ = nullptr;
    double fadeT0_;
    double fadeT1_;
};

} // End namespace Foam

#endif
