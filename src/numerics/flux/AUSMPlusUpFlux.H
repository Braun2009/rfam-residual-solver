/*---------------------------------------------------------------------------*\
  realFluidAllMach - AUSM+-up Flux Scheme

  All-speed flux scheme with low-Mach correction.

  The AUSM+-up low-Mach terms are part of the numerical flux. They are not a
  Turkel/Weiss transformation of the residual, time derivative, or Jacobian.
  Do not report this flux, by itself, as full all-Mach preconditioning.
  Reference: Liou (2006)
\*---------------------------------------------------------------------------*/

#ifndef AUSMPlusUpFlux_H
#define AUSMPlusUpFlux_H

#include "ConservedVars.H"
#include "RealFluidState.H"
#include "AUSMMachSplitting.H"
#include <atomic>
#include <functional>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <limits>

namespace Foam
{

// Shared linear pressure-diffusion coefficient used by AUSM+-up and the
// opt-in pressure-Schur operator. The first helper returns the positive
// pressure-Mach coefficient C_p in M_p=-C_p(p_R-p_L), with units 1/Pa;
// the mass-mobility helper includes a_half and the frozen AUSM upwind density.
// Keeping this expression here prevents the Schur operator from silently
// drifting from the production AUSM flux formula.
//
// Deliberate deviations from Liou (2006) eq. (33), which a reader comparing
// against the paper needs to know about:
//
//   Liou:   C_p = (Kp/fa) * max(1 - sigma*Mbar^2, 0) / (rho_ref * a^2)
//           with M0^2 = min(1, max(Mbar^2, Kp*M_min_length^2)), M_min_length
//           a length-scaled reference Mach, and rho_ref a density-UPWIND value.
//
//   1. M0 uses a GLOBAL floor of 0.3 (Kp * 0.09 with Kp folded in) instead of
//      the wall-distance limiter. Effect: the pressure diffusion is the same
//      magnitude in a 0.22 um injector cell and in a 5 mm free-shear cell,
//      where Liou's limiter would differentiate them. This is why the near-
//      field residual hotspot is where it is, and it is the first thing to
//      revisit when the injector run is unblocked.
//   2. rho_ref is the arithmetic mean (rho_L+rho_R)/2. At a cryogenic
//      interface, where the density can jump by ~3x, this UNDER-estimates the
//      coefficient -- i.e. the pressure coupling is weakest exactly at the
//      propellant interface where the flux is hardest.
//
// Neither is a correctness bug in isolation; both are accuracy items with a
// measurable footprint. Changing either alters the flux matrix, so do it
// behind the acoustic-decay and low-Mach-order tests rather than by eye.
//
// Units: 1/Pa.
inline double ausmPressureDiffusionPressureCoefficient(
    const PrimitiveVars& left,const PrimitiveVars& right,
    const std::vector<double>& normal,
    double rhoLeft,double soundLeft,double rhoRight,double soundRight,
    double Kp=0.25,double sigma=1.0)
{
    if(normal.size()!=3||left.U.size()!=3||right.U.size()!=3
        ||!std::isfinite(rhoLeft)||rhoLeft<=0.0
        ||!std::isfinite(rhoRight)||rhoRight<=0.0
        ||!std::isfinite(soundLeft)||soundLeft<=0.0
        ||!std::isfinite(soundRight)||soundRight<=0.0
        ||!std::isfinite(Kp)||Kp<0.0||!std::isfinite(sigma))
        return std::numeric_limits<double>::quiet_NaN();
    double unLeft=0.0,unRight=0.0;
    for(std::size_t direction=0;direction<3;++direction)
    {
        if(!std::isfinite(normal[direction])
            ||!std::isfinite(left.U[direction])
            ||!std::isfinite(right.U[direction]))
            return std::numeric_limits<double>::quiet_NaN();
        unLeft+=left.U[direction]*normal[direction];
        unRight+=right.U[direction]*normal[direction];
    }
    const double aHalf=0.5*(soundLeft+soundRight);
    const double rhoHalf=0.5*(rhoLeft+rhoRight);
    if(!std::isfinite(aHalf)||aHalf<=0.0
        ||!std::isfinite(rhoHalf)||rhoHalf<=0.0)
        return std::numeric_limits<double>::quiet_NaN();
    const double mLeft=unLeft/aHalf;
    const double mRight=unRight/aHalf;
    // mbar^2 is the face-referenced Mach scale, then floored at 0.3 (see the
    // deviation note above). Both sides are weighted equally by construction.
    const double mbarSquared=0.5*(mLeft*mLeft+mRight*mRight);
    const double m0Squared=std::min(1.0,std::max(mbarSquared,0.3*0.3));
    const double m0=std::sqrt(m0Squared);
    // fa = M0(2-M0) -> 1 as M0 -> 1, 0.51 at the 0.3 floor. It only ever
    // weakens the diffusion, never flips its sign, so the term stays
    // dissipative for any admissible (Kp, sigma).
    const double fa=m0*(2.0-m0);
    if(!std::isfinite(fa)||fa<=0.0) return std::numeric_limits<double>::quiet_NaN();
    return Kp/fa*std::max(1.0-sigma*m0Squared,0.0)
        /(rhoHalf*aHalf*aHalf);
}

inline double ausmPressureDiffusionMassMobility(
    const PrimitiveVars& left,const PrimitiveVars& right,
    const std::vector<double>& normal,
    double rhoLeft,double soundLeft,double rhoRight,double soundRight,
    double Kp=0.25,double sigma=1.0)
{
    const double coefficient=ausmPressureDiffusionPressureCoefficient(
        left,right,normal,rhoLeft,soundLeft,rhoRight,soundRight,Kp,sigma);
    if(!std::isfinite(coefficient))
        return std::numeric_limits<double>::quiet_NaN();
    const double aHalf=0.5*(soundLeft+soundRight);
    double unLeft=0.0,unRight=0.0;
    for(std::size_t direction=0;direction<3;++direction)
    {
        unLeft+=left.U[direction]*normal[direction];
        unRight+=right.U[direction]*normal[direction];
    }
    const double mLeft=unLeft/aHalf;
    const double mRight=unRight/aHalf;
    // Face mobility d(mdot)/d(p_R - p_L), used by the opt-in pressure-Schur
    // low-frequency operator. It re-uses the SAME coefficient and the SAME
    // interface-Mach branch as computeFluxCached, so the linear operator and
    // the nonlinear flux cannot drift apart. The splitPlus/splitMinus lambdas
    // here duplicate AUSMPlusUpFlux::M_plus/M_minus with beta = 1/8 -- if the
    // M^4 form is restored in the flux (see the disabled-term note there),
    // this helper must be switched in the same commit or the Schur operator
    // will linearise a different flux than the one that runs.
    auto splitPlus=[](double value)
    {
        if(std::abs(value)>=1.0) return 0.5*(value+std::abs(value));
        return 0.25*(value+1.0)*(value+1.0)
            +0.125*(value*value-1.0)*(value*value-1.0);
    };
    auto splitMinus=[](double value)
    {
        if(std::abs(value)>=1.0) return 0.5*(value-std::abs(value));
        return -0.25*(value-1.0)*(value-1.0)
            -0.125*(value*value-1.0)*(value*value-1.0);
    };
    // Freeze the committed AUSM upwind branch. Reversing face orientation
    // swaps left/right and changes M_half's sign, so the selected physical
    // upwind density is unchanged and the face mobility remains reciprocal.
    const double pressureTerm=-coefficient*(right.p-left.p);
    const double interfaceMach=splitPlus(mLeft)+splitMinus(mRight)
        +pressureTerm;
    const double upwindDensity=interfaceMach>0.0?rhoLeft:
        (interfaceMach<0.0?rhoRight:0.5*(rhoLeft+rhoRight));
    return aHalf*upwindDensity*coefficient;
}

class AUSMPlusUpFlux
{
private:
    // Interface-consistency hook (see setConsistentDensityHook below).
    using ConsistentDensityHook =
        std::function<bool(double yO2, double p, double h, double& rho)>;
    ConsistentDensityHook consistentDensity_;
    mutable std::atomic<std::size_t> consistentFallbacks_{0};
    // Low-Mach correction parameters
    double Ku_;  // Velocity diffusion parameter
    double Kp_;  // Pressure diffusion parameter
    double sigma_; // Scaling parameter

    // NOTE (2026-10-09): Ku_ is ASSIGNED BUT NEVER READ. There is no velocity-
    // diffusion term anywhere in computeFluxCached, so this flux implements
    // AUSM+-up with the pressure-diffusion term (Kp, sigma) and WITHOUT the
    // velocity-diffusion term that Liou (2006) eq. (33) specifies alongside it.
    //
    // Measured consequence (tests/test_acoustic_decay_1d): the flux is still
    // dissipative -- sigma*/(k*dx) -> 0.5577, entirely from the pressure term --
    // so this is not an instability. It is a deviation from the published scheme
    // that was never recorded, and it means setLowMachParams() silently has no
    // effect on Ku. Anyone tuning Ku and observing no change now knows why.
    //
    // Do not "fix" this by wiring Ku_ into the flux without re-running the
    // acoustic measurement AND the 1D validation: adding dissipation changes the
    // O(dx) coefficient, and the +65 % wave-speed history shows what happens when
    // a splitting term is changed without a measurement to decide it.
    std::vector<std::string> speciesNames_;
    // Which Mach splitting computeFluxCached() runs. Default is the shipped
    // branch; see AUSMMachSplitting.H for why both exist and what is wrong
    // with the M^4 one as currently implemented.
    AUSMMachSplitting splitting_ = DEFAULT_AUSM_MACH_SPLITTING;

    // Mach splitting functions
    double M_plus(double M, double beta) const;
    double M_minus(double M, double beta) const;
    double P_plus(double M, double alpha) const;
    double P_minus(double M, double alpha) const;
    // M^4 splitting (Liou 2006 eq. 26-28): the fourth-order Mach
    // polynomials provide the low-Mach upwind dissipation that the bare
    // M+ + M- splitting lacks.  They are the acoustic-stability backbone
    // of AUSM+-up at Ma << 1 (the 2026-08-04 acoustic-instability gap).
    double M4_plus(double M, double beta) const;
    double M4_minus(double M, double beta) const;

    // Low-Mach correction
    double f_a(double M0) const;

public:
    AUSMPlusUpFlux();
    explicit AUSMPlusUpFlux(const std::vector<std::string>& speciesNames);
    void setSpeciesNames(const std::vector<std::string>& speciesNames)
    { speciesNames_ = speciesNames; }

    // Compute flux from left and right states
    void computeFlux(
        const PrimitiveVars& left,
        const PrimitiveVars& right,
        const std::vector<double>& normal,
        ConservedVars& flux
    ) const;
    void computeFluxCached(
        const PrimitiveVars& left,const PrimitiveVars& right,
        const std::vector<double>& normal,ConservedVars& flux,
        double rhoLeft,double soundLeft,double enthalpyLeft,
        double rhoRight,double soundRight,double enthalpyRight) const;

    // Set low-Mach parameters
    //
    // Ku is currently INERT: no velocity-diffusion term consumes it. See the
    // note on the Ku_ member. Kept in the signature so existing callers and the
    // regression matrix do not have to change, but a caller that sets Ku and
    // expects a different flux is expecting something this implementation does
    // not do; the acoustic measurement in tests/test_acoustic_decay_1d will
    // confirm the invariance.
    void setLowMachParams(double Ku, double Kp, double sigma);

    // Select the interface Mach splitting. Exposed so the acoustic-decay
    // measurement can run both branches through the SAME production flux code
    // path, rather than through a test-only reimplementation (which is how the
    // two would silently disagree).
    void setMachSplitting(AUSMMachSplitting s) { splitting_ = s; }
    AUSMMachSplitting machSplitting() const { return splitting_; }

    // Interface-consistent density hook. When set, the mass flux uses the
    // equilibrium density of the interface state (p_half, h_upwind, Y_upwind)
    // instead of the upwind cell density; the AUSM pressure/velocity splitting
    // is unchanged, so accuracy/stability are preserved while the mass, energy
    // and species fluxes now all correspond to one equilibrium state (the
    // interface projection error drops accordingly). Returning false falls back
    // to the upwind density and is counted.
    void setConsistentDensityHook(ConsistentDensityHook hook)
    { consistentDensity_ = std::move(hook); }
    std::size_t consistentDensityFallbacks() const
    { return consistentFallbacks_.load(std::memory_order_relaxed); }
};

} // End namespace Foam

#endif
