/*---------------------------------------------------------------------------*\
  realFluidAllMach - AUSM+-up Flux Implementation
\*---------------------------------------------------------------------------*/

#include "AUSMPlusUpFlux.H"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Foam
{

AUSMPlusUpFlux::AUSMPlusUpFlux()
:
    Ku_(0.75),
    Kp_(0.25),
    sigma_(1.0)
{}

AUSMPlusUpFlux::AUSMPlusUpFlux(const std::vector<std::string>& speciesNames)
:
    Ku_(0.75),Kp_(0.25),sigma_(1.0),speciesNames_(speciesNames)
{}

void AUSMPlusUpFlux::setLowMachParams(double Ku, double Kp, double sigma)
{
    Ku_ = Ku;
    Kp_ = Kp;
    sigma_ = sigma;
}

// ---------------------------------------------------------------------------
// Mach splitting, Liou (2006) eqs. (24)-(28).
//
//   M+/-(M)   = +/- 1/4 (M +/- 1)^2 + beta (M^2 - 1)^2     for |M| < 1
//             = 1/2 (M +/- |M|)                            for |M| >= 1
//
// beta is a free parameter. AUSM+-up uses beta = 1/8.
//
// The M+ / M- pair ALONE is not a complete low-Mach splitting: for Ma << 1 it
// is essentially centred and contributes no first-order dissipation between
// the two states. The M^4 form is what restores it (see M4_plus below).
// ---------------------------------------------------------------------------
double AUSMPlusUpFlux::M_plus(double M, double beta) const
{
    if (std::abs(M) >= 1.0)
        return 0.5 * (M + std::abs(M));
    else
        return 0.25 * std::pow(M + 1.0, 2) + beta * std::pow(M*M - 1.0, 2);
}

double AUSMPlusUpFlux::M_minus(double M, double beta) const
{
    if (std::abs(M) >= 1.0)
        return 0.5 * (M - std::abs(M));
    else
        return -0.25 * std::pow(M - 1.0, 2) - beta * std::pow(M*M - 1.0, 2);
}

// M^4 splitting, Liou (2006) eqs. (29)-(30):
//
//   M4+/-(M) = M+/-(M) * (1 -/+ 16 M+ M-)
//
// At Ma << 1 the pair sums to
//
//   M4+ + M4- ~= Mbar + 1/4 (M_R - M_L) + O(Ma^2)
//
// The (1/4)(M_R - M_L) term is the low-Mach mass-flux dissipation. Without it
// the interface mass flux reduces to rho * M_p, i.e. it is driven only by the
// pressure difference: no acoustic damping, no upwinding on a uniform field,
// and an implicit diagonal that no longer matches the real flux Jacobian.
// See the use site in computeFluxCached().
double AUSMPlusUpFlux::M4_plus(double M, double beta) const
{
    if (std::abs(M) >= 1.0)
        return 0.5 * (M + std::abs(M));
    else
    {
        const double Mp = M_plus(M, beta);
        const double Mm = M_minus(M, beta);
        return Mp * (1.0 - 16.0 * Mp * Mm);
    }
}

double AUSMPlusUpFlux::M4_minus(double M, double beta) const
{
    if (std::abs(M) >= 1.0)
        return 0.5 * (M - std::abs(M));
    else
    {
        const double Mp = M_plus(M, beta);
        const double Mm = M_minus(M, beta);
        return Mm * (1.0 - 16.0 * Mp * Mm);
    }
}

// Pressure splitting, Liou (2006) eqs. (31)-(32), with alpha = 3/16 for
// AUSM+-up. P+ + P- = 1 identically, so a uniform-pressure field produces no
// spurious momentum flux; that identity is what the flux regression relies on.
double AUSMPlusUpFlux::P_plus(double M, double alpha) const
{
    if (std::abs(M) >= 1.0)
        return 0.5 * (1.0 + std::copysign(1.0, M));
    else
        return 0.25 * std::pow(M + 1.0, 2) * (2.0 - M) + alpha * M * std::pow(M*M - 1.0, 2);
}

double AUSMPlusUpFlux::P_minus(double M, double alpha) const
{
    if (std::abs(M) >= 1.0)
        return 0.5 * (1.0 - std::copysign(1.0, M));
    else
        return 0.25 * std::pow(M - 1.0, 2) * (2.0 + M) - alpha * M * std::pow(M*M - 1.0, 2);
}

// fa = M0 (2 - M0), the low-Mach scaling of the pressure-diffusion term.
// M0 is NOT the local Mach number: it is a reference Mach floored away from
// zero, sqrt(min(1, max(Mbar^2, 0.09))) with the 0.3 floor. It carries the
// length-scale dependence of the limiter (Liou & Edwards 2005); the global
// 0.3 floor used here is a simplification -- see the note in the header
// helper ausmPressureDiffusionPressureCoefficient.
double AUSMPlusUpFlux::f_a(double M0) const
{
    return M0 * (2.0 - M0);
}

void AUSMPlusUpFlux::computeFlux(
    const PrimitiveVars& left,
    const PrimitiveVars& right,
    const std::vector<double>& normal,
    ConservedVars& flux
) const
{
    // Un-cached entry point: builds a fresh RealFluidState on both sides for
    // every call. That is two EOS updates per face per residual evaluation, so
    // it is only appropriate for 1D/unit tests and standalone tooling. The
    // production tetra kernel goes through computeFluxCached(), which receives
    // rho/a/h from the committed cell Q payload and never re-flashes.
    // Get thermodynamic states
    RealFluidState stateL, stateR;
    if(!speciesNames_.empty())
    {
        stateL.setSpecies(speciesNames_);
        stateR.setSpecies(speciesNames_);
    }
    const auto statusL = stateL.updatePTY(left.p, left.T, left.Y);
    const auto statusR = stateR.updatePTY(right.p, right.T, right.Y);
    if (statusL != RealFluidState::Status::SUCCESS ||
        statusR != RealFluidState::Status::SUCCESS)
        throw std::runtime_error(
            "AUSMPlusUpFlux: invalid thermodynamic face state");
    computeFluxCached(left,right,normal,flux,stateL.rho(),stateL.a(),stateL.h(),
                      stateR.rho(),stateR.a(),stateR.h());
}

void AUSMPlusUpFlux::computeFluxCached(
    const PrimitiveVars& left,const PrimitiveVars& right,
    const std::vector<double>& normal,ConservedVars& flux,
    double rhoLeft,double soundLeft,double enthalpyLeft,
    double rhoRight,double soundRight,double enthalpyRight) const
{
    // Cached entry point. The rho/a/h arguments are the committed equilibrium
    // values for left/right; they are validated for finiteness and positivity
    // below and then used as-is. This is deliberate: it guarantees that the
    // flux consumes exactly the state the solver published (one owner for
    // U/W/Q), so a face cannot silently use a different density than the one
    // the cell recovery committed. The cost is that a stale or inconsistent Q
    // payload surfaces as a thrown error here instead of a quiet bias.
    if (normal.size() != 3 || left.U.size() != 3 || right.U.size() != 3
        || left.Y.size() != right.Y.size()
        || flux.rhoY.size() != left.Y.size()
        || !std::isfinite(left.p) || left.p <= 0.0
        || !std::isfinite(right.p) || right.p <= 0.0
        || !std::isfinite(rhoLeft) || rhoLeft <= 0.0
        || !std::isfinite(rhoRight) || rhoRight <= 0.0
        || !std::isfinite(soundLeft) || soundLeft <= 0.0
        || !std::isfinite(soundRight) || soundRight <= 0.0
        || !std::isfinite(enthalpyLeft)
        || !std::isfinite(enthalpyRight))
        throw std::runtime_error(
            "AUSMPlusUpFlux: invalid cached thermodynamic payload");
    for (std::size_t direction = 0; direction < 3; ++direction)
        if (!std::isfinite(normal[direction])
            || !std::isfinite(left.U[direction])
            || !std::isfinite(right.U[direction]))
            throw std::runtime_error(
                "AUSMPlusUpFlux: invalid cached face geometry or velocity");

    // Normal velocity
    double Un_L = 0.0, Un_R = 0.0;
    for (size_t i = 0; i < 3; i++)
    {
        Un_L += left.U[i] * normal[i];
        Un_R += right.U[i] * normal[i];
    }

    // Sound speeds
    double a_L = soundLeft;
    double a_R = soundRight;
    double a_half = 0.5 * (a_L + a_R);

    // Mach numbers
    double M_L = Un_L / a_half;
    double M_R = Un_R / a_half;

    // Pressure diffusion. The shared coefficient owns the all-speed M0/fa
    // floor, so production flux and structured Schur have one formula.
    // The coefficient is shared with the structured
    // pressure-Schur linearisation so production flux and preconditioner
    // cannot silently drift to different AUSM+-up formulas.
    const double pressureCoefficient=
        ausmPressureDiffusionPressureCoefficient(
            left,right,normal,rhoLeft,a_L,rhoRight,a_R,Kp_,sigma_);
    if(!std::isfinite(pressureCoefficient))
        throw std::runtime_error(
            "AUSMPlusUpFlux: invalid pressure-diffusion coefficient");
    // Two deviations from Liou (2006) are folded inside the helper, both
    // documented at their definition in the header:
    //   (1) M0 uses a global 0.3 floor instead of the wall-distance limiter;
    //   (2) the reference density is the arithmetic mean, not density-upwind.
    // Neither is a bug on its own, but together they make the pressure
    // diffusion uniform where the injector needs it differentiated, so the
    // near-field residual hotspot is partly this, not the mesh.
    double Mp=-pressureCoefficient*(right.p-left.p);

    // Interface Mach number.
    //
    // The branch is selectable so the acoustic-decay measurement can compare
    // both splittings through this exact code path. Default is
    // VELOCITY_SPLIT_ONLY (see AUSMMachSplitting.H).
    //
    // WHY THE M^4 BRANCH IS OFF BY DEFAULT (updated 2026-10-09):
    // The recorded reason (2026-08-04) was a +65 % shift of the measured
    // acoustic wave speed. The cause is now identified, and it is not a
    // coefficient detail -- the implemented M4 pair violates the uniform-flow
    // identity that any Mach splitting must satisfy:
    //
    //     required:  M4+(M) + M4-(M) == M       for all |M| < 1
    //     observed:  M4+(M) + M4-(M) == 3.25 M  as M -> 0
    //
    // i.e. +225 % on the convective part at low Mach, which is the same defect
    // the +65 % wave-speed shift was a symptom of. The two M4_plus/M4_minus
    // functions multiply by (1 - 16 M+ M-) with the SAME sign on both halves,
    // so the product term does not cancel between them. Liou (2006) eq. (29)
    // uses opposite signs, and even that form does not restore the identity
    // (verified numerically: it gives 1.69 M at M -> 0), so the correct
    // expression is still to be established from the reference before this
    // branch is enabled.
    //
    // Consequence for a reader: with the default branch there is NO low-Mach
    // mass-flux dissipation on the acoustic mode -- M+(+e) + M-(-e) is exactly
    // 0 for an antisymmetric acoustic pair. That is a real gap, but filling it
    // with a splitting that breaks uniform flow is worse than leaving it open.
    // Use tests/test_acoustic_decay_1d.cpp to measure both branches.
    double M_half;
    if (splitting_ == AUSMMachSplitting::FOURTH_ORDER_LOW_MACH)
        M_half = M4_plus(M_L, 0.125) + M4_minus(M_R, 0.125) + Mp;
    else
        M_half = M_plus(M_L, 0.125) + M_minus(M_R, 0.125) + Mp;

    // Interface pressure
    double p_half = P_plus(M_L, 0.1875) * left.p + P_minus(M_R, 0.1875) * right.p;

    // Mass flux. The density is the upwind cell's, replaced by the equilibrium
    // density of the interface state when the consistency hook is installed.
    double rhoUp = M_half > 0.0 ? rhoLeft : rhoRight;
    if (consistentDensity_)
    {
        const bool upLeft = M_half > 0.0;
        const double hUp = upLeft ? enthalpyLeft : enthalpyRight;
        const double yO2 = upLeft ? left.Y[0] : right.Y[0];
        double rhoConsistent = 0.0;
        if (consistentDensity_(yO2, p_half, hUp, rhoConsistent))
            rhoUp = rhoConsistent;
        else
            consistentFallbacks_.fetch_add(1, std::memory_order_relaxed);
    }
    double mdot = a_half * M_half * rhoUp;

    // Convective fluxes
    flux.rho = mdot;

    for (size_t i = 0; i < 3; i++)
    {
        if (M_half > 0.0)
            flux.rhoU[i] = mdot * left.U[i] + p_half * normal[i];
        else
            flux.rhoU[i] = mdot * right.U[i] + p_half * normal[i];
    }

    if (M_half > 0.0)
    {
        double H_L = enthalpyLeft + 0.5 * (left.U[0]*left.U[0] + left.U[1]*left.U[1] + left.U[2]*left.U[2]);
        flux.rhoE = mdot * H_L;
    }
    else
    {
        double H_R = enthalpyRight + 0.5 * (right.U[0]*right.U[0] + right.U[1]*right.U[1] + right.U[2]*right.U[2]);
        flux.rhoE = mdot * H_R;
    }

    // Species fluxes. Note this is a plain upwind of Y with the same mdot as
    // the mass/energy flux, so species and mass stay exactly consistent under
    // the same interface Mach number. Any future switch to a species-
    // diffusion (Fick) term must keep that consistency or the recovered
    // composition will drift from the transported one.
    for (size_t i = 0; i < left.Y.size(); i++)
    {
        if (M_half > 0.0)
            flux.rhoY[i] = mdot * left.Y[i];
        else
            flux.rhoY[i] = mdot * right.Y[i];
    }
}

} // End namespace Foam
