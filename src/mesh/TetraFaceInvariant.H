#ifndef TetraFaceInvariant_H
#define TetraFaceInvariant_H

/*---------------------------------------------------------------------------*\
    TetraFaceInvariant - the single owner of "how fast information leaves a face".

    WHY THIS FILE EXISTS
    --------------------
    The tetrahedral solver has to answer one question in five different places:

        "What is the fastest wave leaving this face?"

    The answer is the face spectral radius

        Lambda = |u_n| + a,        u_n = U . n_hat,   a = equilibrium sound speed

    and it is consumed by
      1. the explicit / local time-step estimate   (TetraFiniteVolume)
      2. the boundary-aware stability step         (TetraFiniteVolume)
      3. the LU-SGS off-diagonal magnitude         (TetraLUSGS)
      4. the physical-step residual normalisation  (TetraLUSGS)
      5. the worst-cell inlet spectral diagnostic  (TetraLUSGS)
      6. the future RANS eddy-viscosity diagonal   (not written yet)

    Before this header, each of those six sites carried its own copy of the
    formula. That is not a style problem: it is the project's central numerical
    defect, because the LU-SGS diagonal is only a valid splitting bound if its
    Lambda is built from the SAME branch the flux actually runs. With six
    copies, six chances exist for one of them to drift, and a drift does not
    announce itself -- it appears as a relaxation-parameter sensitivity and a
    stalled physical step (`--lusgs-relax 0.7` fails, `1.0` accepts), which is
    what docs/tetra_cold_flow_status.md records.

    THE CONTRACT
    ------------
    Every consumer must obtain Lambda from this header. It is deliberately a
    single function with no overloads and no tunable: a caller cannot
    accidentally ask for a different wave speed, because there is only one.

    TWO CORRECTNESS NOTES, both about things that were previously implicit:

    1. The sound speed is NOT floored. Three previous copies wrote
       `std::max(soundSpeed(), 1.0)`; three wrote it bare. The floor changed the
       value in exactly the cells where it mattered and nowhere else, and since
       a sound speed below 1 m/s is physically impossible for the fluids in this
       project (cryogenic O2/C3H8 above 80 K gives 100-600 m/s), the floor could
       only ever mask a broken state rather than rescue a valid one. This header
       uses the bare value and requires the caller to have admitted the state
       first (see tetraAdmitsAcoustics in TetraAcousticAdmission.H). A NaN or
       non-positive sound speed is therefore a caller error and is reported as
       such instead of being silently absorbed.

    2. The normal must be the outward unit normal of the face, and u_n is
       signed. Only its magnitude enters Lambda, but the sign is retained in the
       signed accessor below because a future characteristic decomposition needs
       it and because it makes the direction-preserving tests possible.

    WHAT THIS FILE DOES NOT DO
    --------------------------
    It does not define the flux, the Jacobian, or the boundary trace. It only
    fixes the scalar that all of them must agree on. The flux Jacobian's
    eigenstructure is a separate (and still open) question -- see the note in
    src/numerics/flux/AUSMPlusUpFlux.C on the disabled M^4 term.

    Pure std C++, header-only. Depends on ConservedVars.H and
    EquilibriumThermoModel.H.
\*---------------------------------------------------------------------------*/

#include "ConservedVars.H"
#include "EquilibriumThermoModel.H"
#include "AUSMMachSplitting.H"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace Foam
{

// A face wave-speed evaluation, kept as a small struct rather than a bare double
// so that a caller which needs the signed normal velocity (characteristic work)
// does not have to recompute it -- recomputing is how the drift started.
struct FaceInvariant
{
    double normalVelocity;   // u_n = U . n_hat  [m/s], signed
    double soundSpeed;       // a  [m/s], as published by the admitted state
    double spectralRadius;   // |u_n| + a  [m/s]

    bool valid() const
    {
        return std::isfinite(normalVelocity)
            && std::isfinite(soundSpeed) && soundSpeed > 0.0
            && std::isfinite(spectralRadius) && spectralRadius > 0.0;
    }
};

// Evaluate the face invariant from an admitted (w, q) pair.
//
// Preconditions, checked: w.U has 3 components, n has 3 components, and the
// thermodynamic payload publishes a finite positive sound speed. An admitted
// state always satisfies the last one; a failure here means the caller skipped
// admission, so it throws rather than returning a plausible-looking number.
//
// This is the ONLY place the formula is written.
inline FaceInvariant tetraFaceInvariant(
    const PrimitiveVars& w,
    const thermo::EquilibriumState& q,
    const double normal[3])
{
    if (w.U.size() != 3 || !normal)
        throw std::runtime_error(
            "tetraFaceInvariant: velocity or normal is not three-dimensional");

    double un = 0.0;
    for (int d = 0; d < 3; ++d)
    {
        if (!std::isfinite(normal[d]) || !std::isfinite(w.U[d]))
            throw std::runtime_error(
                "tetraFaceInvariant: non-finite velocity or normal component");
        un += w.U[d] * normal[d];
    }

    const double a = q.acoustic.soundSpeed();
    if (!std::isfinite(a) || a <= 0.0)
        throw std::runtime_error(
            "tetraFaceInvariant: unadmitted acoustic state (sound speed " +
            std::to_string(a) + "); admit the state before evaluating the face "
            "invariant");

    FaceInvariant out;
    out.normalVelocity = un;
    out.soundSpeed = a;
    out.spectralRadius = std::abs(un) + a;
    return out;
}

// Convenience overload for the vector form used by the face tables.
inline FaceInvariant tetraFaceInvariant(
    const PrimitiveVars& w,
    const thermo::EquilibriumState& q,
    const std::vector<double>& normal)
{
    if (normal.size() != 3)
        throw std::runtime_error(
            "tetraFaceInvariant: normal must have three components");
    return tetraFaceInvariant(w, q, &normal[0]);
}

// ---------------------------------------------------------------------------
// The flux-branch tag.
//
// The M^4 low-Mach term is implemented but rejected (see AUSMMachSplitting.H
// and the note at the use site in AUSMPlusUpFlux.C). That decision changes the
// flux but NOT the way the wave speed is evaluated, so Lambda above is correct
// either way -- which is why the branch is not a parameter of
// tetraFaceInvariant.
//
// The tag is re-exported here so a consumer that already includes this header
// (the tetra kernel, its tests) can assert which flux branch its LU-SGS bound
// was built for, without pulling in the flux implementation. The enum itself
// lives in AUSMMachSplitting.H, shared with the flux, precisely so the two
// cannot disagree about what is running.
using AUSMMachSplitting = Foam::AUSMMachSplitting;

// The single declaration of which branch is live. Delegates to the flux's own
// default; change it there, not here.
inline constexpr AUSMMachSplitting activeAUSMMachSplitting()
{
    return DEFAULT_AUSM_MACH_SPLITTING;
}

} // namespace Foam

#endif
