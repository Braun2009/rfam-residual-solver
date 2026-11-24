#ifndef AUSMMachSplitting_H
#define AUSMMachSplitting_H

/*---------------------------------------------------------------------------*\
    AUSMMachSplitting - which Mach-number splitting the AUSM+-up flux runs.

    WHY THIS IS ITS OWN HEADER
    --------------------------
    The flux lives in src/numerics/, the face-invariant owner lives in
    src/mesh/. Both need to name the splitting -- the flux to select it, the
    invariant to tag which branch its LU-SGS bound was built for. Putting the
    enum in either of those headers would create a dependency in the wrong
    direction (numerics must not include mesh). It is small and dependency-free
    so it lives here.

    THE TWO BRANCHES
    ----------------
    AUSM+-up evaluates the interface Mach number as

        M_half = <splitting>(M_L) + <splitting-minus>(M_R) + M_p

    and <splitting> is what this enum selects.

    VELOCITY_SPLIT_ONLY        M_half = M+(M_L) + M-(M_R) + M_p
        The bare Liou (2006) eq. (24) polynomials. This is what runs in
        production today. Its failing at low Mach is that for an acoustic pair
        M_L = +e, M_R = -e it returns exactly zero:

            M+(+e) + M-(-e) = 0   for all e

        so the mass flux carries no dissipation on the acoustic mode. Checked
        numerically: the expression is exactly 0.0 for e in [1e-6, 0.1].

    FOURTH_ORDER_LOW_MACH      M_half = M4+(M_L) + M4-(M_R) + M_p
        Liou (2006) eqs. (29)-(30), implemented in AUSMPlusUpFlux::M4_plus /
        M4_minus.

        STATUS: REJECTED, and now with a measured reason rather than a hunch.
        tests/test_acoustic_decay_1d.cpp reads two properties straight off the
        production flux (no mesh, no EOS, no time integration):

          P1 uniform flow:   M4+(M) + M4-(M) == M     must hold exactly
          P2 acoustic pair:  M4+(+e) + M4-(-e) != 0   must be non-zero

        Measured 2026-10-09, ideal gas a0 = 340 m/s:

          branch                  P1                       P2
          ----------------------  -----------------------  ------------
          VELOCITY_SPLIT_ONLY     PASS (exact, ~1e-16)     0  (gap)
          FOURTH_ORDER_LOW_MACH   FAIL, 3.25x -> +225 %    0  (gap)

        So the M^4 branch breaks the wave speed AND adds no acoustic damping --
        strictly worse than the default. The +65 % wave-speed shift recorded on
        2026-08-04 is this P1 violation, seen through a time-marching measurement.

        The mechanism: both M4 functions multiply by (1 - 16 M+ M-) with the
        SAME sign, so the product term does not cancel between the two halves.
        Liou's eq. (29) uses opposite signs, but even that form does not restore
        the identity (checked numerically: 1.69 M as M -> 0), so the correct
        expression still has to be established from the reference. P1 and P2 in
        the test above are the gate for it.

    HOW TO ENABLE A BRANCH
    ----------------------
    AUSMPlusUpFlux::setMachSplitting(...). Default is VELOCITY_SPLIT_ONLY, so a
    caller that does nothing gets the shipped behaviour.
\*---------------------------------------------------------------------------*/

namespace Foam
{

enum class AUSMMachSplitting
{
    VELOCITY_SPLIT_ONLY,      // M+ + M-        (production default)
    FOURTH_ORDER_LOW_MACH     // M4+ + M4-      (implemented, currently rejected)
};

// The branch the solver runs unless a caller overrides it. Change this only in
// the same commit that changes AUSMPlusUpFlux::computeFluxCached.
inline constexpr AUSMMachSplitting DEFAULT_AUSM_MACH_SPLITTING =
    AUSMMachSplitting::VELOCITY_SPLIT_ONLY;

inline const char* ausmMachSplittingName(AUSMMachSplitting s)
{
    switch (s)
    {
        case AUSMMachSplitting::VELOCITY_SPLIT_ONLY:   return "M+ + M- + Mp";
        case AUSMMachSplitting::FOURTH_ORDER_LOW_MACH: return "M4+ + M4- + Mp";
    }
    return "unknown";
}

} // namespace Foam

#endif
