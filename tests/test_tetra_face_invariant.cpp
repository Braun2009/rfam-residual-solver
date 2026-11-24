// Guard for the single-owner face invariant (TetraFaceInvariant.H).
//
// Purpose: this test exists to FAIL when the wave-speed formula is duplicated
// again. It does not test the flux. It tests that every consumer gets its
// spectral radius from one place, and that the one place behaves as documented.
//
// The specific regressions it catches:
//   1. A future caller reintroduces `std::max(soundSpeed(), 1.0)`. The
//      no-silent-floor case below fails, because the invariant must report the
//      admitted sound speed exactly and refuse a non-positive one.
//   2. A future caller writes its own |u_n| + a and its value drifts from this
//      one. The consistency case pins the arithmetic for a fixed input, so a
//      drifting copy fails wherever a comparison is made.
//   3. The normal is treated as non-unit, or the sign convention is flipped.
//      Both are pinned below.
//   4. The flux branch tag changes without the consumer being told. The branch
//      assertion at the end fails if the live tag is neither of the two known
//      values, forcing an explicit decision rather than a silent default.
#include "TetraFaceInvariant.H"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace Foam;

namespace
{
// The sound speed is sqrt(dp/drho_T) whenever the acoustic payload is
// structurally valid (EquilibriumAcousticProperties::soundSpeed). Building the
// payload directly keeps this test independent of CoolProp so it runs in every
// build configuration -- important, because the whole point is that it must
// FAIL in CI rather than silently skip if the formula is duplicated.
//
// Note the negative-speeds-are-unrepresentable property: because the accessor
// takes a square root, an acoustic payload can never publish a negative sound
// speed. A negative dp/drho_T is reported as a NaN sound speed. The invalid-input
// cases below therefore exercise (a) zero and (b) non-positive dp/drho_T, which
// are the only two ways a caller can actually hand this header a bad state.
thermo::EquilibriumState stateWithSoundSpeed(double a)
{
    thermo::EquilibriumState q;
    q.acoustic.closure = thermo::AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE;
    q.acoustic.evidence = thermo::AcousticEvidence::ANALYTIC_FROZEN_SINGLE_PHASE;
    q.acoustic.pressureDensityDerivative = a * a;
    return q;
}

// Build a state whose published sound speed is NaN (or non-finite) by giving a
// non-positive dp/drho_T, mirroring what the EOS does for an unstable root.
thermo::EquilibriumState stateWithDensityDerivative(double dpdrho)
{
    thermo::EquilibriumState q;
    q.acoustic.closure = thermo::AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE;
    q.acoustic.pressureDensityDerivative = dpdrho;
    return q;
}

int failures = 0;
void check(bool condition, const std::string& what)
{
    if (!condition) { std::cout << "FAIL: " << what << "\n"; ++failures; }
    else std::cout << "ok  : " << what << "\n";
}
} // namespace

int main()
{
    try
    {
        // ---- 1. arithmetic and sign convention --------------------------
        {
            PrimitiveVars w(2);
            w.U = {3.0, 4.0, 0.0};          // |U| = 5
            const double n[3] = {0.0, 0.0, 1.0};
            auto q = stateWithSoundSpeed(300.0);

            const auto inv = tetraFaceInvariant(w, q, n);
            check(inv.valid(), "invariant is valid for a normal state");
            check(std::abs(inv.normalVelocity - 0.0) < 1e-15,
                  "normal velocity uses the normal component only");
            check(std::abs(inv.soundSpeed - 300.0) < 1e-15,
                  "sound speed is passed through unmodified");
            check(std::abs(inv.spectralRadius - 300.0) < 1e-12,
                  "spectral radius is |u_n| + a");

            // Oblique normal: u_n = 3*cos(45) + 4*sin(45)
            const double diag[3] = {0.7071067811865476, 0.7071067811865476, 0.0};
            const auto inv2 = tetraFaceInvariant(w, q, diag);
            const double expected = std::abs(3.0 * diag[0] + 4.0 * diag[1]) + 300.0;
            check(std::abs(inv2.spectralRadius - expected) < 1e-10,
                  "oblique normal gives the projected velocity");
        }

        // ---- 2. the sign is retained (a future characteristic needs it) --
        {
            PrimitiveVars w(2);
            w.U = {0.0, 0.0, -12.5};
            const double n[3] = {0.0, 0.0, 1.0};
            auto q = stateWithSoundSpeed(180.0);
            const auto inv = tetraFaceInvariant(w, q, n);
            check(std::abs(inv.normalVelocity + 12.5) < 1e-14,
                  "negative normal velocity keeps its sign");
            check(std::abs(inv.spectralRadius - 192.5) < 1e-12,
                  "spectral radius uses the magnitude");
        }

        // ---- 3. no silent sound-speed floor -----------------------------
        // This is the regression that matters. Three former copies of the
        // formula floored a at 1 m/s; that floor made a broken state look like
        // a valid slow-sound-speed one. The invariant must expose it instead.
        {
            PrimitiveVars w(2);
            w.U = {0.0, 0.0, 0.0};
            const double n[3] = {0.0, 0.0, 1.0};

            bool threwOnZero = false;
            try { (void)tetraFaceInvariant(w, stateWithSoundSpeed(0.0), n); }
            catch (const std::exception&) { threwOnZero = true; }
            check(threwOnZero, "zero sound speed is rejected, not floored to 1");

            // An unstable root: dp/drho_T <= 0 publishes a NaN sound speed.
            // This is the case the former `max(a, 1.0)` floors turned from a
            // loud failure into a plausible-looking 1 m/s.
            bool threwOnUnstable = false;
            try { (void)tetraFaceInvariant(w, stateWithDensityDerivative(-1.0e6), n); }
            catch (const std::exception&) { threwOnUnstable = true; }
            check(threwOnUnstable,
                  "negative dp/drho_T (NaN sound speed) is rejected, not floored");

            bool threwOnZeroDerivative = false;
            try { (void)tetraFaceInvariant(w, stateWithDensityDerivative(0.0), n); }
            catch (const std::exception&) { threwOnZeroDerivative = true; }
            check(threwOnZeroDerivative, "zero dp/drho_T is rejected");

            // A small but positive sound speed must pass through EXACTLY.
            // If someone re-adds a 1.0 floor, this fails.
            const auto tiny = tetraFaceInvariant(w, stateWithSoundSpeed(0.25), n);
            check(std::abs(tiny.soundSpeed - 0.25) < 1e-15,
                  "a small positive sound speed is NOT raised to a floor");
            check(std::abs(tiny.spectralRadius - 0.25) < 1e-15,
                  "spectral radius follows the unfloored value");
        }

        // ---- 4. rejected inputs -----------------------------------------
        {
            PrimitiveVars bad(2); bad.U = {1.0, 2.0}; // only two components
            const double n[3] = {0.0, 0.0, 1.0};
            bool threw = false;
            try { (void)tetraFaceInvariant(bad, stateWithSoundSpeed(100.0), n); }
            catch (const std::exception&) { threw = true; }
            check(threw, "a 2-component velocity is rejected");

            PrimitiveVars w(2); w.U = {0.0, 0.0, 0.0};
            bool threwNormal = false;
            try {
                const std::vector<double> shortNormal{0.0, 1.0};
                (void)tetraFaceInvariant(w, stateWithSoundSpeed(100.0), shortNormal);
            } catch (const std::exception&) { threwNormal = true; }
            check(threwNormal, "a 2-component normal is rejected");
        }

        // ---- 5. the flux branch tag is a decision, never a default ------
        // If the M^4 term is restored, activeAUSMMachSplitting() must change in
        // the same commit as AUSMPlusUpFlux::computeFluxCached. This assertion
        // cannot detect that by itself, but it forces the value to be read and
        // named, so the change appears in the diff of this test as well.
        {
            const auto branch = activeAUSMMachSplitting();
            const std::string name = ausmMachSplittingName(branch);
            check(name != std::string("unknown"),
                  "the active flux branch has a known name");
            check(branch == AUSMMachSplitting::VELOCITY_SPLIT_ONLY
                  || branch == AUSMMachSplitting::FOURTH_ORDER_LOW_MACH,
                  "the active flux branch is one of the two implemented splittings");
            std::cout << "info: active AUSM splitting = " << name
                      << " ; face invariant follows it from one owner\n";
        }

        if (failures)
        {
            std::cout << "\n" << failures << " check(s) FAILED\n";
            return 1;
        }
        std::cout << "\nTetraFaceInvariant: all checks passed\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "unexpected exception: " << e.what() << "\n";
        return 1;
    }
}
