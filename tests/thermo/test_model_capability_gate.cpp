#include <iostream>

#include "ModelCapabilities.H"

using namespace Foam;

int main()
{
    static_assert(!ModelCapabilities::twoPhaseFlash,
        "enable only after TP/PH/rho-E flash validation");
    static_assert(!ModelCapabilities::phaseAwareEnergyInversion,
        "enable only after two-phase conserved-state validation");
    static_assert(!ModelCapabilities::latentHeatClosure,
        "enable only after latent-heat validation");
    static_assert(!ModelCapabilities::phaseAwareSoundSpeed,
        "two-phase acoustics stay diagnostic (UNVERIFIED per 2026-08-03 audit); "
        "production AUSM+/CFL must reject HEM acoustics; diagnostics use "
        "RFAM_ALLOW_DIAGNOSTIC_HEM_ACOUSTICS override");
    static_assert(!ModelCapabilities::phaseAwareBoundaryStates,
        "enable only after every production boundary publishes complete states");
    static_assert(!ModelCapabilities::phaseAwareCharacteristicBoundary,
        "enable only after phase-aware characteristic validation");
    static_assert(!ModelCapabilities::denseCryogenicTransport,
        "enable only after cryogenic transport validation");
    static_assert(!ModelCapabilities::phaseAwareChemistry,
        "enable only after phase-restricted chemistry validation");
    static_assert(!ModelCapabilities::condensedCarbonTransport,
        "enable only after dilute soot transport validation");
    static_assert(!ModelCapabilities::condensedCarbonKinetics,
        "enable only after sourced soot kinetics validation");
    static_assert(!ModelCapabilities::condensedCarbonWallDeposition,
        "enable only after wall deposition validation");
    static_assert(!ModelCapabilities::nonReactingSectorProductionReady(),
        "the production nonreacting sector must remain blocked");
    static_assert(!ModelCapabilities::reactingSectorProductionReady(),
        "the production reacting sector must remain blocked");
    static_assert(!ModelCapabilities::condensedCarbonProductionReady(),
        "condensed-carbon production must remain blocked");

    std::cout << "MODEL_CAPABILITY_GATE_BLOCKED\n";
    return 0;
}
