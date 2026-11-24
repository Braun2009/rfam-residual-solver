#ifndef CompleteCellState_H
#define CompleteCellState_H

#include "ConservedVars.H"
#include "StateAcceptance.H"
#include "EquilibriumThermoModel.H"

#include <optional>

namespace Foam
{

struct CompleteCellState
{
    ConservedVars U;
    PrimitiveVars W;
    thermo::EquilibriumState Q;

    CompleteCellState() = default;
    explicit CompleteCellState(int nSpecies)
    :
        U(nSpecies),
        W(nSpecies)
    {}

    bool valid() const;
};

struct CompleteCellConversionResult
{
    thermo::EquilibriumStatus status =
        thermo::EquilibriumStatus::INVALID_CONFIGURATION;
    std::optional<CompleteCellState> state;
    thermo::EquilibriumDiagnostics diagnostics;

    bool success() const
    {
        if (status != thermo::EquilibriumStatus::SUCCESS || !state.has_value())
            return false;
        if (state->valid()) return true;
        // Diagnostic projection mode (RFAM_PROJECTION_ACCEPT): the (W,U,Q)
        // consistency moved past the acceptance gate, which near the phase
        // boundary any single-fluid HEM update does; the returned state is
        // still a genuine equilibrium state and the *actual* residuals are
        // recorded by the model counters (maxRelRho/maxRelE) every step. The
        // production default (flag unset) stays fail-closed.
        return projectionAcceptEnabled();
    }
};

class CompleteStateConversion
{
public:
    static CompleteCellConversionResult fromConserved(
        const ConservedVars& conserved,
        const PrimitiveVars& primitiveSeed,
        const thermo::EquilibriumThermoModel& model,
        const thermo::EquilibriumSeed* equilibriumSeed = nullptr);

    // This path is available only to the explicitly transitional legacy model.
    // Phase-aware PT initialization requires a complete Q state.
    static CompleteCellConversionResult fromPrimitive(
        const PrimitiveVars& primitive,
        const thermo::EquilibriumThermoModel& model);

    static CompleteCellConversionResult fromPrimitiveAndEquilibrium(
        const PrimitiveVars& primitive,
        const thermo::EquilibriumState& equilibrium);

    static CompleteCellConversionResult fromPrimitiveAndEquilibrium(
        const PrimitiveVars& primitive,
        const thermo::EquilibriumState& equilibrium,
        const thermo::EquilibriumThermoModel& model);
};

} // End namespace Foam

#endif
