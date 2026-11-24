#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "ConservedVars.H"
#include "Domain.H"

using namespace Foam;

namespace
{

int failures = 0;
constexpr double universalGasConstant = 8.314462618;

void check(const std::string& name, bool ok)
{
    std::cout << "  [" << (ok ? "PASS" : "FAIL") << "] " << name << '\n';
    if (!ok) ++failures;
}

bool near(double lhs, double rhs, double relativeTolerance = 1.0e-8)
{
    return std::abs(lhs - rhs)
        <= relativeTolerance*std::max({std::abs(lhs), std::abs(rhs), 1.0});
}

PrimitiveVars oxygenState()
{
    PrimitiveVars state(1);
    state.p = 4.0e6;
    state.T = 300.0;
    state.U = {12.0, -3.0, 1.0};
    state.Y = {1.0};
    state.Z = 0.2;
    state.C = 0.1;
    return state;
}

thermo::EquilibriumState singlePhaseState(
    const std::vector<double>& composition = {1.0})
{
    thermo::EquilibriumState state;
    state.region = thermo::EquilibriumRegion::SINGLE_VAPOR;
    state.overallMassFractions = composition;
    state.p = 1.0e6;
    state.T = 350.0;
    state.rho = 12.0;
    state.h = 4.2e5;
    state.e = state.h - state.p/state.rho;
    state.acoustic.closure =
        thermo::AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE;
    state.acoustic.pressureDensityDerivative = 9.0e4;

    thermo::EquilibriumPhaseProperties phase;
    phase.massFractions = composition;
    phase.moleFractions = composition;
    phase.molecularWeight = composition.size() == 1 ? 0.032 : 0.044;
    phase.rho = state.rho;
    phase.h = state.h;
    phase.e = state.e;
    phase.compressibility = state.p*phase.molecularWeight
        /(phase.rho*universalGasConstant*state.T);
    state.singlePhase = phase;
    return state;
}

thermo::EquilibriumState twoPhaseState(double vaporMoleFraction = 0.25)
{
    thermo::EquilibriumState state;
    state.region = thermo::EquilibriumRegion::TWO_PHASE;
    state.p = 1.0e6;
    state.T = 300.0;

    thermo::EquilibriumPhaseProperties liquid;
    liquid.massFractions = {0.95, 0.05};
    liquid.moleFractions = {0.90, 0.10};
    liquid.molecularWeight = 0.032;
    liquid.rho = 1000.0;
    liquid.h = 1.0e5;
    liquid.e = liquid.h - state.p/liquid.rho;
    liquid.s = 1200.0;
    liquid.compressibility = state.p*liquid.molecularWeight
        /(liquid.rho*universalGasConstant*state.T);
    liquid.entropyAvailable = true;

    thermo::EquilibriumPhaseProperties vapor;
    vapor.massFractions = {0.20, 0.80};
    vapor.moleFractions = {0.25, 0.75};
    vapor.molecularWeight = 0.028;
    vapor.rho = 10.0;
    vapor.h = 5.0e5;
    vapor.e = vapor.h - state.p/vapor.rho;
    vapor.s = 4200.0;
    vapor.compressibility = state.p*vapor.molecularWeight
        /(vapor.rho*universalGasConstant*state.T);
    vapor.entropyAvailable = true;

    const double totalMolarMass =
        (1.0 - vaporMoleFraction)*liquid.molecularWeight
        + vaporMoleFraction*vapor.molecularWeight;
    const double vaporMassFraction =
        vaporMoleFraction*vapor.molecularWeight/totalMolarMass;
    const double liquidMassFraction = 1.0 - vaporMassFraction;
    const double specificVolume = liquidMassFraction/liquid.rho
        + vaporMassFraction/vapor.rho;

    state.overallMassFractions = {
        liquidMassFraction*liquid.massFractions[0]
            + vaporMassFraction*vapor.massFractions[0],
        liquidMassFraction*liquid.massFractions[1]
            + vaporMassFraction*vapor.massFractions[1]
    };
    state.rho = 1.0/specificVolume;
    state.h = liquidMassFraction*liquid.h + vaporMassFraction*vapor.h;
    state.e = liquidMassFraction*liquid.e + vaporMassFraction*vapor.e;
    state.s = liquidMassFraction*liquid.s + vaporMassFraction*vapor.s;
    state.acoustic.closure =
        thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM;
    state.acoustic.pressureDensityDerivative = 4.0e4;
    state.vaporMoleFraction = vaporMoleFraction;
    state.vaporMassFraction = vaporMassFraction;
    state.vaporVolumeFraction =
        (vaporMassFraction/vapor.rho)/specificVolume;
    state.entropyAvailable = true;
    state.liquid = liquid;
    state.vapor = vapor;
    return state;
}

class InvalidSuccessModel final : public thermo::EquilibriumThermoModel
{
    std::vector<std::string> species_ = {"O2"};
    std::string id_ = "invalid-success-test-model";

public:
    thermo::EquilibriumResult recoverRhoE(
        double,
        double,
        const std::vector<double>&,
        const thermo::EquilibriumSeed*) const override
    {
        thermo::EquilibriumResult result;
        result.status = thermo::EquilibriumStatus::SUCCESS;
        result.state = singlePhaseState();
        result.state->h += 1000.0;
        return result;
    }

    const std::string& modelId() const override { return id_; }
    const std::vector<std::string>& speciesNames() const override
    { return species_; }
    EnergyConvention energyConvention() const override
    {
        return EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT;
    }
};

} // namespace

int main()
{
    std::cout << "=== StateConversion and Domain fail-close contract ===\n";
    const std::vector<std::string> species = {"O2"};

    PrimitiveVars primitive = oxygenState();
    ConservedVars conserved(1);
    check("normal P2C succeeds",
        StateConversion::primitiveToConserved(
            primitive, conserved, &species));

    PrimitiveVars recovered(1);
    recovered.p = 12345.0;
    recovered.T = 350.0;
    check("normal C2P succeeds",
        StateConversion::conservedToPrimitive(
            conserved, recovered, &species));
    check("normal pressure round trip", near(recovered.p, primitive.p, 1.0e-7));
    check("normal temperature round trip", near(recovered.T, primitive.T, 1.0e-8));
    check("normal velocity round trip", near(recovered.U[0], primitive.U[0]));

    {
        ConservedVars invalid(1);
        PrimitiveVars output(1);
        output.p = 22222.0;
        output.T = 444.0;
        check("zero-density C2P fails",
            !StateConversion::conservedToPrimitive(invalid, output, &species));
        check("zero-density C2P does not overwrite output",
            output.p == 22222.0 && output.T == 444.0);
    }

    {
        ConservedVars invalid = conserved;
        invalid.rhoE = std::numeric_limits<double>::quiet_NaN();
        PrimitiveVars output(1);
        output.p = 33333.0;
        output.T = 555.0;
        check("NaN-energy C2P fails",
            !StateConversion::conservedToPrimitive(invalid, output, &species));
        check("NaN-energy C2P does not overwrite output",
            output.p == 33333.0 && output.T == 555.0);
    }

    {
        PrimitiveVars invalid = primitive;
        invalid.p = -1.0;
        ConservedVars output(1);
        output.rho = 42.0;
        output.rhoE = 43.0;
        check("invalid P2C fails",
            !StateConversion::primitiveToConserved(invalid, output, &species));
        check("invalid P2C does not overwrite output",
            output.rho == 42.0 && output.rhoE == 43.0);
    }

    {
        PrimitiveVars invalid = primitive;
        invalid.Z = -1.0e-12;
        ConservedVars output(1);
        output.rho = 52.0;
        output.rhoE = 53.0;
        check("out-of-range primitive scalar fails",
            !StateConversion::primitiveToConserved(invalid, output, &species));
        check("invalid primitive scalar does not overwrite output",
            output.rho == 52.0 && output.rhoE == 53.0);
    }

    {
        ConservedVars invalid = conserved;
        invalid.rhoZ = std::nextafter(invalid.rho,
            std::numeric_limits<double>::infinity());
        PrimitiveVars output(1);
        output.p = 44444.0;
        output.T = 666.0;
        check("out-of-range conserved scalar fails",
            !StateConversion::conservedToPrimitive(invalid, output, &species));
        check("invalid conserved scalar does not overwrite output",
            output.p == 44444.0 && output.T == 666.0);
    }

    {
        thermo::LegacySinglePhaseThermoModel unknown({"NOT_A_SPECIES"});
        thermo::EquilibriumResult unknownInitialize;
        thermo::EquilibriumResult unknownRecover;
        bool unknownDidNotThrow = true;
        try
        {
            unknownInitialize = unknown.initializePTY(1.0e5, 300.0, {1.0});
            unknownRecover = unknown.recoverRhoE(1.0, 1.0e5, {1.0});
        }
        catch (...)
        {
            unknownDidNotThrow = false;
        }
        check("unknown legacy species fails without throwing",
            unknownDidNotThrow
                && unknownInitialize.status
                    == thermo::EquilibriumStatus::INVALID_CONFIGURATION
                && unknownRecover.status
                    == thermo::EquilibriumStatus::INVALID_CONFIGURATION
                && !unknownInitialize.state.has_value()
                && !unknownRecover.state.has_value()
                && !unknownInitialize.diagnostics.modelId.empty()
                && !unknownRecover.diagnostics.modelId.empty());

        thermo::LegacySinglePhaseThermoModel duplicate({"O2", "O2"});
        thermo::EquilibriumResult duplicateInitialize;
        thermo::EquilibriumResult duplicateRecover;
        bool duplicateDidNotThrow = true;
        try
        {
            duplicateInitialize = duplicate.initializePTY(
                1.0e5, 300.0, {0.5, 0.5});
            duplicateRecover = duplicate.recoverRhoE(
                1.0, 1.0e5, {0.5, 0.5});
        }
        catch (...)
        {
            duplicateDidNotThrow = false;
        }
        check("duplicate legacy species fail without throwing",
            duplicateDidNotThrow
                && duplicateInitialize.status
                    == thermo::EquilibriumStatus::INVALID_CONFIGURATION
                && duplicateRecover.status
                    == thermo::EquilibriumStatus::INVALID_CONFIGURATION
                && !duplicateInitialize.state.has_value()
                && !duplicateRecover.state.has_value());

        thermo::LegacySinglePhaseThermoModel oxygenModel({"O2"});
        const thermo::EquilibriumResult wrongSize =
            oxygenModel.initializePTY(1.0e5, 300.0, {0.5, 0.5});
        check("valid legacy model reports composition size mismatch",
            wrongSize.status == thermo::EquilibriumStatus::SPECIES_MISMATCH
                && !wrongSize.state.has_value());
    }

    {
        const thermo::EquilibriumState single = singlePhaseState();
        check("closed generic single-phase Q is valid", single.valid());
        auto inconsistent = single;
        inconsistent.singlePhase->rho *= 1.01;
        check("single-phase top-level/phase mismatch is rejected",
            !inconsistent.valid());
        inconsistent = single;
        inconsistent.e += 1000.0;
        check("single-phase energy identity mismatch is rejected",
            !inconsistent.valid());
        inconsistent = single;
        inconsistent.vaporMassFraction = 0.0;
        check("single-phase stale vapor fraction is rejected",
            !inconsistent.valid());
        inconsistent = single;
        inconsistent.vaporMoleFraction =
            std::numeric_limits<double>::infinity();
        check("single-phase noncanonical vapor fraction is rejected",
            !inconsistent.valid());
        inconsistent = single;
        inconsistent.s = std::numeric_limits<double>::infinity();
        check("unavailable entropy requires a canonical NaN payload",
            !inconsistent.valid());
        inconsistent = single;
        inconsistent.acoustic.pressureDensityDerivative =
            std::numeric_limits<double>::infinity();
        check("non-finite acoustic derivative is rejected", !inconsistent.valid());
        inconsistent = single;
        inconsistent.acoustic.pressureDensityDerivative = 0.0;
        check("non-positive acoustic derivative is rejected", !inconsistent.valid());
        inconsistent = single;
        inconsistent.region = static_cast<thermo::EquilibriumRegion>(999);
        check("unknown equilibrium-region enum is rejected",
            !inconsistent.valid());
        inconsistent = single;
        inconsistent.acoustic.closure =
            static_cast<thermo::AcousticClosure>(999);
        check("unknown acoustic-closure enum is rejected",
            !inconsistent.valid());
    }

    {
        const thermo::EquilibriumState twoPhase = twoPhaseState();
        check("closed generic two-phase Q is valid", twoPhase.valid());
        check("zero-quality two-phase Q is valid", twoPhaseState(0.0).valid());
        check("unit-quality two-phase Q is valid", twoPhaseState(1.0).valid());

        auto inconsistent = twoPhase;
        inconsistent.vaporMassFraction += 1.0e-3;
        check("mole/mass vapor-fraction mismatch is rejected",
            !inconsistent.valid());
        inconsistent = twoPhase;
        inconsistent.vaporVolumeFraction += 1.0e-3;
        check("mass/volume vapor-fraction mismatch is rejected",
            !inconsistent.valid());
        inconsistent = twoPhase;
        inconsistent.overallMassFractions[0] += 1.0e-3;
        inconsistent.overallMassFractions[1] -= 1.0e-3;
        check("two-phase material mismatch is rejected",
            !inconsistent.valid());
        inconsistent = twoPhase;
        inconsistent.liquid->e += 1000.0;
        check("two-phase phase-energy mismatch is rejected",
            !inconsistent.valid());
        inconsistent = twoPhase;
        inconsistent.s += 10.0;
        check("two-phase entropy lever-rule mismatch is rejected",
            !inconsistent.valid());
        inconsistent = twoPhase;
        inconsistent.vapor->compressibility *= 1.01;
        check("phase compressibility mismatch is rejected",
            !inconsistent.valid());
        inconsistent = twoPhase;
        inconsistent.acoustic.closure =
            thermo::AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE;
        check("two-phase Q requires homogeneous-equilibrium acoustics",
            !inconsistent.valid());

        PrimitiveVars phasePrimitive(2);
        phasePrimitive.p = twoPhase.p;
        phasePrimitive.T = twoPhase.T;
        phasePrimitive.U = {3.0, -2.0, 1.0};
        phasePrimitive.Y = twoPhase.overallMassFractions;
        phasePrimitive.Z = 0.3;
        phasePrimitive.C = 0.4;
        const CompleteCellConversionResult converted =
            CompleteStateConversion::fromPrimitiveAndEquilibrium(
                phasePrimitive, twoPhase);
        check("closed primitive/Q pair produces a complete state",
            converted.success() && converted.state->valid());

        auto invalidScalarState = *converted.state;
        invalidScalarState.U.rhoZ =
            std::nextafter(0.0, -std::numeric_limits<double>::infinity());
        check("complete state rejects negative conserved mixture fraction",
            !invalidScalarState.valid());
        invalidScalarState = *converted.state;
        invalidScalarState.W.C = 1.0;
        invalidScalarState.U.rhoC = std::nextafter(
            invalidScalarState.U.rho,
            std::numeric_limits<double>::infinity());
        check("complete state rejects over-unity conserved progress variable",
            !invalidScalarState.valid());

        const thermo::EquilibriumState zeroSpecies =
            singlePhaseState({0.0, 1.0});
        phasePrimitive.p = zeroSpecies.p;
        phasePrimitive.T = zeroSpecies.T;
        phasePrimitive.Y = {-5.0e-9, 1.0 + 5.0e-9};
        const CompleteCellConversionResult negativeSpecies =
            CompleteStateConversion::fromPrimitiveAndEquilibrium(
                phasePrimitive, zeroSpecies);
        check("tolerance-sized negative species is rejected",
            !negativeSpecies.success() && !negativeSpecies.state.has_value());
    }

    {
        InvalidSuccessModel invalidModel;
        const CompleteCellConversionResult rejected =
            CompleteStateConversion::fromConserved(
                conserved, primitive, invalidModel);
        check("model SUCCESS with invalid Q maps to residual failure",
            rejected.status
                == thermo::EquilibriumStatus::RESIDUAL_CONTRACT_FAILED);
        check("invalid model SUCCESS exposes no complete state",
            !rejected.state.has_value());
    }

    {
        Domain domain(3, 0.3, 1);
        domain.setSpecies(species);
        for (int i = 0; i < domain.N(); ++i) domain.stagedW(i) = primitive;
        check("domain initialization succeeds", domain.primitiveToConserved());
        check("repaired legacy model has a new checkpoint identity",
            domain.thermoModelId()
                == "legacy-real-fluid-single-phase-v3:"
                   "nasa7;species=1;2:O2;");
        check("legacy Q publishes frozen single-phase acoustics",
            domain.Q(0).acoustic.closure
                == thermo::AcousticClosure::FROZEN_COMPOSITION_SINGLE_PHASE
                && domain.Q(0).acoustic.valid());
        check("legacy Q publishes closed caloric and derivative state",
            domain.Q(0).thermodynamics.has_value()
                && domain.Q(0).thermodynamics->valid());
        thermo::EquilibriumState invalidPayload = domain.Q(0);
        invalidPayload.thermodynamics->dedrho_T += 1.0e6;
        check("Q rejects a derivative payload off the committed surface",
            !invalidPayload.valid());
        invalidPayload = domain.Q(0);
        invalidPayload.thermodynamics->version += 1;
        check("Q rejects an unknown thermodynamic payload version",
            !invalidPayload.valid());
        invalidPayload = domain.Q(0);
        invalidPayload.thermodynamics->closure =
            thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM;
        check("Q rejects a payload/acoustic closure mismatch",
            !invalidPayload.valid());
        invalidPayload = domain.Q(0);
        invalidPayload.thermodynamics->dedT_rho *= 1.01;
        check("Q rejects dedT inconsistent with Cv",
            !invalidPayload.valid());
        invalidPayload = domain.Q(0);
        invalidPayload.thermodynamics->Cp *= 1.01;
        check("Q rejects a Cp-Cv identity mismatch",
            !invalidPayload.valid());
        invalidPayload = domain.Q(0);
        invalidPayload.acoustic.pressureDensityDerivative *= 1.01;
        check("Q rejects payload/acoustic derivative mismatch",
            !invalidPayload.valid());
        check("legacy PTY initialization preserves the prescribed cache",
            domain.W(0).p == primitive.p
                && domain.W(0).T == primitive.T
                && domain.W(0).U == primitive.U
                && domain.W(0).Y == primitive.Y);
        const double oldEnergy = domain.U(1).rhoE;
        const double oldPressure = domain.W(1).p;
        const double oldTemperature = domain.W(1).T;

        domain.beginStep();
        domain.Unew(1).rhoE = 1.0e20;
        domain.commit();

        check("domain trips on un-invertible cell", domain.tripped());
        check("failed commit preserves conserved field",
            domain.U(1).rhoE == oldEnergy);
        check("failed commit preserves primitive field",
            domain.W(1).p == oldPressure &&
            domain.W(1).T == oldTemperature);
    }

    {
        Domain domain(2, 0.2, 1);
        domain.setSpecies(species);
        for (int i = 0; i < domain.N(); ++i) domain.stagedW(i) = primitive;
        check("second domain initialization succeeds",
            domain.primitiveToConserved());
        const double oldDensity = domain.U(0).rho;
        domain.stagedW(1).Y[0] = 0.5;
        check("domain P2C rejects bad composition",
            !domain.primitiveToConserved());
        check("failed domain P2C preserves prior conserved field",
            domain.U(0).rho == oldDensity);
    }

    {
        Domain domain(3, 0.3, 1);
        domain.setSpecies(species);
        for (int i = 0; i < domain.N(); ++i) domain.stagedW(i) = primitive;
        check("complete-field rejection setup succeeds",
            domain.primitiveToConserved());

        std::vector<PrimitiveVars> primitives;
        std::vector<thermo::EquilibriumState> equilibrium;
        std::vector<double> oldRhoE;
        std::vector<double> oldCandidateRhoE;
        std::vector<double> oldPressure;
        std::vector<double> oldEnthalpy;
        for (int i = 0; i < domain.N(); ++i)
        {
            primitives.push_back(domain.W(i));
            equilibrium.push_back(domain.Q(i));
            oldRhoE.push_back(domain.U(i).rhoE);
            oldCandidateRhoE.push_back(domain.Unew(i).rhoE);
            oldPressure.push_back(domain.W(i).p);
            oldEnthalpy.push_back(domain.Q(i).h);
        }
        equilibrium[1].h += 1000.0;
        check("Domain rejects one internally inconsistent Q",
            !domain.replaceCompleteField(primitives, equilibrium));
        bool unchanged = true;
        for (int i = 0; i < domain.N(); ++i)
        {
            unchanged = unchanged
                && domain.U(i).rhoE == oldRhoE[i]
                && domain.Unew(i).rhoE == oldCandidateRhoE[i]
                && domain.W(i).p == oldPressure[i]
                && domain.Q(i).h == oldEnthalpy[i];
        }
        check("failed Q replacement preserves the complete field", unchanged);
        check("failed Q replacement reports the first invalid cell",
            domain.failure().cell == 1
                && domain.failure().equilibriumStatus
                    == thermo::EquilibriumStatus::RESIDUAL_CONTRACT_FAILED);
    }

    std::cout << (failures == 0
        ? "STATE_CONVERSION_CONTRACT_PASS\n"
        : "STATE_CONVERSION_CONTRACT_FAIL\n");
    return failures == 0 ? 0 : 1;
}
