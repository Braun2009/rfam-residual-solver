#include <iostream>
#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "PhaseEquilibrium/HEMStateOwnershipContract.H"
#include "PhaseEquilibrium/PREquilibriumThermoModel.H"
#include "PhaseEquilibrium/PRPHFlash.H"
#include "PhaseEquilibrium/PRPhaseCaloric.H"
#include "PhaseEquilibrium/PRTPFlash.H"
#include "MixingRules/SpeciesCriticalDB.H"

using namespace Foam;
using namespace Foam::phaseEquilibrium;

namespace
{

int checks = 0;
int failures = 0;

void check(bool condition, const std::string& message)
{
    ++checks;
    if (!condition)
    {
        ++failures;
        std::cerr << "FAIL: " << message << "\n";
    }
}

bool close(
    double lhs, double rhs, double relative = 1.0e-8,
    double absolute = 1.0e-8)
{
    return std::isfinite(lhs) && std::isfinite(rhs)
        && std::abs(lhs - rhs) <= absolute + relative
            *std::max(std::abs(lhs), std::abs(rhs));
}

ConservedPhaseProperties conserved(
    const PRPhaseCaloricProperties& phase)
{
    ConservedPhaseProperties result;
    result.molecularWeight = phase.molecularWeight;
    result.rho = phase.rho;
    result.h = phase.h;
    result.e = phase.e;
    result.s = phase.s;
    result.massFractions = phase.massFractions;
    return result;
}

std::optional<double> equilibriumEnthalpy(
    const PRTPFlash& flash,
    const PRPhaseCaloric& caloric,
    double pressure,
    double temperature,
    const std::vector<double>& massFractions)
{
    const TPFlashResult tp = flash.flashMassFractions(
        pressure, temperature, massFractions);
    if (!tp.valid()) return std::nullopt;

    if (tp.region == EquilibriumRegion::TWO_PHASE)
    {
        const PRPhaseCaloricResult liquid = caloric.evaluate(
            pressure, temperature, tp.liquidMoleFractions,
            PhaseRoot::LIQUID);
        const PRPhaseCaloricResult vapor = caloric.evaluate(
            pressure, temperature, tp.vaporMoleFractions,
            PhaseRoot::VAPOR);
        if (!liquid.success() || !vapor.success()) return std::nullopt;
        const ConservedMixtureResult mixture = assembleConservedPhaseMixture(
            pressure, tp.vaporMoleFraction,
            conserved(*liquid.state), conserved(*vapor.state),
            tp.feedMassFractions);
        return mixture.success()
            ? std::optional<double>(mixture.state->h) : std::nullopt;
    }

    if (!tp.singlePhaseState.valid()
        || tp.singlePhaseMoleFractions.empty())
        return std::nullopt;
    const PRPhaseCaloricResult single = caloric.evaluate(
        pressure, temperature, tp.singlePhaseMoleFractions,
        tp.singlePhaseState.requestedRoot);
    return single.success()
        ? std::optional<double>(single.state->h) : std::nullopt;
}

RhoEFlashOptions mixingRhoEOptions()
{
    RhoEFlashOptions options;
    options.minimumPressure = 1.0e6;
    options.maximumPressure = 8.0e6;
    options.minimumTemperature = 80.0;
    options.maximumTemperature = 700.0;
    options.initialPressure = 4.0e6;
    options.initialTemperature = 240.0;
    options.pressureScanIntervals = 128;
    options.maximumIterations = 160;
    options.finiteDifferenceStep = 1.0e-4;
    options.logStepLimit = 0.75;
    options.densityAbsoluteTolerance = 1.0e-7;
    options.densityRelativeTolerance = 2.0e-8;
    options.energyAbsoluteTolerance = 1.0e-4;
    options.energyRelativeTolerance = 2.0e-10;
    return options;
}

BinaryInteractionApplicability vortexScreenEnvelope()
{
    BinaryInteractionApplicability envelope;
    envelope.minimumTemperature = 110.0;
    envelope.maximumTemperature = 600.0;
    envelope.minimumPressure = 1.0e5;
    envelope.maximumPressure = 8.0e6;
    return envelope;
}

} // namespace

int main()
{
    const BinaryInteractionApplicability envelope = vortexScreenEnvelope();
    check(envelope.valid(), "Box B envelope is structurally valid");

    // The legacy defaults are a numerical convenience only.  Their O2-C3H8
    // entry must not turn a mixing-box test into a validated VLE claim.
    const BinaryInteractionDB defaults = BinaryInteractionDB::createDefault();
    const BinaryInteractionReadinessAudit defaultAudit =
        defaults.auditProductionReadiness({"O2", "C3H8"}, envelope);
    check(!defaultAudit.productionReady(),
        "legacy/default O2-C3H8 data cannot satisfy Box B");

    HEMBoxBPrerequisites prerequisites;
    prerequisites.o2C3H8KijReadiness = defaultAudit;
    check(evaluateHEMBoxBPrerequisites(prerequisites)
            == HEMBoxPreconditionStatus::OWNERSHIP_TRANSACTION_UNVERIFIED,
        "Box B first requires the rho-e-Y ownership transaction");

    prerequisites.rhoEYOwnershipTransactionValidated = true;
    check(evaluateHEMBoxBPrerequisites(prerequisites)
            == HEMBoxPreconditionStatus::O2_C3H8_KIJ_UNVALIDATED,
        "Box B rejects an unvalidated O2-C3H8 pair");

    BinaryInteractionMetadata fixtureMetadata;
    fixtureMetadata.evidence = BinaryInteractionEvidenceClass::TEST_ONLY;
    fixtureMetadata.zeroClass = BinaryInteractionZeroClass::NOT_ZERO;
    fixtureMetadata.sourceId = "unit-fixture";
    fixtureMetadata.datasetId = "not-physical";
    fixtureMetadata.applicability = envelope;
    fixtureMetadata.absoluteUncertainty = 0.01;
    BinaryInteractionDB testOnly;
    testOnly.add("O2", "C3H8", 0.055, fixtureMetadata);
    const BinaryInteractionReadinessAudit testOnlyAudit =
        testOnly.auditProductionReadiness({"O2", "C3H8"}, envelope);
    check(!testOnlyAudit.productionReady()
            && testOnlyAudit.count(
                BinaryInteractionReadinessStatus::TEST_ONLY) == 1,
        "test-only kij cannot unblock Box B");

    BinaryInteractionMetadata validatedMetadata;
    validatedMetadata.evidence = BinaryInteractionEvidenceClass::VALIDATED;
    validatedMetadata.zeroClass = BinaryInteractionZeroClass::NOT_ZERO;
    validatedMetadata.sourceId = "unit-fixture-not-physical";
    validatedMetadata.datasetId = "unit-fixture-not-physical";
    validatedMetadata.applicability = envelope;
    validatedMetadata.absoluteUncertainty = 0.01;
    BinaryInteractionDB structurallyReadyKij;
    structurallyReadyKij.add("O2", "C3H8", 0.055, validatedMetadata);
    prerequisites.o2C3H8KijReadiness =
        structurallyReadyKij.auditProductionReadiness(
            {"O2", "C3H8"}, envelope);
    prerequisites.requiredApplicability = envelope;
    prerequisites.requiredMinimumO2MoleFraction = 0.0;
    prerequisites.requiredMaximumO2MoleFraction = 1.0;
    check(evaluateHEMBoxBPrerequisites(prerequisites)
            == HEMBoxPreconditionStatus::PHASE_TOPOLOGY_UNVALIDATED,
        "a ready kij cannot substitute for phase-topology evidence");

    BinaryInteractionDB wrongPairKij;
    wrongPairKij.add("N2", "O2", 0.01, validatedMetadata);
    prerequisites.o2C3H8KijReadiness =
        wrongPairKij.auditProductionReadiness({"N2", "O2"}, envelope);
    check(evaluateHEMBoxBPrerequisites(prerequisites)
            == HEMBoxPreconditionStatus::O2_C3H8_KIJ_UNVALIDATED,
        "a ready interaction audit for another species pair is rejected");

    BinaryInteractionApplicability narrowEnvelope = envelope;
    narrowEnvelope.maximumPressure = 4.0e6;
    prerequisites.o2C3H8KijReadiness =
        structurallyReadyKij.auditProductionReadiness(
            {"O2", "C3H8"}, narrowEnvelope);
    check(evaluateHEMBoxBPrerequisites(prerequisites)
            == HEMBoxPreconditionStatus::O2_C3H8_KIJ_UNVALIDATED,
        "a ready interaction audit must cover the complete Box B envelope");

    prerequisites.o2C3H8KijReadiness =
        structurallyReadyKij.auditProductionReadiness(
            {"O2", "C3H8"}, envelope);

    prerequisites.phaseTopology.evidence =
        HEMPhaseTopologyEvidence::TEST_ONLY;
    prerequisites.phaseTopology.firstSpecies = "C3H8";
    prerequisites.phaseTopology.secondSpecies = "O2";
    prerequisites.phaseTopology.applicability = envelope;
    prerequisites.phaseTopology.minimumO2MoleFraction = 0.0;
    prerequisites.phaseTopology.maximumO2MoleFraction = 1.0;
    prerequisites.phaseTopology.maximumEquilibriumPhaseCount = 2;
    prerequisites.phaseTopology.sourceId = "unit-fixture-not-physical";
    prerequisites.phaseTopology.datasetId = "unit-fixture-not-physical";
    check(evaluateHEMBoxBPrerequisites(prerequisites)
            == HEMBoxPreconditionStatus::PHASE_TOPOLOGY_UNVALIDATED,
        "test-only topology evidence cannot unblock Box B");

    prerequisites.phaseTopology.evidence =
        HEMPhaseTopologyEvidence::VALIDATED_MAXIMUM_TWO_PHASES;
    prerequisites.phaseTopology.maximumEquilibriumPhaseCount = 3;
    prerequisites.phaseTopology.vlleObservedOrPredicted = true;
    check(evaluateHEMBoxBPrerequisites(prerequisites)
            == HEMBoxPreconditionStatus::PHASE_TOPOLOGY_UNVALIDATED,
        "observed or predicted VLLE blocks the two-phase-only Box B");

    prerequisites.phaseTopology.maximumEquilibriumPhaseCount = 2;
    prerequisites.phaseTopology.vlleObservedOrPredicted = false;
    prerequisites.phaseTopology.firstSpecies = "N2";
    check(evaluateHEMBoxBPrerequisites(prerequisites)
            == HEMBoxPreconditionStatus::PHASE_TOPOLOGY_UNVALIDATED,
        "phase-topology evidence for another species pair is rejected");

    prerequisites.phaseTopology.firstSpecies = "C3H8";
    prerequisites.phaseTopology.maximumO2MoleFraction = 0.9;
    check(evaluateHEMBoxBPrerequisites(prerequisites)
            == HEMBoxPreconditionStatus::PHASE_TOPOLOGY_UNVALIDATED,
        "phase-topology evidence must cover the required composition range");

    prerequisites.phaseTopology.maximumO2MoleFraction = 1.0;
    check(evaluateHEMBoxBPrerequisites(prerequisites)
            == HEMBoxPreconditionStatus::ELIGIBLE_FOR_NONPRODUCTION_BOX,
        "complete structural prerequisites permit only a nonproduction Box B");

    HEMBoxCPrerequisites chemistryPrerequisites;
    chemistryPrerequisites.boxBPrerequisites =
        HEMBoxPreconditionStatus::PHASE_TOPOLOGY_UNVALIDATED;
    chemistryPrerequisites.boxBConservativeMixingValidated = true;
    chemistryPrerequisites.r2PhysicsValidated = true;
    check(evaluateHEMBoxCPrerequisites(chemistryPrerequisites)
            == HEMBoxPreconditionStatus::BOX_B_UNVALIDATED,
        "Box C cannot bypass Box B through an R2 flag");

    chemistryPrerequisites.boxBPrerequisites =
        HEMBoxPreconditionStatus::ELIGIBLE_FOR_NONPRODUCTION_BOX;
    chemistryPrerequisites.boxBConservativeMixingValidated = false;
    check(evaluateHEMBoxCPrerequisites(chemistryPrerequisites)
            == HEMBoxPreconditionStatus::BOX_B_UNVALIDATED,
        "Box C requires validated conservative Box B closure");

    chemistryPrerequisites.boxBConservativeMixingValidated = true;
    chemistryPrerequisites.r2PhysicsValidated = false;
    check(evaluateHEMBoxCPrerequisites(chemistryPrerequisites)
            == HEMBoxPreconditionStatus::R2_PHYSICS_UNVERIFIED,
        "Box C requires physical R2 validation after Box B");

    chemistryPrerequisites.r2PhysicsValidated = true;
    check(evaluateHEMBoxCPrerequisites(chemistryPrerequisites)
            == HEMBoxPreconditionStatus::ELIGIBLE_FOR_NONPRODUCTION_BOX,
        "all prerequisites permit only a nonproduction Box C experiment");

    // Numerical wiring fixture only. The 0.055 interaction and this binary
    // topology are deliberately not treated as validated physical evidence.
    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const std::vector<std::string> species = {"O2", "C3H8"};
    BinaryInteractionDB numericalKij;
    numericalKij.add("O2", "C3H8", 0.055);
    const double pressure = 4.0e6;
    PRTPFlash tp(critical, numericalKij, species);
    PRPhaseCaloric caloric(critical, numericalKij, species);
    const std::optional<double> loxH = equilibriumEnthalpy(
        tp, caloric, pressure, 130.0, {1.0, 0.0});
    const std::optional<double> fuelH = equilibriumEnthalpy(
        tp, caloric, pressure, 450.0, {0.0, 1.0});
    check(loxH.has_value() && fuelH.has_value(),
        "Box B fixture endpoint enthalpies are available");

    if (loxH && fuelH)
    {
        const double fuelMassFraction = 1.05/3.85;
        const std::vector<double> overallY = {
            1.0 - fuelMassFraction, fuelMassFraction};
        const double targetH = (1.0 - fuelMassFraction)*(*loxH)
            + fuelMassFraction*(*fuelH);

        PHFlashOptions phOptions;
        phOptions.minimumTemperature = 80.0;
        phOptions.maximumTemperature = 700.0;
        phOptions.scanIntervals = 512;
        phOptions.maximumIterations = 240;
        PRPHFlash ph(critical, numericalKij, species, phOptions);
        const PHFlashResult mixed = ph.flashMassFractions(
            pressure, targetH, overallY);
        check(mixed.success() && mixed.conservedMixture.has_value()
                && mixed.tpFlash.region == EquilibriumRegion::TWO_PHASE,
            "test-only Box B reaches a conservative two-phase PH state");

        if (mixed.success() && mixed.conservedMixture)
        {
            const ConservedPhaseMixture& mixture = *mixed.conservedMixture;
            check(close(mixture.h, targetH, 1.0e-9, 1.0e-4),
                "Box B conserves the mixed total enthalpy");
            check(mixture.overallMassFractions.size() == overallY.size()
                    && close(
                        mixture.overallMassFractions[0], overallY[0],
                        1.0e-9, 1.0e-10)
                    && close(
                        mixture.overallMassFractions[1], overallY[1],
                        1.0e-9, 1.0e-10),
                "Box B conserves the overall species mass fractions");

            const thermo::PREquilibriumThermoModel model(
                critical, numericalKij, species, mixingRhoEOptions());
            thermo::EquilibriumSeed seed;
            seed.pressure = pressure;
            seed.temperature = mixed.temperature;
            const thermo::EquilibriumResult recovered = model.recoverRhoE(
                mixture.rho, mixture.e, overallY, &seed);
            check(recovered.success() && recovered.state.has_value()
                    && recovered.state->region
                        == thermo::EquilibriumRegion::TWO_PHASE,
                "the authoritative rho-e-Y model recovers the Box B state");
            if (recovered.success() && recovered.state)
            {
                check(close(recovered.state->rho, mixture.rho)
                        && close(recovered.state->e, mixture.e, 1.0e-8, 1.0e-4)
                        && close(recovered.state->p, pressure, 1.0e-7, 1.0),
                    "Box B PH and rho-e-Y states close on conserved fields");
                check(recovered.state->overallMassFractions == overallY,
                    "quality remains derived while rho-e-Y owns composition");
            }
        }
    }

    if (failures != 0)
    {
        std::cerr << "PR_TWO_PHASE_MIXING_CONTRACT_FAIL checks=" << checks
                  << " failures=" << failures << "\n";
        return 1;
    }
    std::cout << "PR_TWO_PHASE_MIXING_CONTRACT_PASS checks=" << checks
              << " scope=nonproduction-preconditions\n";
    return 0;
}
