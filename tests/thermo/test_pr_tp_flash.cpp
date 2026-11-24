#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "BinaryInteractionDB.H"
#include "PhaseEquilibrium/PRMixtureKernel.H"
#include "PhaseEquilibrium/PRTPFlash.H"
#include "PengRobinsonAlpha.H"
#include "SpeciesCriticalDB.H"
#include "SpeciesRegistry.H"

using namespace Foam;
using namespace Foam::phaseEquilibrium;

namespace
{

int failures = 0;

void check(
    const std::string& name,
    bool condition,
    double got = 0.0,
    double expected = 0.0)
{
    std::cout << "  [" << (condition ? "PASS" : "FAIL") << "] " << name;
    if (!condition)
    {
        std::cout << " got=" << std::setprecision(14) << got
                  << " expected=" << expected;
        ++failures;
    }
    std::cout << '\n';
}

double sum(const std::vector<double>& values)
{
    return std::accumulate(values.begin(), values.end(), 0.0);
}

double relativeError(double lhs, double rhs)
{
    return std::abs(lhs - rhs)
        /std::max({1.0e-30, std::abs(lhs), std::abs(rhs)});
}

BinaryInteractionMetadata interactionMetadata(
    BinaryInteractionEvidenceClass evidence,
    BinaryInteractionZeroClass zeroClass,
    const std::string& datasetId)
{
    BinaryInteractionMetadata metadata;
    metadata.evidence = evidence;
    metadata.zeroClass = zeroClass;
    metadata.sourceId = "rfam-test-source";
    metadata.datasetId = datasetId;
    metadata.applicability.minimumTemperature = 80.0;
    metadata.applicability.maximumTemperature = 3000.0;
    metadata.applicability.minimumPressure = 1.0e5;
    metadata.applicability.maximumPressure = 3.0e7;
    metadata.absoluteUncertainty = 1.0e-4;
    return metadata;
}

void checkTwoPhaseClosure(
    const TPFlashResult& result,
    const std::vector<double>& feed)
{
    check("two-phase flash status", result.valid());
    if (!result.valid())
    {
        std::cout << "    status=" << statusName(result.status) << '\n';
        return;
    }

    check("two-phase region",
        result.region == EquilibriumRegion::TWO_PHASE);
    check("phase fraction is determined", result.phaseFractionDetermined);
    check("phase-split payload predicate passes", result.phaseSplitReady());
    check("vapor fraction lower bound", result.vaporMoleFraction > 0.0,
        result.vaporMoleFraction, 0.0);
    check("vapor fraction upper bound", result.vaporMoleFraction < 1.0,
        result.vaporMoleFraction, 1.0);
    check("liquid composition normalized",
        std::abs(sum(result.liquidMoleFractions) - 1.0) < 1.0e-10,
        sum(result.liquidMoleFractions), 1.0);
    check("vapor composition normalized",
        std::abs(sum(result.vaporMoleFractions) - 1.0) < 1.0e-10,
        sum(result.vaporMoleFractions), 1.0);

    double maximumMaterialResidual = 0.0;
    double maximumFugacityResidual = 0.0;
    for (std::size_t i = 0; i < feed.size(); ++i)
    {
        const double reconstructed =
            (1.0 - result.vaporMoleFraction)
                *result.liquidMoleFractions[i]
            + result.vaporMoleFraction*result.vaporMoleFractions[i];
        maximumMaterialResidual = std::max(
            maximumMaterialResidual,
            std::abs(feed[i] - reconstructed));

        if (feed[i] > 0.0)
        {
            const double residual =
                std::log(result.liquidMoleFractions[i])
                + result.liquidState.lnFugacityCoefficients[i]
                - std::log(result.vaporMoleFractions[i])
                - result.vaporState.lnFugacityCoefficients[i];
            maximumFugacityResidual = std::max(
                maximumFugacityResidual, std::abs(residual));
        }
    }
    check("component material closure",
        maximumMaterialResidual < 1.0e-9,
        maximumMaterialResidual, 0.0);
    check("reported material residual agrees",
        std::abs(maximumMaterialResidual
            - result.maximumMoleMaterialResidual) < 1.0e-12,
        result.maximumMoleMaterialResidual, maximumMaterialResidual);
    check("phase fugacity equality",
        maximumFugacityResidual < 2.0e-8,
        maximumFugacityResidual, 0.0);
    check("reported fugacity residual agrees",
        std::abs(maximumFugacityResidual
            - result.maximumLogFugacityResidual) < 1.0e-12,
        result.maximumLogFugacityResidual, maximumFugacityResidual);
    check("Rachford-Rice residual",
        std::abs(result.rachfordRiceResidual) < 1.0e-9,
        result.rachfordRiceResidual, 0.0);
    check("feed TPD detects instability", result.feedMinimumTpd < -1.0e-8,
        result.feedMinimumTpd, -1.0e-8);
    check("liquid phase passes TPD stability",
        result.liquidMinimumTpd >= -5.0e-8,
        result.liquidMinimumTpd, -5.0e-8);
    check("vapor phase passes TPD stability",
        result.vaporMinimumTpd >= -5.0e-8,
        result.vaporMinimumTpd, -5.0e-8);
}

} // End anonymous namespace

int main()
{
    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB kij = BinaryInteractionDB::createDefault();

    std::cout << "=== Auditable Peng-Robinson TP flash ===\n";

    const SpeciesRegistry registry;
    const BinaryInteractionReadinessAudit globalKijAudit =
        kij.auditProductionReadiness(registry.names());
    const std::vector<std::pair<std::string, std::string>>
        expectedGlobalMissing = {
            {"C2H6", "CO"},
            {"C2H6", "CO2"},
            {"C2H6", "H2"},
            {"C2H6", "H2O"},
            {"C3H8", "CO"},
            {"C3H8", "H2"},
            {"C3H8", "H2O"},
            {"CH4", "CO"},
            {"CO", "H2O"},
            {"CO", "N2"},
            {"CO", "O2"},
            {"CO2", "H2"},
            {"H2", "H2O"},
            {"H2", "N2"},
            {"H2O", "N2"},
            {"H2O", "O2"},
        };
    check("nine-species kij audit is structurally valid",
        globalKijAudit.requestValid && globalKijAudit.pairs.size() == 36);
    check("nine-species kij audit reports 16 missing pairs",
        globalKijAudit.count(BinaryInteractionReadinessStatus::MISSING) == 16);
    check("nine-species kij audit reports 20 legacy pairs",
        globalKijAudit.count(
            BinaryInteractionReadinessStatus::LEGACY_UNPROVENANCED) == 20);
    check("legacy default kij data cannot authorize production",
        !globalKijAudit.productionReady()
            && globalKijAudit.count(
                BinaryInteractionReadinessStatus::READY) == 0);
    const BinaryInteractionEntry* legacyZero = kij.find("N2", "O2");
    check("hard-coded zero remains explicitly unvalidated",
        legacyZero
            && legacyZero->metadata.evidence
                == BinaryInteractionEvidenceClass::LEGACY_UNPROVENANCED
            && legacyZero->metadata.zeroClass
                == BinaryInteractionZeroClass::UNVALIDATED_ZERO);
    check("nine-species missing-pair report is exact and sorted",
        globalKijAudit.pairsWithStatus(
            BinaryInteractionReadinessStatus::MISSING)
            == expectedGlobalMissing);
    std::vector<std::string> reversedSpecies = registry.names();
    std::reverse(reversedSpecies.begin(), reversedSpecies.end());
    check("kij readiness report is independent of input ordering",
        globalKijAudit.report()
            == kij.auditProductionReadiness(reversedSpecies).report());
    check("kij readiness report is LF-only",
        globalKijAudit.report().find('\r') == std::string::npos);

    const BinaryInteractionReadinessAudit reactionKijAudit =
        kij.auditProductionReadiness(
            {"C3H8", "O2", "CO", "CO2", "H2O"});
    const std::vector<std::pair<std::string, std::string>>
        expectedReactionMissing = {
            {"C3H8", "CO"},
            {"C3H8", "H2O"},
            {"CO", "H2O"},
            {"CO", "O2"},
            {"H2O", "O2"},
        };
    check("reacting species set reports its five missing kij pairs",
        reactionKijAudit.pairs.size() == 10
            && reactionKijAudit.pairsWithStatus(
                BinaryInteractionReadinessStatus::MISSING)
                == expectedReactionMissing);

    BinaryInteractionDB readinessDB;
    BinaryInteractionApplicability productionEnvelope;
    productionEnvelope.minimumTemperature = 100.0;
    productionEnvelope.maximumTemperature = 2500.0;
    productionEnvelope.minimumPressure = 2.0e5;
    productionEnvelope.maximumPressure = 2.0e7;
    check("empty kij readiness request fails closed",
        !readinessDB.auditProductionReadiness({}).requestValid);
    readinessDB.add("A", "B", 0.1);
    check("legacy add overload remains numerical but unprovenanced",
        readinessDB.get("B", "A") == 0.1
            && readinessDB.auditProductionReadiness({"A", "B"}).count(
                BinaryInteractionReadinessStatus::LEGACY_UNPROVENANCED) == 1);
    readinessDB.add(
        "A", "B", 0.1,
        interactionMetadata(
            BinaryInteractionEvidenceClass::TEST_ONLY,
            BinaryInteractionZeroClass::NOT_ZERO,
            "test-only"));
    check("test-only kij cannot authorize production",
        readinessDB.auditProductionReadiness({"A", "B"}).count(
            BinaryInteractionReadinessStatus::TEST_ONLY) == 1);
    readinessDB.add(
        "A", "B", 0.1,
        interactionMetadata(
            BinaryInteractionEvidenceClass::PROVISIONAL,
            BinaryInteractionZeroClass::NOT_ZERO,
            "provisional"));
    check("provisional kij cannot authorize production",
        readinessDB.auditProductionReadiness({"A", "B"}).count(
            BinaryInteractionReadinessStatus::PROVISIONAL) == 1);
    readinessDB.add(
        "A", "B", 0.1,
        interactionMetadata(
            BinaryInteractionEvidenceClass::VALIDATED,
            BinaryInteractionZeroClass::NOT_ZERO,
            "validated-nonzero"));
    check("complete validated nonzero kij passes readiness",
        !readinessDB.auditProductionReadiness({"A", "B"})
            .productionReady()
            && readinessDB.auditProductionReadiness(
                {"A", "B"}, productionEnvelope)
            .productionReady());
    BinaryInteractionMetadata blankSource = interactionMetadata(
        BinaryInteractionEvidenceClass::VALIDATED,
        BinaryInteractionZeroClass::NOT_ZERO,
        "validated-nonzero");
    blankSource.sourceId = "   ";
    readinessDB.add("A", "B", 0.1, blankSource);
    check("blank kij source identifier fails closed",
        readinessDB.auditProductionReadiness(
            {"A", "B"}, productionEnvelope).count(
                BinaryInteractionReadinessStatus::INVALID_METADATA) == 1);
    BinaryInteractionMetadata controlDataset = interactionMetadata(
        BinaryInteractionEvidenceClass::VALIDATED,
        BinaryInteractionZeroClass::NOT_ZERO,
        "\t  ");
    readinessDB.add("A", "B", 0.1, controlDataset);
    check("control-character kij dataset identifier fails closed",
        readinessDB.auditProductionReadiness(
            {"A", "B"}, productionEnvelope).count(
                BinaryInteractionReadinessStatus::INVALID_METADATA) == 1);
    readinessDB.add(
        "A", "B", 0.0,
        interactionMetadata(
            BinaryInteractionEvidenceClass::VALIDATED,
            BinaryInteractionZeroClass::UNVALIDATED_ZERO,
            "unvalidated-zero"));
    check("zero kij requires explicit validated-zero evidence",
        readinessDB.auditProductionReadiness({"A", "B"}).count(
            BinaryInteractionReadinessStatus::UNVALIDATED_ZERO) == 1);
    readinessDB.add(
        "A", "B", 0.0,
        interactionMetadata(
            BinaryInteractionEvidenceClass::VALIDATED,
            BinaryInteractionZeroClass::VALIDATED_EXPLICIT_ZERO,
            "validated-zero"));
    check("validated explicit zero can pass readiness",
        readinessDB.auditProductionReadiness(
            {"A", "B"}, productionEnvelope)
            .productionReady());
    BinaryInteractionApplicability outsideEnvelope;
    outsideEnvelope.minimumTemperature = 70.0;
    outsideEnvelope.maximumTemperature = 3000.0;
    outsideEnvelope.minimumPressure = 1.0e5;
    outsideEnvelope.maximumPressure = 3.0e7;
    check("validated kij outside its applicability fails closed",
        readinessDB.auditProductionReadiness(
            {"A", "B"}, outsideEnvelope).count(
                BinaryInteractionReadinessStatus::OUTSIDE_APPLICABILITY)
            == 1);
    check("missing get fallback remains diagnostic-only zero",
        readinessDB.get("A", "C") == 0.0
            && readinessDB.find("A", "C") == nullptr);

    PRMixtureKernel oxygenKernel(critical, kij, {"O2"});
    check("pure O2 kernel configured", oxygenKernel.configured());
    const PRMixtureParameters oxygenParameters =
        oxygenKernel.mixtureParameters(4.0e6, 130.0, {1.0});
    check("pure O2 mixture parameters valid", oxygenParameters.valid());
    if (oxygenParameters.valid())
    {
        const CriticalProperties& oxygen = critical.get("O2");
        const double expectedA = evaluatePRAlpha(
            oxygen.Tc, oxygen.Pc, oxygen.omega, 130.0).a;
        check("phase kernel shares exact PR alpha",
            relativeError(oxygenParameters.speciesA[0], expectedA) < 1.0e-14,
            oxygenParameters.speciesA[0], expectedA);
    }

    PRMixtureKernel sourcedZeroKernel(critical, kij, {"N2", "O2"});
    const PRMixtureParameters sourcedZero = sourcedZeroKernel.mixtureParameters(
        2.0e6, 120.0, {0.5, 0.5});
    check("explicitly sourced zero kij is accepted", sourcedZero.valid());

    PRMixtureKernel missingPairKernel(critical, kij, {"O2", "H2O"});
    const PRMixtureParameters missingPair = missingPairKernel.mixtureParameters(
        2.0e6, 250.0, {0.5, 0.5});
    check("active missing kij fails closed",
        missingPair.status
            == PhaseEquilibriumStatus::MISSING_BINARY_INTERACTION);
    const PRMixtureParameters traceMissingPair =
        missingPairKernel.mixtureParameters(
            2.0e6, 250.0, {1.0 - 1.0e-14, 1.0e-14});
    check("trace nonzero missing kij also fails closed",
        traceMissingPair.status
            == PhaseEquilibriumStatus::MISSING_BINARY_INTERACTION);
    const PRMixtureParameters inactiveMissingPair =
        missingPairKernel.mixtureParameters(2.0e6, 250.0, {1.0, 0.0});
    check("inactive missing kij does not poison pure subspace",
        inactiveMissingPair.valid());

    BinaryInteractionDB mutableKij = BinaryInteractionDB::createDefault();
    PRMixtureKernel frozenKijKernel(
        critical, mutableKij, {"O2", "C3H8"});
    const PRMixtureParameters frozenBefore = frozenKijKernel.mixtureParameters(
        2.0e6, 200.0, {0.5, 0.5});
    mutableKij.add("O2", "C3H8", 0.9);
    const PRMixtureParameters frozenAfter = frozenKijKernel.mixtureParameters(
        2.0e6, 200.0, {0.5, 0.5});
    check("configured kij matrix is immutable",
        frozenBefore.valid() && frozenAfter.valid()
            && frozenBefore.a == frozenAfter.a,
        frozenAfter.a, frozenBefore.a);

    PRTPFlash oxygenFlash(critical, kij, {"O2"});
    const TPFlashResult subcooled =
        oxygenFlash.flashMoleFractions(4.0e6, 130.0, {1.0});
    check("subcooled O2 flash succeeds", subcooled.valid());
    check("subcooled O2 is single liquid",
        subcooled.region == EquilibriumRegion::SINGLE_LIQUID);
    check("single-liquid beta is zero",
        subcooled.phaseFractionDetermined
            && subcooled.vaporMoleFraction == 0.0,
        subcooled.vaporMoleFraction, 0.0);
    check("single-liquid state payload is complete",
        subcooled.singlePhaseState.valid()
            && subcooled.singlePhaseMoleFractions == std::vector<double>{1.0}
            && subcooled.phaseSplitReady());

    const TPFlashResult gaseous =
        oxygenFlash.flashMoleFractions(1.0e5, 300.0, {1.0});
    check("gaseous O2 flash succeeds", gaseous.valid());
    check("gaseous O2 is single vapor",
        gaseous.region == EquilibriumRegion::SINGLE_VAPOR);
    check("single-vapor beta is one",
        gaseous.phaseFractionDetermined
            && gaseous.vaporMoleFraction == 1.0,
        gaseous.vaporMoleFraction, 1.0);
    check("single-vapor state payload is complete",
        gaseous.singlePhaseState.valid()
            && gaseous.singlePhaseMoleFractions == std::vector<double>{1.0}
            && gaseous.phaseSplitReady());

    const TPFlashResult supercritical =
        oxygenFlash.flashMoleFractions(6.0e6, 200.0, {1.0});
    check("supercritical O2 flash succeeds", supercritical.valid());
    check("supercritical region remains distinct",
        supercritical.region == EquilibriumRegion::SINGLE_SUPERCRITICAL);
    check("supercritical state payload is complete",
        supercritical.singlePhaseState.valid()
            && supercritical.singlePhaseMoleFractions
                == std::vector<double>{1.0}
            && supercritical.phaseSplitReady());
    check("supercritical state does not invent vapor fraction",
        !supercritical.phaseFractionDetermined
            && std::isnan(supercritical.vaporMoleFraction));

    const PureSaturationResult saturation =
        oxygenFlash.solvePureSaturationTemperature(4.0e6, 120.0, 154.5);
    check("O2 saturation solve succeeds", saturation.valid());
    if (saturation.valid())
    {
        std::cout << "  O2_SATURATION p=" << saturation.pressure
                  << " T=" << std::setprecision(12)
                  << saturation.temperature
                  << " Zl=" << saturation.liquidCompressibility
                  << " Zv=" << saturation.vaporCompressibility
                  << " ln_f_ratio=" << saturation.logFugacityResidual
                  << '\n';
        check("O2 PR saturation temperature is plausible",
            saturation.temperature > 145.0 && saturation.temperature < 152.0,
            saturation.temperature, 148.5);
        check("O2 4 MPa frozen PR saturation value",
            std::abs(saturation.temperature - 148.52289611) < 2.0e-6,
            saturation.temperature, 148.52289611);
        check("O2 4 MPa CoolProp 7.2 saturation smoke",
            std::abs(saturation.temperature - 148.65925898236446) < 0.25,
            saturation.temperature, 148.65925898236446);
        check("pure saturation fugacity equality",
            std::abs(saturation.logFugacityResidual) < 2.0e-9,
            saturation.logFugacityResidual, 0.0);
        check("saturation has distinct liquid/vapor roots",
            saturation.liquidCompressibility
                < saturation.vaporCompressibility,
            saturation.liquidCompressibility,
            saturation.vaporCompressibility);

        const TPFlashResult coexistence =
            oxygenFlash.flashMoleFractions(
                4.0e6, saturation.temperature, {1.0});
        check("pure coexistence TP solve succeeds", coexistence.valid());
        check("pure coexistence has explicit region",
            coexistence.region
                == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE);
        check("pure coexistence beta is underdetermined",
            !coexistence.phaseFractionDetermined
                && std::isnan(coexistence.vaporMoleFraction));
        check("pure coexistence is not phase-split ready",
            !coexistence.phaseSplitReady());
        check("pure coexistence reports fugacity equality",
            coexistence.maximumLogFugacityResidual < 2.0e-8,
            coexistence.maximumLogFugacityResidual, 0.0);
    }

    const PureSaturationResult saturation2MPa =
        oxygenFlash.solvePureSaturationTemperature(2.0e6, 90.0, 154.5);
    check("O2 2 MPa saturation solve succeeds", saturation2MPa.valid());
    if (saturation2MPa.valid())
    {
        std::cout << "  O2_SATURATION_2MPA T=" << std::setprecision(12)
                  << saturation2MPa.temperature << " ln_f_ratio="
                  << saturation2MPa.logFugacityResidual << '\n';
        check("O2 2 MPa saturation temperature is plausible",
            saturation2MPa.temperature > 115.0
                && saturation2MPa.temperature < 135.0,
            saturation2MPa.temperature, 125.0);
        check("O2 2 MPa frozen PR saturation value",
            std::abs(saturation2MPa.temperature - 132.617012773) < 2.0e-6,
            saturation2MPa.temperature, 132.617012773);
        check("O2 2 MPa CoolProp 7.2 saturation smoke",
            std::abs(saturation2MPa.temperature - 132.74292718240093) < 0.25,
            saturation2MPa.temperature, 132.74292718240093);
    }
    const PureSaturationResult saturation6MPa =
        oxygenFlash.solvePureSaturationTemperature(6.0e6, 120.0, 180.0);
    check("O2 6 MPa has no subcritical saturation bracket",
        saturation6MPa.status
            == PhaseEquilibriumStatus::PURE_SATURATION_NO_BRACKET);

    PRTPFlash mixtureFlash(critical, kij, {"O2", "C3H8"});
    const std::vector<double> feed = {0.5, 0.5};
    const StabilityResult splitStability = mixtureFlash.analyzeStability(
        2.0e6, 200.0, feed);
    check("deterministic multi-start stability succeeds",
        splitStability.valid());
    check("six TPD seeds were audited",
        splitStability.attemptedTrials == 6,
        static_cast<double>(splitStability.attemptedTrials), 6.0);
    check("all TPD seeds converged",
        splitStability.allTrialsConverged(),
        static_cast<double>(splitStability.failedTrials), 0.0);
    const TPFlashResult split = mixtureFlash.flashMoleFractions(
        2.0e6, 200.0, feed);
    if (split.valid())
    {
        std::cout << "  O2_C3H8_FLASH beta_n=" << std::setprecision(12)
                  << split.vaporMoleFraction
                  << " x=[" << split.liquidMoleFractions[0] << ','
                  << split.liquidMoleFractions[1] << "] y=["
                  << split.vaporMoleFractions[0] << ','
                  << split.vaporMoleFractions[1] << "] feed_tpd="
                  << split.feedMinimumTpd << " max_ln_f_res="
                  << split.maximumLogFugacityResidual << " max_z_res="
                  << split.maximumMoleMaterialResidual << '\n';
    }
    if (!split.valid())
    {
        std::cout << "  mixture failure status=" << statusName(split.status)
                  << " iterations=" << split.flashIterations
                  << " fugacityResidual="
                  << split.maximumLogFugacityResidual
                  << " materialResidual="
                  << split.maximumMoleMaterialResidual << '\n';
    }
    checkTwoPhaseClosure(split, feed);

    if (split.valid())
    {
        // CoolProp 7.2 HEOS has no fitted binary departure term for this pair,
        // so these are broad implementation cross-checks, not experimental
        // accuracy claims. The frozen source is tests/reference/.
        check("O2/C3H8 beta agrees with broad CoolProp implementation smoke",
            std::abs(split.vaporMoleFraction - 0.5033307144792941) < 0.10,
            split.vaporMoleFraction, 0.5033307144792941);
        check("O2/C3H8 vapor composition agrees with CoolProp smoke",
            std::abs(split.vaporMoleFractions[0]
                - 0.9852546388372834) < 0.02,
            split.vaporMoleFractions[0], 0.9852546388372834);
        check("O2/C3H8 liquid composition remains within broad smoke bound",
            std::abs(split.liquidMoleFractions[0]
                - 0.008237027754198914) < 0.15,
            split.liquidMoleFractions[0], 0.008237027754198914);

        const TPFlashResult massBasis = mixtureFlash.flashMassFractions(
            2.0e6, 200.0, split.feedMassFractions);
        check("mass-basis flash succeeds", massBasis.valid());
        check("mass-basis and mole-basis beta agree",
            massBasis.valid()
                && std::abs(massBasis.vaporMoleFraction
                    - split.vaporMoleFraction) < 1.0e-9,
            massBasis.vaporMoleFraction, split.vaporMoleFraction);
    }

    BinaryInteractionDB zeroKij;
    zeroKij.add("O2", "C3H8", 0.0);
    PRTPFlash wilsonFirstFlash(critical, zeroKij, {"O2", "C3H8"});
    const std::vector<double> gapMassFractions = {0.9, 0.1};
    for (double temperature : {163.2, 163.5})
    {
        const TPFlashResult gapState =
            wilsonFirstFlash.flashMassFractions(
                6.0e6, temperature, gapMassFractions);
        std::cout << "  WILSON_FIRST_GAP T=" << temperature
                  << " status=" << statusName(gapState.status)
                  << " iterations=" << gapState.flashIterations
                  << " fugacityResidual="
                  << gapState.maximumLogFugacityResidual << '\n';
        check("Wilson-first binary gap state succeeds", gapState.valid());
        check("Wilson-first binary gap state is two phase",
            gapState.valid()
                && gapState.region == EquilibriumRegion::TWO_PHASE);
        if (gapState.valid())
            checkTwoPhaseClosure(gapState, gapState.feedMoleFractions);
    }

    BinaryInteractionDB lowFuelKij;
    lowFuelKij.add("O2", "C3H8", 0.1);
    PRTPFlash lowFuelFlash(critical, lowFuelKij, {"O2", "C3H8"});
    std::vector<double> lowFuelMoleFractions;
    const PhaseEquilibriumStatus lowFuelConversion =
        lowFuelFlash.kernel().massToMoleFractions(
            {0.99, 0.01}, lowFuelMoleFractions);
    check("low-fuel stability feed conversion succeeds",
        lowFuelConversion == PhaseEquilibriumStatus::SUCCESS);
    if (lowFuelConversion == PhaseEquilibriumStatus::SUCCESS)
    {
        for (double temperature : {100.0, 130.0})
        {
            const StabilityResult stability = lowFuelFlash.analyzeStability(
                4.0e6, temperature, lowFuelMoleFractions);
            check("neutral retry closes all low-fuel stability trials",
                stability.valid() && stability.allTrialsConverged(),
                static_cast<double>(stability.failedTrials), 0.0);
            const TPFlashResult state = lowFuelFlash.flashMoleFractions(
                4.0e6, temperature, lowFuelMoleFractions);
            check("low-fuel stable state commits after neutral retry",
                state.valid()
                    && state.region == EquilibriumRegion::SINGLE_LIQUID);
        }

        const StabilityResult unsupportedStability =
            lowFuelFlash.analyzeStability(
                4.0e6, 150.0, lowFuelMoleFractions);
        check("neutral retry separates stability convergence from topology",
            unsupportedStability.valid()
                && unsupportedStability.allTrialsConverged()
                && unsupportedStability.unstable);
        const TPFlashResult unsupportedSplit =
            lowFuelFlash.flashMoleFractions(
                4.0e6, 150.0, lowFuelMoleFractions);
        check("post-split instability remains fail-closed",
            unsupportedSplit.status
                == PhaseEquilibriumStatus::FLASH_NOT_CONVERGED
                && !unsupportedSplit.valid()
                && (unsupportedSplit.liquidMinimumTpd < 0.0
                    || unsupportedSplit.vaporMinimumTpd < 0.0));
    }

    struct BinarySeedFamilyRegressionState
    {
        double pressure;
        double temperature;
        double kij;
        double leftOxygen;
        double rightOxygen;
        const char* label;
    };
    const BinarySeedFamilyRegressionState seedFamilyStates[] = {
        {
            4.0e6, 130.0, 0.10,
            0.42519730138846235, 0.9640176587984397,
            "4 MPa 130 K kij=0.10"
        },
        {
            4.0e6, 150.0, 0.10,
            0.4574313610198662, 0.99928995607587301,
            "4 MPa 150 K kij=0.10"
        },
        {
            4.0e6, 150.0, 0.15,
            0.26848275851064196, 0.99931268630950054,
            "4 MPa 150 K kij=0.15"
        },
        {
            5.0e6, 130.0, 0.10,
            0.42908321561730761, 0.96273807986165549,
            "5 MPa 130 K kij=0.10"
        },
        {
            5.0e6, 150.0, 0.10,
            0.47481840455877045, 0.96633841610613491,
            "5 MPa 150 K kij=0.10"
        }
    };
    for (const BinarySeedFamilyRegressionState& state : seedFamilyStates)
    {
        BinaryInteractionDB seedFamilyKij;
        seedFamilyKij.add("O2", "C3H8", state.kij);
        const PRTPFlash seedFamilyFlash(
            critical, seedFamilyKij, {"O2", "C3H8"});
        const double feedOxygen = state.leftOxygen
            + 0.90*(state.rightOxygen - state.leftOxygen);
        const std::vector<double> seedFamilyFeed = {
            feedOxygen, 1.0 - feedOxygen};
        const TPFlashResult seedFamilySplit =
            seedFamilyFlash.flashMoleFractions(
                state.pressure, state.temperature, seedFamilyFeed);
        const std::string prefix = std::string(state.label) + " ";
        check(prefix + "oxygen-side tie flash succeeds",
            seedFamilySplit.valid());
        check(prefix + "oxygen-side tie is two phase",
            seedFamilySplit.valid()
                && seedFamilySplit.region == EquilibriumRegion::TWO_PHASE);
        if (seedFamilySplit.valid()
            && seedFamilySplit.region == EquilibriumRegion::TWO_PHASE)
        {
            const double splitLeft = std::min(
                seedFamilySplit.liquidMoleFractions[0],
                seedFamilySplit.vaporMoleFractions[0]);
            const double splitRight = std::max(
                seedFamilySplit.liquidMoleFractions[0],
                seedFamilySplit.vaporMoleFractions[0]);
            const double endpointError = std::max(
                std::abs(splitLeft - state.leftOxygen),
                std::abs(splitRight - state.rightOxygen));
            check(prefix + "oxygen-side tie matches global endpoints",
                endpointError <= 1.0e-7, endpointError, 1.0e-7);
            checkTwoPhaseClosure(seedFamilySplit, seedFamilyFeed);
        }
    }

    BinaryInteractionDB topologyKij;
    topologyKij.add("O2", "C3H8", 0.055);
    PRTPFlash topologyFlash(critical, topologyKij, {"O2", "C3H8"});
    struct TopologyRegressionState
    {
        double pressure;
        double leftOxygen;
        double rightOxygen;
        const char* label;
    };
    const TopologyRegressionState topologyStates[] = {
        {
            4.0e6,
            0.22340198679699921,
            0.98810386533919004,
            "4 MPa"
        },
        {
            6.0e6,
            0.32388461471797048,
            0.98722176445920951,
            "6 MPa"
        }
    };
    for (const TopologyRegressionState& state : topologyStates)
    {
        const double midpoint = 0.5
            *(state.leftOxygen + state.rightOxygen);
        const std::vector<double> topologyFeed = {
            midpoint, 1.0 - midpoint};
        const StabilityResult stability = topologyFlash.analyzeStability(
            state.pressure, 200.0, topologyFeed);
        const std::string prefix = std::string(state.label)
            + " 200 K kij=0.055 ";
        check(prefix + "stability succeeds", stability.valid());
        check(prefix + "audits six stability trials",
            stability.attemptedTrials == 6,
            static_cast<double>(stability.attemptedTrials), 6.0);
        check(prefix + "converges all stability trials",
            stability.allTrialsConverged(),
            static_cast<double>(stability.failedTrials), 0.0);

        const TPFlashResult topologySplit =
            topologyFlash.flashMoleFractions(
                state.pressure, 200.0, topologyFeed);
        check(prefix + "TP flash succeeds", topologySplit.valid());
        check(prefix + "TP flash returns two phases",
            topologySplit.valid()
                && topologySplit.region == EquilibriumRegion::TWO_PHASE);
        if (topologySplit.valid()
            && topologySplit.region == EquilibriumRegion::TWO_PHASE)
        {
            const double splitLeft = std::min(
                topologySplit.liquidMoleFractions[0],
                topologySplit.vaporMoleFractions[0]);
            const double splitRight = std::max(
                topologySplit.liquidMoleFractions[0],
                topologySplit.vaporMoleFractions[0]);
            const double endpointError = std::max(
                std::abs(splitLeft - state.leftOxygen),
                std::abs(splitRight - state.rightOxygen));
            check(prefix + "matches global-tie endpoints",
                endpointError <= 1.0e-7, endpointError, 1.0e-7);
            checkTwoPhaseClosure(topologySplit, topologyFeed);
        }
    }

    TPFlashOptions constrainedStabilityOptions;
    constrainedStabilityOptions.maxStabilityIterations = 1;
    constrainedStabilityOptions.iterationTolerance = 1.0e-30;
    PRTPFlash constrainedStabilityFlash(
        critical,
        topologyKij,
        {"O2", "C3H8"},
        constrainedStabilityOptions);
    const double constrainedMidpoint = 0.5
        *(topologyStates[0].leftOxygen + topologyStates[0].rightOxygen);
    const std::vector<double> constrainedFeed = {
        constrainedMidpoint, 1.0 - constrainedMidpoint};
    const StabilityResult constrainedStability =
        constrainedStabilityFlash.analyzeStability(
            topologyStates[0].pressure, 200.0, constrainedFeed);
    check("stability iteration exhaustion remains fail-closed",
        constrainedStability.status
            == PhaseEquilibriumStatus::STABILITY_NOT_CONVERGED
            && constrainedStability.failedTrials > 0);
    const TPFlashResult constrainedSplit =
        constrainedStabilityFlash.flashMoleFractions(
            topologyStates[0].pressure, 200.0, constrainedFeed);
    check("TP flash does not ignore failed stability trials",
        constrainedSplit.status
            == PhaseEquilibriumStatus::STABILITY_NOT_CONVERGED);

    const TPFlashResult invalid = mixtureFlash.flashMoleFractions(
        2.0e6, 200.0, {0.6, -0.1});
    check("negative feed fails with explicit status",
        invalid.status == PhaseEquilibriumStatus::NEGATIVE_COMPOSITION);
    const TPFlashResult invalidPressure = mixtureFlash.flashMoleFractions(
        -2.0e6, 200.0, feed);
    check("negative pressure fails with explicit status",
        invalidPressure.status == PhaseEquilibriumStatus::INVALID_PRESSURE);
    const RachfordRiceResult noBracket =
        mixtureFlash.solveRachfordRice(feed, {2.0, 3.0});
    check("RR no-bracket is explicit and not a phase verdict",
        noBracket.status
            == PhaseEquilibriumStatus::RACHFORD_RICE_NO_BRACKET);

    std::cout << (failures == 0
        ? "PR_TP_FLASH_PASS\n"
        : "PR_TP_FLASH_FAIL\n");
    return failures == 0 ? 0 : 1;
}
