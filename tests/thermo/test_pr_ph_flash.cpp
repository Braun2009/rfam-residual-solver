#include "PhaseEquilibrium/PRPHFlash.H"
#include "SpeciesRegistry.H"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

using namespace Foam;
using namespace Foam::phaseEquilibrium;

namespace
{

int failures = 0;

void check(bool condition, const std::string& label)
{
    std::cout << (condition ? "  PASS " : "  FAIL ") << label << '\n';
    if (!condition) ++failures;
}

bool close(double lhs, double rhs, double relative = 1.0e-9)
{
    const double scale = std::max({std::abs(lhs), std::abs(rhs), 1.0});
    return std::abs(lhs - rhs) <= relative*scale;
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

std::vector<double> composition(
    const SpeciesRegistry& registry,
    std::initializer_list<std::pair<const char*, double>> entries)
{
    std::vector<double> result(
        static_cast<std::size_t>(registry.size()), 0.0);
    for (const auto& entry : entries)
    {
        const int index = registry.index(entry.first);
        if (index >= 0)
            result[static_cast<std::size_t>(index)] = entry.second;
    }
    return result;
}

std::vector<double> blend(
    const std::vector<double>& left,
    const std::vector<double>& right,
    double fraction)
{
    std::vector<double> result(left.size(), 0.0);
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = (1.0 - fraction)*left[i] + fraction*right[i];
    return result;
}

BinaryInteractionDB zeroCompleted(
    const BinaryInteractionDB& source,
    const std::vector<std::string>& names)
{
    BinaryInteractionDB result = source;
    for (std::size_t i = 0; i < names.size(); ++i)
    {
        for (std::size_t j = i + 1; j < names.size(); ++j)
        {
            if (!result.has(names[i], names[j]))
                result.add(names[i], names[j], 0.0);
        }
    }
    return result;
}

std::optional<double> equilibriumEnthalpy(
    const PRTPFlash& tp,
    const PRPhaseCaloric& caloric,
    double pressure,
    double temperature,
    const std::vector<double>& massFractions)
{
    const TPFlashResult flash = tp.flashMassFractions(
        pressure, temperature, massFractions);
    if (!flash.valid()) return std::nullopt;

    if (flash.region == EquilibriumRegion::TWO_PHASE)
    {
        const PRPhaseCaloricResult liquid = caloric.evaluate(
            pressure, temperature, flash.liquidMoleFractions,
            PhaseRoot::LIQUID);
        const PRPhaseCaloricResult vapor = caloric.evaluate(
            pressure, temperature, flash.vaporMoleFractions,
            PhaseRoot::VAPOR);
        if (!liquid.success() || !vapor.success()) return std::nullopt;
        const ConservedMixtureResult mixture = assembleConservedPhaseMixture(
            pressure,
            flash.vaporMoleFraction,
            conserved(*liquid.state),
            conserved(*vapor.state),
            flash.feedMassFractions);
        return mixture.success()
            ? std::optional<double>(mixture.state->h) : std::nullopt;
    }

    PhaseRoot root = PhaseRoot::VAPOR;
    if (flash.region == EquilibriumRegion::SINGLE_LIQUID)
        root = PhaseRoot::LIQUID;
    else if (flash.region != EquilibriumRegion::SINGLE_VAPOR
        && flash.region != EquilibriumRegion::SINGLE_SUPERCRITICAL)
        return std::nullopt;

    const PRPhaseCaloricResult single = caloric.evaluate(
        pressure, temperature, flash.singlePhaseMoleFractions, root);
    return single.success()
        ? std::optional<double>(single.state->h) : std::nullopt;
}

} // End anonymous namespace

int main()
{
    const SpeciesCriticalDB critical = SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB kij = BinaryInteractionDB::createDefault();
    const std::vector<double> pure = {1.0};

    std::cout << "=== Strict Peng-Robinson PH flash ===\n";

    PHFlashOptions propaneOptions;
    propaneOptions.minimumTemperature = 300.0;
    propaneOptions.maximumTemperature = 340.0;
    propaneOptions.scanIntervals = 8;
    PRPHFlash propanePH(
        critical, kij, {"C3H8"}, propaneOptions);
    check(propanePH.configurationStatus() == PHFlashStatus::SUCCESS,
        "C3H8 PH solver configures with an explicit interval");

    PRTPFlash propaneTP(critical, kij, {"C3H8"});
    const PureSaturationResult propaneSaturation =
        propaneTP.solvePureSaturationTemperature(
            1.0e6,
            propaneOptions.minimumTemperature,
            propaneOptions.maximumTemperature);
    check(propaneSaturation.valid(),
        "C3H8 1 MPa saturation lies in the supported caloric interval");

    PRPhaseCaloric propaneCaloric(critical, kij, {"C3H8"});
    PRPhaseCaloricResult saturatedLiquid;
    PRPhaseCaloricResult saturatedVapor;
    if (propaneSaturation.valid())
    {
        saturatedLiquid = propaneCaloric.evaluate(
            1.0e6,
            propaneSaturation.temperature,
            pure,
            PhaseRoot::LIQUID);
        saturatedVapor = propaneCaloric.evaluate(
            1.0e6,
            propaneSaturation.temperature,
            pure,
            PhaseRoot::VAPOR);
    }
    check(saturatedLiquid.success() && saturatedVapor.success(),
        "both C3H8 saturation roots have supported calorics");
    check(saturatedLiquid.success() && saturatedVapor.success()
            && saturatedVapor.state->h > saturatedLiquid.state->h
                + propaneOptions.enthalpyAbsoluteTolerance
                + propaneOptions.enthalpyRelativeTolerance
                    *std::max({
                        std::abs(saturatedLiquid.state->h),
                        std::abs(saturatedVapor.state->h),
                        1.0}),
        "pure PH lever requires positive numerically resolved latent enthalpy");

    constexpr double requestedBeta = 0.37;
    if (saturatedLiquid.success() && saturatedVapor.success())
    {
        const double targetH = (1.0 - requestedBeta)*saturatedLiquid.state->h
            + requestedBeta*saturatedVapor.state->h;
        const PHFlashResult result = propanePH.flashMassFractions(
            1.0e6, targetH, pure);
        if (!result.success())
        {
            std::cout << "    status=" << phFlashStatusName(result.status)
                      << " phase=" << statusName(result.phaseStatus)
                      << " caloric="
                      << phaseCaloricStatusName(result.caloricStatus)
                      << " hResidual=" << result.enthalpyResidual << '\n';
        }
        check(result.success(),
            "pure C3H8 PH coexistence flash succeeds");
        check(result.pureComponentEnthalpyLeverRule,
            "pure coexistence explicitly records the enthalpy lever rule");
        check(close(result.temperature, propaneSaturation.temperature, 1.0e-10),
            "PH temperature equals the independently solved saturation point");
        check(result.tpFlash.region
                == EquilibriumRegion::PURE_COMPONENT_COEXISTENCE
            && !result.tpFlash.phaseFractionDetermined
            && std::isnan(result.tpFlash.vaporMoleFraction),
            "embedded TP payload preserves underdetermined pure beta");
        check(close(result.vaporMoleFraction, requestedBeta, 1.0e-10),
            "PH enthalpy constraint determines pure mass/mole beta");
        check(result.liquidPhase.has_value()
                && result.vaporPhase.has_value()
                && result.conservedMixture.has_value()
                && !result.singlePhase.has_value(),
            "successful coexistence exposes complete mass-basis phase state");
        check(close(result.computedEnthalpy, targetH, 1.0e-12)
                && std::abs(result.enthalpyResidual)
                    <= propaneOptions.enthalpyAbsoluteTolerance
                        + propaneOptions.enthalpyRelativeTolerance
                            *std::max(std::abs(targetH), 1.0),
            "pure PH result satisfies the strict enthalpy residual");
        check(result.tpFugacityResidual
                <= propanePH.tpFlash().options().fugacityTolerance
                && result.tpMaterialResidual
                    <= propanePH.tpFlash().options().materialTolerance,
            "pure PH success rechecks TP fugacity and material closure");
        if (result.conservedMixture)
        {
            check(close(
                    result.conservedMixture->vaporMassFraction,
                    requestedBeta,
                    1.0e-10),
                "pure conserved state carries the PH vapor fraction");
            check(close(
                    result.conservedMixture->e,
                    result.conservedMixture->h
                        - 1.0e6/result.conservedMixture->rho,
                    1.0e-12),
                "pure conserved state closes e=h-p/rho");
        }

        const PHFlashResult noEndpointFallback =
            propanePH.flashMassFractions(1.0e6, 1.0e12, pure);
        check(noEndpointFallback.status
                == PHFlashStatus::NO_ENTHALPY_BRACKET
            && !noEndpointFallback.success()
            && !noEndpointFallback.singlePhase.has_value()
            && !noEndpointFallback.conservedMixture.has_value(),
            "unbracketed target never succeeds at a scan endpoint");

        const double latentEnthalpy =
            saturatedVapor.state->h - saturatedLiquid.state->h;
        PHFlashOptions unresolvedLatentOptions = propaneOptions;
        unresolvedLatentOptions.enthalpyAbsoluteTolerance =
            2.0*latentEnthalpy;
        unresolvedLatentOptions.enthalpyRelativeTolerance = 0.0;
        PRPHFlash unresolvedLatentPH(
            critical, kij, {"C3H8"}, unresolvedLatentOptions);
        const PHFlashResult unresolvedLatent =
            unresolvedLatentPH.flashMassFractions(
                1.0e6,
                0.5*(saturatedLiquid.state->h + saturatedVapor.state->h),
                pure);
        check(unresolvedLatent.status
                == PHFlashStatus::NONPOSITIVE_LATENT_ENTHALPY
            && !unresolvedLatent.success(),
            "pure PH rejects latent enthalpy below its numerical resolution");
        check(!unresolvedLatent.liquidPhase.has_value()
                && !unresolvedLatent.vaporPhase.has_value()
                && !unresolvedLatent.singlePhase.has_value()
                && !unresolvedLatent.conservedMixture.has_value(),
            "unresolved latent enthalpy exposes no committable state");
    }

    PHFlashOptions oxygenOptions;
    oxygenOptions.minimumTemperature = 120.0;
    oxygenOptions.maximumTemperature = 180.0;
    oxygenOptions.scanIntervals = 4;
    PRPHFlash oxygenPH(critical, kij, {"O2"}, oxygenOptions);
    PRPhaseCaloric oxygenReferenceCaloric(critical, kij, {"O2"});
    const PRPhaseCaloricResult oxygenReference =
        oxygenReferenceCaloric.evaluate(
            4.0e6, 160.0, pure, PhaseRoot::VAPOR);
    check(oxygenReference.success(),
        "low-temperature O2 PH target has an independently built state");
    const double oxygenTargetEnthalpy = oxygenReference.state
        ? oxygenReference.state->h : 0.0;
    const PHFlashResult oxygenLowTemperature =
        oxygenPH.flashMassFractions(4.0e6, oxygenTargetEnthalpy, pure);
    check(oxygenLowTemperature.success(),
        "O2 4 MPa PH uses the frozen low-temperature caloric table");
    check(oxygenLowTemperature.caloricStatus ==
            PRPhaseCaloricStatus::SUCCESS
        && oxygenLowTemperature.phaseStatus ==
            PhaseEquilibriumStatus::SUCCESS
        && (oxygenLowTemperature.singlePhase.has_value()
            || oxygenLowTemperature.conservedMixture.has_value()),
        "low-temperature O2 PH returns a complete state payload");
    check(oxygenLowTemperature.enthalpyRelativeResidual < 1.0e-10,
        "low-temperature O2 PH closes the requested enthalpy");

    PHFlashOptions oxygenOutOptions;
    oxygenOutOptions.minimumTemperature = 70.0;
    oxygenOutOptions.maximumTemperature = 75.0;
    oxygenOutOptions.scanIntervals = 4;
    PRPHFlash oxygenOutOfTable(critical, kij, {"O2"}, oxygenOutOptions);
    const PHFlashResult oxygenOutResult =
        oxygenOutOfTable.flashMassFractions(
            4.0e6, oxygenTargetEnthalpy, pure);
    check(!oxygenOutResult.success()
            && !oxygenOutResult.singlePhase.has_value()
            && !oxygenOutResult.conservedMixture.has_value()
            && oxygenOutResult.status != PHFlashStatus::SUCCESS,
        "O2 PH below the frozen table fails closed without an endpoint state");

    PHFlashOptions singleOptions;
    singleOptions.minimumTemperature = 400.0;
    singleOptions.maximumTemperature = 460.0;
    singleOptions.scanIntervals = 6;
    PRPHFlash singlePH(critical, kij, {"C3H8"}, singleOptions);
    constexpr double singleReferenceTemperature = 433.7;
    const auto referenceSingle = propaneCaloric.evaluate(
        1.0e6, singleReferenceTemperature, pure, PhaseRoot::VAPOR);
    check(referenceSingle.success(),
        "single-vapor PH reference state is available");
    if (referenceSingle.success())
    {
        const PHFlashResult recovered = singlePH.flashMassFractions(
            1.0e6, referenceSingle.state->h, pure);
        check(recovered.success() && recovered.singlePhase.has_value(),
            "single-vapor PH inversion returns mass-basis properties");
        check(!recovered.conservedMixture.has_value()
                && !recovered.liquidPhase.has_value()
                && !recovered.vaporPhase.has_value(),
            "single-phase result does not invent a two-phase mixture");
        check(close(
                recovered.temperature,
                singleReferenceTemperature,
                1.0e-10),
            "single-vapor PH bisection recovers an off-grid temperature");
        if (recovered.singlePhase)
        {
            check(close(
                    recovered.singlePhase->e,
                    recovered.singlePhase->h
                        - 1.0e6/recovered.singlePhase->rho,
                    1.0e-12),
                "single-vapor phase closes e=h-p/rho");
        }
    }

    // Frozen supported O2/C3H8 two-phase state. Its temperature is an exact
    // PH scan sample so this test isolates TP-to-caloric-to-conserved assembly
    // from the independent off-grid bisection test above.
    PRTPFlash mixtureTP(critical, kij, {"O2", "C3H8"});
    PRPhaseCaloric mixtureCaloric(critical, kij, {"O2", "C3H8"});
    constexpr double mixturePressure = 2.0e6;
    constexpr double mixtureTemperature = 300.0;
    const TPFlashResult mixtureReference = mixtureTP.flashMoleFractions(
        mixturePressure, mixtureTemperature, {0.2, 0.8});
    check(mixtureReference.valid()
            && mixtureReference.region == EquilibriumRegion::TWO_PHASE,
        "frozen multicomponent reference remains two phase");
    check(mixtureReference.valid()
            && close(
                mixtureReference.vaporMoleFraction,
                0.41961363166,
                2.0e-8),
        "frozen multicomponent TP beta remains stable");
    PRPhaseCaloricResult mixtureLiquid;
    PRPhaseCaloricResult mixtureVapor;
    if (mixtureReference.valid())
    {
        mixtureLiquid = mixtureCaloric.evaluate(
            mixturePressure,
            mixtureTemperature,
            mixtureReference.liquidMoleFractions,
            PhaseRoot::LIQUID);
        mixtureVapor = mixtureCaloric.evaluate(
            mixturePressure,
            mixtureTemperature,
            mixtureReference.vaporMoleFractions,
            PhaseRoot::VAPOR);
    }
    if (mixtureReference.valid()
        && mixtureLiquid.success() && mixtureVapor.success())
    {
        const auto referenceMixture = assembleConservedPhaseMixture(
            mixturePressure,
            mixtureReference.vaporMoleFraction,
            conserved(*mixtureLiquid.state),
            conserved(*mixtureVapor.state),
            mixtureReference.feedMassFractions);
        check(referenceMixture.success(),
            "multicomponent reference conserved state closes");
        if (referenceMixture.success())
        {
            PHFlashOptions mixtureOptions;
            mixtureOptions.minimumTemperature = mixtureTemperature - 5.0;
            mixtureOptions.maximumTemperature = mixtureTemperature + 5.0;
            mixtureOptions.scanIntervals = 2;
            PRPHFlash mixturePH(
                critical,
                kij,
                {"O2", "C3H8"},
                mixtureOptions);
            const PHFlashResult recovered = mixturePH.flashMassFractions(
                mixturePressure,
                referenceMixture.state->h,
                mixtureReference.feedMassFractions);
            if (!recovered.success())
            {
                std::cout << "    mixture status="
                          << phFlashStatusName(recovered.status)
                          << " phase=" << statusName(recovered.phaseStatus)
                          << " hResidual=" << recovered.enthalpyResidual
                          << '\n';
            }
            check(recovered.success()
                    && recovered.tpFlash.region
                        == EquilibriumRegion::TWO_PHASE,
                "multicomponent PH uses the equilibrium TP split");
            check(!recovered.pureComponentEnthalpyLeverRule
                    && close(
                        recovered.vaporMoleFraction,
                        mixtureReference.vaporMoleFraction,
                        1.0e-8),
                "multicomponent PH preserves TP molar beta");
            check(recovered.conservedMixture.has_value()
                    && close(
                        recovered.computedEnthalpy,
                        referenceMixture.state->h,
                        1.0e-10),
            "multicomponent PH returns the complete conserved mixture");
        }
    }

    BinaryInteractionDB zeroKij;
    zeroKij.add("O2", "C3H8", 0.0);
    PRTPFlash wilsonFirstTP(critical, zeroKij, {"O2", "C3H8"});
    PRPhaseCaloric wilsonFirstCaloric(
        critical, zeroKij, {"O2", "C3H8"});
    const std::vector<double> oxygenMassFractions = {1.0, 0.0};
    const std::vector<double> propaneMassFractions = {0.0, 1.0};
    const std::vector<double> gapMassFractions = {0.9, 0.1};
    const std::optional<double> oxygenInletEnthalpy = equilibriumEnthalpy(
        wilsonFirstTP,
        wilsonFirstCaloric,
        6.0e6,
        130.0,
        oxygenMassFractions);
    const std::optional<double> propaneInletEnthalpy = equilibriumEnthalpy(
        wilsonFirstTP,
        wilsonFirstCaloric,
        6.0e6,
        450.0,
        propaneMassFractions);
    check(oxygenInletEnthalpy.has_value()
            && propaneInletEnthalpy.has_value(),
        "6 MPa O2/C3H8 inlet enthalpies are constructible");
    if (oxygenInletEnthalpy && propaneInletEnthalpy)
    {
        const double targetEnthalpy = 0.9*(*oxygenInletEnthalpy)
            + 0.1*(*propaneInletEnthalpy);
        std::vector<double> recoveredTemperatures;
        for (std::size_t intervals : {512u, 1024u, 2048u})
        {
            PHFlashOptions gapOptions;
            gapOptions.minimumTemperature = 80.0;
            gapOptions.maximumTemperature = 700.0;
            gapOptions.scanIntervals = intervals;
            gapOptions.maximumIterations = 240;
            PRPHFlash gapPH(
                critical,
                zeroKij,
                {"O2", "C3H8"},
                gapOptions);
            const PHFlashResult recovered = gapPH.flashMassFractions(
                6.0e6, targetEnthalpy, gapMassFractions);
            std::cout << "    WILSON_FIRST_PH intervals=" << intervals
                      << " status=" << phFlashStatusName(recovered.status)
                      << " phase=" << statusName(recovered.phaseStatus)
                      << " T=" << std::setprecision(17)
                      << recovered.temperature << '\n';
            check(recovered.success()
                    && recovered.tpFlash.region
                        == EquilibriumRegion::TWO_PHASE,
                "Wilson-first PH grid recovers the two-phase root");
            check(recovered.success()
                    && std::abs(recovered.enthalpyResidual)
                        <= gapOptions.enthalpyAbsoluteTolerance
                            + gapOptions.enthalpyRelativeTolerance
                                *std::max(std::abs(targetEnthalpy), 1.0),
                "Wilson-first PH grid closes the inlet enthalpy");
            if (recovered.success())
                recoveredTemperatures.push_back(recovered.temperature);
        }
        check(recoveredTemperatures.size() == 3,
            "all PH scan resolutions return a committable state");
        if (recoveredTemperatures.size() == 3)
        {
            const auto limits = std::minmax_element(
                recoveredTemperatures.begin(),
                recoveredTemperatures.end());
            check(std::abs(*limits.second - *limits.first) < 1.0e-8,
                "PH recovery is invariant across 512/1024/2048 intervals");
            check(std::abs(
                    recoveredTemperatures.front()
                    - 163.05364687505062) < 1.0e-7,
                "Wilson-first PH recovery retains the frozen root");
        }
    }

    struct BoundaryRefinementState
    {
        double pressure;
        double kij;
        double fuelMassFraction;
        double expectedTemperature;
        bool shouldSucceed;
        const char* label;
    };
    const BoundaryRefinementState boundaryRefinementStates[] = {
        {
            4.0e6, 0.15, 0.05, 147.89476370962802, true,
            "4 MPa kij=0.15 Yfuel=0.05"
        },
        {
            4.0e6, 0.055, 0.10, 150.59809070531628, true,
            "4 MPa kij=0.055 Yfuel=0.10"
        },
        {
            5.0e6, 0.055, 0.10, 156.22545782696614, false,
            "5 MPa kij=0.055 Yfuel=0.10"
        },
        {
            5.0e6, 0.15, 0.10, 155.1625361408079, false,
            "5 MPa kij=0.15 Yfuel=0.10"
        }
    };
    for (const BoundaryRefinementState& state : boundaryRefinementStates)
    {
        BinaryInteractionDB stateKij;
        stateKij.add("O2", "C3H8", state.kij);
        const PRTPFlash stateTP(critical, stateKij, {"O2", "C3H8"});
        const PRPhaseCaloric stateCaloric(
            critical, stateKij, {"O2", "C3H8"});
        const std::optional<double> stateOxygenH = equilibriumEnthalpy(
            stateTP,
            stateCaloric,
            state.pressure,
            130.0,
            oxygenMassFractions);
        const std::optional<double> statePropaneH = equilibriumEnthalpy(
            stateTP,
            stateCaloric,
            state.pressure,
            450.0,
            propaneMassFractions);
        const std::string prefix = std::string(state.label) + " ";
        check(stateOxygenH.has_value() && statePropaneH.has_value(),
            prefix + "inlet enthalpies are constructible");
        if (!stateOxygenH || !statePropaneH) continue;

        const std::vector<double> stateFeed = {
            1.0 - state.fuelMassFraction, state.fuelMassFraction};
        const double stateTarget =
            (1.0 - state.fuelMassFraction)*(*stateOxygenH)
            + state.fuelMassFraction*(*statePropaneH);
        std::vector<double> stateTemperatures;
        for (std::size_t intervals : {512u, 1024u, 2048u})
        {
            PHFlashOptions stateOptions;
            stateOptions.minimumTemperature = 80.0;
            stateOptions.maximumTemperature = 700.0;
            stateOptions.scanIntervals = intervals;
            stateOptions.maximumIterations = 240;
            const PRPHFlash statePH(
                critical,
                stateKij,
                {"O2", "C3H8"},
                stateOptions);
            const PHFlashResult recovered = statePH.flashMassFractions(
                state.pressure, stateTarget, stateFeed);
            stateTemperatures.push_back(recovered.temperature);
            std::cout << "    BOUNDARY_REFINEMENT label=" << state.label
                      << " intervals=" << intervals
                      << " status=" << phFlashStatusName(recovered.status)
                      << " phase=" << statusName(recovered.phaseStatus)
                      << " T=" << std::setprecision(17)
                      << recovered.temperature
                      << " hResidual=" << recovered.enthalpyResidual << '\n';
            if (state.shouldSucceed)
            {
                const double enthalpyTolerance =
                    stateOptions.enthalpyAbsoluteTolerance
                    + stateOptions.enthalpyRelativeTolerance
                        *std::max(std::abs(stateTarget), 1.0);
                check(recovered.success()
                        && recovered.tpFlash.region
                            == EquilibriumRegion::TWO_PHASE,
                    prefix + "boundary refinement recovers the PH root");
                check(recovered.success()
                        && std::abs(recovered.enthalpyResidual)
                            <= enthalpyTolerance,
                    prefix + "boundary refinement closes enthalpy");
            }
            else
            {
                check(recovered.status
                            == PHFlashStatus::ENTHALPY_RESIDUAL_FAILURE
                        && !recovered.success()
                        && !recovered.singlePhase.has_value()
                        && !recovered.liquidPhase.has_value()
                        && !recovered.vaporPhase.has_value()
                        && !recovered.conservedMixture.has_value(),
                    prefix + "topology discontinuity remains fail-closed");
            }
        }
        const auto stateLimits = std::minmax_element(
            stateTemperatures.begin(), stateTemperatures.end());
        check(stateTemperatures.size() == 3
                && std::abs(*stateLimits.second - *stateLimits.first) < 1.0e-12,
            prefix + "512/1024/2048 diagnostic temperature is invariant");
        check(stateTemperatures.size() == 3
                && std::abs(
                    stateTemperatures.front()
                    - state.expectedTemperature) < 1.0e-8,
            prefix + "retains the frozen diagnostic temperature");
    }

    // Diagnostic zero completion isolates the PH scan algorithm from the
    // separate production kij-data blocker. The first active H2O/CO2 fraction
    // raises the common caloric lower bound to 200 K.
    const SpeciesRegistry registry;
    const std::vector<std::string> productionSpecies = registry.names();
    const BinaryInteractionDB completedKij = zeroCompleted(
        kij, productionSpecies);
    TPFlashOptions boundaryTPOptions;
    boundaryTPOptions.maxFlashIterations = 2400;
    boundaryTPOptions.maxStabilityIterations = 2400;
    PRTPFlash boundaryTP(
        critical, completedKij, productionSpecies, boundaryTPOptions);
    PRPhaseCaloric boundaryCaloric(
        critical, completedKij, productionSpecies);
    const std::vector<double> lox = composition(
        registry, {{"O2", 1.0}});
    const std::vector<double> background = composition(
        registry, {{"CO2", 0.6468}, {"H2O", 0.3532}});
    const std::optional<double> loxEnthalpy = equilibriumEnthalpy(
        boundaryTP, boundaryCaloric, 4.0e6, 130.0, lox);
    const std::optional<double> backgroundEnthalpy = equilibriumEnthalpy(
        boundaryTP, boundaryCaloric, 4.0e6, 700.0, background);
    check(loxEnthalpy.has_value() && backgroundEnthalpy.has_value(),
        "LOX/background PH boundary endpoints are constructible");
    if (loxEnthalpy && backgroundEnthalpy)
    {
        PHFlashOptions boundaryOptions;
        boundaryOptions.minimumTemperature = 80.0;
        boundaryOptions.maximumTemperature = 2600.0;
        boundaryOptions.scanIntervals = 4096;
        boundaryOptions.maximumIterations = 400;
        PRPHFlash boundaryPH(
            critical,
            completedKij,
            productionSpecies,
            boundaryOptions,
            boundaryTPOptions);

        constexpr double belowBoundaryFraction = 0.1149;
        constexpr double aboveBoundaryFraction = 0.1150;
        const std::vector<double> belowComposition = blend(
            lox, background, belowBoundaryFraction);
        const std::vector<double> aboveComposition = blend(
            lox, background, aboveBoundaryFraction);
        const double belowTarget =
            (1.0 - belowBoundaryFraction)*(*loxEnthalpy)
            + belowBoundaryFraction*(*backgroundEnthalpy);
        const double aboveTarget =
            (1.0 - aboveBoundaryFraction)*(*loxEnthalpy)
            + aboveBoundaryFraction*(*backgroundEnthalpy);

        const PRPhaseCaloricTemperatureRange activeRange =
            boundaryCaloric.supportedTemperatureRange(aboveComposition);
        check(activeRange.success() && close(activeRange.minimum, 200.0),
            "active LOX/background calorics expose the 200 K lower bound");

        const std::optional<double> belowBoundaryEnthalpy =
            equilibriumEnthalpy(
                boundaryTP,
                boundaryCaloric,
                4.0e6,
                activeRange.minimum,
                belowComposition);
        const std::optional<double> aboveBoundaryEnthalpy =
            equilibriumEnthalpy(
                boundaryTP,
                boundaryCaloric,
                4.0e6,
                activeRange.minimum,
                aboveComposition);
        check(belowBoundaryEnthalpy
                && *belowBoundaryEnthalpy > belowTarget,
            "f=0.1149 target remains below the lowest supported state");
        check(aboveBoundaryEnthalpy
                && *aboveBoundaryEnthalpy < aboveTarget,
            "f=0.115 target lies above the 200 K boundary state");

        const PHFlashResult belowBoundary = boundaryPH.flashMassFractions(
            4.0e6, belowTarget, belowComposition);
        check(belowBoundary.status == PHFlashStatus::NO_ENTHALPY_BRACKET
                && !belowBoundary.success()
                && !belowBoundary.singlePhase.has_value()
                && !belowBoundary.liquidPhase.has_value()
                && !belowBoundary.vaporPhase.has_value()
                && !belowBoundary.conservedMixture.has_value(),
            "target below the common caloric domain remains unbracketed");

        const PHFlashResult aboveBoundary = boundaryPH.flashMassFractions(
            4.0e6, aboveTarget, aboveComposition);
        std::cout << "    CALORIC_BOUNDARY_PH status="
                  << phFlashStatusName(aboveBoundary.status)
                  << " phase=" << statusName(aboveBoundary.phaseStatus)
                  << " T=" << std::setprecision(17)
                  << aboveBoundary.temperature
                  << " hResidual=" << aboveBoundary.enthalpyResidual << '\n';
        check(aboveBoundary.success()
                && aboveBoundary.tpFlash.region
                    == EquilibriumRegion::TWO_PHASE
                && aboveBoundary.liquidPhase.has_value()
                && aboveBoundary.vaporPhase.has_value()
                && aboveBoundary.conservedMixture.has_value()
                && close(
                    aboveBoundary.temperature,
                    200.10104199172929,
                    1.0e-8)
                && std::abs(aboveBoundary.enthalpyResidual)
                    <= boundaryOptions.enthalpyAbsoluteTolerance
                        + boundaryOptions.enthalpyRelativeTolerance
                            *std::max(std::abs(aboveTarget), 1.0),
            "PH scan explicitly brackets the root above the caloric boundary");
    }

    PHFlashOptions omittedInterval;
    PRPHFlash invalidOptions(
        critical, kij, {"C3H8"}, omittedInterval);
    check(invalidOptions.configurationStatus()
            == PHFlashStatus::INVALID_TEMPERATURE_INTERVAL,
        "omitting the PH temperature domain fails configuration");

    if (failures != 0)
    {
        std::cerr << "PR_PH_FLASH_FAIL failures=" << failures << '\n';
        return 1;
    }

    std::cout << "PR_PH_FLASH_PASS\n"
              << "NOTE: this is a strict PH foundation; phase derivatives, "
              << "two-phase sound speed, transport, chemistry coupling, and "
              << "rho-E flash remain outside this gate; the first ascending "
              << "T bracket is not a global uniqueness proof.\n";
    return 0;
}
