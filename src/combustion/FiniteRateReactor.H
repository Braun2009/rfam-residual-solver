#ifndef FiniteRateReactor_H
#define FiniteRateReactor_H

/*---------------------------------------------------------------------------*\
    FiniteRateReactor - 0-D finite-rate chemistry sub-step bridging
    chem-space (ChemistryModel + CVODE) and global-space (SpeciesRegistry).

    Stage 1.5 / Phase 7 core. Used by ChemistrySourceOperator (CFD) but kept
    free of the framework (no Domain/EOS) so it can be validated in isolation.

    reactConservative():
      - integrate [Y,e] with real-fluid rho-e-Y temperature recovery
      - scatter the conservative species and sensible-energy increment
      - validate positivity, mass, elements, and the final rho-e state.

    ENERGY CONVENTION (sensible-energy, see [[rfam-species-energy-foundation]]):
    the caller adds dE to rhoE EXPLICITLY. Formation enthalpy is NOT in e.
    The old fixed-Cv [Y,T] reactor remains only behind an explicitly named
    LegacyFixedCvDiagnostic constructor and reactLegacyFixedCv() method. It is
    not constructed or called by CFD production operators.
\*---------------------------------------------------------------------------*/

#include "SpeciesRegistry.H"
#include "SpeciesMapper.H"
#include "FiniteRate/ChemistryModel.H"
#include "FiniteRate/CVODEWrapper.H"
#include "FiniteRate/RealFluidCVODEWrapper.H"
#include "NASA7Thermo.H"
#include "RealFluidState.H"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Foam
{

class FiniteRateReactor
{
public:
    struct LegacyFixedCvDiagnostic
    {
        explicit LegacyFixedCvDiagnostic(double value) : Cv(value) {}
        double Cv;
    };

private:
    ChemistryModel*       chem_;
    bool explicitMode_ = false;
    int explicitSubSteps_ = 256;
    mutable int lastExplicitSubSteps_ = 0;
    const SpeciesRegistry& reg_;
    SpeciesMapper          mapper_;
    double                 legacyCv_;
    bool                   legacyUseNasa7Cv_;
    mutable std::unique_ptr<CVODEWrapper> legacyCvode_;
    mutable RealFluidCVODEWrapper realFluidCvode_;
    mutable std::vector<double> yChem_,yChemNew_,omega_;
    mutable int lastSubcycles_=0;
    mutable bool lastPositivityLimited_=false;
    mutable int lastFailureCode_=0;
    mutable int validationFailureCode_=0;
    mutable double lastElementRelativeError_=0.0;
    mutable int lastRejectedChemIndex_=-1;
    mutable double lastRejectedBefore_=0.0,lastRejectedAfter_=0.0;

    double positivityScale(const std::vector<double>& before,
                           const std::vector<double>& after) const
    {
        double scale=1.0;
        for(size_t c=0;c<before.size();++c)
            if(after[c]<0.0&&after[c]<before[c])
                scale=std::min(scale,before[c]/(before[c]-after[c]));
        return std::clamp(scale,0.0,1.0);
    }

    void limitReactionIncrement(const std::vector<double>& before,
                                std::vector<double>& after,double& temperature,
                                double initialTemperature) const
    {
        const double scale=positivityScale(before,after);
        if(scale>=1.0) return;
        lastPositivityLimited_=true;
        // Leave a small transported-state buffer at depletion. A boundary
        // value at machine epsilon can be driven negative by the next explicit
        // convection/diffusion update even when the chemistry step is valid.
        const double interiorScale=scale*(1.0-1e-6);
        for(size_t c=0;c<before.size();++c)
            after[c]=before[c]+interiorScale*(after[c]-before[c]);
        temperature=initialTemperature
            +interiorScale*(temperature-initialTemperature);
    }

    bool admissibleChemState(const std::vector<double>& before,
                             const std::vector<double>& after) const
    {
        validationFailureCode_=0;
        lastElementRelativeError_=0.0;
        if(before.size()!=after.size()) { validationFailureCode_=2; return false; }
        double sumBefore=0.0,sumAfter=0.0;
        double elementsBefore[4]={0.0,0.0,0.0,0.0};
        double elementsAfter[4]={0.0,0.0,0.0,0.0};
        for(int c=0;c<mapper_.nChem();++c)
        {
            if(!std::isfinite(after[c])||after[c]<-1e-12)
            { validationFailureCode_=2; return false; }
            sumBefore+=before[c]; sumAfter+=after[c];
            const int g=mapper_.globalIndex(c);
            const double invMw=1.0/reg_.molecularWeight(g);
            const int atoms[4]={reg_.nC(g),reg_.nH(g),reg_.nO(g),reg_.nN(g)};
            for(int e=0;e<4;++e)
            {
                elementsBefore[e]+=before[c]*atoms[e]*invMw;
                elementsAfter[e]+=after[c]*atoms[e]*invMw;
            }
        }
        // The registry molecular weights are rounded independently, so exact
        // atom balance leaves an O(1e-7) mass-sum residue for this mechanism.
        if(std::abs(sumAfter-sumBefore)>1e-6*std::max(1.0,std::abs(sumBefore)))
        { validationFailureCode_=4; return false; }
        for(int e=0;e<4;++e)
        {
            const double relativeError=std::abs(elementsAfter[e]-elementsBefore[e])
                /std::max(1.0,std::abs(elementsBefore[e]));
            lastElementRelativeError_=std::max(lastElementRelativeError_,relativeError);
            // Keep invariant acceptance five times tighter than CVODE's 1e-6
            // relative integration tolerance. Hot reversible cells near a
            // depletion boundary can reach just over 1e-7 truncation error.
            if(relativeError>2e-7)
            { validationFailureCode_=5; return false; }
        }
        return true;
    }

public:
    FiniteRateReactor(ChemistryModel* chem, const SpeciesRegistry& reg)
        : chem_(chem), reg_(reg), mapper_(reg, chem->speciesNames()),
          legacyCv_(std::numeric_limits<double>::quiet_NaN()),
          legacyUseNasa7Cv_(false),legacyCvode_(nullptr),
          realFluidCvode_(chem,reg,mapper_) {}

    FiniteRateReactor(
        ChemistryModel* chem,
        const SpeciesRegistry& reg,
        std::shared_ptr<const thermo::EquilibriumThermoModel> thermoModel)
        : chem_(chem), reg_(reg), mapper_(reg, chem->speciesNames()),
          legacyCv_(std::numeric_limits<double>::quiet_NaN()),
          legacyUseNasa7Cv_(false), legacyCvode_(nullptr),
          realFluidCvode_(chem, reg, mapper_, std::move(thermoModel))
    {}

    FiniteRateReactor(
        ChemistryModel* chem,
        const SpeciesRegistry& reg,
        std::shared_ptr<const thermo::EquilibriumThermoModel> thermoModel,
        const phaseEquilibrium::PRPhaseRhoERecovery& phaseRecovery)
        : chem_(chem), reg_(reg), mapper_(reg, chem->speciesNames()),
          legacyCv_(std::numeric_limits<double>::quiet_NaN()),
          legacyUseNasa7Cv_(false), legacyCvode_(nullptr),
          realFluidCvode_(
              chem, reg, mapper_, std::move(thermoModel), phaseRecovery)
    {}

    FiniteRateReactor(
        ChemistryModel* chem,
        const SpeciesRegistry& reg,
        std::shared_ptr<const thermo::EquilibriumThermoModel> thermoModel,
        const phaseEquilibrium::PRPhaseRhoERecovery& phaseRecovery,
        const phaseEquilibrium::PRPhaseDisplacementAuthority&
            phaseDisplacement)
        : chem_(chem), reg_(reg), mapper_(reg, chem->speciesNames()),
          legacyCv_(std::numeric_limits<double>::quiet_NaN()),
          legacyUseNasa7Cv_(false), legacyCvode_(nullptr),
          realFluidCvode_(
              chem, reg, mapper_, std::move(thermoModel), phaseRecovery,
              phaseDisplacement)
    {}

    FiniteRateReactor(
        ChemistryModel* chem,const SpeciesRegistry& reg,
        LegacyFixedCvDiagnostic diagnostic)
        : FiniteRateReactor(chem,reg)
    {
        if(!std::isfinite(diagnostic.Cv)||diagnostic.Cv<=0.0)
            throw std::invalid_argument(
                "FiniteRateReactor: invalid legacy diagnostic Cv");
        legacyCv_=diagnostic.Cv;
        legacyCvode_=std::make_unique<CVODEWrapper>(chem,legacyCv_);
    }

    const SpeciesMapper& mapper() const { return mapper_; }
    int lastSubcycles() const { return lastSubcycles_; }
    bool lastPositivityLimited() const { return lastPositivityLimited_; }

    // Explicit chemistry mode: freeze (rho, T) per sub-step and advance the
    // species/energy with the reaction source terms directly (no CVODE).
    // DIAGNOSTIC: sub-step count is fixed; use for profiling/benchmarking
    // against CVODE, not as a production integrator.
    void setExplicitMode(bool on, int subSteps = 256)
    {
        explicitMode_ = on;
        explicitSubSteps_ = subSteps;
    }
    bool explicitMode() const { return explicitMode_; }
    int lastExplicitSubSteps() const { return lastExplicitSubSteps_; }
    int lastFailureCode() const { return lastFailureCode_; }
    double lastElementRelativeError() const { return lastElementRelativeError_; }
    int lastRejectedChemIndex() const { return lastRejectedChemIndex_; }
    double lastRejectedBefore() const { return lastRejectedBefore_; }
    double lastRejectedAfter() const { return lastRejectedAfter_; }
    bool hasInjectedThermoModel() const
    {
        return realFluidCvode_.hasInjectedThermoModel();
    }
    bool hasPhaseThermoOwner() const
    {
        return realFluidCvode_.hasPhaseThermoOwner();
    }
    phaseEquilibrium::PRPhaseThermoOwnerStatus
    lastPhaseOwnerStatus() const
    {
        return realFluidCvode_.lastPhaseOwnerStatus();
    }
    phaseEquilibrium::PRPhaseDisplacementStatus
    lastPhaseDisplacementStatus() const
    {
        return realFluidCvode_.lastPhaseDisplacementStatus();
    }
    const std::optional<phaseEquilibrium::PRPhaseDisplacementDiagnostics>&
    lastPhaseDisplacementDiagnostics() const
    {
        return realFluidCvode_.lastPhaseDisplacementDiagnostics();
    }
    std::uint64_t lastPhaseProbeCount() const
    {
        return realFluidCvode_.phaseProbeCount();
    }
    std::uint64_t lastPhaseAcceptCount() const
    {
        return realFluidCvode_.phaseAcceptCount();
    }
    const std::string& boundThermoModelId() const
    {
        return realFluidCvode_.boundThermoModelId();
    }

    void enableLegacyNasa7CvDiagnostic()
    {
        if(!legacyCvode_)
            throw std::logic_error(
                "legacy NASA-7 Cv requested on a production reactor");
        legacyUseNasa7Cv_=true;
    }

    double estimateConservativeMaxDt(
        double rho, double p, double T,
        const std::vector<double>& Yglobal,
        double safetyFactor = 0.1) const
    {
        if (!std::isfinite(rho) || !std::isfinite(p) || !std::isfinite(T) ||
            rho <= 0.0 || p <= 0.0 || T <= 0.0 ||
            Yglobal.size() != static_cast<size_t>(reg_.size()))
            return 0.0;

        double thermoDensity = std::numeric_limits<double>::quiet_NaN();
        double thermoEnergy = std::numeric_limits<double>::quiet_NaN();
        double thermoCv = std::numeric_limits<double>::quiet_NaN();
        if (!realFluidCvode_.initializeAdmittedPTY(
                p, T, Yglobal, thermoDensity, thermoEnergy, thermoCv))
            return 0.0;
        const double densityScale = std::max(
            {std::abs(rho), std::abs(thermoDensity), 1.0});
        if (realFluidCvode_.hasInjectedThermoModel()
            && std::abs(thermoDensity-rho) > 1.0e-7*densityScale)
            return 0.0;

        mapper_.mapGlobalToChem(Yglobal, yChem_);
        omega_.assign(mapper_.nChem(), 0.0);
        double qdot = 0.0;
        chem_->computeSourceTerms(T,p,rho,yChem_,omega_,qdot);
        if (!std::isfinite(qdot)) return 0.0;

        double timeScale = 1.0e30;
        if (std::abs(qdot) > 1.0)
        {
            const double energyScale = std::max(
                {std::abs(thermoEnergy), thermoCv*T, 1.0});
            timeScale = std::min(
                timeScale, rho*energyScale/std::abs(qdot));
        }

        for (int c = 0; c < mapper_.nChem(); ++c)
        {
            if (omega_[c] < 0.0 && yChem_[c] > 0.0)
                timeScale = std::min(
                    timeScale, rho*yChem_[c]/(-omega_[c]));
        }

        return timeScale >= 1.0e29
            ? 1.0e30 : safetyFactor*timeScale;
    }

    bool reactConservative(
        double rho, double specificEnergy,
        std::vector<double>& Yglobal, double dt,
        double& deltaRhoE, double& temperature, double& pressure) const
    {
        deltaRhoE = 0.0;
        lastPositivityLimited_ = false;
        lastFailureCode_ = 0;
        lastSubcycles_ = 0;
        lastRejectedChemIndex_ = -1;
        lastRejectedBefore_ = lastRejectedAfter_ = 0.0;

        if (!std::isfinite(rho) || !std::isfinite(specificEnergy) ||
            !std::isfinite(dt) || rho <= 0.0 || dt < 0.0 ||
            Yglobal.size() != static_cast<size_t>(reg_.size()))
        {
            lastFailureCode_ = 2;
            lastSubcycles_ = -1;
            return false;
        }

        struct ConservativeState
        {
            std::vector<double> Y;
            double energy;
            double temperature;
            double pressure;
        };

        const ConservativeState initial{
            Yglobal,specificEnergy,temperature,pressure};

        double initialThermoTemperature =
            std::numeric_limits<double>::quiet_NaN();
        double initialThermoPressure =
            std::numeric_limits<double>::quiet_NaN();
        if (!realFluidCvode_.recoverAdmittedRhoEY(
                rho, initial.energy, initial.Y, initial.temperature,
                initialThermoTemperature, initialThermoPressure))
        {
            lastFailureCode_ = 2;
            lastSubcycles_ = -1;
            return false;
        }

        ConservativeState current = initial;
        current.temperature = initialThermoTemperature;
        current.pressure = initialThermoPressure;

        if (dt == 0.0)
        {
            temperature = current.temperature;
            pressure = current.pressure;
            return true;
        }

        auto validateSegment = [&](
            const std::vector<double>& beforeGlobal,
            double beforeEnergy,
            std::vector<double>& afterGlobal,
            double& afterEnergy,
            double& afterTemperature,
            double& afterPressure) -> bool
        {
            constexpr double roundoffNegative = 1.0e-12;
            for (int c = 0; c < mapper_.nChem(); ++c)
            {
                const int g = mapper_.globalIndex(c);
                if (afterGlobal[g] >= 0.0) continue;
                if (afterGlobal[g] >= -roundoffNegative)
                {
                    afterGlobal[g] = 0.0;
                    continue;
                }
                if (beforeGlobal[g] <= roundoffNegative)
                {
                    // A species absent at the beginning of the segment cannot
                    // provide a meaningful progress-variable scale. Reject
                    // the step so the outer controller refines it instead of
                    // multiplying the complete reaction increment by zero.
                    lastRejectedChemIndex_ = c;
                    lastRejectedBefore_ = beforeGlobal[g];
                    lastRejectedAfter_ = afterGlobal[g];
                    lastFailureCode_ = 1;
                    return false;
                }
            }

            mapper_.mapGlobalToChem(beforeGlobal,yChem_);
            mapper_.mapGlobalToChem(afterGlobal,yChemNew_);
            const double scale = positivityScale(yChem_,yChemNew_);
            if (scale < 1.0)
            {
                if (!std::isfinite(scale) || scale <= 1.0e-8)
                {
                    for (int c = 0; c < mapper_.nChem(); ++c)
                    {
                        if (yChemNew_[c] < 0.0)
                        {
                            lastRejectedChemIndex_ = c;
                            lastRejectedBefore_ = yChem_[c];
                            lastRejectedAfter_ = yChemNew_[c];
                            break;
                        }
                    }
                    lastFailureCode_ = 1;
                    return false;
                }
                lastPositivityLimited_ = true;
                const double interiorScale = scale*(1.0-1.0e-8);
                for (int c = 0; c < mapper_.nChem(); ++c)
                {
                    const int g = mapper_.globalIndex(c);
                    afterGlobal[g] = beforeGlobal[g]
                        + interiorScale*(afterGlobal[g]-beforeGlobal[g]);
                }
                afterEnergy = beforeEnergy
                    + interiorScale*(afterEnergy-beforeEnergy);
            }

            mapper_.mapGlobalToChem(afterGlobal,yChemNew_);
            if (!admissibleChemState(yChem_,yChemNew_))
            {
                lastFailureCode_ = validationFailureCode_;
                return false;
            }

            double recoveredTemperature =
                std::numeric_limits<double>::quiet_NaN();
            double recoveredPressure =
                std::numeric_limits<double>::quiet_NaN();
            if (!realFluidCvode_.recoverAdmittedRhoEY(
                    rho, afterEnergy, afterGlobal, afterTemperature,
                    recoveredTemperature, recoveredPressure))
            {
                lastFailureCode_ = 2;
                return false;
            }
            afterTemperature = recoveredTemperature;
            afterPressure = recoveredPressure;
            return true;
        };

        auto advanceSegment = [&]
        (
            const ConservativeState& before,
            double segmentDt,
            ConservativeState& after
        ) -> bool
        {
            if (explicitMode_)
            {
                // Explicit frozen-(rho,T) substeps: Y += dt*omega/rho,
                // e += dt*qdot/rho, with positivity clamping.  No CVODE.
                const int nSub = std::max(1, explicitSubSteps_);
                const double dtSub = segmentDt / nSub;
                ConservativeState s = before;
                const int nsp = static_cast<int>(s.Y.size());
                std::vector<double> omega(static_cast<std::size_t>(nsp), 0.0);
                double qdot = 0.0;
                for (int k = 0; k < nSub; ++k)
                {
                    const double pEst = rho * before.temperature * 259.8;
                    omega.assign(static_cast<std::size_t>(nsp), 0.0);
                    qdot = 0.0;
                    try
                    {
                        chem_->computeSourceTerms(
                            before.temperature, pEst, rho, s.Y, omega, qdot);
                    }
                    catch (...)
                    {
                        lastFailureCode_ = 3;
                        return false;
                    }
                    bool negative = false;
                    for (int i = 0; i < nsp; ++i)
                    {
                        s.Y[i] += dtSub * omega[static_cast<std::size_t>(i)] / rho;
                        if (s.Y[i] < 0.0)
                        {
                            s.Y[i] = std::max(0.0, s.Y[i]);
                            negative = true;
                        }
                    }
                    s.energy += dtSub * qdot / rho;
                    if (negative) lastPositivityLimited_ = true;
                    if (!std::isfinite(s.energy))
                    {
                        lastFailureCode_ = 4;
                        return false;
                    }
                }
                lastExplicitSubSteps_ = nSub;
                after = s;
                return validateSegment(
                    before.Y,before.energy,after.Y,after.energy,
                    after.temperature,after.pressure);
            }
            if (!realFluidCvode_.hasInjectedThermoModel())
            {
                realFluidCvode_.init(
                    0.0, rho, before.Y, before.energy,
                    before.temperature);
            }
            else
            {
                try
                {
                    realFluidCvode_.init(
                        0.0, rho, before.Y, before.energy,
                        before.temperature);
                }
                catch (...)
                {
                    lastFailureCode_ = 2;
                    return false;
                }
            }
            after = before;
            if (!realFluidCvode_.advanceTo(
                segmentDt,after.Y,after.energy,
                after.temperature,after.pressure))
            {
                lastFailureCode_ = 1;
                return false;
            }
            return validateSegment(
                before.Y,before.energy,after.Y,after.energy,
                after.temperature,after.pressure);
        };

        constexpr int maxAcceptedSubsteps = 16384;
        constexpr int maxRefinements = 40;
        const double timeTolerance = 16.0*std::numeric_limits<double>::epsilon()
            *std::max(dt,1.0e-30);
        double elapsed = 0.0;
        int acceptedSubsteps = 0;
        bool refined = false;

        while (dt-elapsed > timeTolerance)
        {
            const double remaining = dt-elapsed;
            double trialStep = remaining;
            bool accepted = false;
            for (int refinement = 0; refinement < maxRefinements; ++refinement)
            {
                ConservativeState candidate = current;
                if (advanceSegment(current,trialStep,candidate))
                {
                    current = std::move(candidate);
                    accepted = true;
                    break;
                }
                refined = true;
                trialStep *= 0.5;
            }

            if (!accepted || ++acceptedSubsteps > maxAcceptedSubsteps)
            {
                lastSubcycles_ = -1;
                return false;
            }
            elapsed += trialStep;
        }

        Yglobal = std::move(current.Y);
        deltaRhoE = rho*(current.energy-initial.energy);
        temperature = current.temperature;
        pressure = current.pressure;
        lastSubcycles_ = (refined || acceptedSubsteps > 1)
            ? acceptedSubsteps : 0;
        if (!std::isfinite(deltaRhoE))
        {
            Yglobal = initial.Y;
            deltaRhoE = 0.0;
            temperature = initial.temperature;
            pressure = initial.pressure;
            lastFailureCode_ = 2;
            lastSubcycles_ = -1;
            return false;
        }
        return true;
    }

    // Explicit token-bound phase-local path. CVODE may probe bounded rho-e-Y
    // candidates without mutating the owner; one final accepted state replaces
    // the token before the outer HEM owner performs whole-mixture reflash.
    bool reactPhaseConservative(
        const phaseEquilibrium::PRPhaseThermoParentView& parent,
        std::vector<double>& phaseMassFractions,
        double dt,
        double& deltaPhaseRhoE,
        double& temperature,
        double& pressure) const
    {
        deltaPhaseRhoE = 0.0;
        lastPositivityLimited_ = false;
        lastFailureCode_ = 0;
        lastSubcycles_ = 0;
        lastRejectedChemIndex_ = -1;
        lastRejectedBefore_ = lastRejectedAfter_ = 0.0;

        if (!realFluidCvode_.hasPhaseThermoOwner()
            || !std::isfinite(dt) || dt < 0.0
            || phaseMassFractions.size()
                != static_cast<std::size_t>(reg_.size()))
        {
            lastFailureCode_ = 6;
            lastSubcycles_ = -1;
            return false;
        }

        const std::vector<double> originalMassFractions =
            phaseMassFractions;
        const double originalTemperature = temperature;
        const double originalPressure = pressure;
        RealFluidCVODEPhaseState initialPhase;
        try
        {
            if (!realFluidCvode_.initPhase(
                    0.0, parent, temperature, initialPhase))
            {
                lastFailureCode_ = 6;
                lastSubcycles_ = -1;
                return false;
            }
        }
        catch (...)
        {
            realFluidCvode_.abortPhase();
            lastFailureCode_ = 6;
            lastSubcycles_ = -1;
            return false;
        }

        if (initialPhase.massFractions != originalMassFractions)
        {
            realFluidCvode_.abortPhase();
            lastFailureCode_ = 6;
            lastSubcycles_ = -1;
            return false;
        }

        std::vector<double> candidateMassFractions =
            initialPhase.massFractions;
        double candidateEnergy = initialPhase.specificEnergy;
        double candidateTemperature = initialPhase.temperature;
        double candidatePressure = initialPhase.pressure;
        if (!realFluidCvode_.advanceTo(
                dt, candidateMassFractions, candidateEnergy,
                candidateTemperature, candidatePressure))
        {
            phaseMassFractions = originalMassFractions;
            temperature = originalTemperature;
            pressure = originalPressure;
            lastFailureCode_ = realFluidCvode_.lastFailureIsPermanent()
                ? 6 : 1;
            lastSubcycles_ = -1;
            return false;
        }

        mapper_.mapGlobalToChem(initialPhase.massFractions, yChem_);
        mapper_.mapGlobalToChem(candidateMassFractions, yChemNew_);
        if (!admissibleChemState(yChem_, yChemNew_))
        {
            phaseMassFractions = originalMassFractions;
            temperature = originalTemperature;
            pressure = originalPressure;
            lastFailureCode_ = validationFailureCode_;
            lastSubcycles_ = -1;
            return false;
        }

        const double candidateDelta = initialPhase.rho
            *(candidateEnergy - initialPhase.specificEnergy);
        if (!std::isfinite(candidateDelta))
        {
            phaseMassFractions = originalMassFractions;
            temperature = originalTemperature;
            pressure = originalPressure;
            lastFailureCode_ = 2;
            lastSubcycles_ = -1;
            return false;
        }

        phaseMassFractions = std::move(candidateMassFractions);
        deltaPhaseRhoE = candidateDelta;
        temperature = candidateTemperature;
        pressure = candidatePressure;
        return true;
    }

    // Constant-(rho,p) 0-D reactor sub-step.
    bool reactLegacyFixedCv(double rho, double p, double& T,
                            std::vector<double>& Yglobal, double dt,
                            double& dE) const
    {
        if(!legacyCvode_)
            throw std::logic_error(
                "reactLegacyFixedCv requires LegacyFixedCvDiagnostic");
        // Stage 4: compute NASA-7 Cv_mix(T,Y) and update the CVODE closure.
        double Cv_eff = legacyCv_;
        if (legacyUseNasa7Cv_)
        {
            const std::vector<std::string>& spNames = reg_.names();
            const NASA7Thermo& th = nasa7Thermo();
            // Cv_mix = Cp_mix - R_mix  (ideal gas, mass basis)
            // R_mix = R_univ / MW_mix = R_univ * sum(Y_i / MW_i)
            static constexpr double Runiv = 8.314462618;
            double Cp = th.cp_mass_mix(spNames, Yglobal, T);
            double R_mix = 0.0;
            for (int g = 0; g < reg_.size(); ++g)
                R_mix += Yglobal[g] / reg_.molecularWeight(g) * Runiv;
            Cv_eff = Cp - R_mix;
            if (Cv_eff < 200.0) Cv_eff = 200.0;   // physical floor
            legacyCvode_->setCv(Cv_eff);
        }

        mapper_.mapGlobalToChem(Yglobal, yChem_);

        const double initialTemperature=T;
        lastPositivityLimited_=false;
        lastFailureCode_=0;
        legacyCvode_->init(0.0, p, rho, yChem_, initialTemperature);
        double Tnew = T;
        bool success=legacyCvode_->advanceTo(dt, yChemNew_, Tnew);
        if(!success) lastFailureCode_=1;
        if(success) limitReactionIncrement(yChem_,yChemNew_,Tnew,initialTemperature);
        if(success&&!admissibleChemState(yChem_,yChemNew_))
        { success=false; lastFailureCode_=validationFailureCode_; }
        if(success&&(!std::isfinite(Tnew)||Tnew<=0.0))
        { success=false; lastFailureCode_=3; }
        lastSubcycles_=0;
        if(!success)
        {
            for(int divisions:{2,4,8,16,32,64,128,256,512,1024})
            {
                std::vector<double> segmentStart=yChem_;
                double segmentTemperature=initialTemperature;
                success=true;
                for(int substep=0;substep<divisions;++substep)
                {
                    legacyCvode_->init(
                        0.0,p,rho,segmentStart,segmentTemperature);
                    double segmentEndTemperature=segmentTemperature;
                    if(!legacyCvode_->advanceTo(
                        dt/divisions,yChemNew_,segmentEndTemperature))
                    { success=false; lastFailureCode_=1; break; }
                    limitReactionIncrement(segmentStart,yChemNew_,
                                           segmentEndTemperature,segmentTemperature);
                    if(!admissibleChemState(segmentStart,yChemNew_))
                    { success=false; lastFailureCode_=validationFailureCode_; break; }
                    if(!std::isfinite(segmentEndTemperature)||segmentEndTemperature<=0.0)
                    { success=false; lastFailureCode_=3; break; }
                    segmentStart=yChemNew_;
                    segmentTemperature=segmentEndTemperature;
                }
                if(success)
                {
                    yChemNew_=segmentStart;
                    Tnew=segmentTemperature;
                    lastSubcycles_=divisions;
                    lastFailureCode_=0;
                    break;
                }
            }
        }
        if (!success) { yChemNew_ = yChem_; Tnew = initialTemperature; lastSubcycles_=-1; }

        dE = rho * Cv_eff * (Tnew - T);   // integral(qdot dt) under Cv_eff

        for (int j = 0; j < mapper_.nChem(); ++j)
            Yglobal[mapper_.globalIndex(j)] = yChemNew_[j];

        T = Tnew;
        return success;
    }
};

} // End namespace Foam

#endif
