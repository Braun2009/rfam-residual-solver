#ifndef SectorChemistryOperator_H
#define SectorChemistryOperator_H

#include "Operator.H"
#include "FiniteRateReactor.H"
#include "PhaseEquilibrium/HEMStateOwnershipContract.H"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <memory>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace Foam
{

class SectorChemistryOperator : public Operator
{
    bool explicitChemistry_ = false;
    int explicitSubSteps_ = 256;
    enum class PhaseChemistryRoute
    {
        INVALID,
        SKIP_LIQUID,
        OVERALL_SINGLE_PHASE,
        VAPOR_ONLY
    };

    struct Workspace
    {
        std::unique_ptr<FiniteRateReactor> reactor;
        std::vector<double> globalY,oldChem,increments,vaporIncrements;
    };
    ChemistryModel* chemistry_ = nullptr;
    std::shared_ptr<const thermo::EquilibriumThermoModel> boundThermoOwner_;
    std::string boundThermoModelId_;
    std::vector<std::string> boundThermoSpeciesNames_;
    EnergyConvention boundEnergyConvention_ =
        EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT;
    const phaseEquilibrium::PRPhaseRhoERecovery* phaseRecovery_ = nullptr;
    const phaseEquilibrium::PRPhaseDisplacementAuthority*
        phaseDisplacement_ = nullptr;
    ConservativeStateContract phaseStateContract_;
    std::uint64_t processRuntimeEpoch_ = 0U;
    std::string phaseOwnerPrefix_;
    mutable std::vector<Workspace> work_;
    double activationTemperature_;
    mutable long lastActive_=0,lastFailures_=0,lastRetried_=0,lastLimited_=0;
    mutable long lastTwoPhaseActive_=0,lastLiquidSkipped_=0,
                 lastPhaseReflashFailures_=0;
    mutable long lastCvodeFailures_=0,lastStateFailures_=0,lastTemperatureFailures_=0,
                 lastMassFailures_=0,lastElementFailures_=0,
                 lastPhaseAdmissionFailures_=0;
    mutable double lastHeatReleaseMinimum_=0.0,lastHeatReleaseAverage_=0.0,
                   lastHeatReleaseMaximum_=0.0,lastFailureElementError_=0.0;
    mutable int lastFailureI_=-1,lastFailureJ_=-1,lastFailureK_=-1;
    mutable int lastPhaseOwnerStatus_=-1,lastPhaseDisplacementStatus_=-1;
    mutable int lastPhaseRecoveryStatus_=-1;
    mutable std::size_t lastPhaseFailedWaypoint_=
        std::numeric_limits<std::size_t>::max();
    mutable double lastPhaseFailedFraction_=
        std::numeric_limits<double>::quiet_NaN();
    mutable double lastPhaseEnergyDisplacement_=
        std::numeric_limits<double>::quiet_NaN();
    mutable double lastPhaseCompositionDisplacement_=
        std::numeric_limits<double>::quiet_NaN();
    mutable std::uint64_t lastPhaseProbeCount_=0U,lastPhaseAcceptCount_=0U;
    mutable double lastFailureDt_=0.0,lastFailureRho_=0.0,
                   lastFailurePressure_=0.0,lastFailureTemperature_=0.0;
    mutable std::vector<double> lastFailureChemY_;

    int threadIndex() const
    {
#ifdef _OPENMP
        return omp_get_thread_num();
#else
        return 0;
#endif
    }

    static int reactorCount()
    {
        int count=1;
#ifdef _OPENMP
        count=omp_get_max_threads();
#endif
        return count;
    }

    void bindThermoOwner(
        const SpeciesRegistry& registry,
        std::shared_ptr<const thermo::EquilibriumThermoModel> owner)
    {
        if(!owner)
            throw std::invalid_argument(
                "SectorChemistryOperator: null thermo owner");
        try
        {
            boundThermoModelId_=owner->modelId();
            boundThermoSpeciesNames_=owner->speciesNames();
            boundEnergyConvention_=owner->energyConvention();
        }
        catch(...)
        {
            throw std::invalid_argument(
                "SectorChemistryOperator: unreadable thermo owner identity");
        }
        if(boundThermoModelId_.empty()
            ||boundThermoModelId_.find_first_not_of(" \t\r\n")
                ==std::string::npos
            ||boundThermoSpeciesNames_!=registry.names()
            ||boundEnergyConvention_!=EnergyConvention::
                SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT)
            throw std::invalid_argument(
                "SectorChemistryOperator: incompatible thermo owner identity");
        boundThermoOwner_=std::move(owner);
    }

    void initializeWorkspaces(
        ChemistryModel* chemistry,
        const SpeciesRegistry& registry)
    {
        work_.resize(reactorCount());
        for(auto& workspace:work_)
        {
            if(boundThermoOwner_&&phaseRecovery_&&phaseDisplacement_)
                workspace.reactor=std::make_unique<FiniteRateReactor>(
                    chemistry,registry,boundThermoOwner_,*phaseRecovery_,
                    *phaseDisplacement_);
            else if(boundThermoOwner_)
                workspace.reactor=std::make_unique<FiniteRateReactor>(
                    chemistry,registry,boundThermoOwner_);
            else
                workspace.reactor=std::make_unique<FiniteRateReactor>(
                    chemistry,registry);
            if (explicitChemistry_)
                workspace.reactor->setExplicitMode(true, explicitSubSteps_);
            workspace.globalY.resize(registry.size());
            workspace.oldChem.resize(workspace.reactor->mapper().nChem());
            workspace.increments.resize(workspace.reactor->mapper().nChem());
            workspace.vaporIncrements.resize(registry.size());
        }
    }

    bool domainUsesBoundThermoOwner(const Domain& domain) const noexcept
    {
        if(!boundThermoOwner_) return true;
        try
        {
            const auto domainOwner=domain.thermoModelOwner();
            if(!domainOwner
                ||domainOwner.get()!=boundThermoOwner_.get()
                ||boundThermoOwner_->modelId()!=boundThermoModelId_
                ||boundThermoOwner_->speciesNames()
                    !=boundThermoSpeciesNames_
                ||boundThermoOwner_->energyConvention()
                    !=boundEnergyConvention_
                ||domain.thermoModelId()!=boundThermoModelId_
                ||domain.thermoModel().speciesNames()
                    !=boundThermoSpeciesNames_
                ||domain.speciesNames()!=boundThermoSpeciesNames_
                ||domain.energyConvention()!=boundEnergyConvention_)
                return false;
            for(const Workspace& workspace:work_)
                if(!workspace.reactor
                    ||!workspace.reactor->hasInjectedThermoModel()
                    ||workspace.reactor->boundThermoModelId()
                        !=boundThermoModelId_)
                    return false;
            return true;
        }
        catch(...)
        {
            return false;
        }
    }

    static double estimateFromCommittedThermo(
        ChemistryModel* chemistry,
        const SpeciesMapper& mapper,
        double density,
        double pressure,
        double temperature,
        double internalEnergy,
        const std::vector<double>& globalMassFractions,
        double safetyFactor)
    {
        if(!chemistry
            ||!std::isfinite(density)||density<=0.0
            ||!std::isfinite(pressure)||pressure<=0.0
            ||!std::isfinite(temperature)||temperature<=0.0
            ||!std::isfinite(internalEnergy)
            ||!std::isfinite(safetyFactor)||safetyFactor<=0.0
            ||globalMassFractions.empty())
            return 0.0;
        for(double massFraction:globalMassFractions)
            if(!std::isfinite(massFraction)||massFraction<0.0)
                return 0.0;

        std::vector<double> chemistryY(
            static_cast<std::size_t>(mapper.nChem()),0.0);
        std::vector<double> omega;
        mapper.mapGlobalToChem(globalMassFractions,chemistryY);
        double heatRelease=0.0;
        try
        {
            chemistry->computeSourceTerms(
                temperature,pressure,density,
                chemistryY,omega,heatRelease);
        }
        catch(...)
        {
            return 0.0;
        }
        if(!std::isfinite(heatRelease)||omega.size()!=chemistryY.size())
            return 0.0;

        double timeScale=1.0e30;
        if(std::abs(heatRelease)>1.0)
            timeScale=std::min(
                timeScale,
                density*std::max(std::abs(internalEnergy),1.0)
                    /std::abs(heatRelease));
        for(std::size_t species=0;species<chemistryY.size();++species)
        {
            if(!std::isfinite(omega[species])) return 0.0;
            if(omega[species]<0.0&&chemistryY[species]>0.0)
                timeScale=std::min(
                    timeScale,
                    density*chemistryY[species]/(-omega[species]));
        }
        return timeScale>=1.0e29
            ?1.0e30:safetyFactor*timeScale;
    }

    static PhaseChemistryRoute phaseChemistryRoute(
        const Domain& domain,
        thermo::EquilibriumRegion region)
    {
        switch(region)
        {
            case thermo::EquilibriumRegion::TWO_PHASE:
                return PhaseChemistryRoute::VAPOR_ONLY;
            case thermo::EquilibriumRegion::SINGLE_LIQUID:
                return PhaseChemistryRoute::SKIP_LIQUID;
            case thermo::EquilibriumRegion::SINGLE_VAPOR:
            case thermo::EquilibriumRegion::SUPERCRITICAL:
                return PhaseChemistryRoute::OVERALL_SINGLE_PHASE;
            case thermo::EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE:
                return domain.thermoModel()
                        .supportsCommittedLegacyGasTransport()
                    ?PhaseChemistryRoute::OVERALL_SINGLE_PHASE
                    :PhaseChemistryRoute::INVALID;
            case thermo::EquilibriumRegion::UNKNOWN:
                return PhaseChemistryRoute::INVALID;
        }
        return PhaseChemistryRoute::INVALID;
    }

    bool phaseOwnerConfigured() const noexcept
    {
        return boundThermoOwner_&&phaseRecovery_&&phaseDisplacement_
            &&processRuntimeEpoch_!=0U&&!phaseOwnerPrefix_.empty()
            &&phaseStateContract_.hasExplicitSpeciesLayout();
    }

    std::string phaseParentOwnerId(int i,int j,int k) const
    {
        return phaseOwnerPrefix_+"/cell-"+std::to_string(i)+"-"
            +std::to_string(j)+"-"+std::to_string(k);
    }

public:
    // Explicit-chemistry mode switch (diagnostic): forward the flag to the
    // per-thread FiniteRateReactor (frozen rho/T substeps, no CVODE).
    void setExplicitChemistry(bool on, int subSteps = 256)
    {
        explicitChemistry_ = on;
        explicitSubSteps_ = subSteps;
    }

    SectorChemistryOperator(
        ChemistryModel* chemistry,
        const SpeciesRegistry& registry,
        double activationTemperature=800.0)
        :chemistry_(chemistry),activationTemperature_(activationTemperature)
    {
        if(!std::isfinite(activationTemperature_)||activationTemperature_<0.0)
            throw std::invalid_argument(
                "SectorChemistryOperator: invalid activation temperature");
        initializeWorkspaces(chemistry,registry);
    }

    SectorChemistryOperator(
        ChemistryModel* chemistry,
        const SpeciesRegistry& registry,
        std::shared_ptr<const thermo::EquilibriumThermoModel> thermoOwner,
        double activationTemperature=800.0)
        :chemistry_(chemistry),activationTemperature_(activationTemperature)
    {
        if(!std::isfinite(activationTemperature_)||activationTemperature_<0.0)
            throw std::invalid_argument(
                "SectorChemistryOperator: invalid activation temperature");
        bindThermoOwner(registry,std::move(thermoOwner));
        initializeWorkspaces(chemistry,registry);
    }

    SectorChemistryOperator(
        ChemistryModel* chemistry,
        const SpeciesRegistry& registry,
        std::shared_ptr<const thermo::EquilibriumThermoModel> thermoOwner,
        const phaseEquilibrium::PRPhaseRhoERecovery& phaseRecovery,
        const phaseEquilibrium::PRPhaseDisplacementAuthority&
            phaseDisplacement,
        std::uint64_t processRuntimeEpoch,
        std::string phaseOwnerPrefix="sector-chemistry",
        double activationTemperature=800.0)
        :chemistry_(chemistry),phaseRecovery_(&phaseRecovery),
         phaseDisplacement_(&phaseDisplacement),
         processRuntimeEpoch_(processRuntimeEpoch),
         phaseOwnerPrefix_(std::move(phaseOwnerPrefix)),
         activationTemperature_(activationTemperature)
    {
        if(!std::isfinite(activationTemperature_)||activationTemperature_<0.0
            ||processRuntimeEpoch_==0U||phaseOwnerPrefix_.empty())
            throw std::invalid_argument(
                "SectorChemistryOperator: invalid phase-owner configuration");
        bindThermoOwner(registry,std::move(thermoOwner));
        phaseStateContract_.setSpeciesLayout(
            registry.names(),static_cast<std::size_t>(registry.size()));
        initializeWorkspaces(chemistry,registry);
        if(!phaseOwnerConfigured())
            throw std::invalid_argument(
                "SectorChemistryOperator: incomplete phase-owner configuration");
    }

    bool hasBoundThermoOwner() const
    { return static_cast<bool>(boundThermoOwner_); }
    const std::string& boundThermoModelId() const
    { return boundThermoModelId_; }

    void apply(Domain& d,double dt) const override
    {
        if(!domainUsesBoundThermoOwner(d))
            throw std::runtime_error(
                "SectorChemistryOperator: Domain thermo owner mismatch");
        lastFailureI_=lastFailureJ_=lastFailureK_=-1;
        lastPhaseOwnerStatus_=lastPhaseDisplacementStatus_=-1;
        lastPhaseRecoveryStatus_=-1;
        lastPhaseFailedWaypoint_=
            std::numeric_limits<std::size_t>::max();
        lastPhaseFailedFraction_=lastPhaseEnergyDisplacement_=
            lastPhaseCompositionDisplacement_=
                std::numeric_limits<double>::quiet_NaN();
        lastPhaseProbeCount_=lastPhaseAcceptCount_=0U;
        lastFailureDt_=lastFailureRho_=lastFailurePressure_=lastFailureTemperature_=0.0;
        lastFailureChemY_.clear();
        long active=0,failures=0,retried=0,limited=0;
        long twoPhaseActive=0,liquidSkipped=0,phaseReflashFailures=0;
        long cvodeFailures=0,stateFailures=0,temperatureFailures=0,massFailures=0,
             elementFailures=0,phaseAdmissionFailures=0;
        double heatReleaseSum=0.0,heatReleaseMinimum=1e300,heatReleaseMaximum=-1e300;
        double failureElementError=0.0;
        auto recordFailure = [&](int i,int j,int k,double rho,
                                 double pressure,double temperature,
                                 const std::vector<double>& chemistryY)
        {
            #pragma omp critical(rfam_sector_chemistry_failure)
            if(lastFailureI_<0)
            {
                lastFailureI_=i; lastFailureJ_=j; lastFailureK_=k;
                lastFailureDt_=dt; lastFailureRho_=rho;
                lastFailurePressure_=pressure;
                lastFailureTemperature_=temperature;
                lastFailureChemY_=chemistryY;
            }
        };
        #pragma omp parallel for collapse(3) schedule(guided) reduction(+:active,failures,retried,limited,twoPhaseActive,liquidSkipped,phaseReflashFailures,cvodeFailures,stateFailures,temperatureFailures,massFailures,elementFailures,phaseAdmissionFailures,heatReleaseSum) reduction(min:heatReleaseMinimum) reduction(max:heatReleaseMaximum,failureElementError)
        for(int i=1;i<d.Nx()-1;++i) for(int j=1;j<d.Ny()-1;++j)
        for(int k=1;k<d.Nz()-1;++k)
        {
            Workspace& ws=work_[threadIndex()];
            const PrimitiveVars& state=d.W(i,j,k);
            const thermo::EquilibriumState& equilibrium=d.Q(i,j,k);
            const PhaseChemistryRoute phaseRoute=
                phaseChemistryRoute(d,equilibrium.region);
            if(phaseRoute==PhaseChemistryRoute::INVALID)
            {
                ++failures; ++stateFailures;
                recordFailure(
                    i,j,k,d.U(i,j,k).rho,
                    equilibrium.p,equilibrium.T,state.Y);
                continue;
            }
            if(phaseRoute==PhaseChemistryRoute::SKIP_LIQUID)
            {
                ++liquidSkipped;
                continue;
            }
            if(state.T<activationTemperature_) continue;
            const double rho=d.U(i,j,k).rho;
            const bool twoPhase=
                phaseRoute==PhaseChemistryRoute::VAPOR_ONLY;
            if(twoPhase&&boundThermoOwner_
                &&!ws.reactor->hasPhaseThermoOwner())
            {
                ++failures; ++stateFailures;
                recordFailure(
                    i,j,k,d.U(i,j,k).rho,
                    equilibrium.p,equilibrium.T,state.Y);
                continue;
            }
            double reactionRho=rho;
            double reactionEnergy=equilibrium.e;
            double sourceVolumeScale=1.0;
            ws.globalY=state.Y;
            if(twoPhase)
            {
                if(!equilibrium.vapor.has_value()
                    ||!std::isfinite(equilibrium.vaporMassFraction)
                    ||!std::isfinite(equilibrium.vaporVolumeFraction)
                    ||equilibrium.vaporMassFraction<0.0
                    ||equilibrium.vaporMassFraction>1.0
                    ||equilibrium.vaporVolumeFraction<0.0
                    ||equilibrium.vaporVolumeFraction>1.0
                    ||!std::isfinite(equilibrium.vapor->rho)
                    ||equilibrium.vapor->rho<=0.0
                    ||!std::isfinite(equilibrium.vapor->e)
                    ||equilibrium.vapor->massFractions.size()
                        !=ws.globalY.size())
                {
                    ++failures; ++stateFailures;
                    recordFailure(
                        i,j,k,rho,state.p,state.T,state.Y);
                    continue;
                }
                if(equilibrium.vaporMassFraction<=1.0e-14
                    ||equilibrium.vaporVolumeFraction<=1.0e-14)
                    continue;
                reactionRho=equilibrium.vapor->rho;
                reactionEnergy=equilibrium.vapor->e;
                sourceVolumeScale=equilibrium.vaporVolumeFraction;
                ws.globalY=equilibrium.vapor->massFractions;
                const double vaporMassDensity=
                    sourceVolumeScale*reactionRho;
                const double expectedVaporMassDensity=
                    rho*equilibrium.vaporMassFraction;
                if(!std::isfinite(vaporMassDensity)
                    ||std::abs(vaporMassDensity-expectedVaporMassDensity)
                        >1.0e-7*std::max(
                            {std::abs(vaporMassDensity),
                             std::abs(expectedVaporMassDensity),1.0}))
                {
                    ++failures; ++stateFailures;
                    recordFailure(
                        i,j,k,rho,state.p,state.T,ws.globalY);
                    continue;
                }
            }
            if(ws.reactor->estimateConservativeMaxDt(
                reactionRho,equilibrium.p,equilibrium.T,
                ws.globalY,1.0)>=1e29) continue;
            ++active;
            const SpeciesMapper& mapper=ws.reactor->mapper();
            for(int c=0;c<mapper.nChem();++c)
                ws.oldChem[c]=ws.globalY[mapper.globalIndex(c)];
            double chemistryTemperature=equilibrium.T;
            double chemistryPressure=equilibrium.p;
            double energy=0.0;
            bool chemistrySuccess=false;
            if(twoPhase&&ws.reactor->hasPhaseThermoOwner())
            {
                const std::uint64_t accepted=d.acceptedSteps();
                if(accepted==std::numeric_limits<std::uint64_t>::max())
                {
                    ++failures; ++stateFailures;
                    recordFailure(
                        i,j,k,rho,equilibrium.p,equilibrium.T,ws.oldChem);
                    continue;
                }
                CompleteCellState committedCell(d.nSpecies());
                committedCell.U=d.U(i,j,k);
                committedCell.W=d.W(i,j,k);
                committedCell.Q=equilibrium;
                const phaseEquilibrium::PRPhaseThermoParentAdmission
                    admission{
                        phaseParentOwnerId(i,j,k),processRuntimeEpoch_,
                        accepted+1U,
                        phaseEquilibrium::PRPhaseIdentity::VAPOR,
                        phaseEquilibrium::PRPhaseRootIdentity::PR_VAPOR_ROOT};
                const phaseEquilibrium::PRPhaseThermoParentView parent(
                    phaseStateContract_,committedCell,d.thermoModel(),
                    admission);
                chemistrySuccess=ws.reactor->reactPhaseConservative(
                    parent,ws.globalY,dt,energy,
                    chemistryTemperature,chemistryPressure);
            }
            else
                chemistrySuccess=ws.reactor->reactConservative(
                    reactionRho,reactionEnergy,ws.globalY,dt,energy,
                    chemistryTemperature,chemistryPressure);
            if(!chemistrySuccess)
            {
                ++failures;
                if(ws.reactor->lastFailureCode()==1) ++cvodeFailures;
                else if(ws.reactor->lastFailureCode()==2) ++stateFailures;
                else if(ws.reactor->lastFailureCode()==3) ++temperatureFailures;
                else if(ws.reactor->lastFailureCode()==4) ++massFailures;
                else if(ws.reactor->lastFailureCode()==5) ++elementFailures;
                else if(ws.reactor->lastFailureCode()==6)
                    ++phaseAdmissionFailures;
                if(twoPhase&&ws.reactor->hasPhaseThermoOwner())
                {
                    #pragma omp critical(rfam_sector_chemistry_phase_failure)
                    if(lastPhaseOwnerStatus_<0)
                    {
                        lastPhaseOwnerStatus_=static_cast<int>(
                            ws.reactor->lastPhaseOwnerStatus());
                        lastPhaseDisplacementStatus_=static_cast<int>(
                            ws.reactor->lastPhaseDisplacementStatus());
                        lastPhaseProbeCount_=
                            ws.reactor->lastPhaseProbeCount();
                        lastPhaseAcceptCount_=
                            ws.reactor->lastPhaseAcceptCount();
                        const auto& diagnostics =
                            ws.reactor->lastPhaseDisplacementDiagnostics();
                        if(diagnostics)
                        {
                            lastPhaseRecoveryStatus_=static_cast<int>(
                                diagnostics->recoveryStatus);
                            lastPhaseFailedWaypoint_=
                                diagnostics->failedWaypoint;
                            lastPhaseFailedFraction_=
                                diagnostics->failedFraction;
                            lastPhaseEnergyDisplacement_=
                                diagnostics->energyRelativeDisplacement;
                            lastPhaseCompositionDisplacement_=
                                diagnostics->compositionL1Displacement;
                        }
                    }
                }
                failureElementError=std::max(failureElementError,
                    ws.reactor->lastElementRelativeError());
                recordFailure(
                    i,j,k,rho,equilibrium.p,equilibrium.T,ws.oldChem);
                // A rejected chemistry substep is not an admissible source
                // state.  Leave Unew untouched for this cell and fail the
                // operator after the parallel loop; coupling a partially
                // updated workspace would silently turn a chemistry error
                // into a nonconservative flow update.
                continue;
            }
            else if(ws.reactor->lastSubcycles()>0) ++retried;
            if(ws.reactor->lastPositivityLimited()) ++limited;
            ConservedVars& next=d.Unew(i,j,k);
            ConservedVars candidate=next;
            double couplingScale=1.0;
            const double sourceMassDensity=
                sourceVolumeScale*reactionRho;
            for(int c=0;c<mapper.nChem();++c)
            {
                const int g=mapper.globalIndex(c);
                ws.increments[c]=sourceMassDensity
                    *(ws.globalY[g]-ws.oldChem[c]);
                if(ws.increments[c]<0.0
                    &&candidate.rhoY[g]+ws.increments[c]<0.0)
                    couplingScale=std::min(couplingScale,
                        candidate.rhoY[g]/(-ws.increments[c]));
            }
            couplingScale=std::clamp(couplingScale,0.0,1.0);
            if(couplingScale<1.0)
            { couplingScale*=1.0-1e-6; if(!ws.reactor->lastPositivityLimited()) ++limited; }

            if(twoPhase)
            {
                std::fill(
                    ws.vaporIncrements.begin(),
                    ws.vaporIncrements.end(),0.0);
                for(int c=0;c<mapper.nChem();++c)
                {
                    const int global=mapper.globalIndex(c);
                    ws.vaporIncrements[global]=couplingScale*reactionRho
                        *(ws.globalY[global]-ws.oldChem[c]);
                }
                phaseEquilibrium::HEMVaporConservativeIncrement
                    vaporIncrement;
                vaporIncrement.deltaRhoYPerVaporVolume=
                    ws.vaporIncrements;
                vaporIncrement.deltaRhoEPerVaporVolume=
                    couplingScale*energy;
                const phaseEquilibrium::HEMRhoEStateOwnershipTransaction
                    ownership(d.thermoModel());
                const phaseEquilibrium::HEMStateOwnershipResult staged=
                    ownership.stageVaporConservativeUpdate(
                        d.U(i,j,k),equilibrium,vaporIncrement);
                if(!staged.success())
                {
                    ++failures; ++stateFailures; ++phaseReflashFailures;
                    recordFailure(
                        i,j,k,rho,equilibrium.p,equilibrium.T,ws.oldChem);
                    continue;
                }
                const ConservedVars& chemistryCandidate=*staged.candidate;
                const ConservedVars& committed=d.U(i,j,k);
                for(std::size_t species=0;
                    species<candidate.rhoY.size();++species)
                    candidate.rhoY[species]+=
                        chemistryCandidate.rhoY[species]
                        -committed.rhoY[species];
                candidate.rhoE+=chemistryCandidate.rhoE-committed.rhoE;
                energy*=sourceVolumeScale*couplingScale;
                ++twoPhaseActive;
            }
            else
            {
                for(int c=0;c<mapper.nChem();++c)
                    candidate.rhoY[mapper.globalIndex(c)]+=
                        couplingScale*ws.increments[c];
                energy*=couplingScale;
                candidate.rhoE+=energy;
            }

            next=std::move(candidate);
            const double heatRelease=energy/std::max(dt,1e-30);
            heatReleaseSum+=heatRelease;
            heatReleaseMinimum=std::min(heatReleaseMinimum,heatRelease);
            heatReleaseMaximum=std::max(heatReleaseMaximum,heatRelease);
        }
        lastActive_=active; lastFailures_=failures; lastRetried_=retried; lastLimited_=limited;
        lastTwoPhaseActive_=twoPhaseActive;
        lastLiquidSkipped_=liquidSkipped;
        lastPhaseReflashFailures_=phaseReflashFailures;
        lastCvodeFailures_=cvodeFailures; lastStateFailures_=stateFailures;
        lastTemperatureFailures_=temperatureFailures;
        lastMassFailures_=massFailures; lastElementFailures_=elementFailures;
        lastPhaseAdmissionFailures_=phaseAdmissionFailures;
        lastFailureElementError_=failureElementError;
        lastHeatReleaseMinimum_=active?heatReleaseMinimum:0.0;
        lastHeatReleaseAverage_=heatReleaseSum/std::max(active,1L);
        lastHeatReleaseMaximum_=active?heatReleaseMaximum:0.0;
        if (failures > 0)
            throw std::runtime_error(
                "SectorChemistryOperator: conservative chemistry integration failed");
    }

    double maxStableDt(const Domain& d,double /*cfl*/) const override
    {
        if(!domainUsesBoundThermoOwner(d))
            throw std::runtime_error(
                "SectorChemistryOperator: Domain thermo owner mismatch");
        // CVODE mode (implicit, default): the ODE integrator controls its
        // own internal sub-steps, so the GLOBAL explicit dt must NOT be
        // limited by the conservative chemistry stability estimate
        // (observed 1.95e-14 s in the flame — 6 orders below the acoustic
        // CFL, throttling physical time by 1e6x).  Only the explicit
        // diagnostic chemistry path keeps the conservative limit
        // (2026-08-07, all-Mach time advancement).
        if (!explicitChemistry_) return 1.0e30;
        double dtMinimum=1.0e30;
        int invalid=0;
        #pragma omp parallel for collapse(3) reduction(min:dtMinimum) reduction(|:invalid) schedule(static)
        for(int i=1;i<d.Nx()-1;++i) for(int j=1;j<d.Ny()-1;++j)
        for(int k=1;k<d.Nz()-1;++k)
        {
            Workspace& ws=work_[threadIndex()];
            const PrimitiveVars& state=d.W(i,j,k);
            const thermo::EquilibriumState& equilibrium=d.Q(i,j,k);
            const PhaseChemistryRoute phaseRoute=
                phaseChemistryRoute(d,equilibrium.region);
            if(phaseRoute==PhaseChemistryRoute::INVALID)
            { invalid=1; continue; }
            if(phaseRoute==PhaseChemistryRoute::SKIP_LIQUID) continue;
            if(state.T<activationTemperature_) continue;
            double reactionRho=d.U(i,j,k).rho;
            double reactionEnergy=equilibrium.e;
            const std::vector<double>* reactionY=&state.Y;
            if(phaseRoute==PhaseChemistryRoute::VAPOR_ONLY)
            {
                if(boundThermoOwner_
                    &&!ws.reactor->hasPhaseThermoOwner())
                { invalid=1; continue; }
                if(!equilibrium.vapor.has_value()
                    ||!std::isfinite(equilibrium.vaporMassFraction)
                    ||!std::isfinite(equilibrium.vaporVolumeFraction)
                    ||equilibrium.vaporMassFraction<0.0
                    ||equilibrium.vaporMassFraction>1.0
                    ||equilibrium.vaporVolumeFraction<0.0
                    ||equilibrium.vaporVolumeFraction>1.0
                    ||!std::isfinite(equilibrium.vapor->rho)
                    ||equilibrium.vapor->rho<=0.0
                    ||equilibrium.vapor->massFractions.size()
                        !=state.Y.size())
                { invalid=1; continue; }
                if(equilibrium.vaporMassFraction<=1.0e-14
                    ||equilibrium.vaporVolumeFraction<=1.0e-14)
                    continue;
                reactionRho=equilibrium.vapor->rho;
                reactionEnergy=equilibrium.vapor->e;
                reactionY=&equilibrium.vapor->massFractions;
            }
            double dt=phaseRoute==PhaseChemistryRoute::VAPOR_ONLY
                    &&ws.reactor->hasPhaseThermoOwner()
                ?estimateFromCommittedThermo(
                    chemistry_,ws.reactor->mapper(),reactionRho,
                    equilibrium.p,equilibrium.T,reactionEnergy,
                    *reactionY,0.1)
                :ws.reactor->estimateConservativeMaxDt(
                    reactionRho,equilibrium.p,equilibrium.T,
                    *reactionY,0.1);
            if(boundThermoOwner_&&dt<=0.0
                &&!equilibrium.thermodynamics.has_value())
            {
                if(phaseRoute==PhaseChemistryRoute::VAPOR_ONLY)
                {
                    const thermo::EquilibriumResult admittedVapor=
                        boundThermoOwner_->recoverRhoE(
                            reactionRho,reactionEnergy,*reactionY);
                    if(!admittedVapor.success()
                        ||admittedVapor.diagnostics.modelId
                            !=boundThermoModelId_
                        ||!admittedVapor.state.has_value()
                        ||(admittedVapor.state->region
                                !=thermo::EquilibriumRegion::SINGLE_VAPOR
                            &&admittedVapor.state->region
                                !=thermo::EquilibriumRegion::SUPERCRITICAL)
                        ||!boundThermoOwner_->stateValidForModel(
                            *admittedVapor.state))
                    {
                        invalid=1;
                        continue;
                    }
                }
                dt=estimateFromCommittedThermo(
                    chemistry_,ws.reactor->mapper(),
                    reactionRho,equilibrium.p,equilibrium.T,
                    reactionEnergy,*reactionY,0.1);
            }
            if(!std::isfinite(dt)||dt<=0.0) invalid=1;
            else dtMinimum=std::min(dtMinimum,dt);
        }
        if(invalid)
            throw std::runtime_error(
                "SectorChemistryOperator: invalid chemistry timestep state");
        // CVODE mode chemistry dt relaxation was measured (2026-08-07):
        //  - factor 1 (conservative): t=1.9e-9/2000 steps, stable, wall 22 min
        //  - factor 100: t=1.05e-6/60 steps (+540x) but flame chemistry
        //    overshoot at 2000 steps on A406 (qchem 2.6e15, 3 failures)
        //  - factor 10: per-step cost grows non-monotonically (CVODE
        //    internal sub-steps dominated by the stiff flame region), no
        //    wall-time win; 300 steps > 25 min on A406
        // Conclusion: the flame-region chemistry time scale (~1e-14 s) is
        // a REAL physical stiffness; global-dt relaxation is not the lever.
        // The actual lever is per-cell chemistry sub-stepping with an
        // adaptive explicit sub-step count (currently fixed 256) or a
        // relaxed CVODE tolerance for the diagnostic path (pending).
        // Keep the conservative estimate as the production default.
        return dtMinimum;
    }

    long lastActiveCells() const { return lastActive_; }
    long lastFailures() const { return lastFailures_; }
    long lastRetriedCells() const { return lastRetried_; }
    long lastPositivityLimitedCells() const { return lastLimited_; }
    long lastTwoPhaseActiveCells() const { return lastTwoPhaseActive_; }
    long lastLiquidSkippedCells() const { return lastLiquidSkipped_; }
    long lastPhaseReflashFailures() const
    { return lastPhaseReflashFailures_; }
    bool hasPhaseThermoOwner() const noexcept
    { return phaseOwnerConfigured(); }
    long lastCvodeFailures() const { return lastCvodeFailures_; }
    long lastStateFailures() const { return lastStateFailures_; }
    long lastTemperatureFailures() const { return lastTemperatureFailures_; }
    long lastMassFailures() const { return lastMassFailures_; }
    long lastElementFailures() const { return lastElementFailures_; }
    long lastPhaseAdmissionFailures() const
    { return lastPhaseAdmissionFailures_; }
    int lastPhaseOwnerStatus() const { return lastPhaseOwnerStatus_; }
    int lastPhaseDisplacementStatus() const
    { return lastPhaseDisplacementStatus_; }
    int lastPhaseRecoveryStatus() const
    { return lastPhaseRecoveryStatus_; }
    std::size_t lastPhaseFailedWaypoint() const
    { return lastPhaseFailedWaypoint_; }
    double lastPhaseFailedFraction() const
    { return lastPhaseFailedFraction_; }
    double lastPhaseEnergyDisplacement() const
    { return lastPhaseEnergyDisplacement_; }
    double lastPhaseCompositionDisplacement() const
    { return lastPhaseCompositionDisplacement_; }
    std::uint64_t lastPhaseProbeCount() const
    { return lastPhaseProbeCount_; }
    std::uint64_t lastPhaseAcceptCount() const
    { return lastPhaseAcceptCount_; }
    double lastFailureElementError() const { return lastFailureElementError_; }
    int lastFailureI() const { return lastFailureI_; }
    int lastFailureJ() const { return lastFailureJ_; }
    int lastFailureK() const { return lastFailureK_; }
    double lastFailureDt() const { return lastFailureDt_; }
    double lastFailureRho() const { return lastFailureRho_; }
    double lastFailurePressure() const { return lastFailurePressure_; }
    double lastFailureTemperature() const { return lastFailureTemperature_; }
    const std::vector<double>& lastFailureChemY() const { return lastFailureChemY_; }
    double lastHeatReleaseMinimum() const { return lastHeatReleaseMinimum_; }
    double lastHeatReleaseAverage() const { return lastHeatReleaseAverage_; }
    double lastHeatReleaseMaximum() const { return lastHeatReleaseMaximum_; }
    std::string name() const override { return "Sector finite-rate chemistry"; }
};

} // End namespace Foam
#endif
