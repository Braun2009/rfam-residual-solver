#ifndef SectorFrozenSteadyResidualContext_H
#define SectorFrozenSteadyResidualContext_H

/*---------------------------------------------------------------------------*\
    SectorFrozenSteadyResidualContext - production-ordered frozen LES residual.

    Owns the exact operator ordering and transactional boundary refresh used by
    the frozen-chemistry sector FAS tools: AUSM convection, WALE, species and
    thermal diffusion, plus injection/axis/wall/periodic/steady-outlet ghosts.
\*---------------------------------------------------------------------------*/

#include "FluxScheme.H"
#include "SectorConvectionOperator.H"
#include "SectorDiffusionOperator.H"
#include "SectorInjectionBoundary.H"
#include "SectorNonReflectingOutlet.H"
#include "SectorTurbulenceOperator.H"
#include "SectorWallBoundary.H"
#include "SteadyResidualEvaluator.H"
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Foam
{

class SectorFrozenSteadyResidualContext
{
    Domain& domain_;
    const SectorMeshMetrics& mesh_;
    SectorInjectionBoundary injection_;
    SectorNonReflectingOutlet outlet_;
    SectorWallBoundary wall_;
    std::vector<std::unique_ptr<Operator>> operators_;
    double boundaryRamp_=1.0;

    static std::vector<double> background(const SpeciesRegistry& registry)
    {
        std::vector<double> result(
            static_cast<std::size_t>(registry.size()),0.0);
        result[static_cast<std::size_t>(registry.index("CO2"))]=0.6468;
        result[static_cast<std::size_t>(registry.index("H2O"))]=0.3532;
        return result;
    }

    bool refreshBoundary(Domain& trial)
    {
        try
        {
            injection_.applyInStepTransaction(trial,boundaryRamp_);
            std::vector<Domain::CompleteCellReplacement> axis;
            axis.reserve(static_cast<std::size_t>(
                trial.Nx()*trial.Nz()));
            for(int i=0;i<trial.Nx();++i)
            for(int k=0;k<trial.Nz();++k)
            {
                PrimitiveVars primitive=trial.W(i,1,k);
                primitive.U[1]=0.0;
                axis.push_back(Domain::CompleteCellReplacement{
                    i,0,k,primitive,trial.Q(i,1,k)});
            }
            if(!trial.replaceCompleteCellsInStepTransaction(axis))
                return false;
            wall_.applyInStepTransaction(trial);
            trial.applyCyclicZInStepTransaction(120.0);
            if(trial.tripped()) return false;
            outlet_.applySteadyInStepTransaction(trial,0.0);
            trial.applyCyclicZInStepTransaction(120.0);
            return !trial.tripped();
        }
        catch(const std::exception& error)
        {
            if(!trial.tripped())
                trial.recordFailure(
                    Domain::FailureCode::OPERATOR_EXCEPTION,-1,
                    "frozen steady residual boundary",error.what());
            return false;
        }
    }

public:
    SectorFrozenSteadyResidualContext(
        Domain& domain,
        const SectorMeshMetrics& mesh,
        SpeciesRegistry& registry,
        const GERGResidualTable* gerg,
        const SectorTurbulenceOperator* checkpointTurbulence,
        double physicalTime)
        :domain_(domain),mesh_(mesh),
         injection_(mesh,4.0e6,registry.names(),background(registry)),
         outlet_(mesh,101325.0),wall_(700.0)
    {
        if(domain_.speciesNames()!=registry.names())
            throw std::runtime_error(
                "frozen steady residual requires the 9-species registry");
        if(!std::isfinite(physicalTime)||physicalTime<0.0)
            throw std::runtime_error(
                "frozen steady residual physical time is invalid");
        injection_.setGERGTable(gerg);
        boundaryRamp_=0.1+0.9*std::min(physicalTime/2.0e-4,1.0);

        operators_.push_back(std::make_unique<SectorConvectionOperator>(
            std::make_unique<AUSMFluxAdapter>(domain_.speciesNames()),mesh_));
        auto turbulence=std::make_unique<SectorTurbulenceOperator>(
            mesh_,registry,0.05,0.001,0.25,true,
            SectorTurbulenceOperator::TurbulenceMode::WALE);
        if(checkpointTurbulence)
            turbulence->restoreFields(
                checkpointTurbulence->kField(),
                checkpointTurbulence->omegaField());
        SectorTurbulenceOperator* turbulencePointer=turbulence.get();
        operators_.push_back(std::move(turbulence));
        operators_.push_back(std::make_unique<SectorDiffusionOperator>(
            mesh_,registry,turbulencePointer,0.3,true));
    }

    bool evaluate(
        double referenceStep,
        SteadyResidualEvaluator::Result& result)
    {
        auto boundary=[&](Domain& trial)
        { return refreshBoundary(trial); };
        return SteadyResidualEvaluator::evaluate(
            domain_,operators_,referenceStep,result,boundary,true);
    }

    bool structuredImplicitStencil(StructuredImplicitStencil& result)
    {
        result=StructuredImplicitStencil{};
        if(!domain_.stateReady()||domain_.tripped()
            ||domain_.stepTransactionActive())
        {
            if(!domain_.tripped())
                domain_.recordFailure(
                    Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                    "frozen steady implicit stencil",
                    "stencil assembly requires a ready Domain outside a transaction");
            return false;
        }
        try
        {
            detail::StepTransaction transaction(domain_,operators_,false);
            if(!transaction.active()) return false;
            if(!refreshBoundary(domain_)||domain_.tripped()) return false;
            result.ensureSize(static_cast<std::size_t>(domain_.N()));
            for(const auto& operation:operators_)
                operation->addStructuredImplicitStencil(domain_,result);
            // transaction destructor restores Domain/operator state without
            // adding a rejected physical-step ledger entry.
        }
        catch(const std::exception& error)
        {
            if(!domain_.tripped())
                domain_.recordFailure(
                    Domain::FailureCode::OPERATOR_EXCEPTION,-1,
                    "frozen steady implicit stencil",error.what());
            return false;
        }
        return !domain_.tripped();
    }

    static double normalizedRms(
        const SteadyResidualEvaluator::Result& residual,
        const Domain& scaleDomain)
    {
        if(residual.field.size()
                !=static_cast<std::size_t>(scaleDomain.N())
            ||residual.active.size()!=residual.field.size())
            throw std::runtime_error(
                "steady residual RMS field dimensions do not match");
        double sum=0.0;
        std::size_t count=0;
        for(int cell=0;cell<scaleDomain.N();++cell)
        {
            if(residual.active[static_cast<std::size_t>(cell)]==0U)
                continue;
            const ConservedVars& state=scaleDomain.U(cell);
            const ConservedVars& value=
                residual.field[static_cast<std::size_t>(cell)];
            const double rhoScale=std::max(
                std::abs(state.rho),1.0e-12);
            const double acousticScale=std::max(
                scaleDomain.Q(cell).acoustic.soundSpeed(),1.0);
            const double momentumScale=rhoScale*acousticScale;
            const double energyScale=std::max(
                std::abs(state.rhoE),
                rhoScale*acousticScale*acousticScale);
            auto add=[&](double rate,double scale)
            {
                const double normalized=rate/std::max(scale,1.0e-300);
                sum+=normalized*normalized;
                ++count;
            };
            add(value.rho,rhoScale);
            add(value.rhoE,energyScale);
            for(double component:value.rhoU)
                add(component,momentumScale);
            for(double component:value.rhoY)
                add(component,rhoScale);
            add(value.rhoZ,rhoScale);
            add(value.rhoC,rhoScale);
        }
        return std::sqrt(sum/std::max<std::size_t>(count,1U));
    }
};

} // End namespace Foam

#endif
