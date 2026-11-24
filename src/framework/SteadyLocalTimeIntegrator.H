#ifndef SteadyLocalTimeIntegrator_H
#define SteadyLocalTimeIntegrator_H

/*---------------------------------------------------------------------------*\
    SteadyLocalTimeIntegrator - residual-based local pseudo-time iteration.

    This integrator targets a steady fixed point R(U)=0. It deliberately does
    NOT advance SolverBase physical time. Operators assemble a conservative
    residual with a common reference dt; each active cell then scales that
    increment by its geometry/physics-aware local stable dt. The steady fixed
    point is unchanged, while intermediate pseudo-time states are not claimed
    to be time accurate.

    Convergence is measured from the undamped dimensionless residual rate.
    Candidate states still pass the complete rho-e-Y recovery contract and an
    admissibility line search. Reaching the driver iteration limit without the
    residual gate is handled by the driver as a failed steady solve.
\*---------------------------------------------------------------------------*/

#include "TimeIntegrator.H"
#include "AcousticStateAdmission.H"
#include "Domain.H"
#include "Operator.H"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

namespace Foam
{

class SteadyLocalTimeIntegrator : public TimeIntegrator
{
    double stabilityCfl_;
    double residualTolerance_;
    double maximumLocalRatio_;
    int maxBacktracks_;
    double minimumRelaxation_;
    double maximumRelaxation_;
    mutable double initialResidual_ = 0.0;
    mutable double relaxationReferenceResidual_ =
        std::numeric_limits<double>::infinity();
    mutable double relaxation_ = 1.0;
    mutable double lastResidual_ = 0.0;
    mutable double lastResidualRatio_ = 1.0;
    mutable double lastRelativeUpdate_ = 0.0;
    mutable double lastPseudoDtMinimum_ = 0.0;
    mutable double lastPseudoDtMaximum_ = 0.0;
    mutable int lastBacktracks_ = 0;
    mutable bool converged_ = false;

public:
    explicit SteadyLocalTimeIntegrator(
        double stabilityCfl = 0.3,
        double residualTolerance = 0.05,
        double maximumLocalRatio = 100.0,
        int maxBacktracks = 24,
        double minimumRelaxation = 1.0e-4,
        double maximumRelaxation = 1.0)
        : stabilityCfl_(stabilityCfl),
          residualTolerance_(residualTolerance),
          maximumLocalRatio_(maximumLocalRatio),
          maxBacktracks_(maxBacktracks),
          minimumRelaxation_(minimumRelaxation),
          maximumRelaxation_(maximumRelaxation)
    {}

    bool advancesPhysicalTime() const override { return false; }
    bool converged() const { return converged_; }
    double initialResidual() const { return initialResidual_; }
    double lastResidual() const { return lastResidual_; }
    double lastResidualRatio() const { return lastResidualRatio_; }
    double lastRelativeUpdate() const { return lastRelativeUpdate_; }
    double lastPseudoDtMinimum() const { return lastPseudoDtMinimum_; }
    double lastPseudoDtMaximum() const { return lastPseudoDtMaximum_; }
    int lastBacktracks() const { return lastBacktracks_; }
    double relaxation() const { return relaxation_; }

    bool step(
        Domain& d,
        const std::vector<std::unique_ptr<Operator>>& ops,
        double dtReference) const override
    {
        if (!std::isfinite(dtReference) || dtReference <= 0.0)
            throw std::invalid_argument(
                "steady local-time reference step must be positive");
        for (int cell=0; cell<d.N(); ++cell)
            if (!solverConsumerAdmitsAcousticState(d.Q(cell)))
            {
                d.recordFailure(Domain::FailureCode::STATE_RECOVERY_FAILED,
                    cell,name(),"solver acoustic state is not admitted");
                return false;
            }
        for (const auto& op:ops)
            if(op->ownsPersistentState())
            {
                d.recordFailure(Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                    name(),"steady local time is not defined for "
                        "operator-owned state: "+op->name());
                return false;
            }

        const int total=d.N();
        std::vector<double> inverseStableDt(
            static_cast<std::size_t>(total),0.0);
        std::vector<double> operatorDt;
        for(const auto& op:ops)
        {
            op->localStableDt(d,stabilityCfl_,operatorDt);
            if(operatorDt.size()!=static_cast<std::size_t>(total))
                throw std::runtime_error(
                    "operator returned invalid local stability field: "
                    +op->name());
            for(int cell=0;cell<total;++cell)
            {
                const double value=operatorDt[static_cast<std::size_t>(cell)];
                if(std::isfinite(value)&&value>0.0&&value<1e29)
                    inverseStableDt[static_cast<std::size_t>(cell)]+=1.0/value;
            }
        }

        std::vector<double> localScale(static_cast<std::size_t>(total),1.0);
        lastPseudoDtMinimum_=1e300;
        lastPseudoDtMaximum_=0.0;
        for(int cell=0;cell<total;++cell)
        {
            const double inverse=inverseStableDt[static_cast<std::size_t>(cell)];
            if(inverse<=0.0) continue;
            const double pseudoDt=1.0/inverse;
            localScale[static_cast<std::size_t>(cell)]=std::clamp(
                pseudoDt/dtReference,0.01,maximumLocalRatio_);
            const double acceptedPseudoDt=dtReference
                *localScale[static_cast<std::size_t>(cell)];
            lastPseudoDtMinimum_=std::min(
                lastPseudoDtMinimum_,acceptedPseudoDt);
            lastPseudoDtMaximum_=std::max(
                lastPseudoDtMaximum_,acceptedPseudoDt);
        }
        if(lastPseudoDtMaximum_==0.0)
        {
            lastPseudoDtMinimum_=dtReference;
            lastPseudoDtMaximum_=dtReference;
        }

        lastBacktracks_=0;
        detail::StepTransaction transaction(d,ops);
        if(!transaction.active()) return false;
        d.beginStep();
        if(!transaction.apply(dtReference)) return false;

        double residualSum=0.0;
        std::size_t residualComponents=0;
        double relativeUpdateMaximum=0.0;
        for(int cell=0;cell<total;++cell)
        {
            const ConservedVars& current=d.U(cell);
            ConservedVars& candidate=d.Unew(cell);
            const double rhoScale=std::max(std::abs(current.rho),1e-12);
            const double acousticScale=std::max(
                d.Q(cell).acoustic.soundSpeed(),1.0);
            const double momentumScale=rhoScale*acousticScale;
            const double energyScale=std::max(std::abs(current.rhoE),
                rhoScale*acousticScale*acousticScale);
            const double scalarScale=rhoScale;
            auto addResidual=[&](double delta,double scale)
            {
                const double value=delta
                    /(dtReference*std::max(scale,1e-300));
                residualSum+=value*value;
                ++residualComponents;
            };
            addResidual(candidate.rho-current.rho,rhoScale);
            addResidual(candidate.rhoE-current.rhoE,energyScale);
            for(std::size_t s=0;s<candidate.rhoU.size();++s)
                addResidual(candidate.rhoU[s]-current.rhoU[s],
                    std::max(std::abs(current.rhoU[s]),momentumScale));
            for(std::size_t s=0;s<candidate.rhoY.size();++s)
                addResidual(candidate.rhoY[s]-current.rhoY[s],scalarScale);
            addResidual(candidate.rhoZ-current.rhoZ,scalarScale);
            addResidual(candidate.rhoC-current.rhoC,scalarScale);

            const double scale=localScale[static_cast<std::size_t>(cell)];
            const double pseudoDt=dtReference*scale;
            const double pointImplicit=1.0/(1.0+pseudoDt
                *inverseStableDt[static_cast<std::size_t>(cell)]);
            const double updateScale=scale*pointImplicit;
            candidate.rho=current.rho
                +updateScale*(candidate.rho-current.rho);
            candidate.rhoE=current.rhoE
                +updateScale*(candidate.rhoE-current.rhoE);
            for(std::size_t s=0;s<candidate.rhoU.size();++s)
                candidate.rhoU[s]=current.rhoU[s]
                    +updateScale*(candidate.rhoU[s]-current.rhoU[s]);
            for(std::size_t s=0;s<candidate.rhoY.size();++s)
                candidate.rhoY[s]=current.rhoY[s]
                    +updateScale*(candidate.rhoY[s]-current.rhoY[s]);
            candidate.rhoZ=current.rhoZ
                +updateScale*(candidate.rhoZ-current.rhoZ);
            candidate.rhoC=current.rhoC
                +updateScale*(candidate.rhoC-current.rhoC);

            auto update=[&](double delta,double stateScale)
            {
                relativeUpdateMaximum=std::max(relativeUpdateMaximum,
                    std::abs(delta)/std::max(stateScale,1e-300));
            };
            update(candidate.rho-current.rho,rhoScale);
            update(candidate.rhoE-current.rhoE,energyScale);
            for(std::size_t s=0;s<candidate.rhoU.size();++s)
                update(candidate.rhoU[s]-current.rhoU[s],momentumScale);
            for(std::size_t s=0;s<candidate.rhoY.size();++s)
                update(candidate.rhoY[s]-current.rhoY[s],scalarScale);
            update(candidate.rhoZ-current.rhoZ,scalarScale);
            update(candidate.rhoC-current.rhoC,scalarScale);
        }

        lastResidual_=std::sqrt(residualSum/std::max<std::size_t>(
            residualComponents,1U));
        if(initialResidual_<=0.0)
            initialResidual_=std::max(lastResidual_,1e-300);
        lastResidualRatio_=lastResidual_/initialResidual_;
        if(!std::isfinite(relaxationReferenceResidual_))
            relaxationReferenceResidual_=lastResidual_;
        else if(lastResidual_>1.05*relaxationReferenceResidual_)
        {
            relaxation_=std::max(
                minimumRelaxation_,0.5*relaxation_);
            relaxationReferenceResidual_=lastResidual_;
        }
        else if(lastResidual_<0.95*relaxationReferenceResidual_)
        {
            relaxation_=std::min(
                maximumRelaxation_,1.02*relaxation_);
            relaxationReferenceResidual_=lastResidual_;
        }
        lastPseudoDtMinimum_*=relaxation_;
        lastPseudoDtMaximum_*=relaxation_;
        if(relaxation_!=1.0)
        {
            for(int cell=0;cell<total;++cell)
            {
                const ConservedVars& current=d.U(cell);
                ConservedVars& candidate=d.Unew(cell);
                candidate.rho=current.rho+relaxation_
                    *(candidate.rho-current.rho);
                candidate.rhoE=current.rhoE+relaxation_
                    *(candidate.rhoE-current.rhoE);
                for(std::size_t s=0;s<candidate.rhoU.size();++s)
                    candidate.rhoU[s]=current.rhoU[s]+relaxation_
                        *(candidate.rhoU[s]-current.rhoU[s]);
                for(std::size_t s=0;s<candidate.rhoY.size();++s)
                    candidate.rhoY[s]=current.rhoY[s]+relaxation_
                        *(candidate.rhoY[s]-current.rhoY[s]);
                candidate.rhoZ=current.rhoZ+relaxation_
                    *(candidate.rhoZ-current.rhoZ);
                candidate.rhoC=current.rhoC+relaxation_
                    *(candidate.rhoC-current.rhoC);
            }
            relativeUpdateMaximum*=relaxation_;
        }
        lastRelativeUpdate_=relativeUpdateMaximum;
        converged_=lastResidual_<residualTolerance_*initialResidual_
            ||relativeUpdateMaximum
                <=64.0*std::numeric_limits<double>::epsilon();
        if(converged_)
        {
            d.beginStep();
            transaction.accept();
            return true;
        }
        bool committed=d.commit();
        for(int backtrack=0;
            !committed&&backtrack<maxBacktracks_;++backtrack)
        {
            d.resetTrip();
            for(int cell=0;cell<total;++cell)
            {
                const ConservedVars& current=d.U(cell);
                ConservedVars& candidate=d.Unew(cell);
                candidate.rho=current.rho
                    +0.5*(candidate.rho-current.rho);
                candidate.rhoE=current.rhoE
                    +0.5*(candidate.rhoE-current.rhoE);
                for(std::size_t s=0;s<candidate.rhoU.size();++s)
                    candidate.rhoU[s]=current.rhoU[s]
                        +0.5*(candidate.rhoU[s]-current.rhoU[s]);
                for(std::size_t s=0;s<candidate.rhoY.size();++s)
                    candidate.rhoY[s]=current.rhoY[s]
                        +0.5*(candidate.rhoY[s]-current.rhoY[s]);
                candidate.rhoZ=current.rhoZ
                    +0.5*(candidate.rhoZ-current.rhoZ);
                candidate.rhoC=current.rhoC
                    +0.5*(candidate.rhoC-current.rhoC);
            }
            ++lastBacktracks_;
            committed=d.commit();
        }
        if(!committed) return false;
        transaction.accept();
        return true;
    }

    std::string name() const override
    {
        return "SteadyLocalTime (residual pseudo-time)";
    }
};

} // End namespace Foam
#endif
