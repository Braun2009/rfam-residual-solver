#ifndef DualTimeIntegrator_H
#define DualTimeIntegrator_H

/*---------------------------------------------------------------------------*\
    DualTimeIntegrator - Jameson dual-time stepping:
    physical large step + pseudo-time point-implicit inner iteration.

    Physical discretisation (first-order backward Euler):
        (U^{n+1} - U^n) / dtPhys = R(U^{n+1})

    Pseudo-time inner iteration m (point-implicit damping):
        U^{m+1} = U^m + dtInner * [ R(U^m) - (U^m - U^n)/dtPhys ]
                        / (1 + dtInner * lambda_i)

    lambda_i is the cell-local spectral radius ((|u|+a) over grid spacing,
    plus the physical-source spectral radius 1/dtPhys so the fixed-point
    iteration stays stable even at dtMult=1).
    dtAcRef = dtPhys(requested) / dtMult
    dtInner = innerCfl * dtAcRef  (pseudo-time step; held fixed when a physical
    step is retried smaller). The point-implicit denominator includes at least
    1/dtAcRef, so pseudo-CFL > 1 is damped against the same global stability
    scale used by the caller to construct dtPhys.

    Convergence (total-residual drop below resTol relative to first inner
    iteration, or a component-scaled update at the floating-point floor)
    yields U^{n+1} = U^{M}. Reaching maxInnerIters without either gate is a
    rejected step.

    Adaptive stepping: if a physical step fails (inner-iteration divergence,
    non-convergence, or committed-state recovery rejection), the domain
    transaction rolls back to
    U^n and the physical step is retried halved, up to maxStepRetries times.
    This makes aggressive dtMult (40+) safe in injection-near-field regions.

    Purpose: long-time marching / steady state at engineering standard
    (same structure as STAR-CCM+/Fluent pseudo-transient steady solvers).
    The acoustic stability restriction is absorbed by the inner iteration, so
    the physical step can be enlarged dtMult-fold. Rigid split operators such
    as CVODE chemistry must be excluded from the inner chain by the caller and
    applied once per physical step afterwards.

    Constraints:
    - The inner operator chain must be incremental and stateless
      (convection / diffusion / algebraic LES). Any operator with
      ownsPersistentState()==true (e.g. SST k/omega) is rejected.
    - dtMult=1 and innerCfl=1 degenerates to explicit Euler (the inner
      residual vanishes after one iteration and damping has no effect).
\*---------------------------------------------------------------------------*/

#include "TimeIntegrator.H"
#include "AcousticStateAdmission.H"
#include "Domain.H"
#include "Operator.H"
#include "TurkelWeissPreconditioner.H"
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>
#include <limits>
#include <functional>

namespace Foam
{

class DualTimeIntegrator : public TimeIntegrator
{
private:
    double dtMult_;         // physical step enlargement: dtPhys = dtMult * dt_ac
    double innerCfl_;       // pseudo-time CFL enlargement vs acoustic limit
    int    maxInnerIters_;  // inner iteration cap
    double resTol_;         // relative residual convergence tolerance
    int    maxStepRetries_; // physical-step halving retries on state failure
    int    maxInnerBacktracks_; // admissibility line-search halvings
    bool   useLocalPseudoTime_;
    double stabilityCfl_;
    double maxLocalPseudoRatio_;
    bool   useAitkenRelaxation_;
    double minimumAitkenRelaxation_;
    double maximumAitkenRelaxation_;
    bool   useAndersonAcceleration_;
    int    andersonDepth_;
    double andersonRegularization_;
    double andersonMaximumUpdateRatio_;
    std::function<bool(Domain&)> innerIterationBoundaryHook_;
    bool excludeStructuredHaloFromResidual_;
    
    // Low-Mach preconditioning
    TurkelWeissPreconditioner preconditioner_;
    bool useLowMachPreconditioning_;
    double lowMachCutoff_ = 0.05;  ///< minimum reference Mach for preconditioning
    mutable std::vector<ConservedVars> un_;
    mutable std::vector<unsigned char> activeCellMask_;  // U^n (physical step start)
    mutable std::vector<double> localPseudoScale_;
    mutable std::vector<double> localSpectralBound_;
    mutable std::vector<double> fixedPointReferenceScale_;
    mutable std::vector<double> previousFixedPointUpdate_;
    mutable std::vector<double> currentFixedPointUpdate_;
    mutable std::vector<double> currentAndersonState_;
    mutable std::vector<double> currentAndersonResidual_;
    mutable std::vector<std::vector<double>> andersonStateHistory_;
    mutable std::vector<std::vector<double>> andersonResidualHistory_;
    mutable int lastInnerIters_ = 0;         // diagnostics
    mutable double lastResidualRatio_ = 1.0;
    mutable double lastInitialResidual_ = 0.0;
    mutable double lastResidual_ = 0.0;
    mutable double lastRelativeUpdate_ = 1.0;
    mutable int lastStepRetries_ = 0;
    mutable int lastInnerBacktracks_ = 0;
    mutable double lastAcceptedDt_ = 0.0;
    mutable double lastPseudoDtMinimum_ = 0.0;
    mutable double lastPseudoDtMaximum_ = 0.0;
    mutable double lastAitkenRelaxation_ = 1.0;
    mutable int lastAndersonDepthUsed_ = 0;

    static bool solveDenseSystem(
        std::vector<std::vector<double>>& matrix,
        std::vector<double>& rightHandSide)
    {
        const int n=static_cast<int>(rightHandSide.size());
        for(int column=0;column<n;++column)
        {
            int pivot=column;
            for(int row=column+1;row<n;++row)
                if(std::abs(matrix[row][column])
                    >std::abs(matrix[pivot][column])) pivot=row;
            if(std::abs(matrix[pivot][column])<1e-24) return false;
            if(pivot!=column)
            {
                std::swap(matrix[pivot],matrix[column]);
                std::swap(rightHandSide[pivot],rightHandSide[column]);
            }
            const double diagonal=matrix[column][column];
            for(int entry=column;entry<n;++entry)
                matrix[column][entry]/=diagonal;
            rightHandSide[column]/=diagonal;
            for(int row=0;row<n;++row)
            {
                if(row==column) continue;
                const double factor=matrix[row][column];
                for(int entry=column;entry<n;++entry)
                    matrix[row][entry]-=factor*matrix[column][entry];
                rightHandSide[row]-=factor*rightHandSide[column];
            }
        }
        return true;
    }

public:
    int lastInnerIters() const { return lastInnerIters_; }
    double lastResidualRatio() const { return lastResidualRatio_; }
    double lastInitialResidual() const { return lastInitialResidual_; }
    double lastResidual() const { return lastResidual_; }
    double lastRelativeUpdate() const { return lastRelativeUpdate_; }
    int lastStepRetries() const { return lastStepRetries_; }
    int lastInnerBacktracks() const { return lastInnerBacktracks_; }
    double lastAcceptedDt() const { return lastAcceptedDt_; }
    double lastPseudoDtMinimum() const { return lastPseudoDtMinimum_; }
    double lastPseudoDtMaximum() const { return lastPseudoDtMaximum_; }
    double lastAitkenRelaxation() const { return lastAitkenRelaxation_; }
    int lastAndersonDepthUsed() const { return lastAndersonDepthUsed_; }
    double acceptedDt(double requestedDt) const override
    {
        (void)requestedDt;
        return lastAcceptedDt_;
    }

    explicit DualTimeIntegrator(
        double dtMult = 20.0,
        double innerCfl = 5.0,
        int maxInnerIters = 80,
        double resTol = 0.05,
        int maxStepRetries = 4,
        int maxInnerBacktracks = 12,
        bool useLocalPseudoTime = false,
        double stabilityCfl = 0.3,
        double maxLocalPseudoRatio = 100.0,
        bool useAitkenRelaxation = false,
        double minimumAitkenRelaxation = 0.05,
        double maximumAitkenRelaxation = 2.0,
        bool useAndersonAcceleration = false,
        int andersonDepth = 4,
        double andersonRegularization = 1.0e-10,
        double andersonMaximumUpdateRatio = 2.0,
        std::function<bool(Domain&)> innerIterationBoundaryHook = {},
        bool excludeStructuredHaloFromResidual = false)
        : dtMult_(dtMult), innerCfl_(innerCfl),
          maxInnerIters_(maxInnerIters), resTol_(resTol),
          maxStepRetries_(maxStepRetries),
          maxInnerBacktracks_(maxInnerBacktracks),
          useLocalPseudoTime_(useLocalPseudoTime),
          stabilityCfl_(stabilityCfl),
          maxLocalPseudoRatio_(maxLocalPseudoRatio),
          useAitkenRelaxation_(useAitkenRelaxation),
          minimumAitkenRelaxation_(minimumAitkenRelaxation),
          maximumAitkenRelaxation_(maximumAitkenRelaxation),
          useAndersonAcceleration_(useAndersonAcceleration),
          andersonDepth_(std::max(1,andersonDepth)),
          andersonRegularization_(andersonRegularization),
          andersonMaximumUpdateRatio_(andersonMaximumUpdateRatio),
          innerIterationBoundaryHook_(
              std::move(innerIterationBoundaryHook)),
          excludeStructuredHaloFromResidual_(
              excludeStructuredHaloFromResidual),
          preconditioner_(),
          useLowMachPreconditioning_(false) {}

    bool step(
        Domain& d,
        const std::vector<std::unique_ptr<Operator>>& ops,
        double dtPhys
    ) const override
    {
        for (int cell = 0; cell < d.N(); ++cell)
            if (!solverConsumerAdmitsAcousticState(d.Q(cell)))
            {
                d.recordFailure(Domain::FailureCode::STATE_RECOVERY_FAILED,
                    cell,name(),
                    "solver acoustic state is not admitted");
                return false;
            }
        if(useAitkenRelaxation_&&useAndersonAcceleration_)
        {
            d.recordFailure(Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                name(),"Aitken and Anderson acceleration are mutually exclusive");
            return false;
        }
        for (const auto& op : ops)
            if (op->ownsPersistentState())
            {
                d.recordFailure(Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                    name(),
                    "dual-time inner iteration is not defined for "
                    "operator-owned state: " + op->name());
                return false;
            }

        // Snapshot U^n (physical step start) and prepare the optional local
        // pseudo-time preconditioner before opening the transaction.
        const int total = d.N();
        if ((int)un_.size() != total) un_.resize(total);
        for (int k = 0; k < total; ++k) un_[k] = d.U(k);

        activeCellMask_.assign(static_cast<std::size_t>(total),1U);
        if(excludeStructuredHaloFromResidual_)
        {
            const int nx=d.Nx(),ny=d.Ny(),nz=d.Nz();
            for(int linear=0;linear<total;++linear)
            {
                const int i=linear%nx;
                const int j=(linear/nx)%ny;
                const int k=linear/(nx*ny);
                const bool active=i>0&&i<nx-1
                    &&(ny==1||(j>0&&j<ny-1))
                    &&(nz==1||(k>0&&k<nz-1));
                activeCellMask_[static_cast<std::size_t>(linear)]=
                    active?1U:0U;
            }
        }

        const double dtAcousticReference = dtPhys / dtMult_;
        double dtInner = innerCfl_ * dtAcousticReference;
        
        // Apply low-Mach preconditioning: scale pseudo-timestep by β
        // This removes the acoustic stiffness at low Mach numbers
        if (useLowMachPreconditioning_) {
            // Compute domain-average Mach number for initial β estimate
            double avgMach = 0.0;
            int validCells = 0;
            for (int k = 0; k < total; ++k) {
                const auto& w = d.W(k);
                const double a = d.Q(k).acoustic.soundSpeed();
                if (a > 1.0) {
                    double velMag = 0.0;
                    for (const auto& v : w.U) velMag += v * v;
                    avgMach += std::sqrt(velMag) / a;
                    ++validCells;
                }
            }
            if (validCells > 0) {
                avgMach /= validCells;
                const double beta = std::max(avgMach, lowMachCutoff_);
                dtInner *= beta;  // Scale pseudo-time by β
            }
        }
        const double operatorBound = static_cast<double>(
            std::max<std::size_t>(ops.size(),1U));
        localPseudoScale_.assign(static_cast<std::size_t>(total),1.0);
        localSpectralBound_.assign(static_cast<std::size_t>(total),
            operatorBound/dtAcousticReference);
        lastPseudoDtMinimum_=dtInner;
        lastPseudoDtMaximum_=dtInner;
        if(useLocalPseudoTime_)
        {
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
                for(int k=0;k<total;++k)
                {
                    const double value=operatorDt[static_cast<std::size_t>(k)];
                    if(std::isfinite(value)&&value>0.0&&value<1e29)
                        inverseStableDt[static_cast<std::size_t>(k)]+=1.0/value;
                }
            }
            lastPseudoDtMinimum_=1e300;
            lastPseudoDtMaximum_=0.0;
            for(int k=0;k<total;++k)
            {
                const double inverse=inverseStableDt[static_cast<std::size_t>(k)];
                if(inverse<=0.0) continue;
                const double combinedDt=1.0/inverse;
                const double ratio=std::clamp(
                    combinedDt/dtAcousticReference,0.01,maxLocalPseudoRatio_);
                localPseudoScale_[static_cast<std::size_t>(k)]=ratio;
                localSpectralBound_[static_cast<std::size_t>(k)]=inverse;
                const double pseudoDt=dtInner*ratio;
                lastPseudoDtMinimum_=std::min(lastPseudoDtMinimum_,pseudoDt);
                lastPseudoDtMaximum_=std::max(lastPseudoDtMaximum_,pseudoDt);
            }
            if(lastPseudoDtMaximum_==0.0)
            {
                lastPseudoDtMinimum_=dtInner;
                lastPseudoDtMaximum_=dtInner;
            }
        }
        const std::size_t componentsPerCell = static_cast<std::size_t>(
            7+d.nSpecies());
        fixedPointReferenceScale_.assign(
            static_cast<std::size_t>(total)*componentsPerCell,1.0);
        for(int k=0;k<total;++k)
        {
            const ConservedVars& u=d.U(k);
            const double rhoScale=std::max(std::abs(u.rho),1e-12);
            const double acousticScale=std::max(
                d.Q(k).acoustic.soundSpeed(),1.0);
            const double momentumScale=rhoScale*acousticScale;
            const double energyScale=std::max(std::abs(u.rhoE),
                rhoScale*acousticScale*acousticScale);
            std::size_t component=static_cast<std::size_t>(k)
                *componentsPerCell;
            fixedPointReferenceScale_[component++]=rhoScale;
            fixedPointReferenceScale_[component++]=energyScale;
            for(std::size_t s=0;s<u.rhoU.size();++s)
                fixedPointReferenceScale_[component++]=momentumScale;
            for(std::size_t s=0;s<u.rhoY.size();++s)
                fixedPointReferenceScale_[component++]=rhoScale;
            fixedPointReferenceScale_[component++]=rhoScale;
            fixedPointReferenceScale_[component++]=rhoScale;
        }
        if (!d.beginStepTransaction()) return false;

        lastStepRetries_ = 0;
        lastInnerBacktracks_ = 0;
        lastAcceptedDt_ = 0.0;
        lastAndersonDepthUsed_=0;
        double dtPhysTry = dtPhys;
        for (int attempt = 0; attempt <= maxStepRetries_; ++attempt)
        {
            if (attempt > 0)
            {
                // Roll back to U^n, clear the trip raised by the failed
                // state recovery, and retry with a halved physical step.
                dtPhysTry *= 0.5;
                d.rollbackStepTransaction();
                d.resetTrip();
                if (!d.beginStepTransaction()) return false;
                ++lastStepRetries_;
            }
            const double physCoeff = dtInner / dtPhysTry;
            const int Nx = d.Nx(), Ny = d.Ny(), Nz = d.Nz();
            const double dx = d.dx();
            const double dy = (Ny > 1) ? d.dy() : dx;
            const double dz = (Nz > 1) ? d.dz() : dx;

            double res0 = 0.0;
            double previousResidual=std::numeric_limits<double>::infinity();
            bool failed = false;
            bool converged = false;
            lastInnerIters_ = 0;
            lastAitkenRelaxation_ = 1.0;
            previousFixedPointUpdate_.clear();
            andersonStateHistory_.clear();
            andersonResidualHistory_.clear();
            lastAndersonDepthUsed_=0;
            for (int m = 0; m < maxInnerIters_; ++m)
            {
                if(innerIterationBoundaryHook_
                    &&!innerIterationBoundaryHook_(d))
                {
                    failed=true;
                    break;
                }
                d.beginStep();   // Unew <- U^m
                for (const auto& op : ops)
                {
                    op->apply(d,dtInner);
                    if (d.tripped()) { failed = true; break; }
                }
                if (failed) break;

                // Add the backward-Euler physical source, then measure the
                // nonlinear fixed-point residual BEFORE point-implicit
                // damping. Each conserved component is scaled locally so
                // energy units or trace species cannot dominate the norm.
                double resSum = 0.0;
                std::size_t residualComponents = 0;
                double relativeUpdateMax = 0.0;
                currentFixedPointUpdate_.assign(
                    static_cast<std::size_t>(total)*componentsPerCell,0.0);
                currentAndersonState_.assign(
                    static_cast<std::size_t>(total)*componentsPerCell,0.0);
                currentAndersonResidual_.assign(
                    static_cast<std::size_t>(total)*componentsPerCell,0.0);
                for (int k = 0; k < total; ++k)
                {
                    const ConservedVars& u = d.U(k);    // U^m
                    const ConservedVars& un = un_[k];   // U^n
                    ConservedVars& unew = d.Unew(k);    // U^m + dtInner*R(U^m)
                    const bool active = activeCellMask_[
                        static_cast<std::size_t>(k)] != 0U;
                    if(!active)
                    {
                        std::size_t component=static_cast<std::size_t>(k)
                            *componentsPerCell;
                        auto storeBoundaryState=[&](double value)
                        {
                            currentAndersonState_[component]=value
                                /fixedPointReferenceScale_[component];
                            ++component;
                        };
                        storeBoundaryState(u.rho);
                        storeBoundaryState(u.rhoE);
                        for(double value:u.rhoU) storeBoundaryState(value);
                        for(double value:u.rhoY) storeBoundaryState(value);
                        storeBoundaryState(u.rhoZ);
                        storeBoundaryState(u.rhoC);
                        continue;
                    }

                    // Raw backward-Euler residual increment at the common
                    // reference dtInner. This is the convergence quantity;
                    // local pseudo-time only preconditions the update below.
                    unew.rho  -= physCoeff * (u.rho  - un.rho);
                    unew.rhoE -= physCoeff * (u.rhoE - un.rhoE);
                    for (std::size_t s = 0; s < unew.rhoU.size(); ++s)
                        unew.rhoU[s] -= physCoeff * (u.rhoU[s] - un.rhoU[s]);
                    for (std::size_t s = 0; s < unew.rhoY.size(); ++s)
                        unew.rhoY[s] -= physCoeff * (u.rhoY[s] - un.rhoY[s]);
                    unew.rhoZ -= physCoeff * (u.rhoZ - un.rhoZ);
                    unew.rhoC -= physCoeff * (u.rhoC - un.rhoC);

                    const PrimitiveVars& w = d.W(k);
                    const double a = d.Q(k).acoustic.soundSpeed();
                    const double rhoScale = std::max(std::abs(u.rho),1e-12);
                    const double acousticScale = std::max(a,1.0);
                    const double momentumScale = rhoScale*acousticScale;
                    const double energyScale = std::max(std::abs(u.rhoE),
                        rhoScale*acousticScale*acousticScale);
                    // Species/scalars use an absolute mass-fraction scale.
                    // Scaling by their local trace value would let harmless
                    // 1e-30 -> 1e-28 changes dominate the nonlinear norm.
                    const double scalarScale = rhoScale;
                    auto addResidual = [&](double delta,double scale)
                    {
                        const double value=delta/std::max(scale,1e-300);
                        resSum += value*value;
                        ++residualComponents;
                    };

                    // Raw nonlinear residual increment, before preconditioning.
                    addResidual(unew.rho-u.rho,rhoScale);
                    addResidual(unew.rhoE-u.rhoE,energyScale);
                    for (std::size_t s=0;s<unew.rhoU.size();++s)
                        addResidual(unew.rhoU[s]-u.rhoU[s],
                            std::max(std::abs(u.rhoU[s]),momentumScale));
                    for (std::size_t s=0;s<unew.rhoY.size();++s)
                        addResidual(unew.rhoY[s]-u.rhoY[s],
                            std::max(std::abs(u.rhoY[s]),scalarScale));
                    addResidual(unew.rhoZ-u.rhoZ,
                        std::max(std::abs(u.rhoZ),scalarScale));
                    addResidual(unew.rhoC-u.rhoC,
                        std::max(std::abs(u.rhoC),scalarScale));

                    // Local pseudo-time preconditioning scales the raw
                    // residual after it has been measured. It does not alter
                    // the requested physical dt or the nonlinear norm.
                    const double localScale=
                        localPseudoScale_[static_cast<std::size_t>(k)];
                    const double pseudoDt=dtInner*localScale;
                    unew.rho=u.rho+localScale*(unew.rho-u.rho);
                    unew.rhoE=u.rhoE+localScale*(unew.rhoE-u.rhoE);
                    for(std::size_t s=0;s<unew.rhoU.size();++s)
                        unew.rhoU[s]=u.rhoU[s]
                            +localScale*(unew.rhoU[s]-u.rhoU[s]);
                    for(std::size_t s=0;s<unew.rhoY.size();++s)
                        unew.rhoY[s]=u.rhoY[s]
                            +localScale*(unew.rhoY[s]-u.rhoY[s]);
                    unew.rhoZ=u.rhoZ+localScale*(unew.rhoZ-u.rhoZ);
                    unew.rhoC=u.rhoC+localScale*(unew.rhoC-u.rhoC);
                    
                    // Low-Mach preconditioning: apply Γ to the update
                    // Γ scales momentum equations by β²/a² to remove acoustic stiffness
                    if (useLowMachPreconditioning_) {
                        // Local (per-cell) preconditioner cache: step() is const,
                        // so the member preconditioner_ cannot hold per-cell state.
                        TurkelWeissPreconditioner cellPrecond;
                        cellPrecond.setMachCutoff(lowMachCutoff_);
                        cellPrecond.updateLocal(w.U, a);
                        ConservedVars update;
                        update.rho = unew.rho - u.rho;
                        update.rhoE = unew.rhoE - u.rhoE;
                        update.rhoU = {unew.rhoU[0] - u.rhoU[0],
                                       unew.rhoU[1] - u.rhoU[1],
                                       unew.rhoU[2] - u.rhoU[2]};
                        update.rhoY = unew.rhoY;
                        for (std::size_t s = 0; s < update.rhoY.size(); ++s)
                            update.rhoY[s] -= u.rhoY[s];
                        update.rhoZ = unew.rhoZ - u.rhoZ;
                        update.rhoC = unew.rhoC - u.rhoC;
                        
                        ConservedVars precondUpdate = cellPrecond.applyToPseudoTimeDerivative(update);
                        
                        unew.rho = u.rho + precondUpdate.rho;
                        unew.rhoE = u.rhoE + precondUpdate.rhoE;
                        for (std::size_t s = 0; s < unew.rhoU.size(); ++s)
                            unew.rhoU[s] = u.rhoU[s] + precondUpdate.rhoU[s];
                        for (std::size_t s = 0; s < unew.rhoY.size(); ++s)
                            unew.rhoY[s] = u.rhoY[s] + precondUpdate.rhoY[s];
                        unew.rhoZ = u.rhoZ + precondUpdate.rhoZ;
                        unew.rhoC = u.rhoC + precondUpdate.rhoC;
                    }

                    // Point-implicit damping: Unew = U + (Unew-U)/(1+dtPseudo*lambda)
                    double lambda = (std::abs(w.U[0]) + a) / dx;
                    if (Ny > 1) lambda += (std::abs(w.U[1]) + a) / dy;
                    if (Nz > 1) lambda += (std::abs(w.U[2]) + a) / dz;
                    // The caller constructs the requested physical step as
                    // dtMult times the minimum operator stability step. For an
                    // additive chain, spectral radii add; N/dtAcRef is a
                    // conservative upper bound when each operator advertises
                    // dtStable >= dtAcRef. This also protects mapped-sector
                    // coordinates where Cartesian dx/dy/dz is only a rough
                    // local estimate. Physical-source stiffness is additional.
                    lambda += localSpectralBound_[static_cast<std::size_t>(k)];
                    lambda += 1.0 / dtPhysTry;
                    const double imp = 1.0 / (1.0 + pseudoDt * lambda);
                    unew.rho  = u.rho  + imp * (unew.rho  - u.rho);
                    unew.rhoE = u.rhoE + imp * (unew.rhoE - u.rhoE);
                    for (std::size_t s = 0; s < unew.rhoU.size(); ++s)
                        unew.rhoU[s] = u.rhoU[s] + imp * (unew.rhoU[s] - u.rhoU[s]);
                    for (std::size_t s = 0; s < unew.rhoY.size(); ++s)
                        unew.rhoY[s] = u.rhoY[s] + imp * (unew.rhoY[s] - u.rhoY[s]);
                    unew.rhoZ = u.rhoZ + imp * (unew.rhoZ - u.rhoZ);
                    unew.rhoC = u.rhoC + imp * (unew.rhoC - u.rhoC);

                    std::size_t component = static_cast<std::size_t>(k)
                        *componentsPerCell;
                    auto storeUpdate = [&](double current,double candidate,
                                           double scale)
                    {
                        const double delta=candidate-current;
                        const double value=delta/std::max(scale,1e-300);
                        currentFixedPointUpdate_[component]=value;
                        const double fixedScale=
                            fixedPointReferenceScale_[component];
                        currentAndersonState_[component]=current/fixedScale;
                        currentAndersonResidual_[component]=delta/fixedScale;
                        ++component;
                        relativeUpdateMax=std::max(
                            relativeUpdateMax,std::abs(value));
                    };
                    storeUpdate(u.rho,unew.rho,rhoScale);
                    storeUpdate(u.rhoE,unew.rhoE,energyScale);
                    for (std::size_t s=0;s<unew.rhoU.size();++s)
                        storeUpdate(u.rhoU[s],unew.rhoU[s],
                            std::max(std::abs(u.rhoU[s]),momentumScale));
                    for (std::size_t s=0;s<unew.rhoY.size();++s)
                        storeUpdate(u.rhoY[s],unew.rhoY[s],
                            std::max(std::abs(u.rhoY[s]),scalarScale));
                    storeUpdate(u.rhoZ,unew.rhoZ,
                        std::max(std::abs(u.rhoZ),scalarScale));
                    storeUpdate(u.rhoC,unew.rhoC,
                        std::max(std::abs(u.rhoC),scalarScale));
                }
                const double res = std::sqrt(resSum/std::max<std::size_t>(
                    residualComponents,1U));

                double aitkenRelaxation=1.0;
                if(useAitkenRelaxation_
                    && previousFixedPointUpdate_.size()
                        ==currentFixedPointUpdate_.size())
                {
                    double numerator=0.0,denominator=0.0;
                    for(std::size_t q=0;q<currentFixedPointUpdate_.size();++q)
                    {
                        const double difference=currentFixedPointUpdate_[q]
                            -previousFixedPointUpdate_[q];
                        numerator+=previousFixedPointUpdate_[q]*difference;
                        denominator+=difference*difference;
                    }
                    if(denominator>1e-300)
                        aitkenRelaxation=std::clamp(
                            -lastAitkenRelaxation_*numerator/denominator,
                            minimumAitkenRelaxation_,
                            maximumAitkenRelaxation_);
                }
                lastAitkenRelaxation_=aitkenRelaxation;
                previousFixedPointUpdate_=currentFixedPointUpdate_;
                if (m == 0)
                {
                    res0 = std::max(res,1e-30);
                    lastInitialResidual_ = res;
                }
                lastInnerIters_ = m + 1;
                lastResidual_ = res;
                lastResidualRatio_ = res / res0;
                lastRelativeUpdate_ = relativeUpdateMax;

                // The residual belongs to the current authoritative state U^m.
                // If U^m has converged, keep it; do not commit the already
                // assembled U^{m+1} candidate and accidentally take one extra
                // pseudo-time step past the convergence gate.
                if (res < resTol_ * res0
                    || relativeUpdateMax
                        <= 64.0*std::numeric_limits<double>::epsilon())
                {
                    converged = true;
                    break;
                }
                if (m >= 4 && res > 10.0 * res0)
                {
                    failed = true;
                    break;
                }

                if(aitkenRelaxation!=1.0)
                {
                    for(int k=0;k<total;++k)
                    {
                        const ConservedVars& u=d.U(k);
                        ConservedVars& candidate=d.Unew(k);
                        candidate.rho=u.rho+aitkenRelaxation
                            *(candidate.rho-u.rho);
                        candidate.rhoE=u.rhoE+aitkenRelaxation
                            *(candidate.rhoE-u.rhoE);
                        for(std::size_t s=0;s<candidate.rhoU.size();++s)
                            candidate.rhoU[s]=u.rhoU[s]+aitkenRelaxation
                                *(candidate.rhoU[s]-u.rhoU[s]);
                        for(std::size_t s=0;s<candidate.rhoY.size();++s)
                            candidate.rhoY[s]=u.rhoY[s]+aitkenRelaxation
                                *(candidate.rhoY[s]-u.rhoY[s]);
                        candidate.rhoZ=u.rhoZ+aitkenRelaxation
                            *(candidate.rhoZ-u.rhoZ);
                        candidate.rhoC=u.rhoC+aitkenRelaxation
                            *(candidate.rhoC-u.rhoC);
                    }
                    relativeUpdateMax*=std::abs(aitkenRelaxation);
                }

                if(useAndersonAcceleration_)
                {
                    // Safeguard: a rising raw residual invalidates the local
                    // secant model. Restart the history from the current valid
                    // state instead of extrapolating through the LOX recovery
                    // boundary.
                    if(andersonResidualHistory_.size()>1
                        &&res>1.05*previousResidual)
                    {
                        andersonStateHistory_.clear();
                        andersonResidualHistory_.clear();
                        lastAndersonDepthUsed_=0;
                    }
                    andersonStateHistory_.push_back(currentAndersonState_);
                    andersonResidualHistory_.push_back(currentAndersonResidual_);
                    const std::size_t maximumHistory=static_cast<std::size_t>(
                        andersonDepth_+1);
                    if(andersonStateHistory_.size()>maximumHistory)
                    {
                        andersonStateHistory_.erase(
                            andersonStateHistory_.begin());
                        andersonResidualHistory_.erase(
                            andersonResidualHistory_.begin());
                    }
                    const int depth=std::min(andersonDepth_,
                        static_cast<int>(andersonResidualHistory_.size())-1);
                    lastAndersonDepthUsed_=0;
                    if(depth>0)
                    {
                        const int start=static_cast<int>(
                            andersonResidualHistory_.size())-depth-1;
                        std::vector<std::vector<double>> normal(
                            depth,std::vector<double>(depth,0.0));
                        std::vector<double> rhs(depth,0.0);
                        const auto& currentResidual=
                            andersonResidualHistory_.back();
                        double diagonalSum=0.0;
                        for(int aIndex=0;aIndex<depth;++aIndex)
                        for(int bIndex=0;bIndex<depth;++bIndex)
                        {
                            double dot=0.0;
                            const auto& fa0=andersonResidualHistory_[
                                static_cast<std::size_t>(start+aIndex)];
                            const auto& fa1=andersonResidualHistory_[
                                static_cast<std::size_t>(start+aIndex+1)];
                            const auto& fb0=andersonResidualHistory_[
                                static_cast<std::size_t>(start+bIndex)];
                            const auto& fb1=andersonResidualHistory_[
                                static_cast<std::size_t>(start+bIndex+1)];
                            for(std::size_t q=0;q<currentResidual.size();++q)
                                dot+=(fa1[q]-fa0[q])*(fb1[q]-fb0[q]);
                            normal[aIndex][bIndex]=dot;
                            if(aIndex==bIndex) diagonalSum+=dot;
                        }
                        for(int aIndex=0;aIndex<depth;++aIndex)
                        {
                            double dot=0.0;
                            const auto& f0=andersonResidualHistory_[
                                static_cast<std::size_t>(start+aIndex)];
                            const auto& f1=andersonResidualHistory_[
                                static_cast<std::size_t>(start+aIndex+1)];
                            for(std::size_t q=0;q<currentResidual.size();++q)
                                dot+=(f1[q]-f0[q])*currentResidual[q];
                            rhs[aIndex]=dot;
                        }
                        const double regularization=andersonRegularization_
                            *std::max(diagonalSum/std::max(depth,1),1e-30);
                        for(int aIndex=0;aIndex<depth;++aIndex)
                            normal[aIndex][aIndex]+=regularization;
                        if(solveDenseSystem(normal,rhs))
                        {
                            std::vector<double> accelerated(
                                currentAndersonState_.size(),0.0);
                            for(std::size_t q=0;q<accelerated.size();++q)
                                accelerated[q]=currentAndersonState_[q]
                                    +currentAndersonResidual_[q];
                            for(int column=0;column<depth;++column)
                            {
                                const auto& x0=andersonStateHistory_[
                                    static_cast<std::size_t>(start+column)];
                                const auto& x1=andersonStateHistory_[
                                    static_cast<std::size_t>(start+column+1)];
                                const auto& f0=andersonResidualHistory_[
                                    static_cast<std::size_t>(start+column)];
                                const auto& f1=andersonResidualHistory_[
                                    static_cast<std::size_t>(start+column+1)];
                                for(std::size_t q=0;q<accelerated.size();++q)
                                    accelerated[q]-=rhs[column]
                                        *((x1[q]+f1[q])-(x0[q]+f0[q]));
                            }
                            double acceleratedUpdateMaximum=0.0;
                            for(std::size_t q=0;q<accelerated.size();++q)
                                acceleratedUpdateMaximum=std::max(
                                    acceleratedUpdateMaximum,
                                    std::abs(accelerated[q]
                                        -currentAndersonState_[q]));
                            const double updateLimit=andersonMaximumUpdateRatio_
                                *std::max(relativeUpdateMax,1e-16);
                            const double accelerationScale=
                                acceleratedUpdateMaximum>updateLimit
                                ? updateLimit/acceleratedUpdateMaximum : 1.0;
                            if(accelerationScale<1.0)
                                for(std::size_t q=0;q<accelerated.size();++q)
                                    accelerated[q]=currentAndersonState_[q]
                                        +accelerationScale*(accelerated[q]
                                            -currentAndersonState_[q]);
                            relativeUpdateMax=0.0;
                            for(int k=0;k<total;++k)
                            {
                                if(activeCellMask_[static_cast<std::size_t>(k)]
                                    ==0U) continue;
                                const ConservedVars& u=d.U(k);
                                ConservedVars& candidate=d.Unew(k);
                                std::size_t component=static_cast<std::size_t>(k)
                                    *componentsPerCell;
                                auto acceleratedValue=[&](double current)
                                {
                                    const double value=accelerated[component]
                                        *fixedPointReferenceScale_[component];
                                    relativeUpdateMax=std::max(relativeUpdateMax,
                                        std::abs(accelerated[component]
                                            -currentAndersonState_[component]));
                                    ++component;
                                    (void)current;
                                    return value;
                                };
                                candidate.rho=acceleratedValue(u.rho);
                                candidate.rhoE=acceleratedValue(u.rhoE);
                                for(std::size_t s=0;s<candidate.rhoU.size();++s)
                                    candidate.rhoU[s]=acceleratedValue(u.rhoU[s]);
                                for(std::size_t s=0;s<candidate.rhoY.size();++s)
                                    candidate.rhoY[s]=acceleratedValue(u.rhoY[s]);
                                candidate.rhoZ=acceleratedValue(u.rhoZ);
                                candidate.rhoC=acceleratedValue(u.rhoC);
                            }
                            lastAndersonDepthUsed_=depth;
                        }
                    }
                }

                // Admissibility line search: a nonlinear rho-e-Y recovery
                // failure rejects this pseudo-time update, not immediately the
                // requested physical step. Uniformly halve U^{m+1}-U^m until
                // the complete state can be recovered or the line search is
                // exhausted. The current U^m remains authoritative until a
                // commit succeeds, so every trial starts from the same valid
                // state and preserves conservation of the scaled increment.
                double acceptedScale = 1.0;
                bool committed = d.commit();
                for (int backtrack = 0;
                     !committed && backtrack < maxInnerBacktracks_;
                     ++backtrack)
                {
                    d.resetTrip();
                    acceptedScale *= 0.5;
                    for (int k = 0; k < total; ++k)
                    {
                        const ConservedVars& u = d.U(k);
                        ConservedVars& candidate = d.Unew(k);
                        candidate.rho = u.rho
                            + 0.5*(candidate.rho-u.rho);
                        candidate.rhoE = u.rhoE
                            + 0.5*(candidate.rhoE-u.rhoE);
                        for (std::size_t s=0;s<candidate.rhoU.size();++s)
                            candidate.rhoU[s] = u.rhoU[s]
                                + 0.5*(candidate.rhoU[s]-u.rhoU[s]);
                        for (std::size_t s=0;s<candidate.rhoY.size();++s)
                            candidate.rhoY[s] = u.rhoY[s]
                                + 0.5*(candidate.rhoY[s]-u.rhoY[s]);
                        candidate.rhoZ = u.rhoZ
                            + 0.5*(candidate.rhoZ-u.rhoZ);
                        candidate.rhoC = u.rhoC
                            + 0.5*(candidate.rhoC-u.rhoC);
                    }
                    ++lastInnerBacktracks_;
                    committed = d.commit();
                }
                if (!committed) { failed = true; break; }
                // Backtracking scales the accepted update, not the nonlinear
                // residual at U^m. Scaling the residual here would allow an
                // inadmissible fixed point to claim convergence merely by
                // taking a tiny line-search step.
                relativeUpdateMax *= acceptedScale;
                lastRelativeUpdate_ = relativeUpdateMax;
                previousResidual=res;
            }
            // Reaching maxInnerIters without the residual gate is a rejected
            // physical step, not a successful large-step advance.
            if (!failed && !converged) failed = true;
            if (!failed)
            {
                d.acceptStepTransaction();
                lastAcceptedDt_ = dtPhysTry;
                return true;
            }
        }
        // All retries exhausted: restore U^n and report failure.
        d.rollbackStepTransaction();
        d.recordFailure(Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
            name(),"dual-time step failed after "
                + std::to_string(maxStepRetries_) + " halving retries");
        return false;
    }

    std::string name() const override
    {
        return "DualTime (physical " + std::to_string(dtMult_)
            + "x, pseudo-CFL " + std::to_string(innerCfl_) + ")";
    }

    // Low-Mach preconditioning configuration
    void setLowMachPreconditioningEnabled(bool enable) {
        useLowMachPreconditioning_ = enable;
        preconditioner_.setEnabled(enable);
    }

    bool isLowMachPreconditioningEnabled() const {
        return useLowMachPreconditioning_;
    }

    void setLowMachCutoff(double cutoff) {
        lowMachCutoff_ = cutoff;
        preconditioner_.setMachCutoff(cutoff);
    }

    double getLowMachCutoff() const {
        return lowMachCutoff_;
    }
};

} // End namespace Foam

#endif
