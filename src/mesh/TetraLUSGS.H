#ifndef TetraLUSGS_H
#define TetraLUSGS_H

#include "TetraFiniteVolume.H"
#include "TetraStateRecovery.H"
#include "TetraFaceInvariant.H"
#include "FrozenEulerFluxJacobian.H"
#include "TetraBinaryEulerJacobian.H"
#include "TurkelWeissPreconditioner.H"
#include <algorithm>
#include <functional>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>
#include <sstream>
#include <iomanip>
#include <utility>
#include <vector>

namespace Foam
{
// Approximate Euler-block LU-SGS. The off-diagonal is (A_n-lambda I)/2;
// the backward sweep consumes the already corrected upper neighbours.
// Every trial is staged before publication: failure leaves the field intact.
//
// The lambda in that expression is tetraFaceInvariant().spectralRadius -- the
// SAME quantity the flux's own stability bound uses. It must be, or the
// splitting is not a bound for the operator being solved. Three separate copies
// of the wave-speed formula used to live in this file; two of them floored the
// sound speed at 1 m/s and one did not, so the residual normalisation, the
// off-diagonal magnitude and the diagnostics disagreed inside the same solve.
// All three now route through TetraFaceInvariant.H.
class TetraLUSGS
{
public:
    struct Report
    {
        double residualNorm=0.,residualMaxNorm=0.,maxDelta=0.;
        double residualSeconds=0.,linearSeconds=0.,recoverySeconds=0.,boundarySeconds=0.;
        double worstCellInletSpectralFraction=0.;
        std::size_t worstResidualCell=0;
        int worstResidualComponent=0;
        double worstComponentValue=0.;
        std::size_t limited=0,frozen=0,stiffened=0;
        double frozenVolume=0.;
        std::string detail;
    };
private:
    TetraFiniteVolume& fv_;
    const TetraMesh& mesh_;
    std::function<void(std::size_t,std::size_t)> progressHook_;
    std::vector<std::int32_t> order_,rank_;
    std::vector<double> faceLambda_,diagonal_;
    std::vector<ConservedVars> residual_,delta_,solution_;
    bool compositionCoupling_=false; // opt-in audit path; extra EOS calls
    std::vector<TetraBinaryEulerJacobian::Coefficient> compositionDerivative_;
    double tMin_=80.,tMax_=4000.,pMin_=1e4,pMax_=2e7;
    double medianVolume_=0.;
    TurkelWeissPreconditioner preconditioner_;
    std::vector<double> cellBeta_;
    static void add(ConservedVars& a,const ConservedVars& b,double s)
    {
        a.rho+=s*b.rho;a.rhoE+=s*b.rhoE;a.rhoZ+=s*b.rhoZ;a.rhoC+=s*b.rhoC;
        for(int d=0;d<3;++d)a.rhoU[d]+=s*b.rhoU[d];
        for(std::size_t i=0;i<a.rhoY.size();++i)a.rhoY[i]+=s*b.rhoY[i];
    }
    static void scale(ConservedVars& a,double s)
    {
        a.rho*=s;a.rhoE*=s;a.rhoZ*=s;a.rhoC*=s;
        for(double& v:a.rhoU)v*=s;for(double& v:a.rhoY)v*=s;
    }
    void neighbourTerm(std::size_t c,std::size_t f,std::size_t other,
                       const ConservedVars& correction,ConservedVars& target) const
    {
        const double sign=mesh_.face(f).owner==int(c)?1.:-1.;
        const auto& normal=fv_.faceNormal(f);
        std::array<double,3> n{sign*normal[0],sign*normal[1],sign*normal[2]};
        ConservedVars product(fv_.species().size());
        const auto& q=fv_.Q(other);
        if(q.region==thermo::EquilibriumRegion::TWO_PHASE)
        {
            // A scalar spectral preconditioner in HEM cells requires only the
            // verified equilibrium sound speed. Do not invent frozen Cv/Cp.
            if(!tetraAdmitsAcoustics(q))throw std::runtime_error("invalid HEM preconditioner state");
            add(target,correction,.5*faceLambda_[f]);return;
        }
        const bool applied=compositionCoupling_
            ?TetraBinaryEulerJacobian::apply(fv_.W(other),q,correction,n,compositionDerivative_[other],product)
            :FrozenEulerFluxJacobian::apply(fv_.W(other),q,correction,n,product);
        if(!applied)throw std::runtime_error("tetra LU-SGS requires valid thermodynamic derivatives");
        add(target,correction,.5*faceLambda_[f]);
        add(target,product,-.5*fv_.faceArea(f));
    }
public:
    // Optional progress sink for the expensive per-cell stages. The kernel does
    // not know what telemetry is; the driver passes a lambda. Called every 256
    // cells, so the recovery sweep of a 1.76M-cell grid issues ~6.9k updates --
    // free compared with the CoolProp recoveries it is reporting on.
    void setProgressHook(std::function<void(std::size_t,std::size_t)> hook)
    {progressHook_=std::move(hook);}

    explicit TetraLUSGS(TetraFiniteVolume& fv):fv_(fv),mesh_(fv.mesh())
    {
        const auto cells=fv.cellCount();order_.resize(cells);rank_.resize(cells);
        std::vector<double> axial(cells);
        for(std::size_t c=0;c<cells;++c){order_[c]=c;axial[c]=mesh_.cellCentroid(c).z;}
        std::stable_sort(order_.begin(),order_.end(),[&](int a,int b){return axial[a]<axial[b];});
        for(std::size_t i=0;i<cells;++i)rank_[order_[i]]=i;
        std::vector<double> volumes(cells);
        for(std::size_t c=0;c<cells;++c)volumes[c]=mesh_.volume(c);
        std::nth_element(volumes.begin(),volumes.begin()+cells/2,volumes.end());
        medianVolume_=volumes[cells/2];
        faceLambda_.resize(mesh_.faceCount());diagonal_.resize(cells);
        residual_.assign(cells,ConservedVars(fv.species().size()));delta_=residual_;solution_=residual_;
        cellBeta_.assign(cells, 1.0);
    }
    void setGuard(double tMin,double tMax,double pMin,double pMax)
    {tMin_=tMin;tMax_=tMax;pMin_=pMin;pMax_=pMax;}
    void setFreezeVolumeLimit(double) {} // compatibility; freezing is disallowed
    void setCompositionCoupling(bool enabled) { compositionCoupling_=enabled; }
    void setLowMachPreconditioning(bool enable, double machCutoff = 0.05) {
        preconditioner_.setEnabled(enable);
        preconditioner_.setMachCutoff(machCutoff);
    }
    bool lowMachPreconditioningEnabled() const {
        return preconditioner_.isEnabled();
    }
    const std::vector<double>& diagonal()const{return diagonal_;}
    const std::vector<double>& faceLambda()const{return faceLambda_;}
    const std::vector<ConservedVars>& residual()const{return residual_;}
    double residualNorm()const
    {
        double sum=0.,vol=0.;
        for(std::size_t c=0;c<fv_.cellCount();++c)
        {sum+=residual_[c].rho*residual_[c].rho/mesh_.volume(c);vol+=mesh_.volume(c);}
        return std::sqrt(sum/vol);
    }
    std::pair<double,double> physicalResidualNorms(double physicalDt,
                                                   std::size_t* worstCell=nullptr,
                                                   int* worstComponent=nullptr,
                                                   double* worstValue=nullptr)const
    {
        if(!(physicalDt>0.)||!std::isfinite(physicalDt))
            throw std::invalid_argument("physical residual norm needs positive finite dt");
        double sum=0.;std::size_t components=0;double maximum=0.;
        for(std::size_t c=0;c<fv_.cellCount();++c)
        {
            const auto& u=fv_.U(c);const auto& q=fv_.Q(c);const auto& r=residual_[c];
            const double rho=std::max(std::abs(u.rho),1e-12);
            // Bare admitted sound speed, no 1 m/s floor. The reference energy
            // scale below uses it directly; a floor here previously made the
            // normalisation of a genuinely broken state look better than a
            // valid but slow-sound-speed one, which inverted the diagnostic.
            const double acoustic=q.acoustic.soundSpeed();
            if(!std::isfinite(acoustic)||acoustic<=0.)
                throw std::runtime_error(
                    "TetraLUSGS: unadmitted acoustic state in residual norm");
            const double energy=std::max({std::abs(u.rhoE),
                rho*std::max({std::abs(q.e),acoustic*acoustic,1.}),1e-12});
            const double factor=physicalDt/mesh_.volume(c);
            double local=0.;std::size_t localComponents=0;
            int dominantComponent=0;double dominantValue=0.;
            const auto accumulate=[&](double residual,double scaleValue)
            {
                const double value=factor*residual/std::max(scaleValue,1e-300);
                if(std::abs(value)>std::abs(dominantValue))
                {dominantValue=value;dominantComponent=int(localComponents);}
                local+=value*value;sum+=value*value;
                ++localComponents;++components;
            };
            accumulate(r.rho,rho);
            for(int d=0;d<3;++d)
                accumulate(r.rhoU[d],std::max(std::abs(u.rhoU[d]),rho*acoustic));
            accumulate(r.rhoE,energy);
            for(std::size_t s=0;s<r.rhoY.size();++s)
                accumulate(r.rhoY[s],std::max(std::abs(u.rhoY[s]),rho));
            accumulate(r.rhoZ,std::max(std::abs(u.rhoZ),rho));
            accumulate(r.rhoC,std::max(std::abs(u.rhoC),rho));
            const double localNorm=std::sqrt(local/std::max<std::size_t>(localComponents,1));
            if(localNorm>maximum)
            {
                maximum=localNorm;
                if(worstCell)*worstCell=c;
                if(worstComponent)*worstComponent=dominantComponent;
                if(worstValue)*worstValue=dominantValue;
            }
        }
        return {std::sqrt(sum/std::max<std::size_t>(components,1)),maximum};
    }
    std::pair<double,double> measurePhysicalResidual(
        const TetraBoundaryProvider& boundary,
        const std::vector<ConservedVars>& temporalReference,double physicalDt)
    {
        if(!(physicalDt>0.)||!std::isfinite(physicalDt))
            throw std::invalid_argument("physical residual measurement needs positive finite dt");
        if(temporalReference.size()!=fv_.cellCount())
            throw std::invalid_argument("physical residual measurement reference size mismatch");
        fv_.computeResidual(boundary,residual_);
        for(std::size_t c=0;c<fv_.cellCount();++c)
        {
            const double scaleValue=mesh_.volume(c)/physicalDt;
            add(residual_[c],fv_.U(c),scaleValue);
            add(residual_[c],temporalReference[c],-scaleValue);
        }
        return physicalResidualNorms(physicalDt);
    }
    Report step(const TetraBoundaryProvider& boundary,double cfl,double dtScale,
                double relax,const thermo::EquilibriumThermoModel& model,
                const std::function<void()>& validateBoundary={},
                double twoPhaseDtScale=1.0,double injectorDtScale=1.0,
                const std::vector<ConservedVars>* temporalReference=nullptr,
                double physicalDt=0.0)
    {
        const bool physicalStep=temporalReference!=nullptr;
        if(physicalStep!=(physicalDt>0.0))
            throw std::invalid_argument("physical LU-SGS needs reference state and positive dt");
        if(physicalStep&&temporalReference->size()!=fv_.cellCount())
            throw std::invalid_argument("physical LU-SGS reference size mismatch");
        if(compositionCoupling_)
        {
            if(fv_.species().size()!=2)
                throw std::runtime_error("binary composition coupling requires two species");
            compositionDerivative_.assign(fv_.cellCount(),{});
            std::size_t bad=fv_.cellCount();std::string detail;
            #pragma omp parallel for schedule(static)
            for(std::size_t c=0;c<fv_.cellCount();++c)
            {
                if(fv_.Q(c).region==thermo::EquilibriumRegion::TWO_PHASE)continue;
                try{compositionDerivative_[c]=TetraBinaryEulerJacobian::prepare(fv_.Q(c),model);}
                catch(const std::exception& e)
                {
                    #pragma omp critical(rfam_tetra_jacobian_failure)
                    {if(c<bad){bad=c;detail=e.what();}}
                }
            }
            if(bad<fv_.cellCount())
                throw std::runtime_error("composition derivative rejected at cell "
                    +std::to_string(bad)+": "+detail);
        }
        if(!std::isfinite(cfl)||cfl<=0.||!std::isfinite(dtScale)||dtScale<=0.
           ||!std::isfinite(relax)||relax<=0.||relax>1.)
            throw std::invalid_argument("invalid tetra LU-SGS controls");
        const auto clock=std::chrono::steady_clock::now;
        const auto begin=clock();
        Report report;fv_.computeResidual(boundary,residual_);
        if(physicalStep)
            for(std::size_t c=0;c<fv_.cellCount();++c)
            {
                const double scale=mesh_.volume(c)/physicalDt;
                add(residual_[c],fv_.U(c),scale);
                add(residual_[c],(*temporalReference)[c],-scale);
            }
        if(physicalStep)
        {
            const auto norms=physicalResidualNorms(physicalDt,&report.worstResidualCell,
                &report.worstResidualComponent,&report.worstComponentValue);
            report.residualNorm=norms.first;
            report.residualMaxNorm=norms.second;
        }
        else report.residualNorm=residualNorm();
        const auto residualDone=clock();
        report.residualSeconds=std::chrono::duration<double>(residualDone-begin).count();
        // Local preconditioning parameters, computed BEFORE the face loop: the
        // face spectral radius below needs the per-cell beta.
        //
        // SCOPE (2026-10-15). The Turkel-Weiss scaling is a low-Mach
        // gas-dynamics statement: it removes the acoustic/convective scale
        // separation (ratio ~1/M) of the compressible branches. It is applied
        // to the gas-like regions only. Dense liquid, two-phase and condensed
        // cells keep beta = 1: their stiffness is the curvature of the
        // equilibrium manifold (owned by twoPhaseDtScale), not the Mach number,
        // and a 1/beta pseudo-step there pushes the EOS recovery out of bounds
        // -- measured on the B grid, a liquid-LOX cell in the injector channel
        // rejected the update with the relaxation already halved 24 times.
        if(preconditioner_.isEnabled())
        {
            #pragma omp parallel for schedule(static)
            for(std::size_t c=0;c<fv_.cellCount();++c)
            {
                const auto region=fv_.Q(c).region;
                const bool gasLike=
                    region==thermo::EquilibriumRegion::SINGLE_VAPOR
                    ||region==thermo::EquilibriumRegion::SUPERCRITICAL
                    ||region==thermo::EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE;
                if(!gasLike){cellBeta_[c]=1.0;continue;}
                const auto& W=fv_.W(c);
                const auto& Q=fv_.Q(c);
                std::vector<double> velocity={W.U[0],W.U[1],W.U[2]};
                preconditioner_.updateLocal(velocity,Q.acoustic.soundSpeed());
                cellBeta_[c]=preconditioner_.beta();
            }
        }
        else
        {
            std::fill(cellBeta_.begin(),cellBeta_.end(),1.0);
        }
        for(std::size_t f=0;f<mesh_.faceCount();++f)
        {
            const auto& face=mesh_.face(f);const auto& n=fv_.faceNormal(f);
            // Face spectral radius via the single owner. This value forms the
            // LU-SGS bound, so it must match the flux branch exactly; see
            // TetraFaceInvariant.H and the M^4 note in AUSMPlusUpFlux.C.
            //
            // Low-Mach preconditioning (Turkel-Weiss, ADR 0001 P1): replace the
            // acoustic speed `a` by beta*a, which is the preconditioned
            // spectral radius |u_n| + beta*a. The SAME lambda enters the
            // diagonal and the off-diagonal scalar term (neighbourTerm), so the
            // diagonal-dominance margin sum*(1/(cfl*localScale)) is preserved
            // exactly as in the unpreconditioned case. Scaling the diagonal
            // alone would break that bound and the update gets rejected on the
            // first liquid-LOX cell (measured 2026-10-15, cell 337).
            auto speed=[&](const PrimitiveVars& w,const thermo::EquilibriumState& q,
                           double beta)
            {
                const FaceInvariant inv=tetraFaceInvariant(w,q,n);
                return std::abs(inv.normalVelocity)+beta*inv.soundSpeed;
            };
            const double betaOwner=cellBeta_[face.owner];
            double right=0.;
            if(face.neighbour>=0)right=speed(fv_.W(face.neighbour),fv_.Q(face.neighbour),
                                             cellBeta_[face.neighbour]);
            else{PrimitiveVars w(fv_.species().size());thermo::EquilibriumState q;boundary(face.owner,f,w,q);right=speed(w,q,betaOwner);}
            faceLambda_[f]=std::max(speed(fv_.W(face.owner),fv_.Q(face.owner),betaOwner),right)
                *fv_.faceArea(f);
        }
        if(physicalStep)
        {
            double all=0.,inlet=0.;
            for(auto f:mesh_.cellFaces(report.worstResidualCell))
            {
                all+=faceLambda_[f];
                const auto& face=mesh_.face(f);
                if(face.neighbour<0&&face.boundaryZone>=10&&face.boundaryZone<=22)
                    inlet+=faceLambda_[f];
            }
            report.worstCellInletSpectralFraction=all>0.?inlet/all:0.;
        }
        // Two-phase cells plus a one-cell halo: the stiff band's immediate
        // single-phase neighbours share the same steep manifold (A-route
        // step-20, cell 420173 at 4.19 MPa / 309.8 K closes at >1e-2 while its
        // neighbours are two-phase), so they get the same reduced step.
        std::vector<char> stiff(fv_.cellCount(),0);
        if(twoPhaseDtScale!=1.0)
        {
            for(std::size_t c=0;c<fv_.cellCount();++c)
                if(fv_.Q(c).region==thermo::EquilibriumRegion::TWO_PHASE)stiff[c]=1;
            for(std::size_t f=0;f<mesh_.faceCount();++f)
            {
                const auto& face=mesh_.face(f);
                if(face.neighbour<0)continue;
                if(stiff[face.owner]||stiff[face.neighbour])
                    stiff[face.owner]=stiff[face.neighbour]=1;
            }
        }
        std::size_t stiffened=0;
        for(std::size_t c=0;c<fv_.cellCount();++c)
        {
            double sum=0.;for(auto f:mesh_.cellFaces(c))sum+=faceLambda_[f];
            // Per-cell pseudo-time stiffening inside the two-phase band: the
            // equilibrium (rho,e) manifold is so curved there that a full
            // pseudo-time step can move a cell 0.5-2 % off the manifold and
            // exceed the acceptance gate. Those cells advance with a reduced
            // step (larger diagonal); the rest of the mesh keeps full speed,
            // unlike a global dt ramp which froze the injector-channel flush.
            // This is a convergence-acceleration parameter only.
            double localScale=dtScale;
            if(twoPhaseDtScale!=1.0&&stiff[c])
            {
                localScale=dtScale*twoPhaseDtScale;
                ++stiffened;
            }
            // Cell-size-compensated pseudo-time. On this mesh the injector
            // face cells are ~20 um while the chamber is ~0.5 mm, so the small
            // cells' CFL dt sets the pace of the whole pseudo-transient (the
            // jet penetration needs O(1000) steps). Boosting the *relative*
            // step of the smaller cells (factor (V_median/V_cell)^(1/3), capped)
            // equalises the information speed across the mesh. This is a
            // convergence-acceleration parameter only: the converged solution
            // is unchanged, and the per-step projection residuals are recorded.
            if(injectorDtScale!=1.0&&medianVolume_>0.)
            {
                const double ratio=medianVolume_/std::max(mesh_.volume(c),1e-30);
                const double boost=std::clamp(std::cbrt(ratio),1.0,injectorDtScale);
                localScale*=boost;
            }
            // Physical time derivative term: UNPRECONDITIONED (protects transient accuracy).
            // The acoustic term already carries the preconditioned lambda through
            // faceLambda_ (beta scaled in the face loop above), and the same
            // lambda is used by neighbourTerm, so no extra factor belongs here.
            const double physicalTerm=physicalStep?mesh_.volume(c)/physicalDt:0.0;
            const double acousticTerm=sum*(.5+1./(cfl*localScale));
            // Molecular spectral: keep unpreconditioned (viscous scale is physical)
            const double molecularTerm=fv_.molecularSpectral(c)*(1.+1./(cfl*localScale));
            diagonal_[c]=physicalTerm+acousticTerm+molecularTerm;
        }
        report.stiffened=stiffened;
        for(auto cell:order_)
        {
            const std::size_t c=cell;delta_[c]=residual_[c];scale(delta_[c],-1.);
            for(auto f:mesh_.cellFaces(c))
            {
                const auto& face=mesh_.face(f);int other=face.owner==cell?face.neighbour:face.owner;
                if(other>=0&&rank_[other]<rank_[c])neighbourTerm(c,f,other,delta_[other],delta_[c]);
            }
            scale(delta_[c],1./diagonal_[c]);
        }
        for(std::size_t i=order_.size();i-->0;)
        {
            const auto c=order_[i];ConservedVars accum(fv_.species().size());
            for(auto f:mesh_.cellFaces(c))
            {
                const auto& face=mesh_.face(f);int other=face.owner==c?face.neighbour:face.owner;
                if(other>=0&&rank_[other]>rank_[c])neighbourTerm(c,f,other,solution_[other],accum);
            }
            solution_[c]=delta_[c];add(solution_[c],accum,1./diagonal_[c]);
        }
        const auto linearDone=clock();
        report.linearSeconds=std::chrono::duration<double>(linearDone-residualDone).count();
        // Limit the correction, not the accepted state. No energy overwrite,
        // old-velocity reuse, species clipping, or frozen-cell continuation.
        const auto cells=fv_.cellCount();
        std::vector<CompleteCellState> staged(cells,CompleteCellState(fv_.species().size()));
        std::size_t bad=cells,limited=0;std::string detail;double maxDelta=0.;
        #pragma omp parallel for schedule(static) reduction(+:limited) reduction(max:maxDelta)
        for(std::size_t c=0;c<cells;++c)
        {
            if(progressHook_&&(c&255u)==0u)progressHook_(c,cells);
            const auto& u=fv_.U(c);const auto& d=solution_[c];double factor=relax;
            const double acoustic=fv_.Q(c).acoustic.soundSpeed();
            factor=std::min(factor,.2*u.rho/(std::abs(d.rho)+1e-300));
            factor=std::min(factor,.2*u.rho*std::max({std::abs(fv_.Q(c).e),acoustic*acoustic,1.})/(std::abs(d.rhoE)+1e-300));
            double dm=0.,momentum=0.;for(int k=0;k<3;++k){dm+=d.rhoU[k]*d.rhoU[k];momentum+=u.rhoU[k]*u.rhoU[k];}
            factor=std::min(factor,.2*(std::sqrt(momentum)+u.rho*acoustic)/(std::sqrt(dm)+1e-300));
            bool accepted=false;std::string why;
            ConservedVars lastTrial=u;
            for(int attempt=0;attempt<24;++attempt)
            {
                ConservedVars trial=u;add(trial,d,factor);lastTrial=trial;
                CompleteCellConversionResult recovered;
                try{recovered=recoverTetraState(trial,fv_.W(c),fv_.Q(c),model);}
                catch(const std::exception& e){why=e.what();factor*=.5;continue;}
                if(recovered.success())
                {
                    const auto& q=recovered.state->Q;
                    if(q.T>=tMin_&&q.T<=tMax_&&q.p>=pMin_&&q.p<=pMax_)
                    {
                        double change=std::abs(factor*d.rho)/u.rho;
                        change=std::max(change,std::abs(q.p-fv_.Q(c).p)/fv_.Q(c).p);
                        change=std::max(change,std::abs(q.T-fv_.Q(c).T)/fv_.Q(c).T);
                        for(std::size_t s=0;s<q.overallMassFractions.size();++s)
                            change=std::max(change,std::abs(q.overallMassFractions[s]-fv_.W(c).Y[s]));
                        for(int dim=0;dim<3;++dim)
                            change=std::max(change,std::abs(recovered.state->W.U[dim]-fv_.W(c).U[dim])/acoustic);
                        staged[c]=std::move(*recovered.state);accepted=true;
                        limited+=factor<relax;maxDelta=std::max(maxDelta,change);break;
                    }
                    why="physical range guard";
                }
                else
                {
                    // Include the status category: the backend text alone can be a
                    // *success* status string (e.g. COOLPROP_HEOS_SINGLE_PHASE)
                    // when a later check rejects an otherwise valid conversion,
                    // which masked the real failure category.
                    why=recovered.diagnostics.backendStatus
                        +" [status="+thermo::equilibriumStatusName(recovered.status)+"]";
                }
                factor*=.5;
            }
            if(!accepted)
            {
                std::ostringstream message;message<<std::setprecision(17)<<why
                    <<" previous_p="<<fv_.W(c).p<<" previous_T="<<fv_.W(c).T
                    <<" previous_rho="<<u.rho<<" previous_e="<<fv_.Q(c).e
                    <<" region="<<thermo::equilibriumRegionName(fv_.Q(c).region)
                    <<" trial_rho="<<lastTrial.rho<<" trial_rhoE="<<lastTrial.rhoE
                    <<" scale="<<2.*factor<<" trial_momentum=";
                for(double value:lastTrial.rhoU)message<<value<<",";
                message<<" trial_rhoY=";for(double value:lastTrial.rhoY)message<<value<<",";
                #pragma omp critical(rfam_tetra_update_failure)
                {if(c<bad){bad=c;detail=message.str();}}
            }
        }
        if(bad<cells)throw std::runtime_error("tetra LU-SGS update rejected at cell "+std::to_string(bad)+": "+detail);
        const auto recoveryDone=clock();
        report.recoverySeconds=std::chrono::duration<double>(recoveryDone-linearDone).count();
        for(std::size_t c=0;c<cells;++c)fv_.exchangeState(c,staged[c].W,staged[c].U,staged[c].Q);
        try{if(validateBoundary)validateBoundary();}
        catch(...)
        {
            for(std::size_t c=0;c<cells;++c)fv_.exchangeState(c,staged[c].W,staged[c].U,staged[c].Q);
            throw;
        }
        report.boundarySeconds=std::chrono::duration<double>(clock()-recoveryDone).count();
        report.limited=limited;report.maxDelta=maxDelta;return report;
    }
};
}
#endif
