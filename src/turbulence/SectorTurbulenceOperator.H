#ifndef SectorTurbulenceOperator_H
#define SectorTurbulenceOperator_H

#include "Operator.H"
#include "KOmegaSST.H"
#include "SectorMeshMetrics.H"
#include "SpeciesRegistry.H"
#include "TransportProps.H"
#include "WallFunction.H"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>
#include "WALE.H"

namespace Foam
{

class SectorTurbulenceOperator : public Operator
{
public:
    enum class TurbulenceMode { SST, WALE };

private:
    using Tensor = double[3][3];
    const SectorMeshMetrics& mesh_;
    const SpeciesRegistry& registry_;
    mutable std::vector<double> k_, omega_;
    mutable bool initialized_;
    mutable std::vector<double> kStepStart_, omegaStepStart_;
    mutable bool initializedStepStart_ = false;
    mutable bool stepTransactionActive_ = false;
    double kInitial_, omegaInitial_, safety_;
    // 显式诊断模式：允许气体闭合外推（T<200 K 或 rho>100 kg/m3），
    // 但绝不静默——每次 apply 都会审计并输出降级 cell 计数。
    bool diagnosticMode_ = false;
    mutable long diagnosticExtrapolatedCells_ = 0;
    mutable bool diagnosticAuditPrinted_ = false;
    TurbulenceMode mode_ = TurbulenceMode::SST;

    int index(int i,int j,int k) const
    { return k*mesh_.Nx()*mesh_.Nr()+j*mesh_.Nx()+i; }

    void initialize() const
    {
        if(initialized_) return;
        const int count=mesh_.Nx()*mesh_.Nr()*mesh_.Ntheta();
        k_.assign(count,kInitial_);
        omega_.assign(count,omegaInitial_);
        initialized_=true;
    }

    double kValue(int id) const
    { return initialized_?k_[static_cast<size_t>(id)]:kInitial_; }

    double omegaValue(int id) const
    { return initialized_?omega_[static_cast<size_t>(id)]:omegaInitial_; }

    double wallDistance(int i,int j,int k) const
    {
        const int wallFace=mesh_.Nr()-1;
        const FaceAreaVector area=mesh_.radialAreaVector(i,wallFace,k,+1);
        const double areaMagnitude=std::sqrt(
            area.x*area.x+area.y*area.y+area.z*area.z);
        if(areaMagnitude<=1e-30) return 1e-7;
        const FaceAreaVector wallCell=mesh_.cellCenter(i,wallFace-1,k);
        const FaceAreaVector wallGhost=mesh_.cellCenter(i,wallFace,k);
        const double firstCellDistance=0.5*std::abs(
            ((wallGhost.x-wallCell.x)*area.x+(wallGhost.y-wallCell.y)*area.y
             +(wallGhost.z-wallCell.z)*area.z)/areaMagnitude);
        const double layerFactor=2.0*(wallFace-j)-1.0;
        return std::max(layerFactor*firstCellDistance,1e-7);
    }

    void validateSpeciesLayout(const Domain& d) const
    {
        if(d.speciesNames().empty() ||
           d.speciesNames().size()!=static_cast<size_t>(d.nSpecies()))
            throw std::runtime_error(
                "sector SST transport cell -1 status=SPECIES_MISMATCH");
    }

    void validateDomainLayout(const Domain& d) const
    {
        if(mesh_.Nx()!=d.Nx() || mesh_.Nr()!=d.Ny()
           ||mesh_.Ntheta()!=d.Nz())
            throw std::runtime_error(
                "sector SST transport cell -1 status=INVALID_INPUT "
                "detail=DOMAIN_MESH_MISMATCH");
        if(d.Nx()<4 || d.Ny()<4 || d.Nz()<4)
            throw std::runtime_error(
                "sector SST transport cell -1 status=INVALID_INPUT "
                "detail=INSUFFICIENT_HALO_DIMENSIONS");
    }

    std::vector<TransportProps::Point> transportField(const Domain& d) const
    {
        validateDomainLayout(d);
        validateSpeciesLayout(d);
        std::vector<TransportProps::Point> points(static_cast<size_t>(d.N()));
        std::vector<TransportProps::Status> statuses(
            static_cast<size_t>(d.N()),TransportProps::Status::INVALID_INPUT);
        int firstInvalid=d.N();
        // All boundary halos are owned by boundary operators and may be
        // outside the NASA-7 transport range (for example an outlet ghost
        // overshoot). Evaluate only active cells, then reuse the nearest
        // active-cell coefficient for every halo entry. Active-cell failures
        // remain fail-closed and are still reported with their real indices.
        #pragma omp parallel for collapse(3) reduction(min:firstInvalid) schedule(static)
        for(int i=1;i<d.Nx()-1;++i)
        for(int j=1;j<d.Ny()-1;++j)
        for(int k=1;k<d.Nz()-1;++k)
        {
            const int id=index(i,j,k);
            const TransportProps::Result result=
                evaluateCellTransport(d,i,j,k);
            statuses[static_cast<size_t>(id)]=result.status;
            if(result.success())
                points[static_cast<size_t>(id)]=result.point;
            else
                firstInvalid=std::min(firstInvalid,id);
        }
        if(firstInvalid!=d.N())
        {
            const int ci = firstInvalid % d.Nx();
            const int cj = (firstInvalid / d.Nx()) % d.Ny();
            const int ck = firstInvalid / (d.Nx() * d.Ny());
            const PrimitiveVars& cw = d.W(ci, cj, ck);
            throw std::runtime_error(
                "sector SST transport cell "+std::to_string(firstInvalid)
                +" status="+TransportProps::statusName(
                    statuses[static_cast<size_t>(firstInvalid)])
                +" T="+std::to_string(cw.T)
                +" rho="+std::to_string(d.Q(ci,cj,ck).rho)
                +" Y0="+std::to_string(cw.Y.empty()?-1.0:cw.Y[0]));
        }
        // A valid physical wall ghost carries the wall-side conductivity used
        // by the wall heat-flux law. Keep it when available; an unsupported
        // wall state (for example a diagnostic HEM fixture) falls back to the
        // adjacent active point just like the other halos.
        for(int i=1;i<d.Nx()-1;++i)
        for(int k=1;k<d.Nz()-1;++k)
        {
            const int wallId=index(i,d.Ny()-1,k);
            const TransportProps::Result result=
                evaluateCellTransport(d,i,d.Ny()-1,k);
            if(result.success())
            {
                points[static_cast<size_t>(wallId)]=result.point;
                statuses[static_cast<size_t>(wallId)]=result.status;
            }
            else
            {
                const int activeId=index(i,d.Ny()-2,k);
                points[static_cast<size_t>(wallId)]=
                    points[static_cast<size_t>(activeId)];
                statuses[static_cast<size_t>(wallId)]=
                    TransportProps::Status::SUCCESS;
            }
        }
        for(int i=0;i<d.Nx();++i)
        for(int j=0;j<d.Ny();++j)
        for(int k=0;k<d.Nz();++k)
        {
            const int activeI=std::clamp(i,1,d.Nx()-2);
            const int activeJ=std::clamp(j,1,d.Ny()-2);
            const int activeK=std::clamp(k,1,d.Nz()-2);
            const int id=index(i,j,k);
            const int activeId=index(activeI,activeJ,activeK);
            const bool physicalWallHalo=
                j==d.Ny()-1 && i>0 && i<d.Nx()-1
                && k>0 && k<d.Nz()-1;
            if(id!=activeId && !physicalWallHalo)
            {
                points[static_cast<size_t>(id)]=
                    points[static_cast<size_t>(activeId)];
                statuses[static_cast<size_t>(id)]=
                    TransportProps::Status::SUCCESS;
            }
        }
        return points;
    }

    // 后端感知的单元输运评估（docs/cryogenic_transport_backend.md 第 5 节）：
    //   1. 验证的单相气相 -> legacy gas；
    //   2. 纯/近纯 O2/C3H8 稠密态（已注册表）-> 稠密纯流体表；
    //   3. 稠密多组分/两相 -> 诊断模式显式外推（审计），否则 fail-close。
    TransportProps::Result evaluateCellTransport(
        const Domain& d,int i,int j,int k) const
    {
        const auto& names=d.speciesNames();
        const PrimitiveVars& w=d.W(i,j,k);
        const thermo::EquilibriumState& q=d.Q(i,j,k);
        const TransportProps::Backend backend=TransportProps::selectBackend(
            names,w.Y,w.T,q.rho);
        if(backend==TransportProps::Backend::CRYOGENIC_PURE_FLUID)
        {
            if(!q.thermodynamics.has_value())
            {
                TransportProps::Result result;
                result.status=TransportProps::Status::MISSING_THERMODYNAMIC_PAYLOAD;
                return result;
            }
            const auto& payload=*q.thermodynamics;
            const TransportProps::Result cryo=TransportProps::evaluateCryogenicPureFluid(
                names,w.Y,w.p,w.T,q.rho,payload.Cp,payload.Cv);
            if(cryo.success() || !diagnosticMode_)
                return cryo;
            // Diagnostic fallback: a diluted/dilute HEM two-phase or
            // near-pure state can select the cryogenic backend and still
            // fail its pure-fluid gate (mixed composition / two-phase
            // region / table two-phase flag).  In explicit diagnostic
            // mode retry through the audited gas extrapolation.
        }
        if(backend==TransportProps::Backend::LEGACY_GAS)
        {
            const TransportProps::Result committed=
                TransportProps::evaluateCommittedLegacyGas(
                    names,w,q,d.thermoModel(),registry_);
            if(committed.success() || !diagnosticMode_)
                return committed;
            // Backend selection uses the broad gas envelope (T/rho), while
            // the committed route also enforces each active species' NASA-7
            // range.  In explicit diagnostic mode, retry that narrower range
            // rejection through the audited extrapolation path.
        }
        // UNSUPPORTED（稠密多组分/两相/表外）or a committed gas-range reject.
        if(diagnosticMode_)
        {
            const TransportProps::DiagnosticGasResult diag=
                TransportProps::evaluateDiagnosticGas(
                    names,w,q,d.thermoModel(),registry_,true);
            if(diag.extrapolation
                !=TransportProps::DiagnosticExtrapolation::NONE
                &&diag.result.success())
                #pragma omp atomic
                diagnosticExtrapolatedCells_ += 1;
            return diag.result;
        }
        TransportProps::Result result;
        result.status=TransportProps::Status::UNSUPPORTED_PHASE;
        return result;
    }

    void applyFieldBoundaries(
        const Domain& d,const std::vector<TransportProps::Point>& points) const
    {
        const int nx=d.Nx(),ny=d.Ny(),nz=d.Nz();
        for(int j=0;j<ny;++j) for(int k=0;k<nz;++k)
        {
            k_[index(0,j,k)]=kInitial_;omega_[index(0,j,k)]=omegaInitial_;
            k_[index(nx-1,j,k)]=k_[index(nx-2,j,k)];
            omega_[index(nx-1,j,k)]=omega_[index(nx-2,j,k)];
        }
        for(int i=0;i<nx;++i) for(int k=0;k<nz;++k)
        {
            k_[index(i,0,k)]=k_[index(i,1,k)];
            omega_[index(i,0,k)]=omega_[index(i,1,k)];
            k_[index(i,ny-1,k)]=KOmegaSST::k_min;
            const double y=wallDistance(i,ny-2,k);
            const double nu=points[index(i,ny-2,k)].nu;
            omega_[index(i,ny-1,k)]=std::max(
                60.0*nu/(KOmegaSST::beta1*y*y),KOmegaSST::om_min);
        }
        for(int i=0;i<nx;++i) for(int j=0;j<ny;++j)
        {
            k_[index(i,j,0)]=k_[index(i,j,nz-2)];
            omega_[index(i,j,0)]=omega_[index(i,j,nz-2)];
            k_[index(i,j,nz-1)]=k_[index(i,j,1)];
            omega_[index(i,j,nz-1)]=omega_[index(i,j,1)];
        }
    }

    static FaceAreaVector subtract(const FaceAreaVector& a,const FaceAreaVector& b)
    { return {a.x-b.x,a.y-b.y,a.z-b.z}; }

    static void solveTranspose(const FaceAreaVector basis[3],const double rhs[3],double result[3])
    {
        const double a00=basis[0].x,a01=basis[0].y,a02=basis[0].z;
        const double a10=basis[1].x,a11=basis[1].y,a12=basis[1].z;
        const double a20=basis[2].x,a21=basis[2].y,a22=basis[2].z;
        const double det=a00*(a11*a22-a12*a21)-a01*(a10*a22-a12*a20)
                        +a02*(a10*a21-a11*a20);
        if(std::abs(det)<1e-30) { result[0]=result[1]=result[2]=0.0; return; }
        result[0]=(rhs[0]*(a11*a22-a12*a21)-a01*(rhs[1]*a22-a12*rhs[2])
                  +a02*(rhs[1]*a21-a11*rhs[2]))/det;
        result[1]=(a00*(rhs[1]*a22-a12*rhs[2])-rhs[0]*(a10*a22-a12*a20)
                  +a02*(a10*rhs[2]-rhs[1]*a20))/det;
        result[2]=(a00*(a11*rhs[2]-rhs[1]*a21)-a01*(a10*rhs[2]-rhs[1]*a20)
                  +rhs[0]*(a10*a21-a11*a20))/det;
    }

    void velocityGradient(const Domain& d,int i,int j,int k,Tensor& gradient) const
    {
        const int im=std::max(i-1,0),ip=std::min(i+1,d.Nx()-1);
        const int jm=std::max(j-1,0),jp=std::min(j+1,d.Ny()-1);
        const int km=std::max(k-1,0),kp=std::min(k+1,d.Nz()-1);
        const FaceAreaVector basis[3]={
            subtract(mesh_.cellCenter(ip,j,k),mesh_.cellCenter(im,j,k)),
            subtract(mesh_.cellCenter(i,jp,k),mesh_.cellCenter(i,jm,k)),
            subtract(mesh_.cellCenter(i,j,kp),mesh_.cellCenter(i,j,km))};
        for(int m=0;m<3;++m)
        {
            const double rhs[3]={d.W(ip,j,k).U[m]-d.W(im,j,k).U[m],
                                 d.W(i,jp,k).U[m]-d.W(i,jm,k).U[m],
                                 d.W(i,j,kp).U[m]-d.W(i,j,km).U[m]};
            solveTranspose(basis,rhs,gradient[m]);
        }
    }

    void scalarGradient(const std::vector<double>& field,
                        int i,int j,int k,double gradient[3]) const
    {
        const int im=std::max(i-1,0),ip=std::min(i+1,mesh_.Nx()-1);
        const int jm=std::max(j-1,0),jp=std::min(j+1,mesh_.Nr()-1);
        const int km=std::max(k-1,0),kp=std::min(k+1,mesh_.Ntheta()-1);
        const FaceAreaVector basis[3]={
            subtract(mesh_.cellCenter(ip,j,k),mesh_.cellCenter(im,j,k)),
            subtract(mesh_.cellCenter(i,jp,k),mesh_.cellCenter(i,jm,k)),
            subtract(mesh_.cellCenter(i,j,kp),mesh_.cellCenter(i,j,km))};
        const double rhs[3]={field[index(ip,j,k)]-field[index(im,j,k)],
                             field[index(i,jp,k)]-field[index(i,jm,k)],
                             field[index(i,j,kp)]-field[index(i,j,km)]};
        solveTranspose(basis,rhs,gradient);
    }

    double crossDiffusion(int i,int j,int k) const
    {
        double gradK[3],gradOmega[3];
        scalarGradient(k_,i,j,k,gradK);
        scalarGradient(omega_,i,j,k,gradOmega);
        return KOmegaSST::crossDiffusion(
            omega_[index(i,j,k)],gradK,gradOmega);
    }

    double strain(const Domain& d,int i,int j,int k) const
    {
        Tensor g; velocityGradient(d,i,j,k,g);
        double sum=0.0;
        for(int a=0;a<3;++a) for(int b=0;b<3;++b)
        { const double s=0.5*(g[a][b]+g[b][a]); sum+=2.0*s*s; }
        return std::sqrt(sum);
    }

    void stress(const Domain& d,int i,int j,int k,double nut,Tensor& tau) const
    {
        Tensor g; velocityGradient(d,i,j,k,g);
        const int id=index(i,j,k);
        const double rho=d.Q(i,j,k).rho;
        const double muT=rho*nut;
        const double divergence=g[0][0]+g[1][1]+g[2][2];
        for(int a=0;a<3;++a) for(int b=0;b<3;++b)
            tau[a][b]=muT*(g[a][b]+g[b][a])
                     -(a==b ? (2.0/3.0)*(muT*divergence+rho*k_[id]) : 0.0);
    }

    void addStressFace(const Domain& d,int i,int j,int k,int ni,int nj,int nk,
                       const FaceAreaVector& area,
                       const std::vector<TransportProps::Point>& points,
                       const std::vector<double>& nutField,
                       double residual[3],double& energyResidual) const
    {
        if(j==d.Ny()-2&&nj==d.Ny()-1)
        {
            const double areaMagnitude=std::sqrt(
                area.x*area.x+area.y*area.y+area.z*area.z);
            if(areaMagnitude<=1e-30) return;
            const double normal[3]={area.x/areaMagnitude,area.y/areaMagnitude,
                                    area.z/areaMagnitude};
            const PrimitiveVars& gas=d.W(i,j,k);
            const double normalVelocity=gas.U[0]*normal[0]+gas.U[1]*normal[1]
                                       +gas.U[2]*normal[2];
            double tangential[3],speedSquared=0.0;
            for(int m=0;m<3;++m)
            {
                tangential[m]=gas.U[m]-normalVelocity*normal[m];
                speedSquared+=tangential[m]*tangential[m];
            }
            const double speed=std::sqrt(speedSquared);
            if(speed<=1e-12) return;
            const double nu=points[index(i,j,k)].nu;
            const double shear=WallFunction::wallShear(
                wallDistance(i,j,k),speed,nu,d.Q(i,j,k).rho);
            for(int m=0;m<3;++m)
                residual[m]-=shear*tangential[m]/speed*areaMagnitude;
            return;
        }
        Tensor left,right;
        stress(d,i,j,k,nutField[index(i,j,k)],left);
        stress(d,ni,nj,nk,nutField[index(ni,nj,nk)],right);
        double traction[3];
        for(int m=0;m<3;++m)
        {
            traction[m]=0.5*((left[m][0]+right[m][0])*area.x
                            +(left[m][1]+right[m][1])*area.y
                            +(left[m][2]+right[m][2])*area.z);
            residual[m]+=traction[m];
            const double faceVelocity=0.5*(
                d.W(i,j,k).U[m]+d.W(ni,nj,nk).U[m]);
            energyResidual+=faceVelocity*traction[m];
        }
    }

public:
    SectorTurbulenceOperator(const SectorMeshMetrics& mesh,
                             const SpeciesRegistry& registry,
                             double intensity=0.05,double lengthScale=0.001,
                             double safety=0.25,
                             bool diagnosticMode=false,
                             TurbulenceMode mode=TurbulenceMode::SST)
        :mesh_(mesh),registry_(registry),initialized_(false),safety_(safety),
         diagnosticMode_(diagnosticMode),mode_(mode)
    {
        const double referenceVelocity=100.0;
        kInitial_=1.5*std::pow(referenceVelocity*intensity,2);
        omegaInitial_=std::sqrt(kInitial_)/(std::pow(KOmegaSST::beta_star,0.25)*lengthScale);
    }

    // WALE is an algebraic SGS model (apply() advances no k/omega transport
    // state), so it is stateless and may sit inside a dual-time inner chain.
    // SST owns k/omega state and is rejected there (physical-time semantics).
    bool ownsPersistentState() const override
    { return mode_ != TurbulenceMode::WALE; }

    // 显式查询：当前 operator 是否在诊断外推模式（供壁面热流/输出层审计）
    bool diagnosticMode() const { return diagnosticMode_; }

    void beginStepTransaction(const Domain&) const override
    {
        if(stepTransactionActive_)
            throw std::logic_error("nested sector turbulence step transaction");
        kStepStart_=k_;
        omegaStepStart_=omega_;
        initializedStepStart_=initialized_;
        stepTransactionActive_=true;
    }

    void acceptStepTransaction() const noexcept override
    { stepTransactionActive_=false; }

    void rollbackStepTransaction() const noexcept override
    {
        if(!stepTransactionActive_) return;
        k_.swap(kStepStart_);
        omega_.swap(omegaStepStart_);
        initialized_=initializedStepStart_;
        stepTransactionActive_=false;
    }

    void averageWithStepStartState() const override
    {
        if(!stepTransactionActive_ || k_.size()!=omega_.size())
            throw std::logic_error(
                "sector turbulence RK2 state transaction is invalid");

        if(initializedStepStart_ &&
           (k_.size()!=kStepStart_.size()
            || omega_.size()!=omegaStepStart_.size()))
            throw std::logic_error(
                "sector turbulence RK2 state size changed");

        for(std::size_t id=0;id<k_.size();++id)
        {
            const double startK=initializedStepStart_ ? kStepStart_[id] : kInitial_;
            const double startOmega=initializedStepStart_
                ? omegaStepStart_[id] : omegaInitial_;
            k_[id]=0.5*(startK+k_[id]);
            omega_[id]=0.5*(startOmega+omega_[id]);
        }
    }

    double turbulentViscosity(
        int i,int j,int k,const Domain& d,double molecularNu) const
    {
        if(mode_==TurbulenceMode::WALE)
        {
            // Wall-Adapting Local Eddy-viscosity (Nicoud & Ducros 1999):
            // nu_sgs = (Cw*delta)^2 * (S^d:S^d)^(3/2)
            //        / ((S:S)^(5/2) + (S^d:S^d)^(5/4))
            Tensor g;
            velocityGradient(d,i,j,k,g);
            const double delta=WALE::filterWidth(mesh_.volume(i,j,k));
            return WALE::nuSgs(g,delta);
        }
        const int id=index(i,j,k);
        const double kval=kValue(id);
        const double omega=omegaValue(id);
        const double s=strain(d,i,j,k);
        const double f2=KOmegaSST::F2(
            kval,omega,wallDistance(i,j,k),molecularNu);
        return KOmegaSST::nu_t(kval,omega,s,f2);
    }

    double molecularKinematicViscosity(int i,int j,int k,const Domain& d) const
    {
        validateDomainLayout(d);
        validateSpeciesLayout(d);
        const TransportProps::Result result=
            evaluateCellTransport(d,i,j,k);
        if(!result.success())
            throw std::runtime_error(
                "sector SST transport cell "+std::to_string(index(i,j,k))
                +" status="+TransportProps::statusName(result.status));
        return result.point.nu;
    }

    double turbulentViscosity(int i,int j,int k,const Domain& d) const
    {
        return turbulentViscosity(
            i,j,k,d,molecularKinematicViscosity(i,j,k,d));
    }

    double crossDiffusionSource(int i,int j,int k) const
    {
        initialize();
        return crossDiffusion(i,j,k);
    }

    double turbulentKValue(int i,int j,int k) const
    { return kValue(index(i,j,k)); }

    double turbulentOmegaValue(int i,int j,int k) const
    { return omegaValue(index(i,j,k)); }

    bool persistentStateInitialized() const noexcept { return initialized_; }

    const std::vector<double>& kField() const { initialize(); return k_; }
    const std::vector<double>& omegaField() const { initialize(); return omega_; }
    std::size_t expectedFieldSize() const noexcept
    {
        return static_cast<std::size_t>(mesh_.Nx())
            *static_cast<std::size_t>(mesh_.Nr())
            *static_cast<std::size_t>(mesh_.Ntheta());
    }
    bool meshMatches(const Domain& domain) const noexcept
    {
        return mesh_.Nx()==domain.Nx()
            &&mesh_.Nr()==domain.Ny()
            &&mesh_.Ntheta()==domain.Nz();
    }
    bool fieldsAreRestorable(
        const std::vector<double>& kField,
        const std::vector<double>& omegaField) const noexcept
    {
        if(kField.size()!=expectedFieldSize()
           ||omegaField.size()!=expectedFieldSize())
            return false;
        for(std::size_t id=0;id<kField.size();++id)
            if(!std::isfinite(kField[id])||kField[id]<KOmegaSST::k_min
               ||!std::isfinite(omegaField[id])
               ||omegaField[id]<KOmegaSST::om_min)
                return false;
        return true;
    }
    bool stepTransactionActive() const noexcept
    { return stepTransactionActive_; }
    void restoreFields(const std::vector<double>& kField,
                       const std::vector<double>& omegaField)
    {
        std::vector<double> candidateK=kField;
        std::vector<double> candidateOmega=omegaField;
        restoreFields(std::move(candidateK),std::move(candidateOmega));
    }
    void restoreFields(std::vector<double>&& kField,
                       std::vector<double>&& omegaField)
    {
        if(stepTransactionActive_)
            throw std::logic_error(
                "cannot restore sector turbulence during a step transaction");
        if(!fieldsAreRestorable(kField,omegaField))
            throw std::invalid_argument("invalid sector turbulence checkpoint fields");
        k_.swap(kField); omega_.swap(omegaField); initialized_=true;
    }

    void apply(Domain& d,double dt) const override
    {
        if(mode_==TurbulenceMode::WALE)
        {
            // WALE 是代数 SGS 模型：不推进 k/omega 输运方程。
            // k/omega 仅保留初始参考值，供壁面函数/输出兼容。
            initialize();
            return;
        }
        diagnosticExtrapolatedCells_ = 0;
        diagnosticAuditPrinted_ = false;
        const std::vector<TransportProps::Point> points=transportField(d);
        if(diagnosticMode_ && diagnosticExtrapolatedCells_ > 0
           && !diagnosticAuditPrinted_)
        {
            // 绝不静默：诊断模式下每次 apply 显式报告气体闭合外推范围。
#ifdef RFAM_THERMO_DIAGNOSTICS
            std::cerr << "[DIAGNOSTIC] sector SST gas-closure extrapolation: "
                      << diagnosticExtrapolatedCells_ << " cells outside "
                      << "T>=200K/rho<=100kg/m3 (T<200K and/or rho>100kg/m3); "
                      << "Sutherland/Wilke/Eucken values are nonphysical for "
                      << "dense cryogenic states, diagnostic-only\n";
#endif
            diagnosticAuditPrinted_ = true;
        }
        initialize();
        applyFieldBoundaries(d,points);
        std::vector<double> nextK=k_,nextOmega=omega_;
        std::vector<double> nutField(k_.size(),0.0);
        #pragma omp parallel for collapse(3) schedule(static)
        for(int i=0;i<d.Nx();++i)
        for(int j=0;j<d.Ny();++j)
        for(int k=0;k<d.Nz();++k)
            nutField[index(i,j,k)]=turbulentViscosity(
                i,j,k,d,points[index(i,j,k)].nu);
        #pragma omp parallel for collapse(3) schedule(static)
        for(int i=1;i<d.Nx()-1;++i)
        for(int j=1;j<d.Ny()-1;++j)
        for(int k=1;k<d.Nz()-1;++k)
        {
            const int id=index(i,j,k);
            const double s=strain(d,i,j,k);
            const double CDkw=crossDiffusion(i,j,k);
            const double f1=KOmegaSST::F1(
                k_[id],omega_[id],wallDistance(i,j,k),points[id].nu,
                std::max(CDkw,1e-20));
            const double f2=KOmegaSST::F2(
                k_[id],omega_[id],wallDistance(i,j,k),points[id].nu);
            const double nut=KOmegaSST::nu_t(k_[id],omega_[id],s,f2);
            const double pk=KOmegaSST::Pk(nut,s,k_[id],omega_[id]);

            auto transport=[&](const std::vector<double>& field,double sigma)
            {
                double rate=0.0;
                auto face=[&](int ni,int nj,int nk,const FaceAreaVector& area)
                {
                    const double areaMag=std::sqrt(area.x*area.x+area.y*area.y+area.z*area.z);
                    if(areaMag<=1e-30) return;
                    const FaceAreaVector a=mesh_.cellCenter(i,j,k),b=mesh_.cellCenter(ni,nj,nk);
                    const double distance=std::max(((b.x-a.x)*area.x+(b.y-a.y)*area.y
                        +(b.z-a.z)*area.z)/areaMag,1e-12);
                    const int neighbor=index(ni,nj,nk);
                    const PrimitiveVars& wc=d.W(i,j,k);const PrimitiveVars& wn=d.W(ni,nj,nk);
                    const double massFlux=0.5*((d.Q(i,j,k).rho*wc.U[0]+d.Q(ni,nj,nk).rho*wn.U[0])*area.x
                        +(d.Q(i,j,k).rho*wc.U[1]+d.Q(ni,nj,nk).rho*wn.U[1])*area.y
                        +(d.Q(i,j,k).rho*wc.U[2]+d.Q(ni,nj,nk).rho*wn.U[2])*area.z);
                    const double faceValue=massFlux>=0.0?field[id]:field[neighbor];
                    const double muFace=0.5*(points[id].mu+points[neighbor].mu);
                    const double muTCell=d.Q(i,j,k).rho*nutField[id];
                    const double muTNeighbor=d.Q(ni,nj,nk).rho*nutField[neighbor];
                    const double diffusion=(muFace+sigma*0.5*(muTCell+muTNeighbor))
                        *(field[neighbor]-field[id])*areaMag/distance;
                    rate+=(-massFlux*(faceValue-field[id])+diffusion)
                        /(d.Q(i,j,k).rho*mesh_.volume(i,j,k));
                };
                face(i-1,j,k,mesh_.axialAreaVector(i,j,k,-1));
                face(i+1,j,k,mesh_.axialAreaVector(i+1,j,k,+1));
                face(i,j-1,k,mesh_.radialAreaVector(i,j,k,-1));
                face(i,j+1,k,mesh_.radialAreaVector(i,j+1,k,+1));
                face(i,j,k-1,mesh_.azimuthalAreaVector(i,j,k,-1));
                face(i,j,k+1,mesh_.azimuthalAreaVector(i,j,k+1,+1));
                return rate;
            };
            nextK[id]=std::max(k_[id]+dt*(KOmegaSST::src_k(pk,k_[id],omega_[id])
                +transport(k_,KOmegaSST::sigma_k(f1))),KOmegaSST::k_min);
            nextOmega[id]=std::max(omega_[id]+dt*(KOmegaSST::src_omega(
                pk,nut,k_[id],omega_[id],f1,CDkw)
                +transport(omega_,KOmegaSST::sigma_w(f1))),KOmegaSST::om_min);

            double stressResidual[3]={0.0,0.0,0.0};
            double stressEnergyResidual=0.0;
            addStressFace(d,i,j,k,i-1,j,k,mesh_.axialAreaVector(i,j,k,-1),
                          points,nutField,stressResidual,stressEnergyResidual);
            addStressFace(d,i,j,k,i+1,j,k,mesh_.axialAreaVector(i+1,j,k,+1),
                          points,nutField,stressResidual,stressEnergyResidual);
            addStressFace(d,i,j,k,i,j-1,k,mesh_.radialAreaVector(i,j,k,-1),
                          points,nutField,stressResidual,stressEnergyResidual);
            addStressFace(d,i,j,k,i,j+1,k,mesh_.radialAreaVector(i,j+1,k,+1),
                          points,nutField,stressResidual,stressEnergyResidual);
            addStressFace(d,i,j,k,i,j,k-1,mesh_.azimuthalAreaVector(i,j,k,-1),
                          points,nutField,stressResidual,stressEnergyResidual);
            addStressFace(d,i,j,k,i,j,k+1,mesh_.azimuthalAreaVector(i,j,k+1,+1),
                          points,nutField,stressResidual,stressEnergyResidual);
            const double factor=dt/mesh_.volume(i,j,k);
            for(int m=0;m<3;++m) d.Unew(i,j,k).rhoU[m]+=factor*stressResidual[m];
            d.Unew(i,j,k).rhoE+=factor*stressEnergyResidual;
        }
        k_.swap(nextK); omega_.swap(nextOmega);
        applyFieldBoundaries(d,points);
    }

    void localStableDt(
        const Domain& d,double /*cfl*/,std::vector<double>& dt) const override
    {
        const std::vector<TransportProps::Point> points=transportField(d);
        dt.assign(static_cast<std::size_t>(d.N()),1e30);
        #pragma omp parallel for collapse(3) schedule(static)
        for(int i=1;i<d.Nx()-1;++i)
        for(int j=1;j<d.Ny()-1;++j)
        for(int k=1;k<d.Nz()-1;++k)
        {
            const int id=index(i,j,k);
            const double kval=kValue(id);
            const double omega=omegaValue(id);
            const double s=strain(d,i,j,k);
            const double f2=KOmegaSST::F2(
                kval,omega,wallDistance(i,j,k),points[id].nu);
            const double nut=KOmegaSST::nu_t(kval,omega,s,f2);
            const double dx=mesh_.xFace(i+1)-mesh_.xFace(i);
            const double dr=mesh_.wallCenter(i)/mesh_.Nr();
            const double ds=std::max(mesh_.radiusCenter(i,j)*mesh_.sectorAngle()
                /mesh_.Ntheta(),1e-7);
            const double h=std::min(dx,std::min(dr,ds));
            dt[static_cast<std::size_t>(id)] = safety_*h*h
                /(2*(points[id].nu+nut)+1e-30);
        }
    }

    double maxStableDt(const Domain& d,double cfl) const override
    {
        std::vector<double> dt;
        localStableDt(d,cfl,dt);
        return *std::min_element(dt.begin(),dt.end());
    }

    std::string name() const override
    {
        return mode_==TurbulenceMode::WALE ? "Sector WALE LES" : "Sector k-omega SST";
    }
};

} // End namespace Foam
#endif
