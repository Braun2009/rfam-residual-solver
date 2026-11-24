#ifndef SectorDiffusionOperator_H
#define SectorDiffusionOperator_H

#include "Operator.H"
#include "SectorMeshMetrics.H"
#include "SectorWallHeatFlux.H"
#include "SectorTurbulenceOperator.H"
#include "SpeciesRegistry.H"
#include "TransportProps.H"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace Foam
{

class SectorDiffusionOperator : public Operator
{
    using Tensor = double[3][3];
    const SectorMeshMetrics& mesh_;
    const SpeciesRegistry& registry_;
    const SectorTurbulenceOperator* turbulence_;
    double cflDiff_;
    bool diagnosticMode_ = false;
    static constexpr double turbulentPrandtl_=0.9;
    static constexpr double turbulentSchmidt_=0.9;

    int index(int i,int j,int k) const
    { return k*mesh_.Nx()*mesh_.Nr()+j*mesh_.Nx()+i; }

    void validateDomainLayout(const Domain& d) const
    {
        if(mesh_.Nx()!=d.Nx() || mesh_.Nr()!=d.Ny()
           ||mesh_.Ntheta()!=d.Nz())
            throw std::runtime_error(
                "sector diffusion transport cell -1 status=INVALID_INPUT "
                "detail=DOMAIN_MESH_MISMATCH");
        if(d.Nx()<4 || d.Ny()<4 || d.Nz()<4)
            throw std::runtime_error(
                "sector diffusion transport cell -1 status=INVALID_INPUT "
                "detail=INSUFFICIENT_HALO_DIMENSIONS");
        if(turbulence_&&!turbulence_->meshMatches(d))
            throw std::runtime_error(
                "sector diffusion transport cell -1 status=INVALID_INPUT "
                "detail=TURBULENCE_DOMAIN_MESH_MISMATCH");
    }

    std::vector<TransportProps::Point> transportField(const Domain& d) const
    {
        validateDomainLayout(d);
        if(d.speciesNames().empty() ||
           d.speciesNames().size()!=static_cast<size_t>(d.nSpecies()))
            throw std::runtime_error(
                "sector diffusion transport cell -1 status=SPECIES_MISMATCH");
        std::vector<TransportProps::Point> points(static_cast<size_t>(d.N()));
        std::vector<TransportProps::Status> statuses(
            static_cast<size_t>(d.N()),TransportProps::Status::INVALID_INPUT);
        int firstInvalid=d.N();
        // Boundary halos are owned by the boundary operators and can carry
        // states outside the gas-correlation range. Preflight active cells
        // only; all halo coefficients reuse the nearest active cell. This
        // mirrors SectorTurbulenceOperator and prevents an outlet ghost from
        // aborting an otherwise valid diffusion step.
        #pragma omp parallel for collapse(3) reduction(min:firstInvalid) schedule(static)
        for(int i=1;i<d.Nx()-1;++i)
        for(int j=1;j<d.Ny()-1;++j)
        for(int k=1;k<d.Nz()-1;++k)
        {
            const int id=index(i,j,k);
            TransportProps::Result result;
            if(diagnosticMode_)
            {
                const TransportProps::DiagnosticGasResult diag=
                    TransportProps::evaluateDiagnosticGas(
                        d.speciesNames(),d.W(i,j,k),d.Q(i,j,k),
                        d.thermoModel(),registry_,true);
                result=diag.result;
            }
            else
            {
                result=TransportProps::evaluateCommittedLegacyGas(
                    d.speciesNames(),d.W(i,j,k),d.Q(i,j,k),
                    d.thermoModel(),registry_);
            }
            statuses[static_cast<size_t>(id)]=result.status;
            if(result.success())
                points[static_cast<size_t>(id)]=result.point;
            else
                firstInvalid=std::min(firstInvalid,id);
        }
        if(firstInvalid!=d.N())
            throw std::runtime_error(
                "sector diffusion transport cell "+std::to_string(firstInvalid)
                +" status="+TransportProps::statusName(
                    statuses[static_cast<size_t>(firstInvalid)]));
        // Preserve a valid physical wall ghost for the wall heat-flux law.
        // If its state is outside the committed gas closure, use the adjacent
        // active point; outlet, axis, and periodic halos always use active
        // coefficients and are never evaluated as transport states.
        for(int i=1;i<d.Nx()-1;++i)
        for(int k=1;k<d.Nz()-1;++k)
        {
            const int wallId=index(i,d.Ny()-1,k);
            TransportProps::Result result;
            if(diagnosticMode_)
            {
                const TransportProps::DiagnosticGasResult diag=
                    TransportProps::evaluateDiagnosticGas(
                        d.speciesNames(),d.W(i,d.Ny()-1,k),
                        d.Q(i,d.Ny()-1,k),d.thermoModel(),registry_,true);
                result=diag.result;
            }
            else
            {
                result=TransportProps::evaluateCommittedLegacyGas(
                    d.speciesNames(),d.W(i,d.Ny()-1,k),
                    d.Q(i,d.Ny()-1,k),d.thermoModel(),registry_);
            }
            if(result.success())
                points[static_cast<size_t>(wallId)]=result.point;
            else
                points[static_cast<size_t>(wallId)]=
                    points[static_cast<size_t>(index(i,d.Ny()-2,k))];
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
                points[static_cast<size_t>(id)]=
                    points[static_cast<size_t>(activeId)];
        }
        return points;
    }

    std::vector<double> turbulentViscosityField(
        const Domain& d,
        const std::vector<TransportProps::Point>& points) const
    {
        std::vector<double> nut(static_cast<size_t>(d.N()),0.0);
        if(!turbulence_) return nut;

        int firstInvalid=d.N();
        #pragma omp parallel for collapse(3) reduction(min:firstInvalid) schedule(static)
        for(int i=0;i<d.Nx();++i)
        for(int j=0;j<d.Ny();++j)
        for(int k=0;k<d.Nz();++k)
        {
            const int id=index(i,j,k);
            const double value=turbulence_->turbulentViscosity(
                i,j,k,d,points[static_cast<size_t>(id)].nu);
            nut[static_cast<size_t>(id)]=value;
            if(!std::isfinite(value)||value<0.0)
                firstInvalid=std::min(firstInvalid,id);
        }
        if(firstInvalid!=d.N())
            throw std::runtime_error(
                "sector diffusion turbulent viscosity cell "
                +std::to_string(firstInvalid)+" is invalid");
        return nut;
    }

    std::vector<SectorWallHeatFluxResult> wallHeatFluxField(
        const Domain& d,
        const std::vector<TransportProps::Point>& points,
        const std::vector<double>& nut) const
    {
        std::vector<SectorWallHeatFluxResult> wall(
            static_cast<size_t>(d.Nx()*d.Nz()));
        const int cellJ=d.Ny()-2,ghostJ=d.Ny()-1;
        for(int k=0;k<d.Nz();++k)
        for(int i=0;i<d.Nx();++i)
        {
            const int cellId=index(i,cellJ,k);
            const int ghostId=index(i,ghostJ,k);
            const double kval=turbulence_
                ?turbulence_->turbulentKValue(i,cellJ,k):0.0;
            const double omega=turbulence_
                ?turbulence_->turbulentOmegaValue(i,cellJ,k):0.0;
            try
            {
                wall[static_cast<size_t>(k*d.Nx()+i)]=
                    SectorWallHeatFlux::evaluateCached(
                        d,mesh_,i,k,
                        points[static_cast<size_t>(cellId)],
                        points[static_cast<size_t>(ghostId)],
                        turbulence_!=nullptr,
                        nut[static_cast<size_t>(cellId)],kval,omega);
            }
            catch(const std::exception& error)
            {
                throw std::runtime_error(
                    "sector wall heat-flux preflight i="+std::to_string(i)
                    +" k="+std::to_string(k)+" detail="+error.what());
            }
        }
        return wall;
    }

    static double magnitude(const FaceAreaVector& a)
    { return std::sqrt(a.x*a.x+a.y*a.y+a.z*a.z); }

    double normalDistance(int i,int j,int k,int ni,int nj,int nk,
                          const FaceAreaVector& area) const
    {
        const FaceAreaVector a=mesh_.cellCenter(i,j,k);
        const FaceAreaVector b=mesh_.cellCenter(ni,nj,nk);
        const double mag=magnitude(area);
        return std::max(((b.x-a.x)*area.x+(b.y-a.y)*area.y+(b.z-a.z)*area.z)/mag,1e-12);
    }

    static FaceAreaVector subtract(
        const FaceAreaVector& a,const FaceAreaVector& b)
    { return {a.x-b.x,a.y-b.y,a.z-b.z}; }

    static void solveTranspose(
        const FaceAreaVector basis[3],const double rhs[3],double result[3])
    {
        const double a00=basis[0].x,a01=basis[0].y,a02=basis[0].z;
        const double a10=basis[1].x,a11=basis[1].y,a12=basis[1].z;
        const double a20=basis[2].x,a21=basis[2].y,a22=basis[2].z;
        const double det=a00*(a11*a22-a12*a21)-a01*(a10*a22-a12*a20)
                        +a02*(a10*a21-a11*a20);
        if(std::abs(det)<1e-30)
        { result[0]=result[1]=result[2]=0.0; return; }
        result[0]=(rhs[0]*(a11*a22-a12*a21)-a01*(rhs[1]*a22-a12*rhs[2])
                  +a02*(rhs[1]*a21-a11*rhs[2]))/det;
        result[1]=(a00*(rhs[1]*a22-a12*rhs[2])-rhs[0]*(a10*a22-a12*a20)
                  +a02*(a10*rhs[2]-rhs[1]*a20))/det;
        result[2]=(a00*(a11*rhs[2]-rhs[1]*a21)-a01*(a10*rhs[2]-rhs[1]*a20)
                  +rhs[0]*(a10*a21-a11*a20))/det;
    }

    void velocityGradient(
        const Domain& d,int i,int j,int k,Tensor& gradient) const
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

    void molecularStress(const Domain& d,int i,int j,int k,double mu,Tensor& tau) const
    {
        Tensor gradient; velocityGradient(d,i,j,k,gradient);
        const double divergence=gradient[0][0]+gradient[1][1]+gradient[2][2];
        for(int a=0;a<3;++a) for(int b=0;b<3;++b)
            tau[a][b]=mu*(gradient[a][b]+gradient[b][a]
                        -(a==b ? 2.0*divergence/3.0 : 0.0));
    }

    void addMolecularStressFace(
        const Domain& d,int i,int j,int k,int ni,int nj,int nk,
        const FaceAreaVector& area,
        const std::vector<TransportProps::Point>& points,
        double momentum[3],double& energy) const
    {
        const bool wallFace=j==d.Ny()-2&&nj==d.Ny()-1;
        if(wallFace&&turbulence_)
            return; // The SST wall function supplies the total wall traction.
        Tensor left,right;
        molecularStress(d,i,j,k,points[index(i,j,k)].mu,left);
        molecularStress(d,ni,nj,nk,points[index(ni,nj,nk)].mu,right);
        for(int m=0;m<3;++m)
        {
            const double traction=0.5*((left[m][0]+right[m][0])*area.x
                                      +(left[m][1]+right[m][1])*area.y
                                      +(left[m][2]+right[m][2])*area.z);
            momentum[m]+=traction;
            const double faceVelocity=wallFace?d.W(ni,nj,nk).U[m]
                :0.5*(d.W(i,j,k).U[m]+d.W(ni,nj,nk).U[m]);
            energy+=faceVelocity*traction;
        }
    }

    void addFace(const Domain& d,int i,int j,int k,int ni,int nj,int nk,
                 const FaceAreaVector& area,std::vector<double>& species,
                 double momentum[3],double& energy,
                 const std::vector<TransportProps::Point>& points,
                 const std::vector<double>& nut,
                 const std::vector<SectorWallHeatFluxResult>& wall) const
    {
        const PrimitiveVars& a=d.W(i,j,k);
        const PrimitiveVars& b=d.W(ni,nj,nk);
        const TransportProps::Point& transportA=points[index(i,j,k)];
        const TransportProps::Point& transportB=points[index(ni,nj,nk)];
        const double areaMagnitude=magnitude(area);
        if(areaMagnitude<=1e-30)
        {
            std::fill(species.begin(),species.end(),0.0);
            return;
        }
        const bool wallFace=j==d.Ny()-2&&nj==d.Ny()-1;
        if(wallFace)
        {
            std::fill(species.begin(),species.end(),0.0);
            const SectorWallHeatFluxResult& q=
                wall[static_cast<size_t>(k*d.Nx()+i)];
            energy-=q.totalFlux*areaMagnitude;
            addMolecularStressFace(
                d,i,j,k,ni,nj,nk,area,points,momentum,energy);
            return;
        }
        const double scale=areaMagnitude/normalDistance(i,j,k,ni,nj,nk,area);
        double rhoD=0.5*(transportA.mu+transportB.mu)/TransportProps::Sc;
        if(turbulence_)
        {
            // SGS species flux (gradient-diffusion closure, Sc_t = 0.9):
            // rhoD_total = mu/Sc + rho*nu_t/Sc_t.
            const double rhoFace=0.5*(
                d.Q(i,j,k).rho+d.Q(ni,nj,nk).rho);
            const double nutFace=0.5*(
                nut[static_cast<size_t>(index(i,j,k))]
               +nut[static_cast<size_t>(index(ni,nj,nk))]);
            rhoD+=rhoFace*nutFace/turbulentSchmidt_;
        }

        double sumRaw=0.0;
        for(size_t s=0;s<species.size();++s)
        {
            species[s]=rhoD*(b.Y[s]-a.Y[s])*scale;
            sumRaw+=species[s];
        }
        double sumCorrected=0.0;
        for(size_t s=0;s<species.size();++s)
        {
            const double yFace=0.5*(a.Y[s]+b.Y[s]);
            species[s]-=yFace*sumRaw;
            sumCorrected+=species[s];
        }
        if(!species.empty()) species.back()-=sumCorrected;
        double conductivity=0.5*(transportA.kappa+transportB.kappa);
        if(turbulence_)
        {
            const double rhoFace=0.5*(
                d.Q(i,j,k).rho+d.Q(ni,nj,nk).rho);
            const double cpFace=0.5*(transportA.Cp+transportB.Cp);
            const double nutFace=0.5*(
                nut[static_cast<size_t>(index(i,j,k))]
               +nut[static_cast<size_t>(index(ni,nj,nk))]);
            conductivity+=rhoFace*cpFace*nutFace/turbulentPrandtl_;
        }
        energy+=conductivity*(b.T-a.T)*scale;
        addMolecularStressFace(
            d,i,j,k,ni,nj,nk,area,points,momentum,energy);
    }

public:
    SectorDiffusionOperator(const SectorMeshMetrics& mesh,
                            const SpeciesRegistry& registry,
                            const SectorTurbulenceOperator* turbulence=nullptr,
                            double cflDiff=0.3,
                            bool diagnosticMode=false)
        :mesh_(mesh),registry_(registry),turbulence_(turbulence),cflDiff_(cflDiff),
         diagnosticMode_(diagnosticMode) {}

    void apply(Domain& d,double dt) const override
    {
        const std::vector<TransportProps::Point> points=transportField(d);
        const std::vector<double> nut=turbulentViscosityField(d,points);
        const std::vector<SectorWallHeatFluxResult> wall=
            wallHeatFluxField(d,points,nut);
        #pragma omp parallel
        {
        std::vector<double> face(d.nSpecies()),sum(d.nSpecies());
        #pragma omp for collapse(3) schedule(static)
        for(int i=1;i<d.Nx()-1;++i)
        for(int j=1;j<d.Ny()-1;++j)
        for(int k=1;k<d.Nz()-1;++k)
        {
            std::fill(sum.begin(),sum.end(),0.0);
            double momentum[3]={0.0,0.0,0.0};
            double energy=0.0;
            auto accumulate=[&](int ni,int nj,int nk,const FaceAreaVector& area)
            {
                double faceEnergy=0.0;
                double faceMomentum[3]={0.0,0.0,0.0};
                addFace(d,i,j,k,ni,nj,nk,area,face,
                        faceMomentum,faceEnergy,points,nut,wall);
                for(int s=0;s<d.nSpecies();++s) sum[s]+=face[s];
                for(int m=0;m<3;++m) momentum[m]+=faceMomentum[m];
                energy+=faceEnergy;
            };
            accumulate(i-1,j,k,mesh_.axialAreaVector(i,j,k,-1));
            accumulate(i+1,j,k,mesh_.axialAreaVector(i+1,j,k,+1));
            accumulate(i,j-1,k,mesh_.radialAreaVector(i,j,k,-1));
            accumulate(i,j+1,k,mesh_.radialAreaVector(i,j+1,k,+1));
            accumulate(i,j,k-1,mesh_.azimuthalAreaVector(i,j,k,-1));
            accumulate(i,j,k+1,mesh_.azimuthalAreaVector(i,j,k+1,+1));
            const double factor=dt/mesh_.volume(i,j,k);
            for(int s=0;s<d.nSpecies();++s) d.Unew(i,j,k).rhoY[s]+=factor*sum[s];
            for(int m=0;m<3;++m) d.Unew(i,j,k).rhoU[m]+=factor*momentum[m];
            d.Unew(i,j,k).rhoE+=factor*energy;
        }
        }
    }

    void localStableDt(
        const Domain& d,double /*cfl*/,std::vector<double>& dt) const override
    {
        const std::vector<TransportProps::Point> points=transportField(d);
        const std::vector<double> nut=turbulentViscosityField(d,points);
        dt.assign(static_cast<std::size_t>(d.N()),1e30);
        #pragma omp parallel for collapse(3) schedule(static)
        for(int i=1;i<d.Nx()-1;++i)
        for(int j=1;j<d.Ny()-1;++j)
        for(int k=1;k<d.Nz()-1;++k)
        {
            const int id=index(i,j,k);
            const double rho=d.Q(i,j,k).rho;
            const TransportProps::Point& transport=points[id];
            double stateDiff=std::max(transport.diffusivity,transport.nu);
            stateDiff=std::max(
                stateDiff,transport.kappa/(rho*transport.Cv));
            if(turbulence_)
            {
                stateDiff=std::max(stateDiff,
                    (transport.Cp/transport.Cv)
                    *nut[static_cast<size_t>(id)]/turbulentPrandtl_);
            }
            double geometry=0.0;
            auto addGeometry=[&](int ni,int nj,int nk,const FaceAreaVector& area)
            {
                const double areaMagnitude=magnitude(area);
                if(areaMagnitude>1e-30)
                    geometry+=areaMagnitude/normalDistance(i,j,k,ni,nj,nk,area);
            };
            addGeometry(i-1,j,k,mesh_.axialAreaVector(i,j,k,-1));
            addGeometry(i+1,j,k,mesh_.axialAreaVector(i+1,j,k,+1));
            addGeometry(i,j-1,k,mesh_.radialAreaVector(i,j,k,-1));
            addGeometry(i,j+1,k,mesh_.radialAreaVector(i,j+1,k,+1));
            addGeometry(i,j,k-1,mesh_.azimuthalAreaVector(i,j,k,-1));
            addGeometry(i,j,k+1,mesh_.azimuthalAreaVector(i,j,k+1,+1));
            dt[static_cast<std::size_t>(id)] = cflDiff_*mesh_.volume(i,j,k)
                /(stateDiff*geometry+1e-30);
        }
    }

    double maxStableDt(const Domain& d,double cfl) const override
    {
        std::vector<double> dt;
        localStableDt(d,cfl,dt);
        return *std::min_element(dt.begin(),dt.end());
    }

    std::string name() const override
    { return "Sector mapped species/thermal diffusion"; }
};

} // End namespace Foam
#endif
