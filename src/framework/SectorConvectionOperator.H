#ifndef SectorConvectionOperator_H
#define SectorConvectionOperator_H

#include "Operator.H"
#include "AcousticStateAdmission.H"
#include "ExperimentalCflAcousticScaler.H"
#include "FluxScheme.H"
#include "FrozenEulerFluxJacobian.H"
#include "SectorMeshMetrics.H"
#include <cmath>
#include <memory>
#include <vector>

namespace Foam
{

class SectorConvectionOperator : public Operator
{
    struct ThermoCache { double rho,sound,enthalpy; };
    std::unique_ptr<FluxScheme> flux_;
    const SectorMeshMetrics& mesh_;
    ExperimentalCflAcousticScaler scaler_;
    bool useFrozenEulerBlock_=false;
    bool useHydroDiagonalBlock_=false;

    static std::vector<double> unitNormal(const FaceAreaVector& area)
    {
        const double magnitude = std::sqrt(area.x*area.x + area.y*area.y + area.z*area.z);
        return {area.x/magnitude, area.y/magnitude, area.z/magnitude};
    }

    static double magnitude(const FaceAreaVector& area)
    {
        return std::sqrt(area.x*area.x + area.y*area.y + area.z*area.z);
    }

    double waveSpeed(
        const PrimitiveVars& w,
        const thermo::EquilibriumState& q) const
    {
        return std::sqrt(w.U[0]*w.U[0]+w.U[1]*w.U[1]+w.U[2]*w.U[2])
            +q.acoustic.soundSpeed();
    }

    static void requireAdmittedAcoustics(
        const Domain& d,
        const char* operation)
    {
        for (int cell=0; cell<d.N(); ++cell)
            if (!solverConsumerAdmitsAcousticState(d.Q(cell)))
                throw std::runtime_error(
                    std::string("SectorConvectionOperator: ") + operation
                    + " rejected acoustic state at cell "
                    + std::to_string(cell));
    }

    void addFace(const PrimitiveVars& inside, const PrimitiveVars& outside,
                 const ThermoCache& insideThermo,const ThermoCache& outsideThermo,
                 const FaceAreaVector& area, ConservedVars& sum) const
    {
        const double scale = magnitude(area);
        if (scale <= 1e-30) return;
        ConservedVars face(sum.rhoY.size());
        flux_->computeFluxCached(inside,outside,unitNormal(area),face,
            insideThermo.rho,insideThermo.sound,insideThermo.enthalpy,
            outsideThermo.rho,outsideThermo.sound,outsideThermo.enthalpy);
        sum.rho += scale * face.rho;
        sum.rhoE += scale * face.rhoE;
        for (int m=0; m<3; ++m) sum.rhoU[m] += scale * face.rhoU[m];
        for (size_t s=0; s<sum.rhoY.size(); ++s) sum.rhoY[s] += scale * face.rhoY[s];
        sum.rhoZ += scale * face.rhoZ;
        sum.rhoC += scale * face.rhoC;
    }

public:
    SectorConvectionOperator(std::unique_ptr<FluxScheme> flux,
                             const SectorMeshMetrics& mesh)
        : flux_(std::move(flux)), mesh_(mesh) {}

    void apply(Domain& d, double dt) const override
    {
        requireAdmittedAcoustics(d,"apply");
        const int nx=d.Nx(), nr=d.Ny(), nt=d.Nz(), nSpecies=d.nSpecies();
        std::vector<ThermoCache> thermo(nx*nr*nt);
        auto id=[&](int i,int j,int k){return k*nx*nr+j*nx+i;};
        #pragma omp parallel for collapse(3) schedule(static)
        for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
        {
            const thermo::EquilibriumState& q=d.Q(i,j,k);
            thermo[id(i,j,k)]={q.rho,q.acoustic.soundSpeed(),q.h};
        }
        #pragma omp parallel for collapse(3) schedule(static)
        for (int i=1; i<nx-1; ++i)
        for (int j=1; j<nr-1; ++j)
        for (int k=1; k<nt-1; ++k)
        {
            ConservedVars residual(nSpecies);
            const ThermoCache& center=thermo[id(i,j,k)];
            addFace(d.W(i,j,k),d.W(i-1,j,k),center,thermo[id(i-1,j,k)],mesh_.axialAreaVector(i,j,k,-1),residual);
            addFace(d.W(i,j,k),d.W(i+1,j,k),center,thermo[id(i+1,j,k)],mesh_.axialAreaVector(i+1,j,k,+1),residual);
            addFace(d.W(i,j,k),d.W(i,j-1,k),center,thermo[id(i,j-1,k)],mesh_.radialAreaVector(i,j,k,-1),residual);
            addFace(d.W(i,j,k),d.W(i,j+1,k),center,thermo[id(i,j+1,k)],mesh_.radialAreaVector(i,j+1,k,+1),residual);
            addFace(d.W(i,j,k),d.W(i,j,k-1),center,thermo[id(i,j,k-1)],mesh_.azimuthalAreaVector(i,j,k,-1),residual);
            addFace(d.W(i,j,k),d.W(i,j,k+1),center,thermo[id(i,j,k+1)],mesh_.azimuthalAreaVector(i,j,k+1,+1),residual);

            const double factor = dt / mesh_.volume(i,j,k);
            const ConservedVars& old = d.U(i,j,k);
            ConservedVars& next = d.Unew(i,j,k);
            next.rho = old.rho - factor*residual.rho;
            next.rhoE = old.rhoE - factor*residual.rhoE;
            for (int m=0;m<3;++m) next.rhoU[m]=old.rhoU[m]-factor*residual.rhoU[m];
            for (int s=0;s<nSpecies;++s) next.rhoY[s]=old.rhoY[s]-factor*residual.rhoY[s];
            next.rhoZ=old.rhoZ-factor*residual.rhoZ;
            next.rhoC=old.rhoC-factor*residual.rhoC;
        }
    }

    void setCflAcousticScaling(bool enable)
    {
        scaler_.setEnabled(enable);
    }

    void setFrozenEulerBlockCoupling(bool enable)
    {
        useFrozenEulerBlock_=enable;
    }

    void setHydroDiagonalBlockCoupling(bool enable)
    {
        useHydroDiagonalBlock_=enable;
    }

    void localStableDt(
        const Domain& d,double cfl,std::vector<double>& dt) const override
    {
        requireAdmittedAcoustics(d,"CFL scan");
        const int nx=d.Nx(),ny=d.Ny(),nz=d.Nz();
        std::vector<double> speeds(nx*ny*nz);
        dt.assign(static_cast<std::size_t>(d.N()),1e30);
        #pragma omp parallel for collapse(3) schedule(static)
        for(int i=0;i<nx;++i) for(int j=0;j<ny;++j) for(int k=0;k<nz;++k)
            speeds[k*nx*ny+j*nx+i]=waveSpeed(d.W(i,j,k),d.Q(i,j,k));
        auto speed=[&](int i,int j,int k){return speeds[k*nx*ny+j*nx+i];};
        #pragma omp parallel for collapse(3) schedule(static)
        for(int i=1;i<d.Nx()-1;++i)
        for(int j=1;j<d.Ny()-1;++j)
        for(int k=1;k<d.Nz()-1;++k)
        {
            const double centerSpeed=speed(i,j,k);
            double spectralRadius=0.0;
            auto addFaceSpeed=[&](int ni,int nj,int nk,const FaceAreaVector& area)
            {
                spectralRadius+=magnitude(area)
                    *std::max(centerSpeed,speed(ni,nj,nk));
            };
            addFaceSpeed(i-1,j,k,mesh_.axialAreaVector(i,j,k,-1));
            addFaceSpeed(i+1,j,k,mesh_.axialAreaVector(i+1,j,k,+1));
            addFaceSpeed(i,j-1,k,mesh_.radialAreaVector(i,j,k,-1));
            addFaceSpeed(i,j+1,k,mesh_.radialAreaVector(i,j+1,k,+1));
            addFaceSpeed(i,j,k-1,mesh_.azimuthalAreaVector(i,j,k,-1));
            addFaceSpeed(i,j,k+1,mesh_.azimuthalAreaVector(i,j,k+1,+1));
            dt[static_cast<std::size_t>(k*nx*ny+j*nx+i)] =
                cfl*mesh_.volume(i,j,k)/(spectralRadius+1e-30);
        }
    }

    void addStructuredImplicitStencil(
        const Domain& d,StructuredImplicitStencil& stencil) const override
    {
        requireAdmittedAcoustics(d,"implicit stencil");
        const int nx=d.Nx(),ny=d.Ny(),nz=d.Nz();
        stencil.ensureSize(static_cast<std::size_t>(d.N()));
        std::vector<double> speeds(static_cast<std::size_t>(d.N()),0.0);
        #pragma omp parallel for collapse(3) schedule(static)
        for(int i=0;i<nx;++i) for(int j=0;j<ny;++j) for(int k=0;k<nz;++k)
            speeds[static_cast<std::size_t>(k*nx*ny+j*nx+i)]=
                waveSpeed(d.W(i,j,k),d.Q(i,j,k));
        auto speed=[&](int i,int j,int k)
        { return speeds[static_cast<std::size_t>(k*nx*ny+j*nx+i)]; };
        #pragma omp parallel for collapse(3) schedule(static)
        for(int i=1;i<nx-1;++i)
        for(int j=1;j<ny-1;++j)
        for(int k=1;k<nz-1;++k)
        {
            const int cell=k*nx*ny+j*nx+i;
            const double inverseVolume=1.0/mesh_.volume(i,j,k);
            const double center=speed(i,j,k);
            auto add=[&](
                int face,int ni,int nj,int nk,const FaceAreaVector& area)
            {
                const double coupling=magnitude(area)
                    *std::max(center,speed(ni,nj,nk))*inverseVolume;
                stencil.face[face][static_cast<std::size_t>(cell)]+=coupling;
                stencil.diagonal[static_cast<std::size_t>(cell)]+=coupling;
            };
            add(StructuredImplicitStencil::I_MINUS,i-1,j,k,
                mesh_.axialAreaVector(i,j,k,-1));
            add(StructuredImplicitStencil::I_PLUS,i+1,j,k,
                mesh_.axialAreaVector(i+1,j,k,+1));
            add(StructuredImplicitStencil::J_MINUS,i,j-1,k,
                mesh_.radialAreaVector(i,j,k,-1));
            add(StructuredImplicitStencil::J_PLUS,i,j+1,k,
                mesh_.radialAreaVector(i,j+1,k,+1));
            add(StructuredImplicitStencil::K_MINUS,i,j,k-1,
                mesh_.azimuthalAreaVector(i,j,k,-1));
            add(StructuredImplicitStencil::K_PLUS,i,j,k+1,
                mesh_.azimuthalAreaVector(i,j,k+1,+1));
        }
    }

    bool supportsStructuredBlockCoupling() const override
    { return useFrozenEulerBlock_; }

    bool supportsStructuredDiagonalBlock() const override
    { return useHydroDiagonalBlock_; }

    bool addStructuredImplicitDiagonalBlock(
        const Domain& d,int cell,
        StructuredImplicitBlock& contribution,
        double& scalarContribution) const override
    {
        if(!useHydroDiagonalBlock_) return false;
        if(cell<0||cell>=d.N()) return false;
        const int nx=d.Nx(),ny=d.Ny();
        const int i=cell%nx;
        const int j=(cell/nx)%ny;
        const int k=cell/(nx*ny);
        if(i<=0||i>=d.Nx()-1||j<=0||j>=d.Ny()-1
            ||k<=0||k>=d.Nz()-1) return false;
        contribution.setScalarIdentity(0.0);
        scalarContribution=0.0;
        const double inverseVolume=1.0/mesh_.volume(i,j,k);
        const double centerWave=waveSpeed(d.W(i,j,k),d.Q(i,j,k));
        auto addFace=[&](int ni,int nj,int nk,const FaceAreaVector& area)
        {
            const double areaMagnitude=magnitude(area);
            if(areaMagnitude<=1.0e-30) return true;
            const double wave=std::max(
                centerWave,waveSpeed(d.W(ni,nj,nk),d.Q(ni,nj,nk)));
            const double scalar=areaMagnitude*wave*inverseVolume;
            if(!std::isfinite(scalar)||scalar<=0.0) return false;
            scalarContribution+=scalar;
            contribution.add(0,0,0.5*scalar);
            contribution.add(1,1,0.5*scalar);
            contribution.add(2,2,0.5*scalar);
            contribution.add(3,3,0.5*scalar);
            contribution.add(4,4,0.5*scalar);
            const std::array<double,3> normal{
                area.x/areaMagnitude,area.y/areaMagnitude,
                area.z/areaMagnitude};
            // Add a bounded frozen-Euler pressure/velocity/energy coupling
            // around the scalar Rusanov diagonal. The scalar wave bound stays
            // on every diagonal, so this block is never the sole positivity
            // safeguard for a noisy real-fluid state.
            for(int column=0;column<StructuredImplicitBlock::HYDRO_COMPONENTS;
                ++column)
            {
                ConservedVars perturbation(d.nSpecies());
                if(column==0) perturbation.rho=1.0;
                else if(column<=3)
                    perturbation.rhoU[static_cast<std::size_t>(column-1)]=1.0;
                else perturbation.rhoE=1.0;
                ConservedVars product(d.nSpecies());
                if(!FrozenEulerFluxJacobian::applyAbsolute(
                        d.W(i,j,k),d.Q(i,j,k),perturbation,
                        normal,product)) return false;
                const double factor=0.5*areaMagnitude*inverseVolume;
                contribution.add(0,column,factor*product.rho);
                contribution.add(1,column,factor*product.rhoU[0]);
                contribution.add(2,column,factor*product.rhoU[1]);
                contribution.add(3,column,factor*product.rhoU[2]);
                contribution.add(4,column,factor*product.rhoE);
            }
            return true;
        };
        if(!addFace(i-1,j,k,mesh_.axialAreaVector(i,j,k,-1))
            ||!addFace(i+1,j,k,mesh_.axialAreaVector(i+1,j,k,+1))
            ||!addFace(i,j-1,k,mesh_.radialAreaVector(i,j,k,-1))
            ||!addFace(i,j+1,k,mesh_.radialAreaVector(i,j+1,k,+1))
            ||!addFace(i,j,k-1,mesh_.azimuthalAreaVector(i,j,k,-1))
            ||!addFace(i,j,k+1,mesh_.azimuthalAreaVector(i,j,k+1,+1)))
            return false;
        return contribution.finite()&&std::isfinite(scalarContribution)
            &&scalarContribution>0.0;
    }

    bool addStructuredImplicitNeighborContribution(
        const Domain& d,int cell,int neighbor,int face,
        const ConservedVars& neighborCorrection,
        ConservedVars& contribution) const override
    {
        const int nx=d.Nx(),ny=d.Ny();
        const int i=cell%nx;
        const int j=(cell/nx)%ny;
        const int k=cell/(nx*ny);
        FaceAreaVector area{0.0,0.0,0.0};
        switch(face)
        {
            case StructuredImplicitStencil::I_MINUS:
                area=mesh_.axialAreaVector(i,j,k,-1); break;
            case StructuredImplicitStencil::I_PLUS:
                area=mesh_.axialAreaVector(i+1,j,k,+1); break;
            case StructuredImplicitStencil::J_MINUS:
                area=mesh_.radialAreaVector(i,j,k,-1); break;
            case StructuredImplicitStencil::J_PLUS:
                area=mesh_.radialAreaVector(i,j+1,k,+1); break;
            case StructuredImplicitStencil::K_MINUS:
                area=mesh_.azimuthalAreaVector(i,j,k,-1); break;
            case StructuredImplicitStencil::K_PLUS:
                area=mesh_.azimuthalAreaVector(i,j,k+1,+1); break;
            default:return false;
        }
        const double areaMagnitude=magnitude(area);
        if(areaMagnitude<=1.0e-30) return true;
        const std::array<double,3> normal{
            area.x/areaMagnitude,area.y/areaMagnitude,area.z/areaMagnitude};
        ConservedVars fluxProduct(d.nSpecies());
        if(!FrozenEulerFluxJacobian::apply(
                d.W(neighbor),d.Q(neighbor),neighborCorrection,
                normal,fluxProduct))
            return false;
        const double factor=0.5*areaMagnitude/mesh_.volume(i,j,k);
        const double wave=std::max(
            waveSpeed(d.W(cell),d.Q(cell)),
            waveSpeed(d.W(neighbor),d.Q(neighbor)));
        contribution.rho+=factor*(-wave*neighborCorrection.rho
            -fluxProduct.rho);
        contribution.rhoE+=factor*(-wave*neighborCorrection.rhoE
            -fluxProduct.rhoE);
        for(std::size_t component=0;
            component<contribution.rhoU.size();++component)
            contribution.rhoU[component]+=factor
                *(-wave*neighborCorrection.rhoU[component]
                  -fluxProduct.rhoU[component]);
        for(std::size_t component=0;
            component<contribution.rhoY.size();++component)
            contribution.rhoY[component]+=factor
                *(-wave*neighborCorrection.rhoY[component]
                  -fluxProduct.rhoY[component]);
        contribution.rhoZ+=factor*(-wave*neighborCorrection.rhoZ
            -fluxProduct.rhoZ);
        contribution.rhoC+=factor*(-wave*neighborCorrection.rhoC
            -fluxProduct.rhoC);
        return true;
    }

    double maxStableDt(const Domain& d, double cfl) const override
    {
        std::vector<double> dt;
        localStableDt(d,cfl,dt);
        return *std::min_element(dt.begin(),dt.end());
    }

    std::string name() const override { return "Sector mapped convection"; }
};

} // End namespace Foam

#endif
