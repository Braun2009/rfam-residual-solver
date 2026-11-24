#ifndef SectorInjectionBoundary_H
#define SectorInjectionBoundary_H

#include "Domain.H"
#include "FluxScheme.H"
#include "SectorMeshMetrics.H"
#include "RealFluidState.H"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Foam
{

struct SectorInjectionReport
{
    double torchMassFlow,loxMassFlow,fuelMassFlow;
    double torchAxialVelocity,loxAxialVelocity,fuelAxialVelocity;
    int torchFaces,loxFaces,fuelFaces;
};

struct SectorInjectionFluxReport
{
    double background=0.0,torch=0.0,lox=0.0,fuel=0.0;
    double total() const { return background+torch+lox+fuel; }
};

class SectorInjectionBoundary
{
    enum class Patch { Background,Torch,Lox,Fuel };
    const SectorMeshMetrics& mesh_;
    double pressure_;
    std::vector<std::string> species_;
    int oxygenIndex_,fuelIndex_;
    std::vector<double> backgroundY_;
    std::vector<Patch> patches_;
    SectorInjectionReport report_{};
    const GERGResidualTable* gergTable_ = nullptr;

    int index(int j,int k) const { return k*mesh_.Nr()+j; }

    PrimitiveVars state(double p,double T,double yO2,double yFuel) const
    {
        PrimitiveVars w(static_cast<int>(species_.size())); w.p=p; w.T=T;
        w.Y[oxygenIndex_]=yO2; w.Y[fuelIndex_]=yFuel;
        w.U={0.0,0.0,0.0}; return w;
    }

    double density(const PrimitiveVars& w) const
    {
        RealFluidState thermo;
        thermo.setSpecies(species_);
        if (gergTable_) thermo.setGERGTable(gergTable_);
        if(thermo.updatePTY(w.p,w.T,w.Y)!=RealFluidState::Status::SUCCESS)
            throw std::runtime_error("invalid sector injection thermodynamic state");
        return thermo.rho();
    }

    // Register the GERG-2008 residual table so dense-cryogenic boundary
    // states (e.g. 130 K LOX) use the GERG backend instead of the PR
    // two-phase-dome fail-close.

    thermo::EquilibriumState initializeBoundaryThermo(
        const Domain& d,
        const PrimitiveVars& w) const
    {
        const thermo::EquilibriumResult initialized =
            d.thermoModel().initializePTY(w.p,w.T,w.Y);
        if(!initialized.success())
            throw std::runtime_error(
                "sector injection requires a complete model-consistent "
                "boundary equilibrium state");
        return *initialized.state;
    }

    static void requireLegacySinglePhaseBoundaryModel(const Domain& d)
    {
        if(!d.thermoModel().supportsCommittedLegacyGasTransport())
            throw std::runtime_error(
                "SectorInjectionBoundary phase-aware boundary states are not "
                "implemented; requires the committed legacy single-phase model");
    }

    PrimitiveVars injectionState(Patch patch,const FaceAreaVector& c,
                                  double velocityScale) const
    {
        if(patch==Patch::Torch)
        {
            PrimitiveVars w=state(pressure_,2400.0,0.784,0.216);
            w.U[0]=velocityScale*report_.torchAxialVelocity;
            return w;
        }

        // nearest pintle center (full annulus: 3 pintles; sector: 1)
        const bool fullAnnulus=mesh_.sectorAngle()>1.5*3.14159265358979323846;
        const double boltR=0.0125;
        double pcY=0.0,pcZ=0.0,best=1e30;
        const int nP=fullAnnulus?3:1;
        for(int p=0;p<nP;++p)
        {
            const double th=(60.0+120.0*p)*3.14159265358979323846/180.0;
            const double py=boltR*std::cos(th),pz=boltR*std::sin(th);
            const double d2=(c.y-py)*(c.y-py)+(c.z-pz)*(c.z-pz);
            if(d2<best){best=d2;pcY=py;pcZ=pz;}
        }
        const double dy=c.y-pcY;
        const double dz=c.z-pcZ;
        const double inv=1.0/std::max(std::hypot(dy,dz),1e-12);
        const double erY=dy*inv,erZ=dz*inv,etY=-erZ,etZ=erY;
        if(patch==Patch::Lox)
        {
            PrimitiveVars w=state(pressure_,130.0,1.0,0.0);
            w.U[0]=velocityScale*report_.loxAxialVelocity;
            w.U[1]=w.U[0]*(0.933*erY+0.4*etY);
            w.U[2]=w.U[0]*(0.933*erZ+0.4*etZ);
            return w;
        }

        PrimitiveVars w=state(pressure_,450.0,0.0,1.0);
        w.U[0]=velocityScale*report_.fuelAxialVelocity;
        w.U[1]=w.U[0]*0.4*etY; w.U[2]=w.U[0]*0.4*etZ;
        return w;
    }

    double numericalMassFlow(
        const Domain& d,
        Patch target,
        double velocityScale,
        const thermo::EquilibriumState& boundaryThermo) const
    {
        AUSMFluxAdapter flux(species_);
        double result=0.0;
        for(int j=1;j<d.Ny()-1;++j) for(int k=1;k<d.Nz()-1;++k)
        {
            if(patches_[index(j,k)]!=target) continue;
            const PrimitiveVars ghost=injectionState(
                target,mesh_.cellCenter(0,j,k),velocityScale);
            ConservedVars face(d.nSpecies());
            flux.computeFluxFromEquilibrium(
                ghost,boundaryThermo,d.W(1,j,k),d.Q(1,j,k),
                {1.0,0.0,0.0},face);
            result+=face.rho*mesh_.axialArea(1,j,k);
        }
        return result;
    }

    double controlledVelocityScale(
        const Domain& d,
        Patch patch,
        double ramp,
        const thermo::EquilibriumState& boundaryThermo) const
    {
        double prescribed=0.0;
        if(patch==Patch::Torch) prescribed=report_.torchMassFlow;
        else if(patch==Patch::Lox) prescribed=report_.loxMassFlow;
        else prescribed=report_.fuelMassFlow;
        const double target=ramp*prescribed;

        double lower=-std::max(ramp,0.1),upper=std::max(ramp,0.1);
        double lowerFlux=numericalMassFlow(
            d,patch,lower,boundaryThermo);
        double upperFlux=numericalMassFlow(
            d,patch,upper,boundaryThermo);
        while(lowerFlux>target&&lower>-64.0)
        {
            upper=lower; upperFlux=lowerFlux;
            lower*=2.0;
            lowerFlux=numericalMassFlow(
                d,patch,lower,boundaryThermo);
        }
        while(upperFlux<target&&upper<64.0)
        {
            lower=upper; lowerFlux=upperFlux;
            upper*=2.0;
            upperFlux=numericalMassFlow(
                d,patch,upper,boundaryThermo);
        }
        if(lowerFlux>target||upperFlux<target)
            throw std::runtime_error("sector injection mass-flow target is not bracketed");
        const double tolerance=1e-4*prescribed;
        double bestScale=std::abs(lowerFlux-target)<std::abs(upperFlux-target)
            ?lower:upper;
        double bestError=std::min(std::abs(lowerFlux-target),
                                  std::abs(upperFlux-target));
        for(int iteration=0;iteration<12;++iteration)
        {
            const double width=upper-lower;
            double middle=lower+(target-lowerFlux)*width/(upperFlux-lowerFlux);
            middle=std::clamp(middle,lower+0.1*width,upper-0.1*width);
            const double middleFlux=numericalMassFlow(
                d,patch,middle,boundaryThermo);
            const double error=std::abs(middleFlux-target);
            if(error<bestError) { bestScale=middle; bestError=error; }
            if(error<=tolerance) return middle;
            if(middleFlux<target) { lower=middle; lowerFlux=middleFlux; }
            else { upper=middle; upperFlux=middleFlux; }
        }
        if(bestError>tolerance)
            throw std::runtime_error("sector injection mass-flow solve did not converge");
        return bestScale;
    }

public:
    void setGERGTable(const GERGResidualTable* table)
    { gergTable_ = table; }

    // Submerged torch tube (design: lengthened torch tube): the torch core
    // state at the tube end for an in-chamber replacement.
    PrimitiveVars torchTubePrimitive(double scale) const
    {
        return injectionState(
            Patch::Torch, mesh_.cellCenter(0,1,1), scale);
    }
    thermo::EquilibriumState torchTubeThermo(
        const Domain& d, double scale) const
    {
        const FaceAreaVector rep=mesh_.cellCenter(0,1,1);
        const PrimitiveVars w=injectionState(Patch::Torch,rep,scale);
        return initializeBoundaryThermo(d,w);
    }

    static constexpr int faceI_ = 1; // planar boundary of the structured baseline
    explicit SectorInjectionBoundary(const SectorMeshMetrics& mesh,double pressure=4e6,
                                      std::vector<std::string> species={"O2","C3H8"},
                                      std::vector<double> backgroundY={},
                                      double faceOffsetMm=0.0)
        :mesh_(mesh),pressure_(pressure),species_(std::move(species)),
         oxygenIndex_(-1),fuelIndex_(-1),
         backgroundY_(std::move(backgroundY)),
         patches_(mesh.Nr()*mesh.Ntheta(),Patch::Background)
    {
        if(!std::isfinite(faceOffsetMm) || faceOffsetMm!=0.0)
            throw std::invalid_argument(
                "interior injector offsets require a boundary-conforming fluid mesh; "
                "the structured baseline supports only faceOffsetMm=0");
        for(size_t s=0;s<species_.size();++s)
        { if(species_[s]=="O2") oxygenIndex_=static_cast<int>(s);
          if(species_[s]=="C3H8") fuelIndex_=static_cast<int>(s); }
        if(oxygenIndex_<0||fuelIndex_<0)
            throw std::invalid_argument("sector injection requires O2 and C3H8 species");
        if(backgroundY_.empty())
        {
            backgroundY_.assign(species_.size(),0.0);
            backgroundY_[oxygenIndex_]=0.784; backgroundY_[fuelIndex_]=0.216;
        }
        if(backgroundY_.size()!=species_.size())
            throw std::invalid_argument("sector injection background composition size mismatch");
        // Full-annulus (360 deg): three pintles on the bolt circle plus the
        // central torch.  Sector (120 deg): single pintle as before.
        const bool fullAnnulus=mesh_.sectorAngle()>1.5*3.14159265358979323846;
        const double boltR=0.0125;
        std::vector<std::pair<double,double>> pintles;
        if(fullAnnulus)
        {
            for(int p=0;p<3;++p)
            {
                const double th=(60.0+120.0*p)*3.14159265358979323846/180.0;
                pintles.emplace_back(boltR*std::cos(th),boltR*std::sin(th));
            }
        }
        else
        {
            const double pintleTheta=0.5*mesh_.sectorAngle();
            pintles.emplace_back(
                boltR*std::cos(pintleTheta),boltR*std::sin(pintleTheta));
        }
        double torchRhoArea=0.0,loxRhoArea=0.0,fuelRhoArea=0.0;
        const PrimitiveVars torch=state(pressure_,2400.0,0.784,0.216);
        const PrimitiveVars lox=state(pressure_,130.0,1.0,0.0);
        const PrimitiveVars fuel=state(pressure_,450.0,0.0,1.0);
        const double rhoTorch=density(torch),rhoLox=density(lox),rhoFuel=density(fuel);
        for(int j=1;j<mesh_.Nr()-1;++j) for(int k=1;k<mesh_.Ntheta()-1;++k)
        {
            const FaceAreaVector c=mesh_.cellCenter(0,j,k);
            const double radius=std::sqrt(c.y*c.y+c.z*c.z);
            double pintleRadius=1e30;
            for(const auto& pc : pintles)
                pintleRadius=std::min(pintleRadius,
                    std::hypot(c.y-pc.first,c.z-pc.second));
            Patch patch=Patch::Background;
            if(radius<0.002) patch=Patch::Torch;
            else if(pintleRadius<0.003) patch=Patch::Lox;
            else if(pintleRadius>=0.006&&pintleRadius<0.0075) patch=Patch::Fuel;
            patches_[index(j,k)]=patch;
            const double area=mesh_.axialArea(1,j,k);
            if(patch==Patch::Torch) { torchRhoArea+=rhoTorch*area; ++report_.torchFaces; }
            if(patch==Patch::Lox) { loxRhoArea+=rhoLox*area; ++report_.loxFaces; }
            if(patch==Patch::Fuel) { fuelRhoArea+=rhoFuel*area; ++report_.fuelFaces; }
        }
        if(torchRhoArea<=0.0||loxRhoArea<=0.0||fuelRhoArea<=0.0)
            throw std::runtime_error("sector mesh does not resolve every injection patch");
        // Design split (vortex.txt): per-pintle LOX 0.9007 / C3H8 0.3378;
        // Core Torch (LOX 0.0980 + C3H8 0.0367 = 0.1348) is INSIDE the
        // 3.85 kg/s total, so the 120-degree sector target is 3.85/3 =
        // 1.2834 kg/s (NOT 2.8/3+1.05/3+0.1348/3 which double-counts).
        if(fullAnnulus)
        {
            report_.torchMassFlow=0.1348;      // central torch (full annulus)
            report_.loxMassFlow=3.0*0.9007;    // three pintles: LOX 2.7021
            report_.fuelMassFlow=3.0*0.3378;   // three pintles: C3H8 1.0134
        }
        else
        {
            report_.torchMassFlow=0.1348/3.0;
            report_.loxMassFlow=0.9007;
            report_.fuelMassFlow=0.3378;
        }
        report_.torchAxialVelocity=report_.torchMassFlow/torchRhoArea;
        report_.loxAxialVelocity=report_.loxMassFlow/loxRhoArea;
        report_.fuelAxialVelocity=report_.fuelMassFlow/fuelRhoArea;
    }

    void applyImpl(
        Domain& d,double massFlowScale,bool inStepTransaction) const
    {
        requireLegacySinglePhaseBoundaryModel(d);
        massFlowScale=std::clamp(massFlowScale,0.0,1.0);
        const FaceAreaVector representative = mesh_.cellCenter(faceI_-1,1,1);
        const thermo::EquilibriumState torchThermo =
            initializeBoundaryThermo(
                d,injectionState(Patch::Torch,representative,0.0));
        const thermo::EquilibriumState loxThermo =
            initializeBoundaryThermo(
                d,injectionState(Patch::Lox,representative,0.0));
        const thermo::EquilibriumState fuelThermo =
            initializeBoundaryThermo(
                d,injectionState(Patch::Fuel,representative,0.0));
        const double torchScale=controlledVelocityScale(
            d,Patch::Torch,massFlowScale,torchThermo);
        const double loxScale=controlledVelocityScale(
            d,Patch::Lox,massFlowScale,loxThermo);
        const double fuelScale=controlledVelocityScale(
            d,Patch::Fuel,massFlowScale,fuelThermo);

        std::vector<Domain::CompleteCellReplacement> replacements;
        replacements.reserve(
            static_cast<std::size_t>(d.Ny())
            *static_cast<std::size_t>(d.Nz()));
        for(int j=0;j<d.Ny();++j) for(int k=0;k<d.Nz();++k)
        {
            const PrimitiveVars& interior=d.W(faceI_,j,k);
            PrimitiveVars w=interior;
            w.T=700.0;
            for(double& velocity:w.U) velocity=-velocity;
            const Patch patch=patches_[index(j,k)];
            const FaceAreaVector c=mesh_.cellCenter(faceI_-1,j,k);
            if(patch==Patch::Torch)
            {
                w=injectionState(patch,c,torchScale);
            }
            else if(patch==Patch::Lox||patch==Patch::Fuel)
            {
                w=injectionState(patch,c,
                    patch==Patch::Lox?loxScale:fuelScale);
            }
            thermo::EquilibriumState q;
            if(patch==Patch::Torch) q=torchThermo;
            else if(patch==Patch::Lox) q=loxThermo;
            else if(patch==Patch::Fuel) q=fuelThermo;
            else if(interior.T==w.T)
                q=d.Q(faceI_,j,k);
            else
                q=initializeBoundaryThermo(d,w);
            replacements.push_back(
                Domain::CompleteCellReplacement{0,j,k,w,std::move(q)});
        }

        // Only boundary-owned halo cells may be replaced. Writing imposed
        // states into interior cells does not create a nozzle wall or inlet.
        const bool published=inStepTransaction
            ?d.replaceCompleteCellsInStepTransaction(replacements)
            :d.replaceCompleteCells(replacements);
        if(!published)
            throw std::runtime_error(
                "sector injection boundary commit failed after preflight");
    }

    void apply(Domain& d,double massFlowScale=1.0) const
    {
        applyImpl(d,massFlowScale,false);
    }

    void applyInStepTransaction(
        Domain& d,double massFlowScale=1.0) const
    {
        applyImpl(d,massFlowScale,true);
    }

    const SectorInjectionReport& report() const { return report_; }

    std::vector<unsigned char> pressureCorrectionDirichletMask() const
    {
        std::vector<unsigned char> mask(patches_.size(),0U);
        for(int j=1;j<mesh_.Nr()-1;++j)
        for(int k=1;k<mesh_.Ntheta()-1;++k)
            if(patches_[index(j,k)]!=Patch::Background)
                mask[static_cast<std::size_t>(index(j,k))]=1U;
        return mask;
    }

    SectorInjectionFluxReport numericalMassFlows(const Domain& d) const
    {
        requireLegacySinglePhaseBoundaryModel(d);
        SectorInjectionFluxReport result;
        AUSMFluxAdapter flux(species_);
        for(int j=1;j<d.Ny()-1;++j) for(int k=1;k<d.Nz()-1;++k)
        {
            ConservedVars face(d.nSpecies());
            flux.computeFluxFromEquilibrium(
                d.W(0,j,k),d.Q(0,j,k),d.W(1,j,k),d.Q(1,j,k),
                {1.0,0.0,0.0},face);
            const double massFlow=face.rho*mesh_.axialArea(1,j,k);
            switch(patches_[index(j,k)])
            {
                case Patch::Background:result.background+=massFlow;break;
                case Patch::Torch:result.torch+=massFlow;break;
                case Patch::Lox:result.lox+=massFlow;break;
                case Patch::Fuel:result.fuel+=massFlow;break;
            }
        }
        return result;
    }
};

} // End namespace Foam
#endif
