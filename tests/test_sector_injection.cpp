#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "src/framework/Domain.H"
#include "src/framework/ModelCapabilities.H"
#include "src/framework/SectorInjectionBoundary.H"
#include "src/mesh/SectorMeshMetrics.H"
#include "src/thermo/MixingRules/BinaryInteractionDB.H"
#include "src/thermo/MixingRules/SpeciesCriticalDB.H"
#include "src/thermo/PhaseEquilibrium/PREquilibriumThermoModel.H"
#include "src/thermo/SpeciesRegistry.H"

using namespace Foam;

namespace
{

bool sameDouble(double lhs,double rhs)
{
    return lhs==rhs||(std::isnan(lhs)&&std::isnan(rhs));
}

bool sameVector(
    const std::vector<double>& lhs,const std::vector<double>& rhs)
{
    if(lhs.size()!=rhs.size()) return false;
    for(std::size_t i=0;i<lhs.size();++i)
        if(!sameDouble(lhs[i],rhs[i])) return false;
    return true;
}

bool samePhase(
    const thermo::EquilibriumPhaseProperties& lhs,
    const thermo::EquilibriumPhaseProperties& rhs)
{
    return sameVector(lhs.massFractions,rhs.massFractions)
        &&sameVector(lhs.moleFractions,rhs.moleFractions)
        &&sameDouble(lhs.molecularWeight,rhs.molecularWeight)
        &&sameDouble(lhs.rho,rhs.rho)
        &&sameDouble(lhs.h,rhs.h)
        &&sameDouble(lhs.e,rhs.e)
        &&sameDouble(lhs.s,rhs.s)
        &&sameDouble(lhs.compressibility,rhs.compressibility)
        &&lhs.entropyAvailable==rhs.entropyAvailable;
}

bool sameOptionalPhase(
    const std::optional<thermo::EquilibriumPhaseProperties>& lhs,
    const std::optional<thermo::EquilibriumPhaseProperties>& rhs)
{
    return lhs.has_value()==rhs.has_value()
        &&(!lhs.has_value()||samePhase(*lhs,*rhs));
}

bool sameThermodynamics(
    const std::optional<thermo::EquilibriumThermodynamicPayload>& lhs,
    const std::optional<thermo::EquilibriumThermodynamicPayload>& rhs)
{
    if(lhs.has_value()!=rhs.has_value()) return false;
    if(!lhs.has_value()) return true;
    return lhs->version==rhs->version
        &&lhs->closure==rhs->closure
        &&sameDouble(lhs->Cp,rhs->Cp)
        &&sameDouble(lhs->Cv,rhs->Cv)
        &&sameDouble(lhs->dpdrho_T,rhs->dpdrho_T)
        &&sameDouble(lhs->dpdT_rho,rhs->dpdT_rho)
        &&sameDouble(lhs->dedrho_T,rhs->dedrho_T)
        &&sameDouble(lhs->dedT_rho,rhs->dedT_rho);
}

bool sameEquilibrium(
    const thermo::EquilibriumState& lhs,
    const thermo::EquilibriumState& rhs)
{
    return lhs.region==rhs.region
        &&sameVector(lhs.overallMassFractions,rhs.overallMassFractions)
        &&sameDouble(lhs.p,rhs.p)
        &&sameDouble(lhs.T,rhs.T)
        &&sameDouble(lhs.rho,rhs.rho)
        &&sameDouble(lhs.h,rhs.h)
        &&sameDouble(lhs.e,rhs.e)
        &&sameDouble(lhs.s,rhs.s)
        &&sameDouble(lhs.vaporMoleFraction,rhs.vaporMoleFraction)
        &&sameDouble(lhs.vaporMassFraction,rhs.vaporMassFraction)
        &&sameDouble(lhs.vaporVolumeFraction,rhs.vaporVolumeFraction)
        &&lhs.entropyAvailable==rhs.entropyAvailable
        &&lhs.acoustic.closure==rhs.acoustic.closure
        &&sameDouble(
            lhs.acoustic.pressureDensityDerivative,
            rhs.acoustic.pressureDensityDerivative)
        &&sameThermodynamics(lhs.thermodynamics,rhs.thermodynamics)
        &&sameOptionalPhase(lhs.liquid,rhs.liquid)
        &&sameOptionalPhase(lhs.vapor,rhs.vapor)
        &&sameOptionalPhase(lhs.singlePhase,rhs.singlePhase);
}

bool sameConserved(const ConservedVars& lhs,const ConservedVars& rhs)
{
    return sameDouble(lhs.rho,rhs.rho)
        &&sameVector(lhs.rhoU,rhs.rhoU)
        &&sameDouble(lhs.rhoE,rhs.rhoE)
        &&sameVector(lhs.rhoY,rhs.rhoY)
        &&sameDouble(lhs.rhoZ,rhs.rhoZ)
        &&sameDouble(lhs.rhoC,rhs.rhoC);
}

bool samePrimitive(const PrimitiveVars& lhs,const PrimitiveVars& rhs)
{
    return sameDouble(lhs.p,rhs.p)
        &&sameDouble(lhs.T,rhs.T)
        &&sameVector(lhs.U,rhs.U)
        &&sameVector(lhs.Y,rhs.Y)
        &&sameDouble(lhs.Z,rhs.Z)
        &&sameDouble(lhs.C,rhs.C);
}

struct DomainSnapshot
{
    std::vector<ConservedVars> current;
    std::vector<ConservedVars> candidate;
    std::vector<PrimitiveVars> primitive;
    std::vector<thermo::EquilibriumState> equilibrium;
};

DomainSnapshot capture(const Domain& domain)
{
    DomainSnapshot result;
    result.current.reserve(domain.N());
    result.candidate.reserve(domain.N());
    result.primitive.reserve(domain.N());
    result.equilibrium.reserve(domain.N());
    for(int cell=0;cell<domain.N();++cell)
    {
        result.current.push_back(domain.U(cell));
        result.candidate.push_back(domain.Unew(cell));
        result.primitive.push_back(domain.W(cell));
        result.equilibrium.push_back(domain.Q(cell));
    }
    return result;
}

bool sameDomain(const Domain& domain,const DomainSnapshot& expected)
{
    if(domain.N()!=static_cast<int>(expected.current.size())) return false;
    for(int cell=0;cell<domain.N();++cell)
    {
        if(!sameConserved(domain.U(cell),expected.current[cell])
            ||!sameConserved(domain.Unew(cell),expected.candidate[cell])
            ||!samePrimitive(domain.W(cell),expected.primitive[cell])
            ||!sameEquilibrium(domain.Q(cell),expected.equilibrium[cell]))
            return false;
    }
    return true;
}

class CommitRejectingLegacyThermoModel final
    : public thermo::EquilibriumThermoModel
{
    thermo::LegacySinglePhaseThermoModel delegate_;
    double rejectedTemperature_;
    std::string modelId_;

public:
    CommitRejectingLegacyThermoModel(
        std::vector<std::string> species,double rejectedTemperature)
        :delegate_(std::move(species)),
         rejectedTemperature_(rejectedTemperature),
         modelId_("legacy-single-phase:commit-validation-rejection")
    {}

    thermo::EquilibriumResult recoverRhoE(
        double density,double internalEnergy,
        const std::vector<double>& overallMassFractions,
        const thermo::EquilibriumSeed* seed=nullptr) const override
    {
        return delegate_.recoverRhoE(
            density,internalEnergy,overallMassFractions,seed);
    }

    thermo::EquilibriumResult initializePTY(
        double pressure,double temperature,
        const std::vector<double>& overallMassFractions) const override
    {
        return delegate_.initializePTY(
            pressure,temperature,overallMassFractions);
    }

    const std::string& modelId() const override { return modelId_; }
    const std::vector<std::string>& speciesNames() const override
    { return delegate_.speciesNames(); }
    EnergyConvention energyConvention() const override
    { return delegate_.energyConvention(); }
    bool stateValidForModel(
        const thermo::EquilibriumState& state) const override
    {
        return std::abs(state.T-rejectedTemperature_)>1.0e-12
            &&delegate_.stateValidForModel(state);
    }
    bool supportsPrimitivePTYInitialization() const override { return true; }
    bool supportsCommittedLegacyGasTransport() const override { return true; }
};

} // namespace

int main()
{
    static_assert(!ModelCapabilities::phaseAwareBoundaryStates,
        "phase-aware injection must stay capability-gated");
    for(int level=0;level<2;++level)
    {
        const int nx=(level?120:60)+2,nr=(level?30:20)+2,nt=(level?20:12)+2;
        SectorMeshMetrics mesh(vortexProfilePath(),nx,nr,nt,120.0);
        bool rejectedOffset=false;
        try { SectorInjectionBoundary unsupported(mesh,4e6,{"O2","C3H8"},{},9.0); }
        catch(const std::invalid_argument&) { rejectedOffset=true; }
        if(!rejectedOffset)
        {
            std::cerr << "interior cell overwrite was accepted as an injector boundary\n";
            return 1;
        }
        Domain d(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),2);
        for(int i=0;i<nx;++i)for(int j=0;j<nr;++j)for(int k=0;k<nt;++k)
        {auto& w=d.stagedW(i,j,k);w.p=4e6;w.T=900;w.U={20,-3,4};w.Y={0.8,0.2};}
        d.primitiveToConserved();
        SectorInjectionBoundary inlet(mesh); inlet.apply(d);
        const SectorInjectionFluxReport fullFlux=inlet.numericalMassFlows(d);
        if(std::abs(fullFlux.torch-0.1348/3.0)>1e-4*(0.1348/3.0)
           ||std::abs(fullFlux.lox-0.9007)>1e-4*0.9007
           ||std::abs(fullFlux.fuel-0.3378)>1e-4*0.3378)
        {
            std::cerr << "full numerical flux=" << fullFlux.torch << "/"
                      << fullFlux.lox << "/" << fullFlux.fuel << "\n";
            return 1;
        }
        const auto& r=inlet.report();
        if(r.torchFaces<=0||r.loxFaces<=0||r.fuelFaces<=0) return 1;
        const std::vector<unsigned char> pressureMask=
            inlet.pressureCorrectionDirichletMask();
        int pressureDirichletFaces=0;
        for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
            pressureDirichletFaces+=pressureMask[
                static_cast<std::size_t>(k*nr+j)]!=0U?1:0;
        const int injectedFaces=r.torchFaces+r.loxFaces+r.fuelFaces;
        const int activeInletFaces=(nr-2)*(nt-2);
        if(pressureMask.size()!=static_cast<std::size_t>(nr*nt)
            ||pressureDirichletFaces!=injectedFaces
            ||pressureDirichletFaces>=activeInletFaces)
        {
            std::cerr << "pressure correction inlet mask="
                      << pressureDirichletFaces << "/"
                      << activeInletFaces << " injected="
                      << injectedFaces << "\n";
            return 1;
        }
        const double backgroundFlux=inlet.numericalMassFlows(d).background;
        if(std::abs(backgroundFlux)>1e-14)
        {
            std::cerr << "background wall mass flux=" << backgroundFlux << "\n";
            return 1;
        }

        inlet.apply(d,0.1);
        const SectorInjectionFluxReport rampFlux=inlet.numericalMassFlows(d);
        if(std::abs(rampFlux.torch-0.1*(0.1348/3.0))>1e-4*(0.1348/3.0)
           ||std::abs(rampFlux.lox-0.1*0.9007)>1e-4*0.9007
           ||std::abs(rampFlux.fuel-0.1*0.3378)>1e-4*0.3378)
        {
            std::cerr << "ramped numerical flux=" << rampFlux.torch << "/"
                      << rampFlux.lox << "/" << rampFlux.fuel << "\n";
            return 1;
        }

        const DomainSnapshot beforeControlledTransaction=capture(d);
        if(!d.beginStepTransaction()) return 1;
        std::vector<Domain::CompleteCellReplacement> perturbedInterior;
        perturbedInterior.reserve(
            static_cast<std::size_t>((nr-2)*(nt-2)));
        for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
        {
            PrimitiveVars perturbed=d.W(1,j,k);
            perturbed.p*=0.92;
            perturbed.T*=1.03;
            perturbed.U[0]+=15.0;
            const thermo::EquilibriumResult initialized=
                d.thermoModel().initializePTY(
                    perturbed.p,perturbed.T,perturbed.Y);
            if(!initialized.success()) return 1;
            perturbedInterior.push_back(
                Domain::CompleteCellReplacement{
                    1,j,k,perturbed,*initialized.state});
        }
        if(!d.replaceCompleteCellsInStepTransaction(perturbedInterior))
            return 1;
        inlet.applyInStepTransaction(d,0.1);
        const SectorInjectionFluxReport controlledFlux=
            inlet.numericalMassFlows(d);
        if(std::abs(controlledFlux.torch-0.1*(0.1348/3.0))
                >1e-4*(0.1348/3.0)
           ||std::abs(controlledFlux.lox-0.1*0.9007)>1e-4*0.9007
           ||std::abs(controlledFlux.fuel-0.1*0.3378)>1e-4*0.3378)
        {
            std::cerr << "transaction-controlled numerical flux="
                      << controlledFlux.torch << "/"
                      << controlledFlux.lox << "/"
                      << controlledFlux.fuel << "\n";
            return 1;
        }
        d.rollbackStepTransaction(false);
        if(d.tripped()||!sameDomain(d,beforeControlledTransaction))
        {
            std::cerr << "controlled injection transaction did not roll back\n";
            return 1;
        }
    }
    {
        SpeciesRegistry registry;
        const int nx=62,nr=22,nt=14;
        SectorMeshMetrics mesh(vortexProfilePath(),nx,nr,nt,120.0);
        Domain d(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),registry.size());
        d.setSpecies(registry.names());
        for(int i=0;i<nx;++i) for(int j=0;j<nr;++j) for(int k=0;k<nt;++k)
        {
            auto& w=d.stagedW(i,j,k);
            w.p=4e6; w.T=700.0; w.U={0.0,0.0,0.0};
            w.Y.assign(registry.size(),0.0);
            w.Y[registry.index("O2")]=0.8;
            w.Y[registry.index("C3H8")]=0.2;
        }
        if(!d.primitiveToConserved()) return 1;
        SectorInjectionBoundary inlet(mesh,4e6,registry.names()); inlet.apply(d);
        const int o2=registry.index("O2"),fuel=registry.index("C3H8");
        bool foundFuel=false;
        for(int j=1;j<nr-1;++j) for(int k=1;k<nt-1;++k)
        {
            const auto& w=d.W(0,j,k);
            if(w.Y[fuel]>0.999) foundFuel=true;
            for(int s=0;s<registry.size();++s)
                if(s!=o2&&s!=fuel&&std::abs(w.Y[s])>1e-15) return 1;
        }
        if(!foundFuel) return 1;
    }
    {
        const std::vector<std::string> species={"O2","C3H8"};
        const auto model=
            std::make_shared<const CommitRejectingLegacyThermoModel>(
                species,450.0);
        const auto background=model->initializePTY(4e6,700.0,{0.8,0.2});
        if(!background.success()) return 1;

        const int nx=4,nr=22,nt=14;
        SectorMeshMetrics mesh(vortexProfilePath(),nx,nr,nt,120.0);
        Domain d(
            nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),2,model);
        PrimitiveVars primitive(2);
        primitive.p=4e6;
        primitive.T=700.0;
        primitive.U={0.0,0.0,0.0};
        primitive.Y={0.8,0.2};
        if(!d.replaceCompleteField(
                std::vector<PrimitiveVars>(d.N(),primitive),
                std::vector<thermo::EquilibriumState>(
                    d.N(),*background.state)))
            return 1;

        SectorInjectionBoundary inlet(mesh,4e6,species,{0.8,0.2});
        const DomainSnapshot beforeCommitRejection=capture(d);
        bool commitStageRejected=false;
        try
        {
            inlet.apply(d,0.1);
        }
        catch(const std::runtime_error&)
        {
            commitStageRejected=true;
        }
        if(!commitStageRejected||!d.tripped()
            ||d.failure().component!="Domain::replaceCompleteCells"
            ||d.failure().cell<=0
            ||!sameDomain(d,beforeCommitRejection))
        {
            std::cerr<<"legacy commit-stage boundary rejection was not atomic\n";
            return 1;
        }
    }
    {
        const SpeciesCriticalDB critical=SpeciesCriticalDB::createDefault();
        const BinaryInteractionDB kij=BinaryInteractionDB::createDefault();
        const std::vector<std::string> species={"O2","C3H8"};
        phaseEquilibrium::RhoEFlashOptions options;
        options.minimumPressure=1e6;
        options.maximumPressure=6e6;
        options.minimumTemperature=100.0;
        options.maximumTemperature=3000.0;
        options.initialPressure=4e6;
        options.initialTemperature=700.0;
        options.pressureScanIntervals=64;
        options.maximumIterations=120;
        options.densityAbsoluteTolerance=1e-7;
        options.densityRelativeTolerance=2e-8;
        options.energyAbsoluteTolerance=1e-4;
        options.energyRelativeTolerance=2e-10;
        const auto model=
            std::make_shared<const thermo::PREquilibriumThermoModel>(
                critical,kij,species,options);
        const auto background=
            model->initializePTY(4e6,700.0,{0.8,0.2});
        if(!background.success())
        {
            std::cerr<<"phase-aware background initialization failed\n";
            return 1;
        }

        const int nx=4,nr=22,nt=14;
        SectorMeshMetrics mesh(vortexProfilePath(),nx,nr,nt,120.0);
        Domain d(
            nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),2,model);
        PrimitiveVars primitive(2);
        primitive.p=4e6;
        primitive.T=700.0;
        primitive.U={0.0,0.0,0.0};
        primitive.Y={0.8,0.2};
        if(!d.replaceCompleteField(
                std::vector<PrimitiveVars>(d.N(),primitive),
                std::vector<thermo::EquilibriumState>(
                    d.N(),*background.state)))
            return 1;

        SectorInjectionBoundary inlet(mesh,4e6,species,{0.8,0.2});
        const DomainSnapshot beforeBoundaryCall=capture(d);
        bool injectionRejected=false;
        std::string injectionError;
        try
        {
            inlet.apply(d,0.1);
        }
        catch(const std::runtime_error& error)
        {
            injectionRejected=true;
            injectionError=error.what();
        }
        if(!injectionRejected
            ||injectionError.find("phase-aware boundary states are not implemented")
                ==std::string::npos
            ||d.tripped()||!sameDomain(d,beforeBoundaryCall))
        {
            std::cerr<<"phase-aware injection was not rejected before mutation\n";
            return 1;
        }
        bool fluxRejected=false;
        std::string fluxError;
        try
        {
            (void)inlet.numericalMassFlows(d);
        }
        catch(const std::runtime_error& error)
        {
            fluxRejected=true;
            fluxError=error.what();
        }
        if(!fluxRejected
            ||fluxError.find("phase-aware boundary states are not implemented")
                ==std::string::npos
            ||d.tripped()||!sameDomain(d,beforeBoundaryCall)) return 1;
    }
    std::cout<<"SECTOR_INJECTION_PASS\n";
    return 0;
}
