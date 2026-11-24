#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "src/combustion/FiniteRate/TwoStepCO.H"
#include "src/combustion/SectorChemistryOperator.H"
#include "src/framework/Domain.H"
#include "src/framework/TimeIntegrator.H"
#include "src/thermo/MixingRules/BinaryInteractionDB.H"
#include "src/thermo/MixingRules/SpeciesCriticalDB.H"
#include "src/thermo/PhaseEquilibrium/PREquilibriumThermoModel.H"
#include "src/thermo/PhaseEquilibrium/PRPhaseDisplacementAuthority.H"
#include "src/thermo/PhaseEquilibrium/PRPhaseRhoERecovery.H"
#include "src/thermo/PhaseEquilibrium/PRTPFlash.H"
#include "src/thermo/RealFluidState.H"
#include "src/thermo/SpeciesRegistry.H"

using namespace Foam;

namespace
{

void elementInventory(
    const std::vector<double>& Y,const SpeciesRegistry& registry,
    double result[4])
{
    for(int e=0;e<4;++e) result[e]=0.0;
    for(int g=0;g<registry.size();++g)
    {
        const double mol=Y[g]/registry.molecularWeight(g);
        result[0]+=mol*registry.nC(g);
        result[1]+=mol*registry.nH(g);
        result[2]+=mol*registry.nO(g);
        result[3]+=mol*registry.nN(g);
    }
}

bool nonnegativeNormalized(const std::vector<double>& Y)
{
    double sum=0.0;
    for(double value:Y)
    {
        if(!std::isfinite(value)||value<-1.0e-12) return false;
        sum+=value;
    }
    return std::abs(sum-1.0)<=1.0e-6;
}

bool close(
    double lhs,double rhs,double relativeTolerance=1.0e-7,
    double absoluteTolerance=1.0e-10)
{
    return std::isfinite(lhs)&&std::isfinite(rhs)
        &&std::abs(lhs-rhs)<=absoluteTolerance+relativeTolerance
            *std::max(std::abs(lhs),std::abs(rhs));
}

struct ExactDomainState
{
    std::vector<ConservedVars> current;
    std::vector<ConservedVars> candidate;
    std::vector<PrimitiveVars> primitive;
    std::vector<thermo::EquilibriumState> equilibrium;
};

bool exactScalar(double lhs,double rhs)
{
    return lhs==rhs||(std::isnan(lhs)&&std::isnan(rhs));
}

bool exactVector(
    const std::vector<double>& lhs,const std::vector<double>& rhs)
{
    if(lhs.size()!=rhs.size()) return false;
    for(std::size_t i=0;i<lhs.size();++i)
        if(!exactScalar(lhs[i],rhs[i])) return false;
    return true;
}

bool exactSame(
    const thermo::EquilibriumPhaseProperties& lhs,
    const thermo::EquilibriumPhaseProperties& rhs)
{
    return exactVector(lhs.massFractions,rhs.massFractions)
        &&exactVector(lhs.moleFractions,rhs.moleFractions)
        &&exactScalar(lhs.molecularWeight,rhs.molecularWeight)
        &&exactScalar(lhs.rho,rhs.rho)
        &&exactScalar(lhs.h,rhs.h)
        &&exactScalar(lhs.e,rhs.e)
        &&exactScalar(lhs.s,rhs.s)
        &&exactScalar(lhs.compressibility,rhs.compressibility)
        &&lhs.entropyAvailable==rhs.entropyAvailable;
}

bool exactSame(
    const thermo::EquilibriumThermodynamicPayload& lhs,
    const thermo::EquilibriumThermodynamicPayload& rhs)
{
    return lhs.version==rhs.version
        &&lhs.closure==rhs.closure
        &&exactScalar(lhs.Cp,rhs.Cp)
        &&exactScalar(lhs.Cv,rhs.Cv)
        &&exactScalar(lhs.dpdrho_T,rhs.dpdrho_T)
        &&exactScalar(lhs.dpdT_rho,rhs.dpdT_rho)
        &&exactScalar(lhs.dedrho_T,rhs.dedrho_T)
        &&exactScalar(lhs.dedT_rho,rhs.dedT_rho);
}

bool exactSame(
    const thermo::EquilibriumState& lhs,
    const thermo::EquilibriumState& rhs)
{
    if(lhs.liquid.has_value()!=rhs.liquid.has_value()
        ||lhs.vapor.has_value()!=rhs.vapor.has_value()
        ||lhs.singlePhase.has_value()!=rhs.singlePhase.has_value()
        ||lhs.thermodynamics.has_value()!=rhs.thermodynamics.has_value())
        return false;
    return lhs.region==rhs.region
        &&exactVector(
            lhs.overallMassFractions,rhs.overallMassFractions)
        &&exactScalar(lhs.p,rhs.p)
        &&exactScalar(lhs.T,rhs.T)
        &&exactScalar(lhs.rho,rhs.rho)
        &&exactScalar(lhs.h,rhs.h)
        &&exactScalar(lhs.e,rhs.e)
        &&exactScalar(lhs.s,rhs.s)
        &&exactScalar(lhs.vaporMoleFraction,rhs.vaporMoleFraction)
        &&exactScalar(lhs.vaporMassFraction,rhs.vaporMassFraction)
        &&exactScalar(lhs.vaporVolumeFraction,rhs.vaporVolumeFraction)
        &&lhs.acoustic.closure==rhs.acoustic.closure
        &&exactScalar(
            lhs.acoustic.pressureDensityDerivative,
            rhs.acoustic.pressureDensityDerivative)
        &&lhs.entropyAvailable==rhs.entropyAvailable
        &&(!lhs.thermodynamics
            ||exactSame(*lhs.thermodynamics,*rhs.thermodynamics))
        &&(!lhs.liquid||exactSame(*lhs.liquid,*rhs.liquid))
        &&(!lhs.vapor||exactSame(*lhs.vapor,*rhs.vapor))
        &&(!lhs.singlePhase
            ||exactSame(*lhs.singlePhase,*rhs.singlePhase));
}

bool exactSame(const ConservedVars& lhs,const ConservedVars& rhs)
{
    return lhs.rho==rhs.rho&&lhs.rhoU==rhs.rhoU
        &&lhs.rhoE==rhs.rhoE&&lhs.rhoY==rhs.rhoY
        &&lhs.rhoZ==rhs.rhoZ&&lhs.rhoC==rhs.rhoC;
}

bool exactSame(const PrimitiveVars& lhs,const PrimitiveVars& rhs)
{
    return lhs.p==rhs.p&&lhs.T==rhs.T&&lhs.U==rhs.U
        &&lhs.Y==rhs.Y&&lhs.Z==rhs.Z&&lhs.C==rhs.C;
}

ExactDomainState captureExactState(const Domain& domain)
{
    ExactDomainState result;
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

bool exactStateMatches(
    const Domain& domain,const ExactDomainState& expected)
{
    if(domain.N()!=static_cast<int>(expected.current.size())) return false;
    for(int cell=0;cell<domain.N();++cell)
        if(!exactSame(domain.U(cell),expected.current[cell])
            ||!exactSame(domain.Unew(cell),expected.candidate[cell])
            ||!exactSame(domain.W(cell),expected.primitive[cell])
            ||!exactSame(domain.Q(cell),expected.equilibrium[cell]))
            return false;
    return true;
}

class VaporProbeChemistry final : public ChemistryModel
{
    static constexpr double fuelMolecularWeight_=0.04409562;
    static constexpr double oxygenMolecularWeight_=0.0319988;
    static constexpr double coMolecularWeight_=0.0280101;
    static constexpr double waterMolecularWeight_=0.01801528;
    static constexpr double reactionHeat_=1.19502e6;
    double rateConstant_;
    mutable bool observed_=false;
    mutable double observedTemperature_=0.0;
    mutable double observedPressure_=0.0;
    mutable double observedDensity_=0.0;
    mutable std::vector<double> observedY_;

public:
    explicit VaporProbeChemistry(
        const SpeciesRegistry& registry,double rateConstant=50.0)
        :rateConstant_(rateConstant)
    {
        speciesNames_={"C3H8","O2","CO","CO2","H2O","N2"};
        nSpecies_=static_cast<int>(speciesNames_.size());
        nReactions_=1;
        molecularWeights_.reserve(speciesNames_.size());
        for(const std::string& name:speciesNames_)
            molecularWeights_.push_back(
                registry.molecularWeight(registry.index(name)));
    }

    void computeSourceTerms(
        double temperature,double pressure,double density,
        const std::vector<double>& Y,
        std::vector<double>& omega,double& heatRelease) const override
    {
        if(!observed_)
        {
            observed_=true;
            observedTemperature_=temperature;
            observedPressure_=pressure;
            observedDensity_=density;
            observedY_=Y;
        }
        omega.assign(static_cast<std::size_t>(nSpecies_),0.0);
        const int fuel=speciesIndex("C3H8");
        const int oxygen=speciesIndex("O2");
        const int co=speciesIndex("CO");
        const int water=speciesIndex("H2O");
        const double fuelMoles=density*std::max(Y[fuel],0.0)
            /fuelMolecularWeight_;
        const double oxygenMoles=density*std::max(Y[oxygen],0.0)
            /oxygenMolecularWeight_;
        const double extentRate=rateConstant_
            *std::min(fuelMoles,oxygenMoles/3.5);
        omega[fuel]=-extentRate*fuelMolecularWeight_;
        omega[oxygen]=-3.5*extentRate*oxygenMolecularWeight_;
        omega[co]=3.0*extentRate*coMolecularWeight_;
        omega[water]=4.0*extentRate*waterMolecularWeight_;
        heatRelease=extentRate*reactionHeat_;
    }

    void resetObservation() const
    {
        observed_=false;
        observedTemperature_=observedPressure_=observedDensity_=0.0;
        observedY_.clear();
    }
    bool observed() const { return observed_; }
    double observedTemperature() const { return observedTemperature_; }
    double observedPressure() const { return observedPressure_; }
    double observedDensity() const { return observedDensity_; }
    const std::vector<double>& observedY() const { return observedY_; }
};

class SelectiveFailChemistry final : public ChemistryModel
{
    VaporProbeChemistry successfulModel_;
    double failureDensityThreshold_;
    mutable std::atomic<int> successfulCalls_{0};
    mutable std::atomic<int> failingCalls_{0};

public:
    SelectiveFailChemistry(
        const SpeciesRegistry& registry,double failureDensityThreshold)
        :successfulModel_(registry),
         failureDensityThreshold_(failureDensityThreshold)
    {
        if(!std::isfinite(failureDensityThreshold_)
            ||failureDensityThreshold_<=0.0)
            throw std::invalid_argument(
                "SelectiveFailChemistry: invalid density threshold");
        speciesNames_={"C3H8","O2","CO","CO2","H2O","N2"};
        nSpecies_=static_cast<int>(speciesNames_.size());
        nReactions_=1;
        molecularWeights_.reserve(speciesNames_.size());
        for(const std::string& name:speciesNames_)
            molecularWeights_.push_back(
                registry.molecularWeight(registry.index(name)));
    }

    void computeSourceTerms(
        double temperature,double pressure,double density,
        const std::vector<double>& Y,
        std::vector<double>& omega,double& heatRelease) const override
    {
        if(density>failureDensityThreshold_)
        {
            failingCalls_.fetch_add(1,std::memory_order_relaxed);
            omega.assign(
                static_cast<std::size_t>(nSpecies_),
                std::numeric_limits<double>::quiet_NaN());
            heatRelease=std::numeric_limits<double>::quiet_NaN();
            return;
        }
        successfulCalls_.fetch_add(1,std::memory_order_relaxed);
        successfulModel_.computeSourceTerms(
            temperature,pressure,density,Y,omega,heatRelease);
    }

    void resetCounts() const
    {
        successfulCalls_.store(0,std::memory_order_relaxed);
        failingCalls_.store(0,std::memory_order_relaxed);
    }

    int successfulCalls() const
    { return successfulCalls_.load(std::memory_order_relaxed); }
    int failingCalls() const
    { return failingCalls_.load(std::memory_order_relaxed); }
};

BinaryInteractionDB completeTestInteractions(
    const SpeciesRegistry& registry)
{
    BinaryInteractionDB result=BinaryInteractionDB::createDefault();
    const std::vector<std::string> names=registry.names();
    for(std::size_t i=0;i<names.size();++i)
        for(std::size_t j=i+1;j<names.size();++j)
            if(!result.has(names[i],names[j]))
                result.add(names[i],names[j],0.0);
    return result;
}

phaseEquilibrium::RhoEFlashOptions twoPhaseOptions()
{
    phaseEquilibrium::RhoEFlashOptions options;
    options.minimumPressure=1.0e6;
    options.maximumPressure=4.0e6;
    options.minimumTemperature=250.0;
    options.maximumTemperature=1500.0;
    options.initialPressure=2.0e6;
    options.initialTemperature=300.0;
    options.pressureScanIntervals=64;
    options.maximumIterations=120;
    options.finiteDifferenceStep=1.0e-4;
    options.logStepLimit=0.75;
    options.densityAbsoluteTolerance=1.0e-7;
    options.densityRelativeTolerance=2.0e-8;
    options.energyAbsoluteTolerance=1.0e-4;
    options.energyRelativeTolerance=2.0e-10;
    return options;
}

phaseEquilibrium::PRPhaseRhoEOptions phaseChemistryRecoveryOptions()
{
    phaseEquilibrium::PRPhaseRhoEOptions options;
    options.minimumTemperature=250.0;
    options.maximumTemperature=1500.0;
    options.temperatureScanIntervals=128;
    options.maximumIterations=80;
    options.derivativeRelativeStep=3.0e-8;
    options.energyAbsoluteTolerance=1.0e-5;
    options.energyRelativeTolerance=2.0e-10;
    options.densityRelativeTolerance=2.0e-9;
    options.rootRelativeTolerance=2.0e-12;
    options.derivativeIdentityRelativeTolerance=2.0e-6;
    options.derivativeIdentityAbsoluteTolerance=1.0e-4;
    return options;
}

phaseEquilibrium::PRPhaseDisplacementBounds
phaseChemistryDisplacementBounds()
{
    return phaseEquilibrium::PRPhaseDisplacementBounds(
        0.02,
        0.05,
        1.0e-3,
        3,
        1.0e-5,
        1.0e-4,
        0.05,
        0.05,
        0.05,
        2.0e-9,
        1.0e-5,
        2.0e-10,
        2.0e-12);
}

bool sameChemistryState(
    const Domain& lhs,const Domain& rhs,int i,int j,int k)
{
    const ConservedVars& lhsU=lhs.U(i,j,k);
    const ConservedVars& rhsU=rhs.U(i,j,k);
    if(!close(lhsU.rho,rhsU.rho,1.0e-8,1.0e-10)
        ||!close(lhsU.rhoE,rhsU.rhoE,5.0e-5,1.0e-4)
        ||lhsU.rhoY.size()!=rhsU.rhoY.size())
        return false;
    for(std::size_t species=0;species<lhsU.rhoY.size();++species)
        if(!close(
                lhsU.rhoY[species],rhsU.rhoY[species],
                5.0e-5,1.0e-9))
            return false;
    const thermo::EquilibriumState& lhsQ=lhs.Q(i,j,k);
    const thermo::EquilibriumState& rhsQ=rhs.Q(i,j,k);
    return lhsQ.region==thermo::EquilibriumRegion::TWO_PHASE
        &&rhsQ.region==thermo::EquilibriumRegion::TWO_PHASE
        &&close(
            lhsQ.vaporMassFraction,rhsQ.vaporMassFraction,
            5.0e-5,1.0e-8)
        &&close(
            lhsQ.vaporVolumeFraction,rhsQ.vaporVolumeFraction,
            5.0e-5,1.0e-8);
}

bool exerciseReactiveSinglePhaseRoute(
    const std::shared_ptr<const thermo::EquilibriumThermoModel>& model,
    const thermo::EquilibriumState& equilibrium,
    const SpeciesRegistry& registry,
    thermo::EquilibriumRegion expectedRegion,
    const char* label)
{
    PrimitiveVars primitive(registry.size());
    primitive.p=equilibrium.p;
    primitive.T=equilibrium.T;
    primitive.U={0.0,0.0,0.0};
    primitive.Y=equilibrium.overallMassFractions;
    Domain domain(3,3,3,0.03,0.03,0.03,registry.size(),model);
    if(!domain.replaceCompleteField(
            std::vector<PrimitiveVars>(domain.N(),primitive),
            std::vector<thermo::EquilibriumState>(domain.N(),equilibrium)))
    {
        std::cerr<<label<<" Domain initialization failed\n";
        return false;
    }

    VaporProbeChemistry chemistry(registry);
    SectorChemistryOperator chemistryOperator(&chemistry,registry,0.0);
    const double stableDt=chemistryOperator.maxStableDt(domain,0.1);
    if(!std::isfinite(stableDt)||stableDt<=0.0||stableDt>=1.0e29)
    {
        std::cerr<<label<<" stable dt failed: "<<stableDt<<"\n";
        return false;
    }
    const double initialEnergy=domain.U(1,1,1).rhoE;
    chemistry.resetObservation();
    domain.beginStep();
    try
    {
        chemistryOperator.apply(domain,1.0e-3*stableDt);
    }
    catch(const std::runtime_error& error)
    {
        std::cerr<<label<<" apply failed: "<<error.what()
                 <<" failures="<<chemistryOperator.lastFailures()
                 <<" cvode="<<chemistryOperator.lastCvodeFailures()
                 <<" state="<<chemistryOperator.lastStateFailures()
                 <<" T="<<equilibrium.T
                 <<" p="<<equilibrium.p
                 <<" rho="<<equilibrium.rho
                 <<" e="<<equilibrium.e<<"\n";
        return false;
    }
    if(chemistryOperator.lastActiveCells()!=1
        ||chemistryOperator.lastFailures()!=0
        ||chemistryOperator.lastTwoPhaseActiveCells()!=0
        ||chemistryOperator.lastLiquidSkippedCells()!=0
        ||!chemistry.observed()
        ||!close(chemistry.observedTemperature(),equilibrium.T)
        ||!close(chemistry.observedPressure(),equilibrium.p)
        ||!close(chemistry.observedDensity(),equilibrium.rho)
        ||domain.Unew(1,1,1).rhoE<=initialEnergy)
    {
        std::cerr<<label<<" route counters/state failed active="
                 <<chemistryOperator.lastActiveCells()
                 <<" failures="<<chemistryOperator.lastFailures()
                 <<" twoPhase="
                 <<chemistryOperator.lastTwoPhaseActiveCells()
                 <<" liquidSkipped="
                 <<chemistryOperator.lastLiquidSkippedCells()
                 <<" observed="<<chemistry.observed()
                 <<" T="<<chemistry.observedTemperature()
                 <<"/"<<equilibrium.T
                 <<" p="<<chemistry.observedPressure()
                 <<"/"<<equilibrium.p
                 <<" rho="<<chemistry.observedDensity()
                 <<"/"<<equilibrium.rho
                 <<" rhoE="<<domain.Unew(1,1,1).rhoE
                 <<"/"<<initialEnergy<<"\n";
        return false;
    }
    for(int chemistrySpecies=0;
        chemistrySpecies<chemistry.nSpecies();++chemistrySpecies)
    {
        const int globalSpecies=registry.index(
            chemistry.speciesNames()[chemistrySpecies]);
        if(!close(
                chemistry.observedY()[chemistrySpecies],
                equilibrium.overallMassFractions[globalSpecies]))
        {
            std::cerr<<label<<" received wrong overall composition\n";
            return false;
        }
    }
    if(!domain.commit()||domain.Q(1,1,1).region!=expectedRegion)
    {
        std::cerr<<label<<" post-source recovery changed phase route\n";
        return false;
    }
    return true;
}

class UnclassifiedPhaseThermoModel final
    : public thermo::EquilibriumThermoModel
{
    std::shared_ptr<const thermo::EquilibriumThermoModel> delegate_;
    std::string modelId_;

    thermo::EquilibriumResult unclassify(
        thermo::EquilibriumResult result) const
    {
        if(result.success())
            result.state->region=
                thermo::EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE;
        return result;
    }

public:
    explicit UnclassifiedPhaseThermoModel(
        std::shared_ptr<const thermo::EquilibriumThermoModel> delegate)
        :delegate_(std::move(delegate)),
         modelId_(delegate_->modelId()+":forced-unclassified-test")
    {}

    thermo::EquilibriumResult recoverRhoE(
        double density,double internalEnergy,
        const std::vector<double>& overallMassFractions,
        const thermo::EquilibriumSeed* seed=nullptr) const override
    {
        return unclassify(delegate_->recoverRhoE(
            density,internalEnergy,overallMassFractions,seed));
    }

    thermo::EquilibriumResult initializePTY(
        double pressure,double temperature,
        const std::vector<double>& overallMassFractions) const override
    {
        return unclassify(delegate_->initializePTY(
            pressure,temperature,overallMassFractions));
    }

    const std::string& modelId() const override { return modelId_; }
    const std::vector<std::string>& speciesNames() const override
    { return delegate_->speciesNames(); }
    EnergyConvention energyConvention() const override
    { return delegate_->energyConvention(); }
    bool stateValidForModel(
        const thermo::EquilibriumState& state) const override
    { return delegate_->stateValidForModel(state); }
    bool supportsPrimitivePTYInitialization() const override
    { return delegate_->supportsPrimitivePTYInitialization(); }
};

class FailingChemistry final : public ChemistryModel
{
public:
    explicit FailingChemistry(const SpeciesRegistry& registry)
    {
        speciesNames_={"C3H8","O2","CO","CO2","H2O","N2"};
        nSpecies_=static_cast<int>(speciesNames_.size());
        nReactions_=1;
        molecularWeights_.reserve(speciesNames_.size());
        for(const std::string& name:speciesNames_)
            molecularWeights_.push_back(
                registry.molecularWeight(registry.index(name)));
    }

    void computeSourceTerms(
        double,double,double,const std::vector<double>&,
        std::vector<double>& omega,double& heatRelease) const override
    {
        omega.assign(
            static_cast<std::size_t>(nSpecies_),
            std::numeric_limits<double>::quiet_NaN());
        heatRelease=std::numeric_limits<double>::quiet_NaN();
    }
};

} // namespace

int main()
{
    SpeciesRegistry registry;
    TwoStepCO chemistry;
    std::vector<double> mixture(registry.size(),0.0);
    mixture[registry.index("O2")]=0.784;
    mixture[registry.index("C3H8")]=0.216;

    RealFluidState initialThermo;
    initialThermo.setSpecies(registry.names());
    if(initialThermo.updatePTY(4.0e6,1200.0,mixture)
       !=RealFluidState::Status::SUCCESS)
    { std::cerr<<"initial chemistry state failed\n"; return 1; }

    // A bounded single-cell depletion regression exercises the same [Y,e]
    // path used by CFD without reviving any A408 legacy checkpoint state.
    {
        FiniteRateReactor reactor(&chemistry,registry);
        const double sourceDt=reactor.estimateConservativeMaxDt(
            initialThermo.rho(),initialThermo.p(),initialThermo.T(),mixture,0.1);
        std::vector<double> products=mixture;
        double temperature=initialThermo.T(),pressure=initialThermo.p();
        double deltaRhoE=0.0;
        if(!reactor.reactConservative(
            initialThermo.rho(),initialThermo.e(),products,4.0*sourceDt,
            deltaRhoE,temperature,pressure))
        { std::cerr<<"bounded depletion state failed\n"; return 1; }
        if(!nonnegativeNormalized(products) ||
           products[registry.index("C3H8")]>=mixture[registry.index("C3H8")] ||
           !std::isfinite(deltaRhoE)||deltaRhoE<=0.0)
        { std::cerr<<"bounded depletion state is non-conservative\n"; return 1; }
    }

    // The reversible CO oxidation source must remain active at high T.
    {
        std::vector<double> products(chemistry.nSpecies(),0.0),omega;
        products[chemistry.speciesIndex("CO2")]=1.0;
        double qdot=0.0;
        chemistry.computeSourceTerms(2500.0,4.0e6,5.0,products,omega,qdot);
        if(omega[chemistry.speciesIndex("CO2")]>=0.0 ||
           omega[chemistry.speciesIndex("CO")]<=0.0 ||
           omega[chemistry.speciesIndex("O2")]<=0.0 || qdot>=0.0)
        { std::cerr<<"high-temperature CO2 reverse reaction is inactive\n"; return 1; }
    }

    const int nx=3,ny=3,nz=3;
    Domain domain(nx,ny,nz,0.01,0.01,0.01,registry.size());
    domain.setSpecies(registry.names());
    for(int i=0;i<nx;++i) for(int j=0;j<ny;++j) for(int k=0;k<nz;++k)
    {
        PrimitiveVars& state=domain.stagedW(i,j,k);
        state.p=4.0e6;
        state.T=1200.0;
        state.U={0.0,0.0,0.0};
        state.Y=mixture;
    }
    domain.primitiveToConserved();
    if(domain.tripped()) return 1;
    domain.beginStep();

    SectorChemistryOperator chemistryOperator(&chemistry,registry);
    const double stableDt=chemistryOperator.maxStableDt(domain,0.1);
    if(!std::isfinite(stableDt)||stableDt<=0.0)
    { std::cerr<<"invalid chemistry dt\n"; return 1; }
    const double dt=0.5*stableDt;

    const int fuel=registry.index("C3H8");
    const double haloFuel=domain.U(0,0,0).rhoY[fuel];
    const double oldFuel=domain.U(1,1,1).rhoY[fuel];
    const double oldEnergy=domain.U(1,1,1).rhoE;
    double oldElements[4];
    elementInventory(mixture,registry,oldElements);

    chemistryOperator.apply(domain,dt);
    if(chemistryOperator.lastActiveCells()!=1 ||
       chemistryOperator.lastFailures()!=0)
    { std::cerr<<"unexpected sector chemistry counts\n"; return 1; }
    if(domain.Unew(1,1,1).rhoY[fuel]>=oldFuel ||
       domain.Unew(1,1,1).rhoE<=oldEnergy)
    { std::cerr<<"sector chemistry source was not coupled\n"; return 1; }
    if(domain.Unew(0,0,0).rhoY[fuel]!=haloFuel)
    { std::cerr<<"sector chemistry modified a halo\n"; return 1; }

    std::vector<double> newY(registry.size(),0.0);
    for(int g=0;g<registry.size();++g)
        newY[g]=domain.Unew(1,1,1).rhoY[g]/domain.Unew(1,1,1).rho;
    if(!nonnegativeNormalized(newY))
    { std::cerr<<"sector chemistry returned invalid mass fractions\n"; return 1; }
    double newElements[4];
    elementInventory(newY,registry,newElements);
    for(int e=0;e<4;++e)
        if(std::abs(newElements[e]-oldElements[e])>
           2.0e-7*std::max(1.0,std::abs(oldElements[e])))
        { std::cerr<<"sector chemistry element residual too large\n"; return 1; }

    // Operator splitting may leave less transported fuel in Unew than the old
    // primitive state seen by chemistry. The coupled increment must scale
    // species and energy together instead of driving Unew negative.
    domain.beginStep();
    ConservedVars& split=domain.Unew(1,1,1);
    const int product=registry.index("CO2");
    const double transferred=0.999999*split.rhoY[fuel];
    split.rhoY[fuel]-=transferred;
    split.rhoY[product]+=transferred;
    chemistryOperator.apply(domain,dt);
    if(chemistryOperator.lastFailures()!=0 || split.rhoY[fuel]<-1.0e-12)
    { std::cerr<<"split chemistry drove transported fuel negative\n"; return 1; }

    // A failed conservative chemistry solve must not couple a partial
    // workspace into Unew. Inject the failure through the chemistry backend;
    // committed U/W/Q must remain a complete admissible source state.
    const std::vector<double> sentinelSpecies = split.rhoY;
    const double sentinelEnergy = split.rhoE;
    FailingChemistry failureModel(registry);
    SectorChemistryOperator failingOperator(&failureModel,registry);
    bool rejected = false;
    try
    {
        failingOperator.apply(domain,dt);
    }
    catch (const std::runtime_error&)
    {
        rejected = true;
    }
    if(!rejected || failingOperator.lastFailures() != 1 ||
       split.rhoY != sentinelSpecies || split.rhoE != sentinelEnergy)
    { std::cerr<<"failed chemistry state was not rejected atomically\n"; return 1; }

    // Phase-aware operator contract: chemistry sees only committed vapor
    // properties, source increments scale by vapor volume, and the resulting
    // overall rho-e-Y state is reflashed before it can be committed. Missing
    // production binary-interaction data are filled explicitly only in this
    // controlled test fixture; no capability flag is implied by this test.
    const SpeciesCriticalDB critical=SpeciesCriticalDB::createDefault();
    const BinaryInteractionDB testKij=completeTestInteractions(registry);
    phaseEquilibrium::PRTPFlash feedFlash(
        critical,testKij,{"O2","C3H8"});
    const phaseEquilibrium::TPFlashResult feedSplit=
        feedFlash.flashMoleFractions(2.0e6,300.0,{0.2,0.8});
    if(!feedSplit.valid()
        ||feedSplit.region!=phaseEquilibrium::EquilibriumRegion::TWO_PHASE)
    { std::cerr<<"two-phase chemistry feed flash failed\n"; return 1; }
    std::vector<double> twoPhaseY(registry.size(),0.0);
    constexpr double productTrace=1.0e-4;
    const double reactantScale=1.0-2.0*productTrace;
    twoPhaseY[registry.index("O2")]=
        reactantScale*feedSplit.feedMassFractions[0];
    twoPhaseY[registry.index("C3H8")]=
        reactantScale*feedSplit.feedMassFractions[1];
    twoPhaseY[registry.index("CO")]=productTrace;
    twoPhaseY[registry.index("H2O")]=productTrace;
    const auto productionInteractionModel=
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical,BinaryInteractionDB::createDefault(),
            registry.names(),twoPhaseOptions());
    const thermo::EquilibriumResult productionInteractionRejection=
        productionInteractionModel->initializePTY(
            2.0e6,300.0,twoPhaseY);
    if(productionInteractionRejection.status
            !=thermo::EquilibriumStatus::INVALID_CONFIGURATION
        ||productionInteractionRejection.state.has_value()
        ||productionInteractionRejection.diagnostics.backendPhaseStatus
            !="MISSING_BINARY_INTERACTION")
    {
        std::cerr<<"production missing-kij state did not fail closed\n";
        return 1;
    }
    const auto equilibriumModel=
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical,testKij,registry.names(),twoPhaseOptions());

    std::vector<double> routeY(registry.size(),0.0);
    routeY[registry.index("O2")]=0.8;
    routeY[registry.index("C3H8")]=0.2;

    phaseEquilibrium::RhoEFlashOptions liquidOptions=twoPhaseOptions();
    liquidOptions.minimumPressure=1.0e7;
    liquidOptions.maximumPressure=3.0e7;
    liquidOptions.minimumTemperature=140.0;
    liquidOptions.maximumTemperature=250.0;
    liquidOptions.initialPressure=2.0e7;
    liquidOptions.initialTemperature=155.0;
    const auto liquidModel=
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical,BinaryInteractionDB::createDefault(),
            registry.names(),liquidOptions);
    const thermo::EquilibriumResult liquidState=
        liquidModel->initializePTY(2.0e7,155.0,routeY);
    if(!liquidState.success()
        ||liquidState.state->region
            !=thermo::EquilibriumRegion::SINGLE_LIQUID)
    { std::cerr<<"single-liquid chemistry route state failed\n"; return 1; }
    PrimitiveVars liquidPrimitive(registry.size());
    liquidPrimitive.p=liquidState.state->p;
    liquidPrimitive.T=liquidState.state->T;
    liquidPrimitive.U={0.0,0.0,0.0};
    liquidPrimitive.Y=routeY;
    Domain liquidDomain(
        3,3,3,0.03,0.03,0.03,registry.size(),liquidModel);
    if(!liquidDomain.replaceCompleteField(
            std::vector<PrimitiveVars>(liquidDomain.N(),liquidPrimitive),
            std::vector<thermo::EquilibriumState>(
                liquidDomain.N(),*liquidState.state)))
    { std::cerr<<"single-liquid Domain initialization failed\n"; return 1; }
    VaporProbeChemistry liquidChemistry(registry);
    SectorChemistryOperator liquidOperator(
        &liquidChemistry,registry,0.0);
    liquidChemistry.resetObservation();
    const double liquidStableDt=liquidOperator.maxStableDt(liquidDomain,0.1);
    const ConservedVars liquidBefore=liquidDomain.U(1,1,1);
    liquidDomain.beginStep();
    liquidOperator.apply(liquidDomain,1.0e-6);
    const ConservedVars& liquidAfter=liquidDomain.Unew(1,1,1);
    if(liquidStableDt<1.0e29||liquidChemistry.observed()
        ||liquidOperator.lastActiveCells()!=0
        ||liquidOperator.lastLiquidSkippedCells()!=1
        ||liquidOperator.lastFailures()!=0
        ||liquidAfter.rho!=liquidBefore.rho
        ||liquidAfter.rhoU!=liquidBefore.rhoU
        ||liquidAfter.rhoE!=liquidBefore.rhoE
        ||liquidAfter.rhoY!=liquidBefore.rhoY
        ||liquidAfter.rhoZ!=liquidBefore.rhoZ
        ||liquidAfter.rhoC!=liquidBefore.rhoC)
    { std::cerr<<"single-liquid cell entered chemistry\n"; return 1; }

    phaseEquilibrium::RhoEFlashOptions vaporOptions=twoPhaseOptions();
    vaporOptions.minimumPressure=5.0e4;
    vaporOptions.maximumPressure=5.0e5;
    vaporOptions.minimumTemperature=280.0;
    vaporOptions.maximumTemperature=500.0;
    vaporOptions.initialPressure=1.0e5;
    vaporOptions.initialTemperature=300.0;
    const auto vaporModel=
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical,testKij,registry.names(),vaporOptions);
    const thermo::EquilibriumResult vaporState=
        vaporModel->initializePTY(1.0e5,300.0,routeY);
    if(!vaporState.success()
        ||vaporState.state->region
            !=thermo::EquilibriumRegion::SINGLE_VAPOR
        ||!exerciseReactiveSinglePhaseRoute(
            vaporModel,*vaporState.state,registry,
            thermo::EquilibriumRegion::SINGLE_VAPOR,"single-vapor"))
    { std::cerr<<"single-vapor chemistry route failed\n"; return 1; }

    // SectorChemistryOperator publishes successful-cell candidates before it
    // reports a different cell's failure. The integrator transaction, not a
    // direct apply(), owns whole-step atomicity across those cells.
    {
        const thermo::EquilibriumResult failingVaporState=
            vaporModel->initializePTY(1.1e5,300.0,routeY);
        if(!failingVaporState.success()
            ||failingVaporState.state->region
                !=thermo::EquilibriumRegion::SINGLE_VAPOR
            ||failingVaporState.state->rho<=vaporState.state->rho)
        {
            std::cerr<<"transaction failure vapor state failed\n";
            return 1;
        }

        constexpr int transactionNx=4;
        constexpr int transactionNy=3;
        constexpr int transactionNz=3;
        const int successfulCell=
            1*transactionNx*transactionNy+1*transactionNx+1;
        const int failingCell=
            1*transactionNx*transactionNy+1*transactionNx+2;
        auto makeTransactionDomain=[&]()
        {
            auto result=std::make_unique<Domain>(
                transactionNx,transactionNy,transactionNz,
                0.04,0.03,0.03,registry.size(),vaporModel);
            PrimitiveVars successfulPrimitive(registry.size());
            successfulPrimitive.p=vaporState.state->p;
            successfulPrimitive.T=vaporState.state->T;
            successfulPrimitive.U={0.0,0.0,0.0};
            successfulPrimitive.Y=routeY;
            PrimitiveVars failingPrimitive=successfulPrimitive;
            failingPrimitive.p=failingVaporState.state->p;
            failingPrimitive.T=failingVaporState.state->T;
            failingPrimitive.Y=routeY;
            std::vector<PrimitiveVars> primitive(
                result->N(),successfulPrimitive);
            std::vector<thermo::EquilibriumState> equilibrium(
                result->N(),*vaporState.state);
            primitive[static_cast<std::size_t>(failingCell)]=
                failingPrimitive;
            equilibrium[static_cast<std::size_t>(failingCell)]=
                *failingVaporState.state;
            if(!result->replaceCompleteField(primitive,equilibrium))
                result.reset();
            return result;
        };

        VaporProbeChemistry timestepChemistry(registry);
        FiniteRateReactor timestepReactor(&timestepChemistry,registry);
        const double transactionStableDt=
            timestepReactor.estimateConservativeMaxDt(
                vaporState.state->rho,vaporState.state->p,
                vaporState.state->T,routeY,0.1);
        if(!std::isfinite(transactionStableDt)
            ||transactionStableDt<=0.0||transactionStableDt>=1.0e29)
        {
            std::cerr<<"transaction chemistry timestep failed\n";
            return 1;
        }
        const double transactionDt=1.0e-3*transactionStableDt;
        const double failureDensityThreshold=
            0.5*(vaporState.state->rho+failingVaporState.state->rho);

#ifdef _OPENMP
        const int previousThreadCount=omp_get_max_threads();
        omp_set_num_threads(1);
        const auto restoreThreadCount=[&]()
        { omp_set_num_threads(previousThreadCount); };
#else
        const auto restoreThreadCount=[]() {};
#endif

        auto directDomain=makeTransactionDomain();
        if(!directDomain)
        {
            restoreThreadCount();
            std::cerr<<"direct multi-cell chemistry Domain failed\n";
            return 1;
        }
        SelectiveFailChemistry directChemistry(
            registry,failureDensityThreshold);
        SectorChemistryOperator directOperator(
            &directChemistry,registry,0.0);
        directDomain->beginStep();
        const ConservedVars successfulCandidateBefore=
            directDomain->Unew(successfulCell);
        const ConservedVars failingCandidateBefore=
            directDomain->Unew(failingCell);
        bool directRejected=false;
        try
        {
            directOperator.apply(*directDomain,transactionDt);
        }
        catch(const std::runtime_error&)
        {
            directRejected=true;
        }
        if(!directRejected||directOperator.lastActiveCells()!=2
            ||directOperator.lastFailures()!=1
            ||directOperator.lastFailureI()!=2
            ||directOperator.lastFailureJ()!=1
            ||directOperator.lastFailureK()!=1
            ||directChemistry.successfulCalls()<=0
            ||directChemistry.failingCalls()<=0
            ||exactSame(
                directDomain->Unew(successfulCell),
                successfulCandidateBefore)
            ||!exactSame(
                directDomain->Unew(failingCell),
                failingCandidateBefore))
        {
            restoreThreadCount();
            std::cerr
                <<"direct multi-cell chemistry failure path was not exposed\n";
            return 1;
        }

        auto transactionDomain=makeTransactionDomain();
        if(!transactionDomain)
        {
            restoreThreadCount();
            std::cerr<<"transactional multi-cell chemistry Domain failed\n";
            return 1;
        }
        transactionDomain->Unew(0).rhoE+=17.0;
        const ExactDomainState transactionStart=
            captureExactState(*transactionDomain);
        SelectiveFailChemistry transactionChemistry(
            registry,failureDensityThreshold);
        std::vector<std::unique_ptr<Operator>> operators;
        auto transactionOperator=
            std::make_unique<SectorChemistryOperator>(
                &transactionChemistry,registry,0.0);
        SectorChemistryOperator* transactionOperatorView=
            transactionOperator.get();
        operators.push_back(std::move(transactionOperator));
        OperatorSplittingIntegrator integrator;
        bool transactionRejected=false;
        bool transactionReturned=false;
        try
        {
            transactionReturned=integrator.step(
                *transactionDomain,operators,transactionDt);
        }
        catch(const std::runtime_error&)
        {
            transactionRejected=true;
        }
        const bool transactionValid=
            transactionRejected&&!transactionReturned
            &&transactionDomain->tripped()
            &&!transactionDomain->stepTransactionActive()
            &&transactionDomain->acceptedSteps()==0
            &&transactionDomain->rejectedSteps()==1
            &&transactionDomain->failure().code
                ==Domain::FailureCode::OPERATOR_EXCEPTION
            &&transactionOperatorView->lastActiveCells()==2
            &&transactionOperatorView->lastFailures()==1
            &&transactionChemistry.successfulCalls()>0
            &&transactionChemistry.failingCalls()>0
            &&exactStateMatches(*transactionDomain,transactionStart);
        restoreThreadCount();
        if(!transactionValid)
        {
            std::cerr
                <<"chemistry-specific step transaction did not roll back\n";
            return 1;
        }
    }

    phaseEquilibrium::RhoEFlashOptions supercriticalOptions=twoPhaseOptions();
    supercriticalOptions.minimumPressure=1.0e7;
    supercriticalOptions.maximumPressure=3.0e7;
    supercriticalOptions.minimumTemperature=250.0;
    supercriticalOptions.maximumTemperature=500.0;
    supercriticalOptions.initialPressure=2.0e7;
    supercriticalOptions.initialTemperature=300.0;
    const auto supercriticalModel=
        std::make_shared<const thermo::PREquilibriumThermoModel>(
            critical,testKij,registry.names(),supercriticalOptions);
    const thermo::EquilibriumResult supercriticalState=
        supercriticalModel->initializePTY(2.0e7,300.0,routeY);
    if(!supercriticalState.success()
        ||supercriticalState.state->region
            !=thermo::EquilibriumRegion::SUPERCRITICAL
        ||!exerciseReactiveSinglePhaseRoute(
            supercriticalModel,*supercriticalState.state,registry,
            thermo::EquilibriumRegion::SUPERCRITICAL,"supercritical"))
    { std::cerr<<"supercritical chemistry route failed\n"; return 1; }

    const auto unclassifiedModel=
        std::make_shared<const UnclassifiedPhaseThermoModel>(vaporModel);
    const thermo::EquilibriumResult unclassifiedState=
        unclassifiedModel->initializePTY(1.0e5,300.0,routeY);
    if(!unclassifiedState.success()
        ||unclassifiedState.state->region
            !=thermo::EquilibriumRegion::UNCLASSIFIED_SINGLE_PHASE)
    { std::cerr<<"unclassified phase test state failed\n"; return 1; }
    PrimitiveVars unclassifiedPrimitive(registry.size());
    unclassifiedPrimitive.p=unclassifiedState.state->p;
    unclassifiedPrimitive.T=unclassifiedState.state->T;
    unclassifiedPrimitive.U={0.0,0.0,0.0};
    unclassifiedPrimitive.Y=routeY;
    Domain unclassifiedDomain(
        3,3,3,0.03,0.03,0.03,registry.size(),unclassifiedModel);
    if(!unclassifiedDomain.replaceCompleteField(
            std::vector<PrimitiveVars>(
                unclassifiedDomain.N(),unclassifiedPrimitive),
            std::vector<thermo::EquilibriumState>(
                unclassifiedDomain.N(),*unclassifiedState.state)))
    { std::cerr<<"unclassified phase Domain failed\n"; return 1; }
    VaporProbeChemistry unclassifiedChemistry(registry);
    SectorChemistryOperator unclassifiedOperator(
        &unclassifiedChemistry,registry,0.0);
    bool unclassifiedDtRejected=false;
    try
    {
        (void)unclassifiedOperator.maxStableDt(unclassifiedDomain,0.1);
    }
    catch(const std::runtime_error&)
    {
        unclassifiedDtRejected=true;
    }
    const ConservedVars unclassifiedBefore=unclassifiedDomain.U(1,1,1);
    unclassifiedDomain.beginStep();
    bool unclassifiedApplyRejected=false;
    try
    {
        unclassifiedOperator.apply(unclassifiedDomain,1.0e-6);
    }
    catch(const std::runtime_error&)
    {
        unclassifiedApplyRejected=true;
    }
    const ConservedVars& unclassifiedAfter=unclassifiedDomain.Unew(1,1,1);
    if(!unclassifiedDtRejected||!unclassifiedApplyRejected
        ||unclassifiedChemistry.observed()
        ||unclassifiedOperator.lastActiveCells()!=0
        ||unclassifiedOperator.lastStateFailures()!=1
        ||unclassifiedAfter.rho!=unclassifiedBefore.rho
        ||unclassifiedAfter.rhoU!=unclassifiedBefore.rhoU
        ||unclassifiedAfter.rhoE!=unclassifiedBefore.rhoE
        ||unclassifiedAfter.rhoY!=unclassifiedBefore.rhoY)
    { std::cerr<<"unclassified phase chemistry did not fail closed\n"; return 1; }

    const thermo::EquilibriumResult initialEquilibrium=
        equilibriumModel->initializePTY(2.0e6,300.0,twoPhaseY);
    if(!initialEquilibrium.success()
        ||initialEquilibrium.state->region
            !=thermo::EquilibriumRegion::TWO_PHASE
        ||!initialEquilibrium.state->liquid.has_value()
        ||!initialEquilibrium.state->vapor.has_value())
    { std::cerr<<"complete two-phase chemistry state failed\n"; return 1; }
    const thermo::EquilibriumState& initialQ=*initialEquilibrium.state;

    PrimitiveVars twoPhasePrimitive(registry.size());
    twoPhasePrimitive.p=initialQ.p;
    twoPhasePrimitive.T=initialQ.T;
    twoPhasePrimitive.U={0.0,0.0,0.0};
    twoPhasePrimitive.Y=twoPhaseY;
    auto makeTwoPhaseDomain=[&]()
    {
        auto result=std::make_unique<Domain>(
            3,3,3,0.03,0.03,0.03,registry.size(),equilibriumModel);
        if(!result->replaceCompleteField(
                std::vector<PrimitiveVars>(result->N(),twoPhasePrimitive),
                std::vector<thermo::EquilibriumState>(result->N(),initialQ)))
            throw std::runtime_error(
                "two-phase chemistry Domain initialization failed");
        return result;
    };

    VaporProbeChemistry referenceChemistry(registry);
    FiniteRateReactor referenceReactor(&referenceChemistry,registry);
    const double stableTwoPhaseDt=referenceReactor.estimateConservativeMaxDt(
        initialQ.vapor->rho,initialQ.p,initialQ.T,
        initialQ.vapor->massFractions,0.1);
    if(!std::isfinite(stableTwoPhaseDt)||stableTwoPhaseDt<=0.0
        ||stableTwoPhaseDt>=1.0e29)
    { std::cerr<<"two-phase vapor chemistry dt failed\n"; return 1; }
    const double twoPhaseDt=0.025*stableTwoPhaseDt;
    std::vector<double> reactedVaporY=initialQ.vapor->massFractions;
    double referenceTemperature=initialQ.T;
    double referencePressure=initialQ.p;
    double referenceDeltaRhoE=0.0;
    if(!referenceReactor.reactConservative(
            initialQ.vapor->rho,initialQ.vapor->e,
            reactedVaporY,twoPhaseDt,referenceDeltaRhoE,
            referenceTemperature,referencePressure))
    { std::cerr<<"reference vapor chemistry integration failed\n"; return 1; }

    // A shared phase-aware thermo owner cannot silently fall back to the
    // overall-state reactor for a vapor-only update. Without the explicit
    // phase recovery/displacement authority, both the timestep and source
    // paths reject before publishing a candidate increment.
    auto missingPhaseAuthority=makeTwoPhaseDomain();
    VaporProbeChemistry missingPhaseChemistry(registry);
    SectorChemistryOperator missingPhaseOperator(
        &missingPhaseChemistry,registry,equilibriumModel,0.0);
    bool missingPhaseDtRejected=false;
    try
    {
        (void)missingPhaseOperator.maxStableDt(*missingPhaseAuthority,0.1);
    }
    catch(const std::runtime_error&)
    {
        missingPhaseDtRejected=true;
    }
    missingPhaseAuthority->beginStep();
    const ConservedVars missingPhaseBefore=
        missingPhaseAuthority->Unew(1,1,1);
    bool missingPhaseApplyRejected=false;
    try
    {
        missingPhaseOperator.apply(*missingPhaseAuthority,twoPhaseDt);
    }
    catch(const std::runtime_error&)
    {
        missingPhaseApplyRejected=true;
    }
    if(!missingPhaseDtRejected||!missingPhaseApplyRejected
        ||missingPhaseOperator.hasPhaseThermoOwner()
        ||missingPhaseOperator.lastActiveCells()!=0
        ||missingPhaseOperator.lastStateFailures()!=1
        ||missingPhaseChemistry.observed()
        ||!exactSame(
            missingPhaseAuthority->Unew(1,1,1),missingPhaseBefore))
    {
        std::cerr<<"missing phase chemistry authority did not fail closed\n";
        return 1;
    }

    const phaseEquilibrium::PRPhaseRhoEOptions phaseRecoveryOptions=
        phaseChemistryRecoveryOptions();
    const phaseEquilibrium::PRPhaseRhoERecovery phaseRecovery(
        critical,testKij,registry.names(),phaseRecoveryOptions);
    const phaseEquilibrium::PRPhaseDisplacementAuthority phaseDisplacement(
        critical,testKij,registry.names(),phaseRecoveryOptions,
        phaseChemistryDisplacementBounds());
    if(!phaseRecovery.configured()||!phaseDisplacement.configured()
        ||phaseRecovery.physicalModelAuthorityFingerprint()
            !=phaseDisplacement.physicalModelAuthorityFingerprint()
        ||phaseRecovery.physicalModelAuthorityFingerprint()
            !=equilibriumModel->phaseLocalPhysicalModelFingerprint())
    {
        std::cerr<<"phase chemistry authorities do not share one model\n";
        return 1;
    }

    auto phaseOwnedStep=makeTwoPhaseDomain();
    const ConservedVars phaseOwnedInitial=phaseOwnedStep->U(1,1,1);
    VaporProbeChemistry phaseOwnedChemistry(registry);
    SectorChemistryOperator phaseOwnedOperator(
        &phaseOwnedChemistry,registry,equilibriumModel,
        phaseRecovery,phaseDisplacement,7U,
        "sector-chemistry-test",0.0);
    const double phaseOwnedStableDt=
        phaseOwnedOperator.maxStableDt(*phaseOwnedStep,0.1);
    if(!phaseOwnedOperator.hasPhaseThermoOwner()
        ||!std::isfinite(phaseOwnedStableDt)||phaseOwnedStableDt<=0.0
        ||phaseOwnedStableDt>=1.0e29)
    {
        std::cerr<<"phase-owned vapor chemistry dt failed\n";
        return 1;
    }
    const double phaseOwnedDt=1.0e-5*phaseOwnedStableDt;
    phaseOwnedStep->beginStep();
    phaseOwnedChemistry.resetObservation();
    try
    {
        phaseOwnedOperator.apply(*phaseOwnedStep,phaseOwnedDt);
    }
    catch(const std::runtime_error& error)
    {
        std::cerr<<"phase-owned vapor chemistry apply failed: "
                 <<error.what()<<" failureCode="
                 <<phaseOwnedOperator.lastStateFailures()
                 <<" cvode="<<phaseOwnedOperator.lastCvodeFailures()
                 <<" temperature="
                 <<phaseOwnedOperator.lastTemperatureFailures()
                 <<" mass="<<phaseOwnedOperator.lastMassFailures()
                 <<" element="<<phaseOwnedOperator.lastElementFailures()
                 <<" phaseAdmission="
                 <<phaseOwnedOperator.lastPhaseAdmissionFailures()
                 <<" ownerStatus="
                 <<phaseOwnedOperator.lastPhaseOwnerStatus()
                 <<" displacementStatus="
                 <<phaseOwnedOperator.lastPhaseDisplacementStatus()
                 <<" recoveryStatus="
                 <<phaseOwnedOperator.lastPhaseRecoveryStatus()
                 <<" failedWaypoint="
                 <<phaseOwnedOperator.lastPhaseFailedWaypoint()
                 <<" failedFraction="
                 <<phaseOwnedOperator.lastPhaseFailedFraction()
                 <<" energyDisplacement="
                 <<phaseOwnedOperator.lastPhaseEnergyDisplacement()
                 <<" compositionDisplacement="
                 <<phaseOwnedOperator.lastPhaseCompositionDisplacement()
                 <<" probes="<<phaseOwnedOperator.lastPhaseProbeCount()
                 <<" accepts="<<phaseOwnedOperator.lastPhaseAcceptCount()
                 <<" dt="<<phaseOwnedDt<<" stable="
                 <<phaseOwnedStableDt<<"\n";
        return 1;
    }
    if(phaseOwnedOperator.lastActiveCells()!=1
        ||phaseOwnedOperator.lastTwoPhaseActiveCells()!=1
        ||phaseOwnedOperator.lastFailures()!=0
        ||phaseOwnedOperator.lastPhaseReflashFailures()!=0
        ||!phaseOwnedChemistry.observed()
        ||phaseOwnedStep->Unew(1,1,1).rhoE<=phaseOwnedInitial.rhoE
        ||!phaseOwnedStep->commit()
        ||phaseOwnedStep->Q(1,1,1).region
            !=thermo::EquilibriumRegion::TWO_PHASE)
    {
        std::cerr<<"token-bound phase chemistry did not publish one reflash\n";
        return 1;
    }

    auto oneStep=makeTwoPhaseDomain();
    const ConservedVars oneStepInitial=oneStep->U(1,1,1);
    VaporProbeChemistry oneStepChemistry(registry);
    SectorChemistryOperator oneStepOperator(
        &oneStepChemistry,registry,0.0);
    oneStep->beginStep();
    oneStepChemistry.resetObservation();
    oneStepOperator.apply(*oneStep,twoPhaseDt);
    if(oneStepOperator.lastActiveCells()!=1
        ||oneStepOperator.lastTwoPhaseActiveCells()!=1
        ||oneStepOperator.lastFailures()!=0
        ||oneStepOperator.lastPhaseReflashFailures()!=0
        ||oneStepOperator.lastPositivityLimitedCells()!=0
        ||!oneStepChemistry.observed()
        ||!close(oneStepChemistry.observedTemperature(),initialQ.T)
        ||!close(oneStepChemistry.observedPressure(),initialQ.p)
        ||!close(oneStepChemistry.observedDensity(),initialQ.vapor->rho))
    { std::cerr<<"vapor-only chemistry routing failed\n"; return 1; }
    for(int chemistrySpecies=0;
        chemistrySpecies<oneStepChemistry.nSpecies();++chemistrySpecies)
    {
        const int globalSpecies=registry.index(
            oneStepChemistry.speciesNames()[chemistrySpecies]);
        if(!close(
                oneStepChemistry.observedY()[chemistrySpecies],
                initialQ.vapor->massFractions[globalSpecies]))
        { std::cerr<<"chemistry received non-vapor composition\n"; return 1; }
    }

    const double vaporMassDensity=
        initialQ.vaporVolumeFraction*initialQ.vapor->rho;
    for(int species=0;species<registry.size();++species)
    {
        const double expected=oneStepInitial.rhoY[species]
            +vaporMassDensity
                *(reactedVaporY[species]
                  -initialQ.vapor->massFractions[species]);
        if(!close(oneStep->Unew(1,1,1).rhoY[species],expected,2.0e-7,1.0e-9))
        { std::cerr<<"vapor species source scaling failed\n"; return 1; }
    }
    const double expectedEnergy=oneStepInitial.rhoE
        +initialQ.vaporVolumeFraction*referenceDeltaRhoE;
    if(!close(oneStep->Unew(1,1,1).rhoE,expectedEnergy,2.0e-7,1.0e-4))
    { std::cerr<<"vapor energy source scaling failed\n"; return 1; }

    // The HEM ownership transaction validates chemistry against the exact
    // committed U/Q pair. Its conservative delta must then accumulate into
    // Unew without consuming an increment already staged by another operator.
    auto additiveStep=makeTwoPhaseDomain();
    const ConservedVars additiveInitial=additiveStep->U(1,1,1);
    additiveStep->beginStep();
    const int additiveDonor=registry.index("CO");
    const int additiveReceiver=registry.index("CO2");
    const double priorMassIncrement=1.0e-8*additiveInitial.rho;
    const double priorEnergyIncrement=1.0e-8
        *std::max(std::abs(additiveInitial.rhoE),1.0);
    additiveStep->Unew(1,1,1).rhoY[additiveDonor]-=
        priorMassIncrement;
    additiveStep->Unew(1,1,1).rhoY[additiveReceiver]+=
        priorMassIncrement;
    additiveStep->Unew(1,1,1).rhoE+=priorEnergyIncrement;
    VaporProbeChemistry additiveChemistry(registry);
    SectorChemistryOperator additiveOperator(
        &additiveChemistry,registry,0.0);
    additiveOperator.apply(*additiveStep,twoPhaseDt);
    for(int species=0;species<registry.size();++species)
    {
        double priorIncrement=0.0;
        if(species==additiveDonor) priorIncrement=-priorMassIncrement;
        if(species==additiveReceiver) priorIncrement=priorMassIncrement;
        if(!close(
                additiveStep->Unew(1,1,1).rhoY[species],
                oneStep->Unew(1,1,1).rhoY[species]+priorIncrement,
                2.0e-7,1.0e-9))
        {
            std::cerr<<"two-phase chemistry replaced an additive species "
                     <<"increment\n";
            return 1;
        }
    }
    if(!close(
            additiveStep->Unew(1,1,1).rhoE,
            oneStep->Unew(1,1,1).rhoE+priorEnergyIncrement,
            2.0e-7,1.0e-4)
        ||!additiveStep->commit())
    {
        std::cerr<<"two-phase chemistry replaced an additive energy "
                 <<"increment or failed final reflash\n";
        return 1;
    }

    if(!oneStep->commit())
    { std::cerr<<"post-chemistry phase reflash failed\n"; return 1; }
    const thermo::EquilibriumState& committedQ=oneStep->Q(1,1,1);
    if(committedQ.region!=thermo::EquilibriumRegion::TWO_PHASE
        ||!committedQ.liquid.has_value()||!committedQ.vapor.has_value()
        ||!std::isfinite(committedQ.vaporMassFraction)
        ||!std::isfinite(committedQ.vaporVolumeFraction)
        ||committedQ.vaporMassFraction<=0.0
        ||committedQ.vaporMassFraction>=1.0
        ||committedQ.vaporVolumeFraction<=0.0
        ||committedQ.vaporVolumeFraction>=1.0)
    { std::cerr<<"post-chemistry phase fractions are invalid\n"; return 1; }
    std::vector<double> initialElements(4,0.0),committedElements(4,0.0);
    elementInventory(twoPhaseY,registry,initialElements.data());
    std::vector<double> committedY(registry.size(),0.0);
    for(int species=0;species<registry.size();++species)
        committedY[species]=oneStep->U(1,1,1).rhoY[species]
            /oneStep->U(1,1,1).rho;
    elementInventory(committedY,registry,committedElements.data());
    if(!nonnegativeNormalized(committedY))
    { std::cerr<<"post-chemistry species are invalid\n"; return 1; }
    for(int element=0;element<4;++element)
        if(!close(
                committedElements[element],initialElements[element],
                2.0e-7,1.0e-10))
        { std::cerr<<"post-chemistry element closure failed\n"; return 1; }
    if(!close(
            committedQ.e,
            oneStep->U(1,1,1).rhoE/oneStep->U(1,1,1).rho,
            2.0e-8,1.0e-5))
    { std::cerr<<"post-chemistry rho-e closure failed\n"; return 1; }

    auto twoHalfSteps=makeTwoPhaseDomain();
    VaporProbeChemistry halfStepChemistry(registry);
    SectorChemistryOperator halfStepOperator(
        &halfStepChemistry,registry,0.0);
    for(int half=0;half<2;++half)
    {
        twoHalfSteps->beginStep();
        halfStepOperator.apply(*twoHalfSteps,0.5*twoPhaseDt);
        if(!twoHalfSteps->commit())
        { std::cerr<<"half-step chemistry reflash failed\n"; return 1; }
    }
    if(!sameChemistryState(*oneStep,*twoHalfSteps,1,1,1))
    {
        const ConservedVars& one=oneStep->U(1,1,1);
        const ConservedVars& half=twoHalfSteps->U(1,1,1);
        std::cerr<<"phase chemistry timestep refinement failed"
                 <<" rhoE="<<one.rhoE<<"/"<<half.rhoE
                 <<" betaMass="<<oneStep->Q(1,1,1).vaporMassFraction
                 <<"/"<<twoHalfSteps->Q(1,1,1).vaporMassFraction
                 <<" betaVolume="<<oneStep->Q(1,1,1).vaporVolumeFraction
                 <<"/"<<twoHalfSteps->Q(1,1,1).vaporVolumeFraction<<"\n";
        for(int species=0;species<registry.size();++species)
            std::cerr<<"  "<<registry.name(species)<<" rhoY="
                     <<one.rhoY[species]<<"/"<<half.rhoY[species]<<"\n";
        return 1;
    }

    std::cout<<"SECTOR_CHEMISTRY_PASS dt="<<dt
             <<" active="<<chemistryOperator.lastActiveCells()<<"\n";
    return 0;
}
