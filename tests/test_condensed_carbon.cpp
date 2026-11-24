#include "src/combustion/CondensedCarbonCheckpoint.H"
#include "src/combustion/CondensedCarbonOutput.H"
#include "src/framework/PointImplicit.H"
#include "src/framework/SolverBase.H"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace Foam;

namespace
{

class TestSuite
{
    int checks_ = 0;

public:
    void require(bool condition,const std::string& message)
    {
        ++checks_;
        if(!condition) throw std::runtime_error(message);
    }

    template<class Callable>
    void requireThrows(Callable&& callable,const std::string& message)
    {
        bool threw=false;
        try { callable(); }
        catch(const std::exception&) { threw=true; }
        require(threw,message);
    }

    int checks() const noexcept { return checks_; }
};

bool near(double a,double b,double absolute=1.0e-12,double relative=1.0e-10)
{
    return std::isfinite(a)&&std::isfinite(b)
        &&std::abs(a-b)<=absolute+relative*std::max(std::abs(a),std::abs(b));
}

class TemporaryWorkspace
{
    std::filesystem::path root_;

public:
    TemporaryWorkspace()
    {
        const auto stamp=std::chrono::high_resolution_clock::now()
            .time_since_epoch().count();
        root_=std::filesystem::temp_directory_path()
            /("rfam_condensed_carbon_"+std::to_string(stamp));
        std::filesystem::create_directories(root_);
    }

    ~TemporaryWorkspace()
    {
        std::error_code error;
        std::filesystem::remove_all(root_,error);
    }

    std::filesystem::path path(const std::string& name) const
    {
        return root_/name;
    }

    std::string writeStraightProfile() const
    {
        return writeProfile("straight_profile.json",65.0,65.0);
    }

    std::string writeProfile(
        const std::string& name,double inletRadiusMm,double outletRadiusMm) const
    {
        const std::filesystem::path profile=path(name);
        std::ofstream out(profile);
        out<<"[{\"x_mm\":0.0,\"r_mm\":"<<inletRadiusMm<<"},"
              "{\"x_mm\":100.0,\"r_mm\":"<<outletRadiusMm<<"}]\n";
        if(!out) throw std::runtime_error("cannot write straight test profile");
        return profile.string();
    }
};

std::string readText(const std::filesystem::path& path)
{
    std::ifstream in(path,std::ios::binary);
    if(!in) throw std::runtime_error("cannot read test output: "+path.string());
    return std::string(
        std::istreambuf_iterator<char>(in),std::istreambuf_iterator<char>());
}

std::vector<FaceAreaVector> readVtkPoints(const std::string& text)
{
    std::istringstream input(text);
    std::string token;
    while(input>>token)
    {
        if(token!="POINTS") continue;
        std::size_t count=0;
        std::string numericType;
        if(!(input>>count>>numericType)||numericType!="double")
            throw std::runtime_error("invalid condensed-carbon VTK point header");
        std::vector<FaceAreaVector> points(count);
        for(FaceAreaVector& point:points)
            if(!(input>>point.x>>point.y>>point.z))
                throw std::runtime_error("truncated condensed-carbon VTK points");
        return points;
    }
    throw std::runtime_error("condensed-carbon VTK lacks points");
}

CondensedCarbonOptions baseOptions(double rhoSoot=0.02)
{
    CondensedCarbonOptions options;
    options.momentModel=CondensedCarbonMomentModel::ONE_MOMENT;
    options.applicabilityPolicy=
        CondensedCarbonApplicabilityPolicy::DIAGNOSTIC_ONLY;
    options.carbonMolarMass=0.0120107;
    options.carbonFormationEnthalpy=0.0;
    options.materialDensity=1800.0;
    options.diagnosticSpecificHeat=710.0;
    options.oneMomentParticleDiameter=50.0e-9;
    options.initialState.rhoSoot=rhoSoot;
    options.inletState=options.initialState;
    options.limits.maximumMassLoading=0.2;
    options.limits.maximumVolumeFraction=1.0e-3;
    options.limits.maximumHeatCapacityRatio=0.1;
    options.limits.maximumStokesNumber=0.1;
    options.productionValidated=false;
    options.validate();
    return options;
}

std::unique_ptr<Domain> makeSectorDomain(
    const SectorMeshMetrics& mesh,
    const SpeciesRegistry& registry,
    double axialVelocity=0.0,
    double temperature=1000.0)
{
    std::unique_ptr<Domain> domain(new Domain(
        mesh.Nx(),mesh.Nr(),mesh.Ntheta(),mesh.length(),0.065,
        mesh.sectorAngle(),registry.size()));
    domain->setSpecies(registry.names());
    std::vector<double> composition(
        static_cast<std::size_t>(registry.size()),0.0);
    composition[static_cast<std::size_t>(registry.index("O2"))]=1.0;
    for(int i=0;i<mesh.Nx();++i)
    for(int j=0;j<mesh.Nr();++j)
    for(int k=0;k<mesh.Ntheta();++k)
    {
        PrimitiveVars& state=domain->stagedW(i,j,k);
        state.p=2.0e5;
        state.T=temperature;
        state.U={axialVelocity,0.0,0.0};
        state.Y=composition;
        state.Z=0.0;
        state.C=0.0;
    }
    if(!domain->primitiveToConserved())
        throw std::runtime_error("cannot initialize condensed-carbon sector");
    return domain;
}

std::unique_ptr<Domain> makeCellDomain(
    const SpeciesRegistry& registry,double axialVelocity=0.0)
{
    std::unique_ptr<Domain> domain(new Domain(1,1.0,registry.size()));
    domain->setSpecies(registry.names());
    PrimitiveVars& state=domain->stagedW(0);
    state.p=2.0e5;
    state.T=1000.0;
    state.U={axialVelocity,0.0,0.0};
    state.Y.assign(static_cast<std::size_t>(registry.size()),0.0);
    state.Y[static_cast<std::size_t>(registry.index("O2"))]=1.0;
    state.Z=0.0;
    state.C=0.0;
    if(!domain->primitiveToConserved())
        throw std::runtime_error("cannot initialize condensed-carbon cell");
    return domain;
}

CondensedCarbonExchange oxidationExchange(
    const PrimitiveVars& fluid,
    const std::vector<std::string>& speciesNames,
    const SpeciesRegistry& registry,
    const CondensedCarbonOptions& options,
    double molarRate)
{
    CondensedCarbonExchange exchange;
    exchange.fluidSpeciesMassRate.assign(speciesNames.size(),0.0);
    const int localO2=static_cast<int>(std::find(
        speciesNames.begin(),speciesNames.end(),"O2")-speciesNames.begin());
    const int localCO2=static_cast<int>(std::find(
        speciesNames.begin(),speciesNames.end(),"CO2")-speciesNames.begin());
    if(localO2>=static_cast<int>(speciesNames.size())
       ||localCO2>=static_cast<int>(speciesNames.size()))
        throw std::runtime_error("test species layout lacks O2 or CO2");
    exchange.fluidSpeciesMassRate[static_cast<std::size_t>(localO2)]=
        -molarRate*registry.molecularWeight(registry.index("O2"));
    exchange.fluidSpeciesMassRate[static_cast<std::size_t>(localCO2)]=
        molarRate*registry.molecularWeight(registry.index("CO2"));
    exchange.sootMassRate=-molarRate*options.carbonMolarMass;
    exchange.sootNumberRate=0.0;
    exchange.fluidEnergyRate=
        CondensedCarbonExchangeContract::expectedFluidEnergyRate(
            fluid,speciesNames,registry,options,exchange);
    return exchange;
}

long double modeledEnergyDensity(
    const ConservedVars& fluid,
    const PrimitiveVars& primitive,
    const CondensedCarbonState& carbon,
    const std::vector<std::string>& speciesNames,
    const SpeciesRegistry& registry,
    const CondensedCarbonOptions& options)
{
    long double value=fluid.rhoE;
    for(std::size_t species=0;species<speciesNames.size();++species)
    {
        const int global=registry.index(speciesNames[species]);
        value+=static_cast<long double>(fluid.rhoY[species])
            /registry.molecularWeight(global)
            *registry.formationEnthalpy(global);
    }
    value+=static_cast<long double>(carbon.rhoSoot)
        /options.carbonMolarMass*options.carbonFormationEnthalpy;
    long double velocitySquared=0.0L;
    for(double velocity:primitive.U)
        velocitySquared+=static_cast<long double>(velocity)*velocity;
    value+=0.5L*carbon.rhoSoot*velocitySquared;
    return value;
}

class ConstantOxidationKinetics : public CondensedCarbonKinetics
{
    const SpeciesRegistry& registry_;
    std::vector<std::string> speciesNames_;
    CondensedCarbonOptions options_;
    double molarRate_;

public:
    ConstantOxidationKinetics(
        const SpeciesRegistry& registry,
        std::vector<std::string> speciesNames,
        CondensedCarbonOptions options,
        double molarRate)
    :
        registry_(registry),
        speciesNames_(std::move(speciesNames)),
        options_(std::move(options)),
        molarRate_(molarRate)
    {}

    std::string modelId() const override
    { return "test.constant-c-oxidation/v1"; }

    CondensedCarbonExchange evaluate(
        const CondensedCarbonCellContext& context) const override
    {
        return oxidationExchange(
            context.fluidPrimitive,speciesNames_,registry_,options_,molarRate_);
    }
};

class ConstantWallDeposition : public WallCarbonDepositionModel
{
    double massFlux_;
    double heatFlux_;

public:
    explicit ConstantWallDeposition(double massFlux,double heatFlux=0.0)
    : massFlux_(massFlux),heatFlux_(heatFlux)
    {}

    std::string modelId() const override
    { return "test.constant-wall-deposition/v1"; }

    WallCarbonTransfer evaluate(
        const WallCarbonCellContext&) const override
    {
        WallCarbonTransfer transfer;
        transfer.sootMassFluxToWall=massFlux_;
        transfer.wallHeatFlux=heatFlux_;
        return transfer;
    }
};

class SwitchableFailure : public Operator
{
    bool enabled_;

public:
    explicit SwitchableFailure(bool enabled):enabled_(enabled) {}
    void setEnabled(bool enabled) { enabled_=enabled; }

    void apply(Domain&,double) const override
    {
        if(enabled_) throw std::runtime_error("injected downstream failure");
    }

    std::string name() const override
    { return "condensed-carbon downstream failure"; }
};

void testExchangeContract(TestSuite& test,const SpeciesRegistry& registry)
{
    const CondensedCarbonOptions options=baseOptions(0.05);
    std::unique_ptr<Domain> domain=makeCellDomain(registry,125.0);
    const CondensedCarbonState carbonBase=options.initialState;
    const CondensedCarbonExchange exchange=oxidationExchange(
        domain->W(0),domain->speciesNames(),registry,options,1.0);
    const double dt=1.0e-4;
    const CondensedCarbonExchangeResult accepted=
        CondensedCarbonExchangeContract::evaluate(
            domain->U(0),domain->W(0),carbonBase,domain->speciesNames(),
            registry,options,exchange,dt);
    test.require(accepted.accepted,"balanced carbon exchange was rejected");
    test.require(near(
        accepted.fluidCandidate.rho+accepted.carbonCandidate.rhoSoot,
        domain->U(0).rho+carbonBase.rhoSoot),
        "balanced carbon exchange changed material mass");
    for(int direction=0;direction<3;++direction)
        test.require(near(
            accepted.fluidCandidate.rhoU[direction]
                /accepted.fluidCandidate.rho,
            domain->W(0).U[direction]),
            "mass exchange changed the zero-slip fluid velocity");
    const long double energyBefore=modeledEnergyDensity(
        domain->U(0),domain->W(0),carbonBase,domain->speciesNames(),
        registry,options);
    const long double energyAfter=modeledEnergyDensity(
        accepted.fluidCandidate,domain->W(0),accepted.carbonCandidate,
        domain->speciesNames(),registry,options);
    test.require(std::abs(energyAfter-energyBefore)
            <1.0e-8L*std::max(1.0L,std::abs(energyBefore)),
        "balanced carbon exchange changed the modeled energy ledger");

    CondensedCarbonExchange badMass=exchange;
    badMass.sootMassRate+=1.0e-3;
    test.require(!CondensedCarbonExchangeContract::evaluate(
        domain->U(0),domain->W(0),carbonBase,domain->speciesNames(),registry,
        options,badMass,dt).accepted,"unbalanced mass exchange was accepted");

    CondensedCarbonExchange badElement=exchange;
    const int o2=registry.index("O2");
    const int co2=registry.index("CO2");
    badElement.fluidSpeciesMassRate[static_cast<std::size_t>(o2)]+=1.0e-3;
    badElement.fluidSpeciesMassRate[static_cast<std::size_t>(co2)]-=1.0e-3;
    test.require(!CondensedCarbonExchangeContract::evaluate(
        domain->U(0),domain->W(0),carbonBase,domain->speciesNames(),registry,
        options,badElement,dt).accepted,"element-unbalanced exchange was accepted");

    CondensedCarbonExchange badEnergy=exchange;
    badEnergy.fluidEnergyRate+=100.0;
    test.require(!CondensedCarbonExchangeContract::evaluate(
        domain->U(0),domain->W(0),carbonBase,domain->speciesNames(),registry,
        options,badEnergy,dt).accepted,"energy-unbalanced exchange was accepted");

    CondensedCarbonExchange badNumber=exchange;
    badNumber.sootNumberRate=1.0;
    test.require(!CondensedCarbonExchangeContract::evaluate(
        domain->U(0),domain->W(0),carbonBase,domain->speciesNames(),registry,
        options,badNumber,dt).accepted,"one-moment number source was accepted");

    CondensedCarbonExchange aggregateOverflow=exchange;
    for(double& rate:aggregateOverflow.fluidSpeciesMassRate) rate*=5.0e307;
    aggregateOverflow.sootMassRate*=5.0e307;
    aggregateOverflow.fluidEnergyRate=0.0;
    test.require(!CondensedCarbonExchangeContract::evaluate(
        domain->U(0),domain->W(0),carbonBase,domain->speciesNames(),registry,
        options,aggregateOverflow,dt).accepted,
        "finite-input aggregate overflow was accepted");
}

void testConfigurationAndApplicability(TestSuite& test)
{
    CondensedCarbonOptions oneMoment=baseOptions();
    oneMoment.initialState.numberDensity=1.0;
    test.requireThrows([&]{ oneMoment.validate(); },
        "one-moment state accepted number density");

    CondensedCarbonOptions twoMoment=baseOptions();
    twoMoment.momentModel=CondensedCarbonMomentModel::TWO_MOMENT;
    twoMoment.initialState.numberDensity=0.0;
    test.requireThrows([&]{ twoMoment.validate(); },
        "two-moment state accepted mass without number density");
    twoMoment.initialState.numberDensity=1.0e12;
    twoMoment.inletState=twoMoment.initialState;
    twoMoment.validate();
    test.require(twoMoment.initialState.valid(twoMoment.momentModel),
        "valid two-moment state was rejected");

    CondensedCarbonOptions unsupported=baseOptions();
    unsupported.thermophoresisEnabled=true;
    test.requireThrows([&]{ unsupported.validate(); },
        "v1 accepted unsupported thermophoresis");

    CondensedCarbonOptions options=baseOptions();
    options.productionValidated=true;
    options.limits.maximumMassLoading=1.0;
    options.limits.maximumVolumeFraction=1.0;
    options.limits.maximumHeatCapacityRatio=0.01;
    options.limits.maximumStokesNumber=1.0;
    CondensedCarbonApplicabilityInputs inputs;
    inputs.rhoFluid=1.0;
    inputs.rhoSoot=1.0e-3;
    inputs.fluidCp=1000.0;
    inputs.dynamicViscosity=2.0e-5;
    inputs.fluidTimeScale=1.0e-3;
    inputs.particleDiameter=50.0e-9;
    test.require(evaluateCondensedCarbonApplicability(options,inputs).status
            ==CondensedCarbonApplicabilityStatus::SUPPORTED,
        "in-range validated soot state was not supported");
    inputs.rhoSoot=2.0;
    test.require(evaluateCondensedCarbonApplicability(options,inputs).status
            ==CondensedCarbonApplicabilityStatus::OUTSIDE_DILUTE_LIMIT,
        "high soot loading did not leave the dilute domain");
    inputs.rhoSoot=1.0e-3;
    inputs.particleDiameter=1.0e-3;
    test.require(evaluateCondensedCarbonApplicability(options,inputs).status
            ==CondensedCarbonApplicabilityStatus::SLIP_MODEL_REQUIRED,
        "large particle response time did not require slip");
    inputs.rhoSoot=0.1;
    inputs.particleDiameter=50.0e-9;
    test.require(evaluateCondensedCarbonApplicability(options,inputs).status
            ==CondensedCarbonApplicabilityStatus::SOOT_ENERGY_MODEL_REQUIRED,
        "large soot heat capacity did not require an energy model");
    options.productionValidated=false;
    inputs.rhoSoot=1.0e-3;
    test.require(evaluateCondensedCarbonApplicability(options,inputs).status
            ==CondensedCarbonApplicabilityStatus::WARNING_RESEARCH_ONLY,
        "unvalidated in-range model did not remain research-only");
}

void testUniformAdvection(
    TestSuite& test,const SectorMeshMetrics& mesh,
    const SpeciesRegistry& registry)
{
    CondensedCarbonOptions options=baseOptions(0.002);
    std::unique_ptr<Domain> domain=makeSectorDomain(mesh,registry,40.0);
    std::unique_ptr<SectorCondensedCarbonOperator> carbon(
        new SectorCondensedCarbonOperator(mesh,registry,options));
    SectorCondensedCarbonOperator* carbonPtr=carbon.get();
    SolverBase solver(*domain);
    solver.addOperator(std::move(carbon));
    test.require(solver.step(1.0e-5),"uniform soot advection step failed");
    for(int i=1;i<domain->Nx()-1;++i)
    for(int j=1;j<domain->Ny()-1;++j)
    for(int k=1;k<domain->Nz()-1;++k)
        test.require(near(carbonPtr->rhoSoot(i,j,k),0.002,1.0e-14,1.0e-12),
            "straight-duct uniform soot field drifted");
    test.require(near(
        carbonPtr->sootBoundaryMassNetIntoDomain(),0.0,1.0e-14,1.0e-12),
        "uniform soot advection produced a net boundary mass imbalance");
    const CondensedCarbonInventoryAudit audit=carbonPtr->audit(*domain);
    test.require(near(audit.sootMass,0.002*mesh.totalVolume(),1.0e-14,1.0e-11),
        "uniform soot inventory audit is incorrect");

    CondensedCarbonOptions twoMoment=options;
    twoMoment.momentModel=CondensedCarbonMomentModel::TWO_MOMENT;
    twoMoment.initialState.numberDensity=1.0e12;
    twoMoment.inletState=twoMoment.initialState;
    std::unique_ptr<Domain> twoMomentDomain=makeSectorDomain(
        mesh,registry,40.0);
    std::unique_ptr<SectorCondensedCarbonOperator> twoMomentCarbon(
        new SectorCondensedCarbonOperator(mesh,registry,twoMoment));
    SectorCondensedCarbonOperator* twoMomentPtr=twoMomentCarbon.get();
    SolverBase twoMomentSolver(*twoMomentDomain);
    twoMomentSolver.addOperator(std::move(twoMomentCarbon));
    test.require(twoMomentSolver.step(1.0e-5),
        "uniform two-moment soot advection step failed");
    test.require(near(twoMomentPtr->rhoSoot(2,2,2),0.002,1.0e-14,1.0e-12),
        "two-moment soot mass drifted under uniform advection");
    test.require(near(
        twoMomentPtr->numberDensity(2,2,2),1.0e12,1.0e-3,1.0e-12),
        "two-moment number density drifted under uniform advection");
    test.require(near(
        twoMomentPtr->sootBoundaryNumberNetIntoDomain(),0.0,1.0e-4,1.0e-12),
        "uniform two-moment advection changed the number boundary ledger");
}

void testUnownedBoundaryFailClose(
    TestSuite& test,const SectorMeshMetrics& mesh,
    const SpeciesRegistry& registry)
{
    CondensedCarbonOptions options=baseOptions();
    {
        std::unique_ptr<Domain> domain=makeSectorDomain(mesh,registry);
        for(int i=0;i<domain->Nx();++i)
        for(int j=0;j<domain->Ny();++j)
        for(int k=0;k<domain->Nz();++k)
        {
            PrimitiveVars state=domain->W(i,j,k);
            const double theta=mesh.thetaCenter(k);
            state.U={0.0,-25.0*std::sin(theta),25.0*std::cos(theta)};
            test.require(domain->setCell(i,j,k,state),
                "cannot construct rotationally periodic swirl state");
        }
        domain->applyCyclicZ(120.0);
        test.require(!domain->tripped(),
            "rotationally periodic swirl boundary setup failed");
        std::unique_ptr<SectorCondensedCarbonOperator> carbon(
            new SectorCondensedCarbonOperator(mesh,registry,options));
        SolverBase solver(*domain);
        solver.addOperator(std::move(carbon));
        test.require(solver.step(1.0e-6),
            "valid rotationally periodic soot transport was rejected");
    }

    {
        std::unique_ptr<Domain> domain=makeSectorDomain(mesh,registry);
        const int i=2;
        const int k=2;
        const int wallGhost=domain->Ny()-1;
        PrimitiveVars ghost=domain->W(i,wallGhost,k);
        const double theta=mesh.thetaCenter(k);
        ghost.U={0.0,std::cos(theta),std::sin(theta)};
        test.require(domain->setCell(i,wallGhost,k,ghost),
            "cannot construct nonzero wall-normal flux counterexample");
        std::unique_ptr<SectorCondensedCarbonOperator> carbon(
            new SectorCondensedCarbonOperator(mesh,registry,options));
        SectorCondensedCarbonOperator* carbonPtr=carbon.get();
        SolverBase solver(*domain);
        solver.addOperator(std::move(carbon));
        test.requireThrows([&]{ (void)solver.step(1.0e-5); },
            "condensed carbon accepted unowned wall-normal flux");
        test.require(domain->failure().detail.find("WALL_NORMAL_VOLUME_FLUX")
                !=std::string::npos,
            "wall-normal rejection did not identify its ownership failure");
        test.require(!carbonPtr->persistentStateInitialized(),
            "wall-normal rejection published condensed-carbon state");
    }

    {
        std::unique_ptr<Domain> domain=makeSectorDomain(mesh,registry);
        const int i=2;
        const int j=2;
        PrimitiveVars lowerGhost=domain->W(i,j,0);
        lowerGhost.U={0.0,0.0,-1.0};
        test.require(domain->setCell(i,j,0,lowerGhost),
            "cannot construct periodic-flux counterexample");
        std::unique_ptr<SectorCondensedCarbonOperator> carbon(
            new SectorCondensedCarbonOperator(mesh,registry,options));
        SectorCondensedCarbonOperator* carbonPtr=carbon.get();
        SolverBase solver(*domain);
        solver.addOperator(std::move(carbon));
        test.requireThrows([&]{ (void)solver.step(1.0e-5); },
            "condensed carbon accepted unmatched periodic volume flux");
        test.require(domain->failure().detail.find(
                "PERIODIC_VOLUME_FLUX_MISMATCH")!=std::string::npos,
            "periodic rejection did not identify its ownership failure");
        test.require(!carbonPtr->persistentStateInitialized(),
            "periodic-flux rejection published condensed-carbon state");
    }
}

struct ReactiveCaseResult
{
    double temperature=0.0;
    double axialVelocity=0.0;
    CondensedCarbonInventoryAudit before;
    CondensedCarbonInventoryAudit after;
};

ReactiveCaseResult runReactiveCase(
    const SectorMeshMetrics& mesh,
    const SpeciesRegistry& registry,
    double axialVelocity)
{
    CondensedCarbonOptions options=baseOptions(0.05);
    std::unique_ptr<Domain> domain=makeSectorDomain(
        mesh,registry,axialVelocity);
    ConstantOxidationKinetics kinetics(
        registry,domain->speciesNames(),options,0.5);
    std::unique_ptr<SectorCondensedCarbonOperator> carbon(
        new SectorCondensedCarbonOperator(
            mesh,registry,options,&kinetics,nullptr));
    SectorCondensedCarbonOperator* carbonPtr=carbon.get();
    ReactiveCaseResult result;
    result.before=carbonPtr->audit(*domain);
    SolverBase solver(*domain);
    solver.addOperator(std::move(carbon));
    if(!solver.step(1.0e-4))
        throw std::runtime_error("reactive condensed-carbon step failed");
    result.after=carbonPtr->audit(*domain);
    result.temperature=domain->W(2,2,2).T;
    result.axialVelocity=domain->W(2,2,2).U[0];
    return result;
}

void testGalileanEnergyAudit(
    TestSuite& test,const SectorMeshMetrics& mesh,
    const SpeciesRegistry& registry)
{
    const ReactiveCaseResult stationary=runReactiveCase(mesh,registry,0.0);
    const ReactiveCaseResult shifted=runReactiveCase(mesh,registry,250.0);
    test.require(near(stationary.temperature,shifted.temperature,1.0e-8,1.0e-11),
        "condensed-carbon energy source is not Galilean temperature invariant");
    test.require(near(shifted.axialVelocity,250.0,1.0e-9,1.0e-11),
        "zero-slip mass transfer changed the fluid velocity");
    for(const ReactiveCaseResult* result:{&stationary,&shifted})
    {
        test.require(near(
            result->before.materialMass,result->after.materialMass,
            1.0e-11,1.0e-10),
            "closed-box condensed-carbon reaction changed material mass");
        test.require(near(
            result->before.modeledEnergyLedger,
            result->after.modeledEnergyLedger,1.0e-7,1.0e-9),
            "closed-box condensed-carbon reaction changed modeled energy");
        for(int element=0;element<4;++element)
            test.require(near(
                result->before.totalElementMoles[element],
                result->after.totalElementMoles[element],1.0e-11,1.0e-10),
                "closed-box condensed-carbon reaction changed element inventory");
    }
}

bool sameCarbonField(
    const std::vector<CondensedCarbonState>& a,
    const std::vector<CondensedCarbonState>& b)
{
    if(a.size()!=b.size()) return false;
    for(std::size_t i=0;i<a.size();++i)
        if(a[i].rhoSoot!=b[i].rhoSoot
           ||a[i].numberDensity!=b[i].numberDensity)
            return false;
    return true;
}

void testLifecycle(
    TestSuite& test,const SectorMeshMetrics& mesh,
    const SpeciesRegistry& registry)
{
    {
        CondensedCarbonOptions options=baseOptions(0.05);
        std::unique_ptr<Domain> domain=makeSectorDomain(mesh,registry);
        ConstantOxidationKinetics kinetics(
            registry,domain->speciesNames(),options,0.25);
        std::unique_ptr<SectorCondensedCarbonOperator> carbon(
            new SectorCondensedCarbonOperator(
                mesh,registry,options,&kinetics,nullptr));
        SectorCondensedCarbonOperator* carbonPtr=carbon.get();
        SolverBase solver(*domain);
        solver.setIntegrator(std::unique_ptr<TimeIntegrator>(new RK2Integrator()));
        solver.addOperator(std::move(carbon));
        const double dt=1.0e-4;
        test.require(solver.step(dt),"RK2 condensed-carbon step failed");
        const double expected=options.initialState.rhoSoot
            -dt*0.25*options.carbonMolarMass;
        test.require(near(carbonPtr->rhoSoot(2,2,2),expected,1.0e-12,1.0e-10),
            "RK2 condensed-carbon source did not average from the lazy start");
    }

    {
        CondensedCarbonOptions options=baseOptions(0.05);
        std::unique_ptr<Domain> domain=makeSectorDomain(mesh,registry);
        ConstantOxidationKinetics kinetics(
            registry,domain->speciesNames(),options,0.25);
        std::unique_ptr<SectorCondensedCarbonOperator> carbon(
            new SectorCondensedCarbonOperator(
                mesh,registry,options,&kinetics,nullptr));
        SectorCondensedCarbonOperator* carbonPtr=carbon.get();
        std::unique_ptr<SwitchableFailure> failure(new SwitchableFailure(false));
        SwitchableFailure* failurePtr=failure.get();
        SolverBase solver(*domain);
        solver.addOperator(std::move(carbon));
        solver.addOperator(std::move(failure));
        test.require(solver.step(1.0e-4),"rollback baseline step failed");
        const std::vector<CondensedCarbonState> before=carbonPtr->stateField();
        failurePtr->setEnabled(true);
        test.requireThrows([&]{ (void)solver.step(1.0e-4); },
            "downstream failure was not propagated");
        test.require(sameCarbonField(before,carbonPtr->stateField()),
            "downstream failure did not roll back condensed-carbon state");
        test.require(!carbonPtr->stepTransactionActive(),
            "condensed-carbon rollback left a transaction active");
    }

    {
        CondensedCarbonOptions options=baseOptions();
        std::unique_ptr<Domain> domain=makeSectorDomain(mesh,registry);
        std::unique_ptr<SectorCondensedCarbonOperator> carbon(
            new SectorCondensedCarbonOperator(mesh,registry,options));
        SectorCondensedCarbonOperator* carbonPtr=carbon.get();
        std::unique_ptr<SwitchableFailure> failure(new SwitchableFailure(true));
        SolverBase solver(*domain);
        solver.addOperator(std::move(carbon));
        solver.addOperator(std::move(failure));
        test.requireThrows([&]{ (void)solver.step(1.0e-4); },
            "lazy-initialization failure was not propagated");
        test.require(!carbonPtr->persistentStateInitialized(),
            "failed first step published lazy condensed-carbon state");
    }

    {
        CondensedCarbonOptions options=baseOptions();
        std::unique_ptr<Domain> domain=makeSectorDomain(mesh,registry);
        std::unique_ptr<SectorCondensedCarbonOperator> carbon(
            new SectorCondensedCarbonOperator(mesh,registry,options));
        SectorCondensedCarbonOperator* carbonPtr=carbon.get();
        SolverBase solver(*domain);
        solver.setIntegrator(std::unique_ptr<TimeIntegrator>(
            new PointImplicitIntegrator()));
        solver.addOperator(std::move(carbon));
        test.require(!solver.step(1.0e-4),
            "point-implicit accepted condensed-carbon persistent state");
        test.require(domain->failure().code
                ==Domain::FailureCode::TRANSACTION_PROTOCOL,
            "point-implicit rejection used the wrong failure class");
        test.require(!carbonPtr->persistentStateInitialized(),
            "point-implicit rejection initialized condensed-carbon state");
    }
}

void testWallTransfer(
    TestSuite& test,const SectorMeshMetrics& mesh,
    const SpeciesRegistry& registry)
{
    CondensedCarbonOptions options=baseOptions(0.02);
    ConstantWallDeposition wallModel(1.0e-5);
    std::unique_ptr<Domain> domain=makeSectorDomain(mesh,registry);
    std::unique_ptr<SectorCondensedCarbonOperator> carbon(
        new SectorCondensedCarbonOperator(
            mesh,registry,options,nullptr,&wallModel));
    SectorCondensedCarbonOperator* carbonPtr=carbon.get();
    const CondensedCarbonInventoryAudit before=carbonPtr->audit(*domain);
    SolverBase solver(*domain);
    solver.addOperator(std::move(carbon));
    test.require(solver.step(1.0e-2),"wall deposition step failed");
    const CondensedCarbonInventoryAudit after=carbonPtr->audit(*domain);
    test.require(after.wallCarbonMass>0.0,
        "wall deposition did not increase wall carbon inventory");
    test.require(after.sootMass<before.sootMass,
        "wall deposition did not decrease suspended soot");
    test.require(near(
        before.sootMass,after.sootMass+after.wallCarbonMass,
        1.0e-12,1.0e-10),
        "soot-to-wall transfer was not an atomic mass exchange");
    test.require(near(before.materialMass,after.materialMass,1.0e-12,1.0e-10),
        "wall deposition changed the system material inventory");
    test.require(!after.wallDepositionEnergyModeled,
        "v1 incorrectly claimed deposition energy closure");

    ConstantWallDeposition unsupportedHeat(0.0,1.0e-13);
    std::unique_ptr<Domain> rejectedDomain=makeSectorDomain(mesh,registry);
    std::unique_ptr<SectorCondensedCarbonOperator> rejectedCarbon(
        new SectorCondensedCarbonOperator(
            mesh,registry,options,nullptr,&unsupportedHeat));
    SectorCondensedCarbonOperator* rejectedCarbonPtr=rejectedCarbon.get();
    SolverBase rejectedSolver(*rejectedDomain);
    rejectedSolver.addOperator(std::move(rejectedCarbon));
    test.require(!rejectedSolver.step(1.0e-2),
        "v1 accepted an unmodeled nonzero wall heat flux");
    test.require(!rejectedCarbonPtr->persistentStateInitialized(),
        "rejected wall heat flux published condensed-carbon state");
}

void testCheckpointAndOutput(
    TestSuite& test,const TemporaryWorkspace& temporary,
    const SectorMeshMetrics& mesh,const SpeciesRegistry& registry)
{
    CondensedCarbonOptions options=baseOptions(0.01);
    std::unique_ptr<Domain> domain=makeSectorDomain(mesh,registry);
    SectorCondensedCarbonOperator source(mesh,registry,options);
    std::vector<CondensedCarbonState> state(
        static_cast<std::size_t>(domain->N()),options.initialState);
    std::vector<WallCarbonInventory> wall(
        static_cast<std::size_t>(domain->Nx())*domain->Nz());
    const int cell=2*domain->Nx()*domain->Ny()+2*domain->Nx()+2;
    state[static_cast<std::size_t>(cell)].rhoSoot=0.0125;
    wall[static_cast<std::size_t>(2*domain->Nx()+2)].massPerArea=0.3;
    source.restoreFields(state,wall,1.25,0.0);

    const std::filesystem::path checkpoint=temporary.path("carbon.rcc");
    CondensedCarbonCheckpointMetadata metadata;
    metadata.time=2.5e-3;
    metadata.iteration=42;
    metadata.coupledFluidCheckpointId="fluid-checkpoint-sha256:test";
    CondensedCarbonCheckpoint::write(
        checkpoint.string(),*domain,source,metadata);

    SectorCondensedCarbonOperator restored(mesh,registry,options);
    const CondensedCarbonCheckpointMetadata loaded=
        CondensedCarbonCheckpoint::read(
            checkpoint.string(),*domain,restored,
            metadata.coupledFluidCheckpointId);
    test.require(loaded.time==metadata.time&&loaded.iteration==metadata.iteration,
        "condensed-carbon checkpoint metadata did not round trip");
    test.require(sameCarbonField(source.stateField(),restored.stateField()),
        "condensed-carbon checkpoint field did not round trip");
    test.require(restored.wallField()[static_cast<std::size_t>(
            2*domain->Nx()+2)].massPerArea==0.3,
        "wall carbon checkpoint field did not round trip");
    test.require(restored.sootBoundaryMassNetIntoDomain()==1.25,
        "condensed-carbon boundary ledger did not round trip");

    SectorCondensedCarbonOperator wrongCoupling(mesh,registry,options);
    test.requireThrows([&]{
        (void)CondensedCarbonCheckpoint::read(
            checkpoint.string(),*domain,wrongCoupling,"wrong-fluid-id");
    },"checkpoint accepted the wrong coupled fluid identity");
    test.require(!wrongCoupling.persistentStateInitialized(),
        "failed coupling check published checkpoint state");

    CondensedCarbonOptions changed=options;
    changed.limits.maximumMassLoading=0.19;
    SectorCondensedCarbonOperator wrongFingerprint(mesh,registry,changed);
    test.requireThrows([&]{
        (void)CondensedCarbonCheckpoint::read(
            checkpoint.string(),*domain,wrongFingerprint,
            metadata.coupledFluidCheckpointId);
    },"checkpoint accepted a stale model fingerprint");
    test.require(!wrongFingerprint.persistentStateInitialized(),
        "failed fingerprint check published checkpoint state");

    CondensedCarbonOptions twoMoment=options;
    twoMoment.momentModel=CondensedCarbonMomentModel::TWO_MOMENT;
    twoMoment.initialState={0.0,0.0};
    twoMoment.inletState={0.0,0.0};
    SectorCondensedCarbonOperator wrongMode(mesh,registry,twoMoment);
    test.requireThrows([&]{
        (void)CondensedCarbonCheckpoint::read(
            checkpoint.string(),*domain,wrongMode,
            metadata.coupledFluidCheckpointId);
    },"checkpoint accepted a one/two-moment restart mismatch");
    test.require(!wrongMode.persistentStateInitialized(),
        "failed moment-mode check published checkpoint state");

    const std::string changedProfile=temporary.writeProfile(
        "changed_profile.json",65.0,70.0);
    SectorMeshMetrics changedMesh(
        changedProfile,mesh.Nx(),mesh.Nr(),mesh.Ntheta(),120.0);
    std::unique_ptr<Domain> changedDomain=makeSectorDomain(
        changedMesh,registry);
    SectorCondensedCarbonOperator wrongGeometry(
        changedMesh,registry,options);
    test.require(source.meshFingerprint()!=wrongGeometry.meshFingerprint(),
        "different sector geometry produced the same mesh fingerprint");
    test.requireThrows([&]{
        (void)CondensedCarbonCheckpoint::read(
            checkpoint.string(),*changedDomain,wrongGeometry,
            metadata.coupledFluidCheckpointId);
    },"checkpoint accepted a same-dimension different sector geometry");
    test.require(!wrongGeometry.persistentStateInitialized(),
        "failed geometry check published checkpoint state");
    test.requireThrows([&]{
        CondensedCarbonOutput::writeSectorFieldCsv(
            temporary.path("wrong_geometry.csv").string(),
            *domain,changedMesh,restored);
    },"condensed-carbon output accepted a mismatched sector geometry");
    test.requireThrows([&]{
        CondensedCarbonOutput::writeSectorVtk(
            temporary.path("wrong_geometry.vtk").string(),
            *domain,changedMesh,restored);
    },"condensed-carbon VTK accepted a mismatched sector geometry");
    test.requireThrows([&]{
        CondensedCarbonOutput::writeSectorVtk(
            temporary.path("wrong_carbon_geometry.vtk").string(),
            *domain,mesh,wrongGeometry);
    },"condensed-carbon VTK accepted fields owned by another sector geometry");

    const std::filesystem::path fieldCsv=temporary.path("carbon_field.csv");
    const std::filesystem::path wallCsv=temporary.path("carbon_wall.csv");
    CondensedCarbonOutput::writeSectorFieldCsv(
        fieldCsv.string(),*domain,mesh,restored);
    CondensedCarbonOutput::writeWallInventoryCsv(
        wallCsv.string(),*domain,mesh,restored);
    const std::string fieldText=readText(fieldCsv);
    const std::string wallText=readText(wallCsv);
    test.require(fieldText.find("rho_soot_kg_m3")!=std::string::npos,
        "condensed-carbon field CSV lacks its unit-bearing header");
    test.require(wallText.find("wall_carbon_kg_m2")!=std::string::npos,
        "wall carbon CSV lacks its unit-bearing header");

    const std::filesystem::path vtk=temporary.path("carbon.vtk");
    CondensedCarbonOutput::writeSectorVtk(
        vtk.string(),*domain,mesh,restored);
    const std::string vtkText=readText(vtk);
    test.require(vtkText.find("DATASET STRUCTURED_GRID")!=std::string::npos,
        "condensed-carbon VTK is not a sector structured grid");
    test.require(vtkText.find("RECTILINEAR_GRID")==std::string::npos,
        "condensed-carbon VTK silently used rectilinear geometry");
    test.require(vtkText.find("SCALARS rho_soot double")!=std::string::npos,
        "VTK output lacks rho_soot");
    test.require(vtkText.find("SCALARS N_soot double")==std::string::npos,
        "one-moment VTK output unexpectedly contains N_soot");
    const std::vector<FaceAreaVector> vtkPoints=readVtkPoints(vtkText);
    const int pointNx=domain->Nx()-1;
    const int pointNr=domain->Ny()-1;
    const int pointI=2,pointJ=2,pointK=2;
    const std::size_t pointIndex=static_cast<std::size_t>(pointK-1)
        *pointNr*pointNx+static_cast<std::size_t>(pointJ-1)*pointNx
        +static_cast<std::size_t>(pointI-1);
    const double pointRadius=mesh.wallFace(pointI)*mesh.etaFace(pointJ);
    const FaceAreaVector expectedPoint={
        mesh.xFace(pointI),
        pointRadius*std::cos(mesh.thetaFace(pointK)),
        pointRadius*std::sin(mesh.thetaFace(pointK))};
    test.require(pointIndex<vtkPoints.size()
            &&near(vtkPoints[pointIndex].x,expectedPoint.x)
            &&near(vtkPoints[pointIndex].y,expectedPoint.y)
            &&near(vtkPoints[pointIndex].z,expectedPoint.z),
        "condensed-carbon VTK points do not follow sector geometry");

    const std::filesystem::path twoMomentVtk=temporary.path("carbon_two_moment.vtk");
    CondensedCarbonOutput::writeSectorVtk(
        twoMomentVtk.string(),*domain,mesh,wrongMode);
    test.require(readText(twoMomentVtk).find("SCALARS N_soot double")
            !=std::string::npos,
        "two-moment VTK output lacks N_soot");
}

} // namespace

int main()
{
    try
    {
        TestSuite test;
        TemporaryWorkspace temporary;
        SpeciesRegistry registry;
        const std::string profile=temporary.writeStraightProfile();
        SectorMeshMetrics mesh(profile,6,6,6,120.0);

        testExchangeContract(test,registry);
        testConfigurationAndApplicability(test);
        testUniformAdvection(test,mesh,registry);
        testUnownedBoundaryFailClose(test,mesh,registry);
        testGalileanEnergyAudit(test,mesh,registry);
        testLifecycle(test,mesh,registry);
        testWallTransfer(test,mesh,registry);
        testCheckpointAndOutput(test,temporary,mesh,registry);

        std::cout<<"CONDENSED_CARBON_V1_PASS checks="<<test.checks()<<'\n';
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<"CONDENSED_CARBON_V1_FAIL: "<<error.what()<<'\n';
        return 1;
    }
}
