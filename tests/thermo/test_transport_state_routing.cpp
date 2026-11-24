#include "DiffusionOperator.H"
#include "SectorDiffusionOperator.H"
#include "SectorMeshMetrics.H"
#include "SectorTurbulenceOperator.H"
#include "SectorWallHeatFlux.H"
#include "TransportProps.H"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Foam;

namespace
{

int failures=0;

void check(bool condition,const std::string& label)
{
    std::cout<<"  ["<<(condition?"PASS":"FAIL")<<"] "<<label<<'\n';
    if(!condition) ++failures;
}

bool same(const ConservedVars& lhs,const ConservedVars& rhs)
{
    return lhs.rho==rhs.rho
        &&lhs.rhoU==rhs.rhoU
        &&lhs.rhoE==rhs.rhoE
        &&lhs.rhoY==rhs.rhoY
        &&lhs.rhoZ==rhs.rhoZ
        &&lhs.rhoC==rhs.rhoC;
}

PrimitiveVars primitive(
    double pressure,double temperature,
    const std::vector<double>& massFractions={1.0})
{
    PrimitiveVars result(static_cast<int>(massFractions.size()));
    result.p=pressure;
    result.T=temperature;
    result.U={0.0,0.0,0.0};
    result.Y=massFractions;
    result.Z=0.0;
    result.C=0.0;
    return result;
}

class SelectiveHemModel final : public thermo::EquilibriumThermoModel
{
    std::vector<std::string> species_={"O2"};
    thermo::LegacySinglePhaseThermoModel delegate_{species_};
    double threshold_;
    std::string id_=delegate_.modelId();

    thermo::EquilibriumResult route(
        thermo::EquilibriumResult result) const
    {
        result.diagnostics.modelId=id_;
        if(!result.success() || result.state->T<threshold_) return result;
        result.state->acoustic.closure=
            thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM;
        result.state->thermodynamics.reset();
        if(!result.state->valid())
        {
            result.status=thermo::EquilibriumStatus::ACOUSTIC_CLOSURE_FAILED;
            result.state.reset();
        }
        return result;
    }

public:
    explicit SelectiveHemModel(double threshold) : threshold_(threshold) {}

    thermo::EquilibriumResult recoverRhoE(
        double density,double internalEnergy,
        const std::vector<double>& massFractions,
        const thermo::EquilibriumSeed* seed) const override
    {
        return route(delegate_.recoverRhoE(
            density,internalEnergy,massFractions,seed));
    }

    thermo::EquilibriumResult initializePTY(
        double pressure,double temperature,
        const std::vector<double>& massFractions) const override
    {
        return route(delegate_.initializePTY(
            pressure,temperature,massFractions));
    }

    const std::string& modelId() const override { return id_; }
    const std::vector<std::string>& speciesNames() const override
    { return species_; }
    EnergyConvention energyConvention() const override
    {
        return EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT;
    }
    bool stateValidForModel(
        const thermo::EquilibriumState& state) const override
    {
        if(state.acoustic.closure
            ==thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM)
        {
            return thermo::EquilibriumThermoModel::stateValidForModel(state);
        }
        return delegate_.stateValidForModel(state);
    }
    bool supportsPrimitivePTYInitialization() const override { return true; }
    bool supportsCommittedLegacyGasTransport() const override { return true; }
};

thermo::EquilibriumState twoPhaseOxygen(
    const SpeciesRegistry& registry)
{
    constexpr double pressure=1.0e6;
    constexpr double temperature=300.0;
    constexpr double vaporMoleFraction=0.5;
    const double molecularWeight=
        registry.molecularWeight(registry.index("O2"));

    thermo::EquilibriumPhaseProperties liquid;
    liquid.massFractions={1.0};
    liquid.moleFractions={1.0};
    liquid.molecularWeight=molecularWeight;
    liquid.rho=1000.0;
    liquid.h=100000.0;
    liquid.e=liquid.h-pressure/liquid.rho;
    liquid.compressibility=pressure*molecularWeight
        /(liquid.rho*8.314462618*temperature);

    thermo::EquilibriumPhaseProperties vapor;
    vapor.massFractions={1.0};
    vapor.moleFractions={1.0};
    vapor.molecularWeight=molecularWeight;
    vapor.rho=10.0;
    vapor.h=500000.0;
    vapor.e=vapor.h-pressure/vapor.rho;
    vapor.compressibility=pressure*molecularWeight
        /(vapor.rho*8.314462618*temperature);

    const double vaporMassFraction=vaporMoleFraction;
    const double liquidMassFraction=1.0-vaporMassFraction;
    const double specificVolume=liquidMassFraction/liquid.rho
        +vaporMassFraction/vapor.rho;

    thermo::EquilibriumState result;
    result.region=thermo::EquilibriumRegion::TWO_PHASE;
    result.overallMassFractions={1.0};
    result.p=pressure;
    result.T=temperature;
    result.rho=1.0/specificVolume;
    result.h=liquidMassFraction*liquid.h+vaporMassFraction*vapor.h;
    result.e=liquidMassFraction*liquid.e+vaporMassFraction*vapor.e;
    result.vaporMoleFraction=vaporMoleFraction;
    result.vaporMassFraction=vaporMassFraction;
    result.vaporVolumeFraction=(vaporMassFraction/vapor.rho)/specificVolume;
    result.acoustic.closure=
        thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM;
    result.acoustic.pressureDensityDerivative=1.0e5;
    result.liquid=liquid;
    result.vapor=vapor;
    return result;
}

bool initialize1D(Domain& domain,const std::vector<double>& temperatures)
{
    for(int cell=0;cell<domain.N();++cell)
        domain.stagedW(cell)=primitive(4.0e6,temperatures[cell]);
    return domain.primitiveToConserved();
}

bool initializeSector(
    Domain& domain,int badI=-1,int badJ=-1,int badK=-1)
{
    for(int i=0;i<domain.Nx();++i)
    for(int j=0;j<domain.Ny();++j)
    for(int k=0;k<domain.Nz();++k)
    {
        const double temperature=(i==badI&&j==badJ&&k==badK)
            ?600.0:300.0;
        domain.stagedW(i,j,k)=primitive(4.0e6,temperature);
    }
    return domain.primitiveToConserved();
}

bool initializeHotOutletGhost(Domain& domain)
{
    for(int i=0;i<domain.Nx();++i)
    for(int j=0;j<domain.Ny();++j)
    for(int k=0;k<domain.Nz();++k)
        domain.stagedW(i,j,k)=primitive(1.0e5,300.0);
    domain.stagedW(domain.Nx()-1,1,1)=primitive(1.0e5,4000.0);
    return domain.primitiveToConserved();
}

} // namespace

int main()
{
    std::cout<<"=== Committed-state transport routing contract ===\n";
    SpeciesRegistry registry;
    const std::vector<std::string> names={"O2"};
    const std::vector<double> oxygen={1.0};
    const auto legacy=thermo::makeLegacySinglePhaseThermoModel(names);

    const auto base=legacy->initializePTY(4.0e6,300.0,oxygen);
    check(base.success(),"legacy O2 committed state initializes");
    if(!base.success()) return 1;
    check(base.state->thermodynamics.has_value()
          &&base.state->thermodynamics->valid(),
        "legacy Q publishes a versioned thermodynamic payload");
    const PrimitiveVars basePrimitive=primitive(4.0e6,300.0);
    const TransportProps::Result accepted=
        TransportProps::evaluateCommittedLegacyGas(
            names,basePrimitive,*base.state,*legacy,registry);
    check(accepted.success(),"legacy frozen single-phase gas is accepted");
    check(accepted.provenance
              ==TransportProps::Provenance::COMMITTED_LEGACY_Q,
        "successful transport records the committed-Q route");
    check(std::isfinite(accepted.point.mu)&&accepted.point.mu>0.0
          &&std::isfinite(accepted.point.kappa)&&accepted.point.kappa>0.0
          &&std::isfinite(accepted.point.Cp)&&accepted.point.Cp>accepted.point.Cv
          &&std::isfinite(accepted.point.Cv)&&accepted.point.Cv>0.0,
        "accepted transport properties are finite and positive");
    check(accepted.point.Cp==base.state->thermodynamics->Cp
          &&accepted.point.Cv==base.state->thermodynamics->Cv,
        "committed transport consumes Q Cp/Cv without a second caloric state");
    thermo::EquilibriumState modifiedPayload=*base.state;
    const double originalCp=modifiedPayload.thermodynamics->Cp;
    const double originalCv=modifiedPayload.thermodynamics->Cv;
    modifiedPayload.thermodynamics->Cv=1.05*originalCv;
    modifiedPayload.thermodynamics->dedT_rho=
        modifiedPayload.thermodynamics->Cv;
    modifiedPayload.thermodynamics->Cp=
        modifiedPayload.thermodynamics->Cv+(originalCp-originalCv);
    modifiedPayload.acoustic.pressureDensityDerivative=
        (modifiedPayload.thermodynamics->Cp
            /modifiedPayload.thermodynamics->Cv)
        *modifiedPayload.thermodynamics->dpdrho_T;
    check(modifiedPayload.valid()
          &&modifiedPayload.thermodynamics->Cp!=originalCp
          &&modifiedPayload.thermodynamics->Cv!=originalCv,
        "synthetic direct-Q payload remains algebraically closed");
    const TransportProps::Result modifiedTransport=
        TransportProps::evaluateCommittedLegacyGas(
            names,basePrimitive,modifiedPayload,*legacy,registry);
    check(modifiedTransport.success()
          &&modifiedTransport.point.Cp==modifiedPayload.thermodynamics->Cp
          &&modifiedTransport.point.Cv==modifiedPayload.thermodynamics->Cv,
        "transport uses committed Cp/Cv rather than a PTY reflash");

    thermo::EquilibriumState forgedSpeciesBasis=*base.state;
    forgedSpeciesBasis.singlePhase->molecularWeight*=1.1;
    forgedSpeciesBasis.singlePhase->compressibility=
        forgedSpeciesBasis.p*forgedSpeciesBasis.singlePhase->molecularWeight
        /(forgedSpeciesBasis.singlePhase->rho*8.314462618
            *forgedSpeciesBasis.T);
    check(forgedSpeciesBasis.valid()
          &&!legacy->stateValidForModel(forgedSpeciesBasis),
        "structural Q validation does not forge a model species basis");
    check(TransportProps::evaluateCommittedLegacyGas(
            names,basePrimitive,forgedSpeciesBasis,*legacy,registry).status
            ==TransportProps::Status::INVALID_EQUILIBRIUM_STATE,
        "transport binds committed Q to the supplied thermo model");
    TransportProps::Result forgedSuccess=accepted;
    forgedSuccess.point.Cv=std::numeric_limits<double>::quiet_NaN();
    check(!forgedSuccess.success(),
        "SUCCESS cannot carry an invalid transport payload");
    TransportProps::Result missingProvenance=accepted;
    missingProvenance.provenance=TransportProps::Provenance::NONE;
    check(!missingProvenance.success(),
        "SUCCESS cannot omit thermodynamic provenance");

    thermo::EquilibriumState missingPayload=*base.state;
    missingPayload.thermodynamics.reset();
    check(missingPayload.valid()
          &&TransportProps::evaluateCommittedLegacyGas(
              names,basePrimitive,missingPayload,*legacy,registry).status
              ==TransportProps::Status::MISSING_THERMODYNAMIC_PAYLOAD,
        "legacy transport rejects a committed Q without caloric derivatives");
    thermo::EquilibriumState hem=*base.state;
    hem.acoustic.closure=thermo::AcousticClosure::HOMOGENEOUS_EQUILIBRIUM;
    hem.thermodynamics.reset();
    check(hem.valid(),"HEM-routed single-phase test state remains structurally valid");
    check(TransportProps::evaluateCommittedLegacyGas(
            names,basePrimitive,hem,*legacy,registry).status
            ==TransportProps::Status::UNSUPPORTED_ACOUSTIC_CLOSURE,
        "HEM acoustics cannot be reinterpreted as legacy gas transport");

    const thermo::EquilibriumState twoPhase=twoPhaseOxygen(registry);
    check(twoPhase.valid()&&twoPhase.rho<TransportProps::gasClosureMaxDensity,
        "two-phase routing state is valid and below the density gate");
    check(TransportProps::evaluateCommittedLegacyGas(
            names,primitive(twoPhase.p,twoPhase.T),twoPhase,
            *legacy,registry).status
            ==TransportProps::Status::UNSUPPORTED_PHASE,
        "two-phase state is rejected by phase before gas correlations");

    const auto cold=legacy->initializePTY(4.0e6,130.0,oxygen);
    check(cold.success(),"cold O2 thermo state initializes for routing audit");
    if(cold.success())
        check(TransportProps::evaluateCommittedLegacyGas(
                names,primitive(4.0e6,130.0),*cold.state,
                *legacy,registry).status
                ==TransportProps::Status::UNSUPPORTED_TEMPERATURE_RANGE,
            "NASA-7 low-temperature extrapolation is transport-gated");
    const auto belowSutherlandFloor=
        legacy->initializePTY(1.0e5,90.0,oxygen);
    check(belowSutherlandFloor.success(),
        "sub-Sutherland-floor O2 thermo state initializes for gate ordering");
    if(belowSutherlandFloor.success())
        check(TransportProps::evaluateCommittedLegacyGas(
                names,primitive(1.0e5,90.0),*belowSutherlandFloor.state,
                *legacy,registry).status
                ==TransportProps::Status::UNSUPPORTED_TEMPERATURE_RANGE,
            "temperature capability gate precedes gas correlations");

    const std::vector<std::string> oxygenNitrogenNames={"O2","N2"};
    const auto oxygenNitrogen=
        thermo::makeLegacySinglePhaseThermoModel(oxygenNitrogenNames);
    const std::vector<double> oxygenOnly={1.0,0.0};
    const std::vector<double> nitrogenOnly={0.0,1.0};
    const auto oxygenAt250=
        oxygenNitrogen->initializePTY(1.0e5,250.0,oxygenOnly);
    check(oxygenAt250.success(),"250 K O2/N2-layout O2 state initializes");
    if(oxygenAt250.success())
        check(TransportProps::evaluateCommittedLegacyGas(
                oxygenNitrogenNames,primitive(1.0e5,250.0,oxygenOnly),
                *oxygenAt250.state,*oxygenNitrogen,registry).success(),
            "inactive N2 does not raise the active-species temperature floor");
    const auto nitrogenAt250=
        oxygenNitrogen->initializePTY(1.0e5,250.0,nitrogenOnly);
    check(nitrogenAt250.success(),"250 K N2 thermo state initializes for routing audit");
    if(nitrogenAt250.success())
        check(TransportProps::evaluateCommittedLegacyGas(
                oxygenNitrogenNames,primitive(1.0e5,250.0,nitrogenOnly),
                *nitrogenAt250.state,*oxygenNitrogen,registry).status
                ==TransportProps::Status::UNSUPPORTED_TEMPERATURE_RANGE,
            "active N2 enforces its 300 K NASA-7 lower bound");
    const std::vector<std::string> propaneNames={"C3H8"};
    const std::vector<double> propaneOnly={1.0};
    const auto propaneModel=
        thermo::makeLegacySinglePhaseThermoModel(propaneNames);
    const auto propaneBelowNasaRange=
        propaneModel->initializePTY(1.0e5,299.8,propaneOnly);
    check(propaneBelowNasaRange.success(),
        "sub-300 K C3H8 state initializes for transport routing audit");
    if(propaneBelowNasaRange.success())
    {
        check(TransportProps::evaluateCommittedLegacyGas(
                propaneNames,primitive(1.0e5,299.8,propaneOnly),
                *propaneBelowNasaRange.state,*propaneModel,registry).status
                ==TransportProps::Status::UNSUPPORTED_TEMPERATURE_RANGE,
            "production C3H8 transport remains fail-closed below 300 K");
        const TransportProps::DiagnosticGasResult propaneDiagnostic=
            TransportProps::evaluateDiagnosticGas(
                propaneNames,primitive(1.0e5,299.8,propaneOnly),
                *propaneBelowNasaRange.state,*propaneModel,registry);
        check(propaneDiagnostic.result.success()
              &&propaneDiagnostic.extrapolation
                  ==TransportProps::DiagnosticExtrapolation::CALORIC_RANGE,
            "diagnostic C3H8 transport explicitly audits the range extrapolation");
    }

    const auto oxygenAt4000=legacy->initializePTY(1.0e5,4000.0,oxygen);
    check(oxygenAt4000.success(),
        "high-temperature O2 thermo state initializes for routing audit");
    if(oxygenAt4000.success())
        check(TransportProps::evaluateCommittedLegacyGas(
                names,primitive(1.0e5,4000.0),*oxygenAt4000.state,
                *legacy,registry).status
                ==TransportProps::Status::UNSUPPORTED_TEMPERATURE_RANGE,
            "active O2 enforces its NASA-7 upper bound beyond slack");

    const auto dense=legacy->initializePTY(2.0e7,200.0,oxygen);
    check(dense.success(),"dense O2 thermo state initializes for routing audit");
    if(dense.success())
    {
        check(dense.state->rho>TransportProps::gasClosureMaxDensity,
            "dense O2 state exceeds the legacy gas density limit");
        check(TransportProps::evaluateCommittedLegacyGas(
                names,primitive(2.0e7,200.0),*dense.state,
                *legacy,registry).status
                ==TransportProps::Status::UNSUPPORTED_DENSITY_RANGE,
            "dense state is rejected before gas transport evaluation");
    }

    PrimitiveVars mismatched=basePrimitive;
    mismatched.T+=1.0;
    check(TransportProps::evaluateCommittedLegacyGas(
            names,mismatched,*base.state,*legacy,registry).status
            ==TransportProps::Status::PRIMITIVE_EQUILIBRIUM_MISMATCH,
        "primitive and committed equilibrium states must match");
    const auto nitrogenModel=thermo::makeLegacySinglePhaseThermoModel({"N2"});
    check(TransportProps::evaluateCommittedLegacyGas(
            names,basePrimitive,*base.state,*nitrogenModel,registry).status
            ==TransportProps::Status::UNSUPPORTED_THERMO_MODEL,
        "transport rejects a different thermo-model provenance");

    const auto selective=std::make_shared<const SelectiveHemModel>(550.0);
    Domain diffusionDomain(5,0.05,1,selective);
    check(initialize1D(diffusionDomain,{300.0,350.0,500.0,600.0,700.0}),
        "selective-HEM diffusion Domain initializes");
    diffusionDomain.beginStep();
    std::vector<ConservedVars> diffusionBefore;
    for(int cell=0;cell<diffusionDomain.N();++cell)
        diffusionBefore.push_back(diffusionDomain.Unew(cell));
    DiffusionOperator diffusion(registry);
    bool diffusionRejected=false;
    std::string diffusionError;
    try
    {
        diffusion.apply(diffusionDomain,1.0e-7);
    }
    catch(const std::runtime_error& error)
    {
        diffusionRejected=true;
        diffusionError=error.what();
    }
    bool diffusionUnchanged=true;
    for(int cell=0;cell<diffusionDomain.N();++cell)
        diffusionUnchanged=diffusionUnchanged
            &&same(diffusionDomain.Unew(cell),diffusionBefore[cell]);
    check(diffusionRejected
          &&diffusionError.find("cell 3")!=std::string::npos
          &&diffusionError.find("UNSUPPORTED_ACOUSTIC_CLOSURE")
              !=std::string::npos,
        "diffusion reports the first unsupported committed transport cell");
    check(diffusionUnchanged,
        "diffusion transport preflight fails before any Unew mutation");

    constexpr int nx=4,nr=4,nt=4;
    SectorMeshMetrics mesh(vortexProfilePath(),nx,nr,nt,120.0);
    SectorTurbulenceOperator turbulence(mesh,registry);
    check(!turbulence.persistentStateInitialized(),
        "sector SST starts without materialized persistent fields");

    Domain diagnosticPropane(
        nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),1,propaneModel);
    for(int i=0;i<diagnosticPropane.Nx();++i)
    for(int j=0;j<diagnosticPropane.Ny();++j)
    for(int k=0;k<diagnosticPropane.Nz();++k)
        diagnosticPropane.stagedW(i,j,k)=
            primitive(1.0e5,299.8,propaneOnly);
    check(diagnosticPropane.primitiveToConserved(),
        "sub-300 K C3H8 diagnostic sector initializes");
    diagnosticPropane.beginStep();
    SectorTurbulenceOperator diagnosticPropaneTurbulence(
        mesh,registry,0.05,0.001,0.25,true);
    diagnosticPropaneTurbulence.beginStepTransaction(diagnosticPropane);
    bool diagnosticPropaneRejected=false;
    try
    {
        diagnosticPropaneTurbulence.apply(diagnosticPropane,1.0e-7);
    }
    catch(const std::runtime_error&)
    {
        diagnosticPropaneRejected=true;
    }
    check(!diagnosticPropaneRejected,
        "diagnostic SST retries a committed NASA-range reject through audited extrapolation");
    diagnosticPropaneTurbulence.rollbackStepTransaction();

    Domain validSector(
        nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),1,selective);
    check(initializeSector(validSector),"valid sector transport field initializes");
    const SectorWallHeatFluxResult validWall=SectorWallHeatFlux::evaluate(
        validSector,mesh,registry,&turbulence,2,2);
    check(std::isfinite(validWall.totalFlux)
          &&!turbulence.persistentStateInitialized(),
        "wall transport reads initial SST state without lazy initialization");
    bool wallIndexRejected=false;
    std::string wallIndexError;
    try
    {
        (void)SectorWallHeatFlux::evaluate(
            validSector,mesh,registry,&turbulence,-1,2);
    }
    catch(const std::runtime_error& error)
    {
        wallIndexRejected=true;
        wallIndexError=error.what();
    }
    check(wallIndexRejected
          &&wallIndexError.find("WALL_INDEX_OUT_OF_RANGE")
              !=std::string::npos
          &&!turbulence.persistentStateInitialized(),
        "wall transport rejects an invalid index before state access");

    Domain insufficientDomain(
        3,3,3,mesh.length(),0.065,mesh.sectorAngle(),1,selective);
    bool insufficientHaloRejected=false;
    std::string insufficientHaloError;
    try
    {
        (void)SectorWallHeatFlux::evaluate(
            insufficientDomain,mesh,registry,nullptr,1,1);
    }
    catch(const std::runtime_error& error)
    {
        insufficientHaloRejected=true;
        insufficientHaloError=error.what();
    }
    check(insufficientHaloRejected
          &&insufficientHaloError.find("INSUFFICIENT_HALO_DIMENSIONS")
              !=std::string::npos,
        "wall transport rejects insufficient halo dimensions before state access");
    SectorDiffusionOperator validSectorDiffusion(mesh,registry,&turbulence);
    const double sectorDt=validSectorDiffusion.maxStableDt(validSector,1.0);
    check(std::isfinite(sectorDt)&&sectorDt>0.0
          &&!turbulence.persistentStateInitialized(),
        "sector diffusion CFL preflights nut without OMP lazy state writes");

    SectorMeshMetrics mismatchedMesh(
        vortexProfilePath(),nx+1,nr,nt,120.0);
    SectorTurbulenceOperator mismatchedTurbulence(
        mismatchedMesh,registry);
    validSector.beginStep();
    std::vector<ConservedVars> mismatchBefore;
    for(int cell=0;cell<validSector.N();++cell)
        mismatchBefore.push_back(validSector.Unew(cell));
    SectorDiffusionOperator mismatchedSectorDiffusion(
        mesh,registry,&mismatchedTurbulence);
    bool diffusionMeshRejected=false;
    std::string diffusionMeshError;
    try
    {
        mismatchedSectorDiffusion.apply(validSector,1.0e-7);
    }
    catch(const std::runtime_error& error)
    {
        diffusionMeshRejected=true;
        diffusionMeshError=error.what();
    }
    bool mismatchUnchanged=true;
    for(int cell=0;cell<validSector.N();++cell)
        mismatchUnchanged=mismatchUnchanged
            &&same(validSector.Unew(cell),mismatchBefore[cell]);
    check(diffusionMeshRejected
          &&diffusionMeshError.find("TURBULENCE_DOMAIN_MESH_MISMATCH")
              !=std::string::npos
          &&mismatchUnchanged
          &&!mismatchedTurbulence.persistentStateInitialized(),
        "sector diffusion rejects mismatched SST mesh before state mutation");

    mismatchedTurbulence.beginStepTransaction(validSector);
    bool turbulenceMeshRejected=false;
    std::string turbulenceMeshError;
    try
    {
        mismatchedTurbulence.apply(validSector,1.0e-7);
    }
    catch(const std::runtime_error& error)
    {
        turbulenceMeshRejected=true;
        turbulenceMeshError=error.what();
    }
    bool turbulenceMismatchUnchanged=true;
    for(int cell=0;cell<validSector.N();++cell)
        turbulenceMismatchUnchanged=turbulenceMismatchUnchanged
            &&same(validSector.Unew(cell),mismatchBefore[cell]);
    check(turbulenceMeshRejected
          &&turbulenceMeshError.find("DOMAIN_MESH_MISMATCH")
              !=std::string::npos
          &&turbulenceMismatchUnchanged
          &&!mismatchedTurbulence.persistentStateInitialized(),
        "sector SST rejects Domain/mesh mismatch before Unew and first-step state");
    mismatchedTurbulence.rollbackStepTransaction();

    Domain sectorDomain(
        nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),1,selective);
    check(initializeSector(sectorDomain,2,nr-1,2),
        "wall-ghost HEM sector Domain initializes");
    sectorDomain.beginStep();
    std::vector<ConservedVars> sectorBefore;
    for(int cell=0;cell<sectorDomain.N();++cell)
        sectorBefore.push_back(sectorDomain.Unew(cell));

    bool wallRejected=false;
    std::string wallError;
    try
    {
        (void)SectorWallHeatFlux::evaluate(
            sectorDomain,mesh,registry,&turbulence,2,2);
    }
    catch(const std::runtime_error& error)
    {
        wallRejected=true;
        wallError=error.what();
    }
    check(wallRejected
          &&wallError.find("cell 46")!=std::string::npos
          &&wallError.find("UNSUPPORTED_ACOUSTIC_CLOSURE")
              !=std::string::npos
          &&!turbulence.persistentStateInitialized(),
        "wall heat flux rejects committed wall Q before SST state mutation");

    SectorDiffusionOperator sectorDiffusion(mesh,registry,&turbulence);
    bool sectorDiffusionRejected=false;
    std::string sectorDiffusionError;
    try
    {
        sectorDiffusion.apply(sectorDomain,1.0e-7);
    }
    catch(const std::runtime_error& error)
    {
        sectorDiffusionRejected=true;
        sectorDiffusionError=error.what();
    }
    bool sectorDiffusionUnchanged=true;
    for(int cell=0;cell<sectorDomain.N();++cell)
        sectorDiffusionUnchanged=sectorDiffusionUnchanged
            &&same(sectorDomain.Unew(cell),sectorBefore[cell]);
    check(!sectorDiffusionRejected&&sectorDiffusionError.empty(),
        "sector diffusion ignores unsupported wall-ghost transport state");
    check(!sectorDiffusionUnchanged&&!turbulence.persistentStateInitialized(),
        "sector diffusion applies cached wall flux without SST initialization");

    turbulence.beginStepTransaction(sectorDomain);
    bool turbulenceRejected=false;
    std::string turbulenceError;
    try
    {
        turbulence.apply(sectorDomain,1.0e-7);
    }
    catch(const std::runtime_error& error)
    {
        turbulenceRejected=true;
        turbulenceError=error.what();
    }
    bool sectorUnchanged=true;
    for(int cell=0;cell<sectorDomain.N();++cell)
        sectorUnchanged=sectorUnchanged
            &&same(sectorDomain.Unew(cell),sectorBefore[cell]);
    check(!turbulenceRejected&&turbulenceError.empty(),
        "sector SST ignores unsupported wall-ghost transport state");
    check(!sectorUnchanged&&turbulence.persistentStateInitialized(),
        "sector SST advances after active-cell transport preflight");
    turbulence.rollbackStepTransaction();

    Domain hotOutlet(nx,nr,nt,mesh.length(),0.065,mesh.sectorAngle(),1);
    hotOutlet.setSpecies(names);
    check(initializeHotOutletGhost(hotOutlet),
        "hot outlet ghost state initializes for halo transport regression");
    hotOutlet.beginStep();
    SectorTurbulenceOperator hotOutletTurbulence(mesh,registry);
    SectorDiffusionOperator hotOutletDiffusion(
        mesh,registry,&hotOutletTurbulence);
    bool hotOutletRejected=false;
    try
    {
        hotOutletDiffusion.apply(hotOutlet,1.0e-9);
        hotOutletTurbulence.apply(hotOutlet,1.0e-9);
    }
    catch(const std::exception&)
    {
        hotOutletRejected=true;
    }
    check(!hotOutletRejected,
        "4000 K outlet ghost does not abort active transport preflight");

    if(failures!=0)
    {
        std::cerr<<"TRANSPORT_STATE_ROUTING_FAIL failures="<<failures<<'\n';
        return 1;
    }
    std::cout<<"TRANSPORT_STATE_ROUTING_PASS\n";
    return 0;
}
