#include "CondensedCarbonOperator.H"

#include "ThermoDataProvenance.H"
#include "TransportProps.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Foam
{
namespace
{

constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double boundaryFluxAbsoluteTolerance = 1.0e-14;
constexpr double boundaryFluxRelativeTolerance = 1.0e-12;

double magnitude(const FaceAreaVector& area)
{
    return std::sqrt(area.x*area.x+area.y*area.y+area.z*area.z);
}

double dot(const std::vector<double>& velocity,const FaceAreaVector& area)
{
    return velocity[0]*area.x+velocity[1]*area.y+velocity[2]*area.z;
}

int applicabilitySeverity(CondensedCarbonApplicabilityStatus status)
{
    switch(status)
    {
        case CondensedCarbonApplicabilityStatus::SUPPORTED: return 0;
        case CondensedCarbonApplicabilityStatus::WARNING_RESEARCH_ONLY: return 1;
        case CondensedCarbonApplicabilityStatus::SOOT_ENERGY_MODEL_REQUIRED: return 2;
        case CondensedCarbonApplicabilityStatus::SLIP_MODEL_REQUIRED: return 3;
        case CondensedCarbonApplicabilityStatus::OUTSIDE_DILUTE_LIMIT: return 4;
        case CondensedCarbonApplicabilityStatus::FAILED: return 5;
    }
    return 5;
}

bool finiteNonnegative(double value)
{
    return std::isfinite(value)&&value>=0.0;
}

bool closedBoundaryFlux(double residual,double scale)
{
    return std::isfinite(residual)&&std::isfinite(scale)
        &&std::abs(residual)<=boundaryFluxAbsoluteTolerance
            +boundaryFluxRelativeTolerance*scale;
}

} // namespace

std::string condensedCarbonSectorMeshFingerprint(
    const SectorMeshMetrics& mesh)
{
    std::string manifest="schema=rfam.sector-mesh/v1;";
    thermoProvenance::appendSigned(manifest,"Nx",mesh.Nx());
    thermoProvenance::appendSigned(manifest,"Nr",mesh.Nr());
    thermoProvenance::appendSigned(manifest,"Ntheta",mesh.Ntheta());
    thermoProvenance::appendDouble(manifest,"length",mesh.length());
    thermoProvenance::appendDouble(
        manifest,"sector-angle",mesh.sectorAngle());
    for(int i=0;i<=mesh.Nx();++i)
    {
        const std::string prefix="axial-face["+std::to_string(i)+"].";
        thermoProvenance::appendDouble(
            manifest,prefix+"x",mesh.xFace(i));
        thermoProvenance::appendDouble(
            manifest,prefix+"wall-radius",mesh.wallFace(i));
    }
    for(int j=0;j<=mesh.Nr();++j)
        thermoProvenance::appendDouble(
            manifest,"eta-face["+std::to_string(j)+"]",mesh.etaFace(j));
    for(int k=0;k<=mesh.Ntheta();++k)
        thermoProvenance::appendDouble(
            manifest,"theta-face["+std::to_string(k)+"]",mesh.thetaFace(k));
    return thermoProvenance::sha256Hex(manifest);
}

int SectorCondensedCarbonOperator::index(int i,int j,int k) const noexcept
{
    return k*mesh_.Nx()*mesh_.Nr()+j*mesh_.Nx()+i;
}

int SectorCondensedCarbonOperator::wallIndex(int i,int k) const noexcept
{
    return k*mesh_.Nx()+i;
}

SectorCondensedCarbonOperator::SectorCondensedCarbonOperator(
    const SectorMeshMetrics& mesh,
    const SpeciesRegistry& registry,
    CondensedCarbonOptions options,
    const CondensedCarbonKinetics* kinetics,
    const WallCarbonDepositionModel* wallModel)
:
    mesh_(mesh),
    registry_(registry),
    options_(std::move(options)),
    kinetics_(kinetics),
    wallModel_(wallModel),
    kineticsModelId_(kinetics?kinetics->modelId():"none"),
    wallModelId_(wallModel?wallModel->modelId():"none"),
    meshFingerprint_(condensedCarbonSectorMeshFingerprint(mesh_))
{
    options_.validate();
    if(kineticsModelId_.empty()||wallModelId_.empty())
        throw std::invalid_argument(
            "condensed-carbon model identifiers must be nonempty");
    modelFingerprint_=condensedCarbonModelFingerprint(
        options_,registry_,kineticsModelId_,wallModelId_);
}

void SectorCondensedCarbonOperator::initialize() const
{
    if(initialized_) return;
    state_.assign(
        static_cast<std::size_t>(mesh_.Nx())*mesh_.Nr()*mesh_.Ntheta(),
        options_.initialState);
    wall_.assign(
        static_cast<std::size_t>(mesh_.Nx())*mesh_.Ntheta(),
        WallCarbonInventory());
    initialized_=true;
}

void SectorCondensedCarbonOperator::validateDomainLayout(
    const Domain& domain) const
{
    if(!meshMatches(domain))
        throw std::runtime_error(
            "condensed-carbon v1 status=INVALID_INPUT detail=DOMAIN_MESH_MISMATCH");
    if(domain.Nx()<4||domain.Ny()<4||domain.Nz()<4)
        throw std::runtime_error(
            "condensed-carbon v1 status=INVALID_INPUT detail=INSUFFICIENT_HALOS");
    if(!domain.stateReady())
        throw std::runtime_error(
            "condensed-carbon v1 status=INVALID_INPUT detail=DOMAIN_NOT_READY");
    if(domain.speciesNames().size()!=static_cast<std::size_t>(domain.nSpecies()))
        throw std::runtime_error(
            "condensed-carbon v1 status=SPECIES_MISMATCH");
    if(domain.energyConvention()
       !=EnergyConvention::SENSIBLE_INTERNAL_ENERGY_EXPLICIT_CHEMICAL_HEAT)
        throw std::runtime_error(
            "condensed-carbon v1 status=UNSUPPORTED_ENERGY_CONVENTION");

    const int wallCell=domain.Ny()-2;
    const int wallGhost=domain.Ny()-1;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const FaceAreaVector area=
            mesh_.radialAreaVector(i,wallGhost,k,+1);
        const double interior=dot(domain.W(i,wallCell,k).U,area);
        const double ghost=dot(domain.W(i,wallGhost,k).U,area);
        const double centeredFlux=0.5*(interior+ghost);
        const double scale=0.5*(std::abs(interior)+std::abs(ghost));
        if(!closedBoundaryFlux(centeredFlux,scale))
            throw std::runtime_error(
                "condensed-carbon v1 status=INVALID_INPUT "
                "detail=WALL_NORMAL_VOLUME_FLUX");
    }

    const int lowerCell=1;
    const int lowerGhost=0;
    const int upperCell=domain.Nz()-2;
    const int upperGhost=domain.Nz()-1;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    {
        const FaceAreaVector lowerArea=
            mesh_.azimuthalAreaVector(i,j,lowerCell,-1);
        const FaceAreaVector upperArea=
            mesh_.azimuthalAreaVector(i,j,upperGhost,+1);
        const double lowerFlux=0.5*(
            dot(domain.W(i,j,lowerCell).U,lowerArea)
            +dot(domain.W(i,j,lowerGhost).U,lowerArea));
        const double upperFlux=0.5*(
            dot(domain.W(i,j,upperCell).U,upperArea)
            +dot(domain.W(i,j,upperGhost).U,upperArea));
        if(!closedBoundaryFlux(
                lowerFlux+upperFlux,
                std::abs(lowerFlux)+std::abs(upperFlux)))
            throw std::runtime_error(
                "condensed-carbon v1 status=INVALID_INPUT "
                "detail=PERIODIC_VOLUME_FLUX_MISMATCH");
    }
}

void SectorCondensedCarbonOperator::applyFieldBoundaries(
    const Domain& domain,
    std::vector<CondensedCarbonState>& field) const
{
    const int nx=domain.Nx(),ny=domain.Ny(),nz=domain.Nz();
    for(int j=0;j<ny;++j) for(int k=0;k<nz;++k)
    {
        field[static_cast<std::size_t>(index(0,j,k))]=options_.inletState;
        field[static_cast<std::size_t>(index(nx-1,j,k))]=
            field[static_cast<std::size_t>(index(nx-2,j,k))];
    }
    for(int i=0;i<nx;++i) for(int k=0;k<nz;++k)
    {
        field[static_cast<std::size_t>(index(i,0,k))]=
            field[static_cast<std::size_t>(index(i,1,k))];
        field[static_cast<std::size_t>(index(i,ny-1,k))]=
            field[static_cast<std::size_t>(index(i,ny-2,k))];
    }
    for(int i=0;i<nx;++i) for(int j=0;j<ny;++j)
    {
        field[static_cast<std::size_t>(index(i,j,0))]=
            field[static_cast<std::size_t>(index(i,j,nz-2))];
        field[static_cast<std::size_t>(index(i,j,nz-1))]=
            field[static_cast<std::size_t>(index(i,j,1))];
    }
}

double SectorCondensedCarbonOperator::scalarOutwardFlux(
    const Domain& domain,
    const std::vector<CondensedCarbonState>& field,
    int i,int j,int k,int ni,int nj,int nk,
    const FaceAreaVector& area,
    bool numberField) const
{
    const double volumeFlux=0.5*(
        dot(domain.W(i,j,k).U,area)+dot(domain.W(ni,nj,nk).U,area));
    const CondensedCarbonState& upwind=volumeFlux>=0.0
        ?field[static_cast<std::size_t>(index(i,j,k))]
        :field[static_cast<std::size_t>(index(ni,nj,nk))];
    return volumeFlux*(numberField?upwind.numberDensity:upwind.rhoSoot);
}

double SectorCondensedCarbonOperator::particleDiameter(
    const CondensedCarbonState& state) const
{
    if(options_.momentModel==CondensedCarbonMomentModel::ONE_MOMENT)
        return options_.oneMomentParticleDiameter;
    if(state.rhoSoot==0.0&&state.numberDensity==0.0)
        return std::numeric_limits<double>::quiet_NaN();
    if(state.rhoSoot<=0.0||state.numberDensity<=0.0)
        return std::numeric_limits<double>::quiet_NaN();
    const double volumePerParticle=
        state.rhoSoot/(options_.materialDensity*state.numberDensity);
    return std::cbrt(6.0*volumePerParticle/pi);
}

SootApplicabilityReport SectorCondensedCarbonOperator::applicability(
    const Domain& domain,
    int i,int j,int k,
    const ConservedVars& fluid,
    const CondensedCarbonState& carbon) const
{
    if(carbon.rhoSoot==0.0)
    {
        SootApplicabilityReport result;
        result.massLoading=0.0;
        result.volumeFraction=0.0;
        result.heatCapacityRatio=0.0;
        result.stokesNumber=0.0;
        result.status=options_.productionValidated
            ?CondensedCarbonApplicabilityStatus::SUPPORTED
            :CondensedCarbonApplicabilityStatus::WARNING_RESEARCH_ONLY;
        return result;
    }
    const thermo::EquilibriumState& equilibrium=domain.Q(i,j,k);
    if(!equilibrium.thermodynamics)
        return SootApplicabilityReport();
    const TransportProps::Result transport=
        TransportProps::evaluateCommittedLegacyGas(
            domain.speciesNames(),domain.W(i,j,k),equilibrium,
            domain.thermoModel(),registry_);
    if(!transport.success()) return SootApplicabilityReport();
    const double dx=mesh_.xFace(i+1)-mesh_.xFace(i);
    const double dr=mesh_.wallCenter(i)
        *(mesh_.etaFace(j+1)-mesh_.etaFace(j));
    const double ds=mesh_.radiusCenter(i,j)
        *(mesh_.thetaFace(k+1)-mesh_.thetaFace(k));
    const double lengthScale=std::max(
        std::min(dx,std::min(dr,ds)),1.0e-12);
    const PrimitiveVars& primitive=domain.W(i,j,k);
    const double speed=std::sqrt(
        primitive.U[0]*primitive.U[0]
        +primitive.U[1]*primitive.U[1]
        +primitive.U[2]*primitive.U[2]);
    const double signalSpeed=speed+equilibrium.acoustic.soundSpeed();
    CondensedCarbonApplicabilityInputs inputs;
    inputs.rhoFluid=fluid.rho;
    inputs.rhoSoot=carbon.rhoSoot;
    inputs.fluidCp=equilibrium.thermodynamics->Cp;
    inputs.dynamicViscosity=transport.point.mu;
    inputs.fluidTimeScale=lengthScale/std::max(signalSpeed,1.0e-12);
    inputs.particleDiameter=particleDiameter(carbon);
    return evaluateCondensedCarbonApplicability(options_,inputs);
}

void SectorCondensedCarbonOperator::beginStepTransaction(
    const Domain&) const
{
    if(stepTransactionActive_)
        throw std::logic_error(
            "nested condensed-carbon step transaction");
    stateStepStart_=state_;
    wallStepStart_=wall_;
    initializedStepStart_=initialized_;
    boundaryMassStepStart_=sootBoundaryMassNetIntoDomain_;
    boundaryNumberStepStart_=sootBoundaryNumberNetIntoDomain_;
    applicabilityStepStart_=lastApplicability_;
    stepTransactionActive_=true;
}

void SectorCondensedCarbonOperator::acceptStepTransaction() const noexcept
{
    stateStepStart_.clear();
    wallStepStart_.clear();
    stepTransactionActive_=false;
}

void SectorCondensedCarbonOperator::rollbackStepTransaction() const noexcept
{
    if(!stepTransactionActive_) return;
    state_.swap(stateStepStart_);
    wall_.swap(wallStepStart_);
    initialized_=initializedStepStart_;
    sootBoundaryMassNetIntoDomain_=boundaryMassStepStart_;
    sootBoundaryNumberNetIntoDomain_=boundaryNumberStepStart_;
    lastApplicability_=applicabilityStepStart_;
    stepTransactionActive_=false;
}

void SectorCondensedCarbonOperator::averageWithStepStartState() const
{
    if(!stepTransactionActive_)
        throw std::logic_error(
            "condensed-carbon RK2 transaction is inactive");
    const std::size_t expectedState=static_cast<std::size_t>(mesh_.Nx())
        *mesh_.Nr()*mesh_.Ntheta();
    const std::size_t expectedWall=static_cast<std::size_t>(mesh_.Nx())
        *mesh_.Ntheta();
    if(state_.size()!=expectedState||wall_.size()!=expectedWall)
        throw std::logic_error(
            "condensed-carbon RK2 state size changed");
    if(initializedStepStart_
       &&(state_.size()!=stateStepStart_.size()
          ||wall_.size()!=wallStepStart_.size()))
        throw std::logic_error(
            "condensed-carbon RK2 start state size changed");
    for(std::size_t id=0;id<state_.size();++id)
    {
        const CondensedCarbonState& start=initializedStepStart_
            ?stateStepStart_[id]:options_.initialState;
        state_[id].rhoSoot=0.5*(start.rhoSoot+state_[id].rhoSoot);
        state_[id].numberDensity=0.5*(
            start.numberDensity+state_[id].numberDensity);
    }
    for(std::size_t id=0;id<wall_.size();++id)
    {
        const double start=initializedStepStart_
            ?wallStepStart_[id].massPerArea:0.0;
        wall_[id].massPerArea=0.5*(start+wall_[id].massPerArea);
    }
    sootBoundaryMassNetIntoDomain_=0.5*(
        boundaryMassStepStart_+sootBoundaryMassNetIntoDomain_);
    sootBoundaryNumberNetIntoDomain_=0.5*(
        boundaryNumberStepStart_+sootBoundaryNumberNetIntoDomain_);
}

void SectorCondensedCarbonOperator::apply(Domain& domain,double dt) const
{
    validateDomainLayout(domain);
    if(!std::isfinite(dt)||dt<0.0)
        throw std::invalid_argument(
            "invalid condensed-carbon timestep");
    initialize();
    std::vector<CondensedCarbonState> transportBase=state_;
    applyFieldBoundaries(domain,transportBase);

    std::vector<CondensedCarbonState> nextState=transportBase;
    std::vector<WallCarbonInventory> nextWall=wall_;
    std::vector<ConservedVars> fluidCandidates;
    fluidCandidates.reserve(static_cast<std::size_t>(domain.N()));
    for(int k=0;k<domain.Nz();++k) for(int j=0;j<domain.Ny();++j)
    for(int i=0;i<domain.Nx();++i)
        fluidCandidates.push_back(domain.Unew(i,j,k));

    double boundaryMassIncrement=0.0;
    double boundaryNumberIncrement=0.0;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        double massResidual=0.0,numberResidual=0.0;
        auto addFace=[&](int ni,int nj,int nk,const FaceAreaVector& area)
        {
            massResidual+=scalarOutwardFlux(
                domain,transportBase,i,j,k,ni,nj,nk,area,false);
            if(numberDensityEnabled())
                numberResidual+=scalarOutwardFlux(
                    domain,transportBase,i,j,k,ni,nj,nk,area,true);
        };
        addFace(i-1,j,k,mesh_.axialAreaVector(i,j,k,-1));
        addFace(i+1,j,k,mesh_.axialAreaVector(i+1,j,k,+1));
        addFace(i,j-1,k,mesh_.radialAreaVector(i,j,k,-1));
        addFace(i,j+1,k,mesh_.radialAreaVector(i,j+1,k,+1));
        addFace(i,j,k-1,mesh_.azimuthalAreaVector(i,j,k,-1));
        addFace(i,j,k+1,mesh_.azimuthalAreaVector(i,j,k+1,+1));
        const double factor=dt/mesh_.volume(i,j,k);
        CondensedCarbonState& candidate=
            nextState[static_cast<std::size_t>(index(i,j,k))];
        candidate.rhoSoot=
            transportBase[static_cast<std::size_t>(index(i,j,k))].rhoSoot
            -factor*massResidual;
        candidate.numberDensity=
            transportBase[static_cast<std::size_t>(index(i,j,k))].numberDensity
            -factor*numberResidual;
        if(!candidate.valid(options_.momentModel))
        {
            domain.recordFailure(
                Domain::FailureCode::TRANSACTION_PROTOCOL,index(i,j,k),name(),
                "CONDENSED_CARBON_TRANSPORT_INVALID_STATE");
            return;
        }
        if(i==1)
        {
            const double massFlux=scalarOutwardFlux(
                domain,transportBase,i,j,k,0,j,k,
                mesh_.axialAreaVector(i,j,k,-1),false);
            boundaryMassIncrement-=dt*massFlux;
            if(numberDensityEnabled())
                boundaryNumberIncrement-=dt*scalarOutwardFlux(
                    domain,transportBase,i,j,k,0,j,k,
                    mesh_.axialAreaVector(i,j,k,-1),true);
        }
        if(i==domain.Nx()-2)
        {
            const double massFlux=scalarOutwardFlux(
                domain,transportBase,i,j,k,i+1,j,k,
                mesh_.axialAreaVector(i+1,j,k,+1),false);
            boundaryMassIncrement-=dt*massFlux;
            if(numberDensityEnabled())
                boundaryNumberIncrement-=dt*scalarOutwardFlux(
                    domain,transportBase,i,j,k,i+1,j,k,
                    mesh_.axialAreaVector(i+1,j,k,+1),true);
        }
    }

    if(kinetics_)
    {
        for(int i=1;i<domain.Nx()-1;++i)
        for(int j=1;j<domain.Ny()-1;++j)
        for(int k=1;k<domain.Nz()-1;++k)
        {
            const int id=index(i,j,k);
            const CondensedCarbonCellContext context={
                i,j,k,domain.U(i,j,k),domain.W(i,j,k),domain.Q(i,j,k),
                transportBase[static_cast<std::size_t>(id)]};
            const CondensedCarbonExchange exchange=kinetics_->evaluate(context);
            const CondensedCarbonExchangeResult result=
                CondensedCarbonExchangeContract::evaluate(
                    fluidCandidates[static_cast<std::size_t>(id)],
                    domain.W(i,j,k),nextState[static_cast<std::size_t>(id)],
                    domain.speciesNames(),registry_,options_,exchange,dt);
            if(!result.accepted)
            {
                domain.recordFailure(
                    Domain::FailureCode::TRANSACTION_PROTOCOL,id,name(),
                    "CONDENSED_CARBON_EXCHANGE_REJECTED: "+result.reason);
                return;
            }
            fluidCandidates[static_cast<std::size_t>(id)]=result.fluidCandidate;
            nextState[static_cast<std::size_t>(id)]=result.carbonCandidate;
        }
    }

    if(wallModel_)
    {
        const int cellJ=domain.Ny()-2;
        const int wallFace=domain.Ny()-1;
        for(int i=1;i<domain.Nx()-1;++i)
        for(int k=1;k<domain.Nz()-1;++k)
        {
            const int id=index(i,cellJ,k);
            const int wid=wallIndex(i,k);
            const WallCarbonCellContext context={
                i,k,domain.U(i,cellJ,k),domain.W(i,cellJ,k),
                domain.Q(i,cellJ,k),nextState[static_cast<std::size_t>(id)],
                nextWall[static_cast<std::size_t>(wid)]};
            const WallCarbonTransfer transfer=wallModel_->evaluate(context);
            if(!finiteNonnegative(transfer.sootMassFluxToWall)
               ||!finiteNonnegative(transfer.sootNumberFluxToWall)
               ||!std::isfinite(transfer.wallHeatFlux)
               ||transfer.wallHeatFlux!=0.0
               ||(!numberDensityEnabled()
                  &&transfer.sootNumberFluxToWall!=0.0))
            {
                domain.recordFailure(
                    Domain::FailureCode::TRANSACTION_PROTOCOL,id,name(),
                    "CONDENSED_CARBON_WALL_TRANSFER_INVALID");
                return;
            }
            const double area=magnitude(
                mesh_.radialAreaVector(i,wallFace,k,+1));
            const double volume=mesh_.volume(i,cellJ,k);
            CondensedCarbonState& carbon=nextState[static_cast<std::size_t>(id)];
            carbon.rhoSoot-=dt*transfer.sootMassFluxToWall*area/volume;
            carbon.numberDensity-=
                dt*transfer.sootNumberFluxToWall*area/volume;
            nextWall[static_cast<std::size_t>(wid)].massPerArea+=
                dt*transfer.sootMassFluxToWall;
            if(!carbon.valid(options_.momentModel)
               ||!nextWall[static_cast<std::size_t>(wid)].valid())
            {
                domain.recordFailure(
                    Domain::FailureCode::TRANSACTION_PROTOCOL,id,name(),
                    "CONDENSED_CARBON_WALL_TRANSFER_NEGATIVE_INVENTORY");
                return;
            }
        }
    }

    applyFieldBoundaries(domain,nextState);
    SootApplicabilityReport worst;
    worst.massLoading=0.0;
    worst.volumeFraction=0.0;
    worst.heatCapacityRatio=0.0;
    worst.stokesNumber=0.0;
    worst.status=options_.productionValidated
        ?CondensedCarbonApplicabilityStatus::SUPPORTED
        :CondensedCarbonApplicabilityStatus::WARNING_RESEARCH_ONLY;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const int id=index(i,j,k);
        const SootApplicabilityReport report=applicability(
            domain,i,j,k,fluidCandidates[static_cast<std::size_t>(id)],
            nextState[static_cast<std::size_t>(id)]);
        worst.massLoading=std::max(worst.massLoading,report.massLoading);
        worst.volumeFraction=std::max(
            worst.volumeFraction,report.volumeFraction);
        worst.heatCapacityRatio=std::max(
            worst.heatCapacityRatio,report.heatCapacityRatio);
        worst.stokesNumber=std::max(worst.stokesNumber,report.stokesNumber);
        if(applicabilitySeverity(report.status)
           >applicabilitySeverity(worst.status))
            worst.status=report.status;
    }
    if(options_.applicabilityPolicy
            ==CondensedCarbonApplicabilityPolicy::FAIL_CLOSED
       &&!worst.productionSupported())
    {
        domain.recordFailure(
            Domain::FailureCode::TRANSACTION_PROTOCOL,-1,name(),
            std::string("CONDENSED_CARBON_APPLICABILITY_")
                +condensedCarbonApplicabilityStatusName(worst.status));
        return;
    }

    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
        domain.Unew(i,j,k)=
            fluidCandidates[static_cast<std::size_t>(index(i,j,k))];
    state_.swap(nextState);
    wall_.swap(nextWall);
    sootBoundaryMassNetIntoDomain_+=boundaryMassIncrement;
    sootBoundaryNumberNetIntoDomain_+=boundaryNumberIncrement;
    lastApplicability_=worst;
}

double SectorCondensedCarbonOperator::maxStableDt(
    const Domain& domain,double cfl) const
{
    validateDomainLayout(domain);
    if(!std::isfinite(cfl)||cfl<=0.0) return 0.0;
    initialize();
    double dtMin=std::numeric_limits<double>::infinity();
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        double spectralVolume=0.0;
        auto addFace=[&](int ni,int nj,int nk,const FaceAreaVector& area)
        {
            spectralVolume+=std::abs(0.5*(
                dot(domain.W(i,j,k).U,area)
                +dot(domain.W(ni,nj,nk).U,area)));
        };
        addFace(i-1,j,k,mesh_.axialAreaVector(i,j,k,-1));
        addFace(i+1,j,k,mesh_.axialAreaVector(i+1,j,k,+1));
        addFace(i,j-1,k,mesh_.radialAreaVector(i,j,k,-1));
        addFace(i,j+1,k,mesh_.radialAreaVector(i,j+1,k,+1));
        addFace(i,j,k-1,mesh_.azimuthalAreaVector(i,j,k,-1));
        addFace(i,j,k+1,mesh_.azimuthalAreaVector(i,j,k+1,+1));
        if(spectralVolume>0.0)
            dtMin=std::min(
                dtMin,cfl*mesh_.volume(i,j,k)/spectralVolume);
        if(kinetics_)
        {
            const int id=index(i,j,k);
            const CondensedCarbonCellContext context={
                i,j,k,domain.U(i,j,k),domain.W(i,j,k),domain.Q(i,j,k),
                state_[static_cast<std::size_t>(id)]};
            const double sourceDt=kinetics_->maxStableDt(context);
            if(std::isnan(sourceDt)||sourceDt<=0.0) return 0.0;
            if(std::isfinite(sourceDt)) dtMin=std::min(dtMin,sourceDt);
        }
    }
    return dtMin;
}

bool SectorCondensedCarbonOperator::meshMatches(
    const Domain& domain) const noexcept
{
    return mesh_.Nx()==domain.Nx()
        &&mesh_.Nr()==domain.Ny()
        &&mesh_.Ntheta()==domain.Nz();
}

double SectorCondensedCarbonOperator::rhoSoot(int i,int j,int k) const
{
    initialize();
    return state_.at(static_cast<std::size_t>(index(i,j,k))).rhoSoot;
}

double SectorCondensedCarbonOperator::numberDensity(int i,int j,int k) const
{
    initialize();
    if(!numberDensityEnabled())
        throw std::logic_error(
            "number density requested from one-moment condensed carbon");
    return state_.at(static_cast<std::size_t>(index(i,j,k))).numberDensity;
}

double SectorCondensedCarbonOperator::wallCarbonMassPerArea(
    int i,int k) const
{
    initialize();
    return wall_.at(static_cast<std::size_t>(wallIndex(i,k))).massPerArea;
}

const std::vector<CondensedCarbonState>&
SectorCondensedCarbonOperator::stateField() const
{
    initialize();
    return state_;
}

const std::vector<WallCarbonInventory>&
SectorCondensedCarbonOperator::wallField() const
{
    initialize();
    return wall_;
}

bool SectorCondensedCarbonOperator::fieldsAreRestorable(
    const std::vector<CondensedCarbonState>& state,
    const std::vector<WallCarbonInventory>& wall,
    double boundaryMass,
    double boundaryNumber) const noexcept
{
    const std::size_t expectedState=static_cast<std::size_t>(mesh_.Nx())
        *mesh_.Nr()*mesh_.Ntheta();
    const std::size_t expectedWall=static_cast<std::size_t>(mesh_.Nx())
        *mesh_.Ntheta();
    if(state.size()!=expectedState||wall.size()!=expectedWall
       ||!std::isfinite(boundaryMass)||!std::isfinite(boundaryNumber)
       ||(!numberDensityEnabled()&&boundaryNumber!=0.0))
        return false;
    for(const CondensedCarbonState& value:state)
        if(!value.valid(options_.momentModel)) return false;
    for(const WallCarbonInventory& value:wall)
        if(!value.valid()) return false;
    return true;
}

void SectorCondensedCarbonOperator::restoreFields(
    std::vector<CondensedCarbonState> state,
    std::vector<WallCarbonInventory> wall,
    double boundaryMass,
    double boundaryNumber)
{
    if(stepTransactionActive_)
        throw std::logic_error(
            "cannot restore condensed carbon during a step transaction");
    if(!fieldsAreRestorable(state,wall,boundaryMass,boundaryNumber))
        throw std::invalid_argument(
            "invalid condensed-carbon checkpoint fields");
    state_.swap(state);
    wall_.swap(wall);
    sootBoundaryMassNetIntoDomain_=boundaryMass;
    sootBoundaryNumberNetIntoDomain_=boundaryNumber;
    initialized_=true;
}

CondensedCarbonInventoryAudit SectorCondensedCarbonOperator::audit(
    const Domain& domain) const
{
    validateDomainLayout(domain);
    initialize();
    CondensedCarbonInventoryAudit result;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const double volume=mesh_.volume(i,j,k);
        const ConservedVars& fluid=domain.U(i,j,k);
        const PrimitiveVars& primitive=domain.W(i,j,k);
        const CondensedCarbonState& carbon=
            state_[static_cast<std::size_t>(index(i,j,k))];
        result.fluidMass+=fluid.rho*volume;
        result.sootMass+=carbon.rhoSoot*volume;
        result.fluidResolvedEnergy+=fluid.rhoE*volume;
        const double velocitySquared=
            primitive.U[0]*primitive.U[0]
            +primitive.U[1]*primitive.U[1]
            +primitive.U[2]*primitive.U[2];
        result.sootDiagnosticKineticEnergy+=
            0.5*carbon.rhoSoot*velocitySquared*volume;
        for(int species=0;species<domain.nSpecies();++species)
        {
            const int global=registry_.index(domain.speciesNames()[species]);
            if(global<0)
                throw std::runtime_error(
                    "condensed-carbon audit species mismatch");
            const double moles=fluid.rhoY[species]*volume
                /registry_.molecularWeight(global);
            result.fluidElementMoles[0]+=moles*registry_.nC(global);
            result.fluidElementMoles[1]+=moles*registry_.nH(global);
            result.fluidElementMoles[2]+=moles*registry_.nO(global);
            result.fluidElementMoles[3]+=moles*registry_.nN(global);
            result.formationEnergyLedger+=
                moles*registry_.formationEnthalpy(global);
        }
    }
    const int wallFace=domain.Ny()-1;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const double area=magnitude(
            mesh_.radialAreaVector(i,wallFace,k,+1));
        result.wallCarbonMass+=
            wall_[static_cast<std::size_t>(wallIndex(i,k))].massPerArea*area;
    }
    const double carbonMass=result.sootMass+result.wallCarbonMass;
    const double carbonMoles=carbonMass/options_.carbonMolarMass;
    result.totalElementMoles=result.fluidElementMoles;
    result.totalElementMoles[0]+=carbonMoles;
    result.formationEnergyLedger+=
        carbonMoles*options_.carbonFormationEnthalpy;
    result.materialMass=
        result.fluidMass+result.sootMass+result.wallCarbonMass;
    result.modeledEnergyLedger=result.fluidResolvedEnergy
        +result.formationEnergyLedger+result.sootDiagnosticKineticEnergy;
    result.sootBoundaryMassNetIntoDomain=
        sootBoundaryMassNetIntoDomain_;
    result.sootBoundaryNumberNetIntoDomain=
        sootBoundaryNumberNetIntoDomain_;
    result.physicalEnergyComplete=false;
    result.wallDepositionEnergyModeled=wallModel_==nullptr;
    return result;
}

} // End namespace Foam
