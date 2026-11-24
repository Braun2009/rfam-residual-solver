#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "CryogenicTransportTable.H"
#include "GERGResidualTable.H"
#include "SectorCheckpoint.H"
#include "SectorConservativeGridTransfer.H"
#include "SectorFASDefect.H"
#include "SectorFrozenSteadyResidualContext.H"
#include "SectorMeshMetrics.H"
#include "SectorResidualField.H"
#include "SectorTurbulenceOperator.H"
#include "SpeciesRegistry.H"
#include "TransportProps.H"

using namespace Foam;

namespace
{

bool allowedFrozenInputMode(const std::string& mode)
{
    return mode=="finite-rate-les"||mode=="frozen-chemistry-les";
}

struct Level
{
    std::string checkpointPath;
    SectorCheckpointInfo info;
    SectorMeshMetrics mesh;
    Domain domain;
    std::unique_ptr<SectorTurbulenceOperator> turbulence;
    SectorCheckpointMetadata metadata;

    Level(
        const std::string& path,const std::string& profile,
        SpeciesRegistry& registry,GERGResidualTable& gerg)
        :checkpointPath(path),
         info(SectorCheckpoint::inspect(path)),
         mesh(profile,info.nx,info.ny,info.nz,120.0),
         domain(info.nx,info.ny,info.nz,
             mesh.length(),0.065,mesh.sectorAngle(),info.nSpecies)
    {
        if(!allowedFrozenInputMode(info.metadata.mode))
            throw std::runtime_error(
                "FAS source builder requires finite-rate-les or "
                "frozen-chemistry-les checkpoints");
        if(info.species!=registry.names())
            throw std::runtime_error(
                "FAS source checkpoint species order is not the frozen "
                "9-species registry");
        domain.setSpecies(info.species);
        domain.thermoModel().setGERGTable(&gerg);
        turbulence=std::make_unique<SectorTurbulenceOperator>(
            mesh,registry,0.05,0.001,0.25,true,
            SectorTurbulenceOperator::TurbulenceMode::WALE);
        metadata=SectorCheckpoint::read(
            checkpointPath,domain,turbulence.get());
    }
};

std::vector<ConservedVars> captureField(const Domain& domain)
{
    std::vector<ConservedVars> result;
    result.reserve(static_cast<std::size_t>(domain.N()));
    for(int cell=0;cell<domain.N();++cell) result.push_back(domain.U(cell));
    return result;
}

double maximumFieldDifference(const Domain& left,const Domain& right)
{
    if(left.Nx()!=right.Nx()||left.Ny()!=right.Ny()
        ||left.Nz()!=right.Nz()||left.nSpecies()!=right.nSpecies())
        return std::numeric_limits<double>::infinity();
    double maximum=0.0;
    auto compare=[&](double a,double b)
    {
        maximum=std::max(maximum,
            std::abs(a-b)/std::max({std::abs(a),std::abs(b),1.0e-300}));
    };
    for(int i=1;i<left.Nx()-1;++i)
    for(int j=1;j<left.Ny()-1;++j)
    for(int k=1;k<left.Nz()-1;++k)
    {
        const ConservedVars& a=left.U(i,j,k);
        const ConservedVars& b=right.U(i,j,k);
        compare(a.rho,b.rho);
        compare(a.rhoE,b.rhoE);
        for(std::size_t component=0;component<a.rhoU.size();++component)
            compare(a.rhoU[component],b.rhoU[component]);
        for(std::size_t component=0;component<a.rhoY.size();++component)
            compare(a.rhoY[component],b.rhoY[component]);
        compare(a.rhoZ,b.rhoZ);
        compare(a.rhoC,b.rhoC);
    }
    return maximum;
}

} // namespace

int main(int argc,char** argv)
{
    std::string finePath,coarsePath,outputPath;
    double referenceStep=1.0e-9;
    for(int argument=1;argument<argc;++argument)
    {
        if(std::strcmp(argv[argument],"--fine")==0&&argument+1<argc)
            finePath=argv[++argument];
        else if(std::strcmp(argv[argument],"--coarse")==0
            &&argument+1<argc) coarsePath=argv[++argument];
        else if(std::strcmp(argv[argument],"--output")==0
            &&argument+1<argc) outputPath=argv[++argument];
        else if(std::strcmp(argv[argument],"--reference-step")==0
            &&argument+1<argc) referenceStep=std::atof(argv[++argument]);
    }
    if(finePath.empty()||coarsePath.empty()||outputPath.empty())
    {
        std::cerr << "usage: rfam_sector_build_fas_source "
            "--fine FINE_CHECKPOINT --coarse COARSE_CHECKPOINT "
            "--output SOURCE [--reference-step DT]\n";
        return 2;
    }
    if(!std::isfinite(referenceStep)||referenceStep<=0.0)
    {
        std::cerr << "invalid FAS residual reference step\n";
        return 2;
    }

    try
    {
        const std::string profile=vortexProfilePath();
        SpeciesRegistry registry;
        GERGResidualTable gerg;
        const std::string gergPath=std::string(RFAM_SOURCE_DIR)
            +"/tests/reference/gerg_residual_table_v3.csv";
        if(!gerg.load(gergPath))
            throw std::runtime_error("cannot load GERG table: "+gergPath);
        CryogenicTransportTable transport;
        const std::string transportPath=std::string(RFAM_SOURCE_DIR)
            +"/tests/reference/cryogenic_transport_runtime_coolprop72_v1.csv";
        if(!transport.load(transportPath))
            throw std::runtime_error(
                "cannot load cryogenic transport table: "+transportPath);
        TransportProps::registerCryogenicTable(&transport);

        Level fine(finePath,profile,registry,gerg);
        Level coarse(coarsePath,profile,registry,gerg);
        const std::uint64_t profileHash=SectorCheckpoint::hashFile(profile);
        if(fine.metadata.profileHash!=profileHash
            ||coarse.metadata.profileHash!=profileHash)
            throw std::runtime_error("FAS checkpoint profile hash mismatch");
        if(fine.metadata.time!=coarse.metadata.time)
            throw std::runtime_error("FAS fine/coarse physical clocks differ");

        Domain restricted(
            coarse.domain.Nx(),coarse.domain.Ny(),coarse.domain.Nz(),
            coarse.mesh.length(),0.065,coarse.mesh.sectorAngle(),
            coarse.domain.nSpecies());
        restricted.setSpecies(coarse.domain.speciesNames());
        restricted.thermoModel().setGERGTable(&gerg);
        SectorConservativeGridTransfer::apply(
            fine.domain,fine.mesh,nullptr,
            restricted,coarse.mesh,nullptr);
        const double restrictionError=maximumFieldDifference(
            restricted,coarse.domain);
        if(restrictionError>1.0e-9)
            throw std::runtime_error(
                "coarse FAS state is not the aligned restriction of fine state");

        SectorFrozenSteadyResidualContext fineContext(
            fine.domain,fine.mesh,registry,&gerg,
            fine.turbulence.get(),fine.metadata.time);
        SectorFrozenSteadyResidualContext coarseContext(
            coarse.domain,coarse.mesh,registry,&gerg,
            coarse.turbulence.get(),coarse.metadata.time);
        SteadyResidualEvaluator::Result fineResidual,coarseResidual;
        if(!fineContext.evaluate(referenceStep,fineResidual))
            throw std::runtime_error(
                "fine FAS residual evaluation failed at cell "
                +std::to_string(fine.domain.failure().cell)+": "
                +fine.domain.failure().detail);
        if(!coarseContext.evaluate(referenceStep,coarseResidual))
            throw std::runtime_error(
                "coarse FAS residual evaluation failed at cell "
                +std::to_string(coarse.domain.failure().cell)+": "
                +coarse.domain.failure().detail);

        const SectorFASDefect::Result defect=SectorFASDefect::build(
            fine.domain,fine.mesh,fineResidual,
            coarse.domain,coarse.mesh,coarseResidual);
        if(defect.maximumInitialConsistencyError>1.0e-10
            ||defect.maximumRestrictionInventoryError>1.0e-10)
            throw std::runtime_error(
                "FAS source construction failed its consistency gate");

        SectorResidualField::Data source;
        source.metadata.kind=SectorResidualField::Kind::FAS_SOURCE;
        source.metadata.nx=coarse.domain.Nx();
        source.metadata.ny=coarse.domain.Ny();
        source.metadata.nz=coarse.domain.Nz();
        source.metadata.nSpecies=coarse.domain.nSpecies();
        source.metadata.referenceStep=referenceStep;
        source.metadata.time=coarse.metadata.time;
        source.metadata.iteration=coarse.metadata.iteration;
        source.metadata.profileHash=profileHash;
        source.metadata.stateCheckpointHash=
            SectorCheckpoint::hashFile(coarsePath);
        source.metadata.fineStateCheckpointHash=
            SectorCheckpoint::hashFile(finePath);
        source.metadata.mode="frozen-chemistry-les";
#ifdef RFAM_GIT_COMMIT
        source.metadata.gitCommit=RFAM_GIT_COMMIT;
#else
        source.metadata.gitCommit="unknown";
#endif
        source.metadata.thermoModelId=coarse.domain.thermoModelId();
        source.metadata.species=coarse.domain.speciesNames();
        source.active=coarseResidual.active;
        source.field=defect.coarseSource;
        SectorResidualField::write(outputPath,source);

        std::cout << std::setprecision(12)
                  << "fine/coarse_active="
                  << (fine.domain.Nx()-2)*(fine.domain.Ny()-2)
                        *(fine.domain.Nz()-2) << "/"
                  << (coarse.domain.Nx()-2)*(coarse.domain.Ny()-2)
                        *(coarse.domain.Nz()-2) << "\n"
                  << "restriction_error=" << restrictionError << "\n"
                  << "fine/coarse_residual_rms="
                  << SectorFrozenSteadyResidualContext::normalizedRms(
                        fineResidual,fine.domain) << "/"
                  << SectorFrozenSteadyResidualContext::normalizedRms(
                        coarseResidual,coarse.domain) << "\n"
                  << "fas_consistency/inventory="
                  << defect.maximumInitialConsistencyError << "/"
                  << defect.maximumRestrictionInventoryError << "\n"
                  << "coarse/fine_checkpoint_hash="
                  << source.metadata.stateCheckpointHash << "/"
                  << source.metadata.fineStateCheckpointHash << "\n"
                  << "output=" << outputPath << " git="
                  << source.metadata.gitCommit << "\n"
                  << "SECTOR_FAS_SOURCE_BUILD_PASS\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << "SECTOR_FAS_SOURCE_BUILD_FAIL: "
                  << error.what() << "\n";
        return 1;
    }
}
