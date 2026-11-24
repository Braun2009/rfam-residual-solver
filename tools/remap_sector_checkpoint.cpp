#include <algorithm>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "Domain.H"
#include "SectorCheckpoint.H"
#include "SectorCheckpointRemap.H"
#include "SectorConservativeGridTransfer.H"
#include "SectorDiagnostics.H"
#include "SectorMeshMetrics.H"
#include "SectorTurbulenceOperator.H"
#include "GERGResidualTable.H"

using namespace Foam;

int main(int argc,char** argv)
{
    std::string input,output,grid;
    for(int a=1;a<argc;++a)
    {
        if(std::strcmp(argv[a],"--input")==0&&a+1<argc) input=argv[++a];
        else if(std::strcmp(argv[a],"--output")==0&&a+1<argc) output=argv[++a];
        else if(std::strcmp(argv[a],"--grid")==0&&a+1<argc) grid=argv[++a];
    }
    if(input.empty()||output.empty()||grid.empty())
    {
        std::cerr << "usage: rfam_sector_remap_checkpoint --input FILE --output FILE "
                     "--grid coarse|smoke|pathfinder|intermediate\n";
        return 2;
    }
    int activeNx=0,activeNr=0,activeNt=0;
    if(grid=="coarse") {activeNx=30;activeNr=10;activeNt=6;}
    else if(grid=="smoke") {activeNx=60;activeNr=20;activeNt=12;}
    else if(grid=="pathfinder") {activeNx=120;activeNr=30;activeNt=20;}
    else if(grid=="intermediate") {activeNx=240;activeNr=50;activeNt=36;}
    else {std::cerr << "unsupported target grid: " << grid << "\n";return 2;}

    const SectorCheckpointInfo info=SectorCheckpoint::inspect(input);
    SpeciesRegistry registry;
    const std::string profile=vortexProfilePath();
    SectorMeshMetrics sourceMesh(profile,info.nx,info.ny,info.nz,120.0);
    Domain source(info.nx,info.ny,info.nz,sourceMesh.length(),0.065,
                  sourceMesh.sectorAngle(),info.nSpecies);
    source.setSpecies(info.species);
    GERGResidualTable gerg;
    const std::string gergPath=std::string(RFAM_SOURCE_DIR)
        +"/tests/reference/gerg_residual_table_v3.csv";
    if(!gerg.load(gergPath))
        throw std::runtime_error(
            "cannot load conservative-remap GERG table: "+gergPath);
    source.thermoModel().setGERGTable(&gerg);
    const bool rans=info.metadata.mode.find("rans")!=std::string::npos;
    const bool les=info.metadata.mode.find("les")!=std::string::npos;
    std::unique_ptr<SectorTurbulenceOperator> sourceTurbulence;
    if(rans||les) sourceTurbulence=
        std::make_unique<SectorTurbulenceOperator>(
            sourceMesh,registry,0.05,0.001,0.25,false,
            les ? SectorTurbulenceOperator::TurbulenceMode::WALE
                : SectorTurbulenceOperator::TurbulenceMode::SST);
    const SectorCheckpointMetadata loaded=SectorCheckpoint::read(
        input,source,sourceTurbulence.get());

    const int nx=activeNx+2,nr=activeNr+2,nt=activeNt+2;
    SectorMeshMetrics targetMesh(profile,nx,nr,nt,120.0);
    Domain target(nx,nr,nt,targetMesh.length(),0.065,targetMesh.sectorAngle(),info.nSpecies);
    target.setSpecies(info.species);
    target.thermoModel().setGERGTable(&gerg);
    std::unique_ptr<SectorTurbulenceOperator> targetTurbulence;
    if(rans||les) targetTurbulence=
        std::make_unique<SectorTurbulenceOperator>(
            targetMesh,registry,0.05,0.001,0.25,false,
            les ? SectorTurbulenceOperator::TurbulenceMode::WALE
                : SectorTurbulenceOperator::TurbulenceMode::SST);
    std::string transferMethod="trilinear";
    try
    {
        const SectorConservativeGridTransfer::Plan transfer=
            SectorConservativeGridTransfer::plan(source,target);
        SectorConservativeGridTransfer::apply(
            source,sourceMesh,sourceTurbulence.get(),
            target,targetMesh,targetTurbulence.get(),
            transfer.direction
                ==SectorConservativeGridTransfer::Direction::PROLONG?4:0);
        transferMethod="aligned-conservative";
    }
    catch(const std::invalid_argument&)
    {
        SectorCheckpointRemap::apply(
            source,sourceTurbulence.get(),target,targetTurbulence.get());
    }

    int threads=1;
#ifdef _OPENMP
    threads=omp_get_max_threads();
#endif
    SectorCheckpointMetadata metadata=loaded;
    metadata.threadCount=threads;
#ifdef RFAM_GIT_COMMIT
    metadata.gitCommit=RFAM_GIT_COMMIT;
#endif
    SectorCheckpoint::write(output,target,metadata,targetTurbulence.get());

    const double sourceMass=SectorDiagnostics(sourceMesh).evaluate(source).massInventory;
    const double targetMass=SectorDiagnostics(targetMesh).evaluate(target).massInventory;
    const double relativeMassDifference=std::abs(targetMass-sourceMass)
        /std::max(std::abs(sourceMass),1e-30);
    std::cout << "source_grid=" << info.nx-2 << "x" << info.ny-2 << "x" << info.nz-2
              << " target_grid=" << activeNx << "x" << activeNr << "x" << activeNt << "\n"
              << "iteration=" << loaded.iteration << " time=" << loaded.time
              << " source_git=" << loaded.gitCommit << " remap_git=" << metadata.gitCommit << "\n"
              << "transfer_method=" << transferMethod << "\n"
              << "mass_source/target/relative_difference=" << sourceMass << "/"
              << targetMass << "/" << relativeMassDifference << "\n"
              << "SECTOR_CHECKPOINT_REMAP_PASS\n";
    return 0;
}
