#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "CryogenicTransportTable.H"
#include "GERGResidualTable.H"
#include "SectorCheckpoint.H"
#include "SectorConservativeGridTransfer.H"
#include "SectorFASDefect.H"
#include "SectorFrozenSteadyResidualContext.H"
#include "SectorGalerkinCoarseSolve.H"
#include "SectorGalerkinCoarseStencil.H"
#include "SectorMeshMetrics.H"
#include "SectorTurbulenceOperator.H"
#include "SpeciesRegistry.H"
#include "TransportProps.H"

using namespace Foam;

namespace
{

struct LoadedSector
{
    std::string path;
    SectorCheckpointInfo info;
    SectorMeshMetrics mesh;
    Domain domain;
    std::unique_ptr<SectorTurbulenceOperator> turbulence;
    SectorCheckpointMetadata metadata;

    LoadedSector(
        const std::string& checkpointPath,
        const std::string& profile,
        SpeciesRegistry& registry,
        GERGResidualTable& gerg)
        :path(checkpointPath),
         info(SectorCheckpoint::inspect(checkpointPath)),
         mesh(profile,info.nx,info.ny,info.nz,120.0),
         domain(info.nx,info.ny,info.nz,
             mesh.length(),0.065,mesh.sectorAngle(),info.nSpecies)
    {
        if(info.metadata.mode!="finite-rate-les"
            &&info.metadata.mode!="frozen-chemistry-les")
            throw std::runtime_error(
                "Galerkin correction requires finite-rate-les or "
                "frozen-chemistry-les checkpoints");
        if(info.species!=registry.names())
            throw std::runtime_error(
                "Galerkin correction requires the frozen 9-species layout");
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

std::vector<double> captureTemperature(const Domain& domain)
{
    std::vector<double> result(
        static_cast<std::size_t>(domain.N()),300.0);
    for(int cell=0;cell<domain.N();++cell)
        result[static_cast<std::size_t>(cell)]=domain.W(cell).T;
    return result;
}

void addScaled(
    ConservedVars& target,const ConservedVars& source,double scale)
{
    target.rho+=scale*source.rho;
    target.rhoE+=scale*source.rhoE;
    for(std::size_t component=0;component<target.rhoU.size();++component)
        target.rhoU[component]+=scale*source.rhoU[component];
    for(std::size_t component=0;component<target.rhoY.size();++component)
        target.rhoY[component]+=scale*source.rhoY[component];
    target.rhoZ+=scale*source.rhoZ;
    target.rhoC+=scale*source.rhoC;
}

bool projectSpeciesComposition(
    const ConservedVars& base,ConservedVars& candidate,
    double& maximumRelativeProjection)
{
    maximumRelativeProjection=0.0;
    if(!std::isfinite(candidate.rho)||candidate.rho<=0.0
        ||candidate.rhoY.size()!=base.rhoY.size()
        ||candidate.rhoY.empty()) return false;
    double baseSpeciesMass=0.0;
    for(double value:base.rhoY)
    {
        if(!std::isfinite(value)||value<0.0) return false;
        baseSpeciesMass+=value;
    }
    if(!std::isfinite(baseSpeciesMass)||baseSpeciesMass<=0.0)
        return false;

    std::vector<double> densityFollowing(candidate.rhoY.size(),0.0);
    std::vector<double> composition(candidate.rhoY.size(),0.0);
    double compositionSum=0.0;
    for(std::size_t species=0;species<candidate.rhoY.size();++species)
    {
        if(!std::isfinite(candidate.rhoY[species])) return false;
        const double massFraction=base.rhoY[species]/baseSpeciesMass;
        densityFollowing[species]=massFraction*candidate.rho;
        composition[species]=candidate.rhoY[species]
            -densityFollowing[species];
        compositionSum+=composition[species];
    }
    for(std::size_t species=0;species<candidate.rhoY.size();++species)
        composition[species]-=
            (base.rhoY[species]/baseSpeciesMass)*compositionSum;

    double beta=1.0;
    for(std::size_t species=0;species<candidate.rhoY.size();++species)
        if(composition[species]<0.0)
            beta=std::min(beta,
                densityFollowing[species]/(-composition[species]));
    beta=std::clamp(beta,0.0,1.0);

    double sum=0.0;
    std::size_t largest=0;
    for(std::size_t species=0;species<candidate.rhoY.size();++species)
    {
        const double previous=candidate.rhoY[species];
        candidate.rhoY[species]=std::max(
            densityFollowing[species]+beta*composition[species],0.0);
        sum+=candidate.rhoY[species];
        if(candidate.rhoY[species]>candidate.rhoY[largest]) largest=species;
        maximumRelativeProjection=std::max(maximumRelativeProjection,
            std::abs(candidate.rhoY[species]-previous)/candidate.rho);
    }
    candidate.rhoY[largest]+=candidate.rho-sum;
    return candidate.rhoY[largest]>=0.0
        &&std::isfinite(candidate.rhoY[largest]);
}

double maximumRelativeCorrection(
    const Domain& domain,const std::vector<ConservedVars>& correction)
{
    double maximum=0.0;
    for(int i=1;i<domain.Nx()-1;++i)
    for(int j=1;j<domain.Ny()-1;++j)
    for(int k=1;k<domain.Nz()-1;++k)
    {
        const int cell=k*domain.Nx()*domain.Ny()+j*domain.Nx()+i;
        const ConservedVars& state=domain.U(cell);
        const ConservedVars& delta=correction[static_cast<std::size_t>(cell)];
        const double rhoScale=std::max(std::abs(state.rho),1.0e-12);
        const double acousticScale=std::max(
            domain.Q(cell).acoustic.soundSpeed(),1.0);
        const double momentumScale=rhoScale*acousticScale;
        const double energyScale=std::max(
            std::abs(state.rhoE),rhoScale*acousticScale*acousticScale);
        auto update=[&](double value,double scale)
        {
            maximum=std::max(maximum,
                std::abs(value)/std::max(scale,1.0e-300));
        };
        update(delta.rho,rhoScale);
        update(delta.rhoE,energyScale);
        for(double value:delta.rhoU) update(value,momentumScale);
        for(double value:delta.rhoY) update(value,rhoScale);
        update(delta.rhoZ,rhoScale);
        update(delta.rhoC,rhoScale);
    }
    return maximum;
}

double maximumFieldDifference(const Domain& left,const Domain& right)
{
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
    double referenceStep=1.0e-8;
    double pseudoCfl=1.0;
    double maximumUpdate=0.05;
    int coarseSweeps=20;
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
        else if(std::strcmp(argv[argument],"--pseudo-cfl")==0
            &&argument+1<argc) pseudoCfl=std::atof(argv[++argument]);
        else if(std::strcmp(argv[argument],"--maximum-update")==0
            &&argument+1<argc) maximumUpdate=std::atof(argv[++argument]);
        else if(std::strcmp(argv[argument],"--coarse-sweeps")==0
            &&argument+1<argc)
            coarseSweeps=std::clamp(std::atoi(argv[++argument]),1,200);
    }
    if(finePath.empty()||coarsePath.empty()||outputPath.empty())
    {
        std::cerr << "usage: rfam_sector_build_galerkin_correction "
            "--fine FINE_CHECKPOINT --coarse COARSE_CHECKPOINT "
            "--output COARSE_CORRECTED_CHECKPOINT "
            "[--reference-step DT] [--pseudo-cfl CFL] "
            "[--maximum-update VALUE] [--coarse-sweeps N]\n";
        return 2;
    }
    if(!std::isfinite(referenceStep)||referenceStep<=0.0
        ||!std::isfinite(pseudoCfl)||pseudoCfl<=0.0
        ||!std::isfinite(maximumUpdate)||maximumUpdate<=0.0
        ||maximumUpdate>0.5)
    {
        std::cerr << "invalid Galerkin correction controls\n";
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

        LoadedSector fine(finePath,profile,registry,gerg);
        LoadedSector coarse(coarsePath,profile,registry,gerg);
        const std::uint64_t profileHash=SectorCheckpoint::hashFile(profile);
        if(fine.metadata.profileHash!=profileHash
            ||coarse.metadata.profileHash!=profileHash
            ||fine.metadata.time!=coarse.metadata.time)
            throw std::runtime_error(
                "Galerkin fine/coarse profile or physical clock mismatch");

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
                "Galerkin coarse checkpoint is not the restriction of fine state");

        SectorFrozenSteadyResidualContext fineContext(
            fine.domain,fine.mesh,registry,&gerg,
            fine.turbulence.get(),fine.metadata.time);
        SteadyResidualEvaluator::Result fineResidual;
        if(!fineContext.evaluate(referenceStep,fineResidual))
            throw std::runtime_error(
                "Galerkin fine residual evaluation failed at cell "
                +std::to_string(fine.domain.failure().cell)+": "
                +fine.domain.failure().detail);
        const SteadyResidualEvaluator::Result restrictedResidual=
            SectorFASDefect::restrictResidual(
                fine.domain,fine.mesh,fineResidual,
                coarse.domain,coarse.mesh);
        StructuredImplicitStencil fineStencil;
        if(!fineContext.structuredImplicitStencil(fineStencil))
            throw std::runtime_error(
                "Galerkin fine stencil assembly failed at cell "
                +std::to_string(fine.domain.failure().cell)+": "
                +fine.domain.failure().detail);
        const SectorGalerkinCoarseStencil::Result coarseStencil=
            SectorGalerkinCoarseStencil::aggregate(
                fine.domain,fine.mesh,fineStencil,
                coarse.domain,coarse.mesh,pseudoCfl);
        const SectorGalerkinCoarseSolve::Result solved=
            SectorGalerkinCoarseSolve::solve(
                coarse.domain,coarseStencil.stencil,pseudoCfl,
                restrictedResidual.field,coarseSweeps,1.0e-8);

        const double unscaledMaximum=maximumRelativeCorrection(
            coarse.domain,solved.correction);
        double acceptedScale=unscaledMaximum>maximumUpdate
            ?maximumUpdate/unscaledMaximum:1.0;
        std::vector<ConservedVars> candidate;
        const std::vector<ConservedVars> base=captureField(coarse.domain);
        const std::vector<double> temperature=captureTemperature(coarse.domain);
        int recoveryBacktracks=0;
        int projectedCells=0;
        double maximumProjection=0.0;
        bool published=false;
        for(;recoveryBacktracks<=24;++recoveryBacktracks)
        {
            candidate=base;
            projectedCells=0;
            maximumProjection=0.0;
            for(int i=1;i<coarse.domain.Nx()-1;++i)
            for(int j=1;j<coarse.domain.Ny()-1;++j)
            for(int k=1;k<coarse.domain.Nz()-1;++k)
            {
                const int cell=k*coarse.domain.Nx()*coarse.domain.Ny()
                    +j*coarse.domain.Nx()+i;
                ConservedVars& state=candidate[static_cast<std::size_t>(cell)];
                addScaled(state,
                    solved.correction[static_cast<std::size_t>(cell)],
                    acceptedScale);
                if(std::isfinite(state.rho)&&state.rho>0.0)
                {
                    state.rhoZ=std::clamp(state.rhoZ,0.0,state.rho);
                    state.rhoC=std::clamp(state.rhoC,0.0,state.rho);
                }
                double projection=0.0;
                if(!projectSpeciesComposition(
                        base[static_cast<std::size_t>(cell)],state,projection))
                    throw std::runtime_error(
                        "Galerkin coarse composition projection failed");
                if(projection>0.0)
                {
                    ++projectedCells;
                    maximumProjection=std::max(maximumProjection,projection);
                }
            }
            if(coarse.domain.tripped()) coarse.domain.resetTrip();
            if(coarse.domain.replaceConservedField(candidate,temperature))
            {
                published=true;
                break;
            }
            acceptedScale*=0.5;
        }
        if(!published)
            throw std::runtime_error(
                "Galerkin coarse correction recovery line search exhausted");

        SectorCheckpointMetadata metadata=coarse.metadata;
        metadata.mode="frozen-chemistry-les";
#ifdef RFAM_GIT_COMMIT
        metadata.gitCommit=RFAM_GIT_COMMIT;
#endif
#ifdef _OPENMP
        metadata.threadCount=omp_get_max_threads();
#endif
        SectorCheckpoint::write(
            outputPath,coarse.domain,metadata,coarse.turbulence.get());

        std::cout << std::setprecision(12)
                  << "fine/coarse_active="
                  << (fine.domain.Nx()-2)*(fine.domain.Ny()-2)
                        *(fine.domain.Nz()-2) << "/"
                  << (coarse.domain.Nx()-2)*(coarse.domain.Ny()-2)
                        *(coarse.domain.Nz()-2) << "\n"
                  << "restriction_error=" << restrictionError << "\n"
                  << "fine/restricted_residual_rms="
                  << SectorFrozenSteadyResidualContext::normalizedRms(
                        fineResidual,fine.domain) << "/"
                  << SectorFrozenSteadyResidualContext::normalizedRms(
                        restrictedResidual,coarse.domain) << "\n"
                  << "galerkin_diagonal_min/max/dominance/face="
                  << coarseStencil.minimumActualDiagonal << "/"
                  << coarseStencil.maximumActualDiagonal << "/"
                  << coarseStencil.minimumDiagonalDominance << "/"
                  << coarseStencil.maximumFaceCoupling << "\n"
                  << "coarse_solve_sweeps/residual_ratio="
                  << solved.sweeps << "/" << solved.residualRatio << "\n"
                  << "correction_unscaled/accepted_scale/update="
                  << unscaledMaximum << "/" << acceptedScale << "/"
                  << acceptedScale*unscaledMaximum << "\n"
                  << "composition_projected/max_projection/recovery_backtracks="
                  << projectedCells << "/" << maximumProjection << "/"
                  << recoveryBacktracks << "\n"
                  << "fine/coarse_hash="
                  << SectorCheckpoint::hashFile(finePath) << "/"
                  << SectorCheckpoint::hashFile(coarsePath) << "\n"
                  << "output=" << outputPath << " git="
                  << metadata.gitCommit << "\n"
                  << "SECTOR_GALERKIN_CORRECTION_BUILD_PASS\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << "SECTOR_GALERKIN_CORRECTION_BUILD_FAIL: "
                  << error.what() << "\n";
        return 1;
    }
}
