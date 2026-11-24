#include <algorithm>
#include <cmath>
#include <cstring>
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
#include "Domain.H"
#include "GERGResidualTable.H"
#include "SectorCheckpoint.H"
#include "SectorCoarseCorrection.H"
#include "SectorConservativeGridTransfer.H"
#include "SectorDiagnostics.H"
#include "SectorFASDefect.H"
#include "SectorFrozenSteadyResidualContext.H"
#include "SectorMeshMetrics.H"
#include "SectorTurbulenceOperator.H"
#include "TransportProps.H"

using namespace Foam;

namespace
{

struct LoadedSector
{
    SectorCheckpointInfo info;
    SectorMeshMetrics mesh;
    Domain domain;
    std::unique_ptr<SectorTurbulenceOperator> turbulence;
    SectorCheckpointMetadata metadata;

    LoadedSector(
        const std::string& path,
        const std::string& profile,
        SpeciesRegistry& registry,
        GERGResidualTable& gerg)
        :info(SectorCheckpoint::inspect(path)),
         mesh(profile,info.nx,info.ny,info.nz,120.0),
         domain(info.nx,info.ny,info.nz,
             mesh.length(),0.065,mesh.sectorAngle(),info.nSpecies)
    {
        domain.setSpecies(info.species);
        domain.thermoModel().setGERGTable(&gerg);
        const bool rans=info.metadata.mode.find("rans")!=std::string::npos;
        const bool les=info.metadata.mode.find("les")!=std::string::npos;
        if(rans)
            throw std::runtime_error(
                "offline coarse correction does not support persistent SST state");
        if(les)
            turbulence=std::make_unique<SectorTurbulenceOperator>(
                mesh,registry,0.05,0.001,0.25,false,
                SectorTurbulenceOperator::TurbulenceMode::WALE);
        metadata=SectorCheckpoint::read(path,domain,turbulence.get());
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
    std::vector<double> result(static_cast<std::size_t>(domain.N()),300.0);
    for(int cell=0;cell<domain.N();++cell)
        result[static_cast<std::size_t>(cell)]=domain.W(cell).T;
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

bool allowedModeTransition(
    const std::string& base,const std::string& improved)
{
    if(base==improved) return true;
    const std::string finitePrefix="finite-rate";
    const std::string frozenPrefix="frozen-chemistry";
    return base.rfind(finitePrefix,0)==0
        &&improved==frozenPrefix+base.substr(finitePrefix.size());
}

} // namespace

int main(int argc,char** argv)
{
    std::string fineBasePath,coarseBasePath,coarseImprovedPath,outputPath;
    int smoothingSweeps=0;
    double maximumRelativeCorrection=0.05;
    bool excludeBoundaryParentLayer=true;
    bool requireFineResidualDecrease=false;
    bool reverseCorrection=false;
    bool freezeFineComposition=false;
    bool conservativeParentComposition=false;
    double residualReferenceStep=1.0e-9;
    int residualBacktracks=12;
    for(int argument=1;argument<argc;++argument)
    {
        if(std::strcmp(argv[argument],"--fine-base")==0&&argument+1<argc)
            fineBasePath=argv[++argument];
        else if(std::strcmp(argv[argument],"--coarse-base")==0
            &&argument+1<argc) coarseBasePath=argv[++argument];
        else if(std::strcmp(argv[argument],"--coarse-improved")==0
            &&argument+1<argc) coarseImprovedPath=argv[++argument];
        else if(std::strcmp(argv[argument],"--output")==0&&argument+1<argc)
            outputPath=argv[++argument];
        else if(std::strcmp(argv[argument],"--smoothing-sweeps")==0
            &&argument+1<argc)
            smoothingSweeps=std::clamp(std::atoi(argv[++argument]),0,20);
        else if(std::strcmp(argv[argument],"--max-relative-correction")==0
            &&argument+1<argc)
            maximumRelativeCorrection=std::atof(argv[++argument]);
        else if(std::strcmp(
            argv[argument],"--include-boundary-parent-layer")==0)
            excludeBoundaryParentLayer=false;
        else if(std::strcmp(
            argv[argument],"--require-fine-residual-decrease")==0)
            requireFineResidualDecrease=true;
        else if(std::strcmp(argv[argument],"--reverse-correction")==0)
            reverseCorrection=true;
        else if(std::strcmp(
            argv[argument],"--freeze-fine-composition")==0)
            freezeFineComposition=true;
        else if(std::strcmp(
            argv[argument],"--conservative-parent-composition")==0)
            conservativeParentComposition=true;
        else if(std::strcmp(argv[argument],"--residual-reference-step")==0
            &&argument+1<argc)
            residualReferenceStep=std::atof(argv[++argument]);
        else if(std::strcmp(argv[argument],"--residual-backtracks")==0
            &&argument+1<argc)
            residualBacktracks=std::clamp(std::atoi(argv[++argument]),0,24);
    }
    if(fineBasePath.empty()||coarseBasePath.empty()
        ||coarseImprovedPath.empty()||outputPath.empty())
    {
        std::cerr << "usage: rfam_sector_apply_coarse_correction "
            "--fine-base FILE --coarse-base FILE --coarse-improved FILE "
            "--output FILE [--smoothing-sweeps N] "
            "[--max-relative-correction VALUE] "
            "[--include-boundary-parent-layer] "
            "[--require-fine-residual-decrease] [--reverse-correction] "
            "[--freeze-fine-composition] "
            "[--conservative-parent-composition] "
            "[--residual-reference-step DT] [--residual-backtracks N]\n";
        return 2;
    }

    try
    {
        SpeciesRegistry registry;
        const std::string profile=vortexProfilePath();
        GERGResidualTable gerg;
        const std::string gergPath=std::string(RFAM_SOURCE_DIR)
            +"/tests/reference/gerg_residual_table_v3.csv";
        if(!gerg.load(gergPath))
            throw std::runtime_error("cannot load GERG table: "+gergPath);

        LoadedSector fineBase(
            fineBasePath,profile,registry,gerg);
        LoadedSector coarseBase(
            coarseBasePath,profile,registry,gerg);
        LoadedSector coarseImproved(
            coarseImprovedPath,profile,registry,gerg);
        if(fineBase.metadata.profileHash!=coarseBase.metadata.profileHash
            ||fineBase.metadata.profileHash
                !=coarseImproved.metadata.profileHash)
            throw std::runtime_error("coarse correction profile hash mismatch");
        if(!allowedModeTransition(
                fineBase.metadata.mode,coarseImproved.metadata.mode)
            ||!allowedModeTransition(
                coarseBase.metadata.mode,coarseImproved.metadata.mode))
            throw std::runtime_error("coarse correction mode mismatch");

        Domain restrictedCheck(
            coarseBase.domain.Nx(),coarseBase.domain.Ny(),
            coarseBase.domain.Nz(),coarseBase.mesh.length(),0.065,
            coarseBase.mesh.sectorAngle(),coarseBase.domain.nSpecies());
        restrictedCheck.setSpecies(coarseBase.domain.speciesNames());
        restrictedCheck.thermoModel().setGERGTable(&gerg);
        SectorConservativeGridTransfer::apply(
            fineBase.domain,fineBase.mesh,nullptr,
            restrictedCheck,coarseBase.mesh,nullptr);
        const double baseRestrictionError=maximumFieldDifference(
            restrictedCheck,coarseBase.domain);
        if(baseRestrictionError>1.0e-9)
            throw std::runtime_error(
                "coarse base is not the aligned restriction of fine base");

        if((freezeFineComposition||conservativeParentComposition)
            &&!requireFineResidualDecrease)
            throw std::runtime_error(
                "fine composition modes require the fine residual safeguard");
        if(requireFineResidualDecrease
            &&(!std::isfinite(residualReferenceStep)
               ||residualReferenceStep<=0.0))
            throw std::runtime_error(
                "fine residual reference step must be positive");
        if(requireFineResidualDecrease)
        {
            const std::string expectedFineMode="finite-rate-les";
            const std::string expectedFrozenMode="frozen-chemistry-les";
            if((fineBase.metadata.mode!=expectedFineMode
                &&fineBase.metadata.mode!=expectedFrozenMode)
                ||coarseImproved.metadata.mode!=expectedFrozenMode)
                throw std::runtime_error(
                    "fine residual safeguard requires frozen LES correction");
            static CryogenicTransportTable transport;
            const std::string transportPath=std::string(RFAM_SOURCE_DIR)
                +"/tests/reference/cryogenic_transport_runtime_coolprop72_v1.csv";
            if(!transport.load(transportPath))
                throw std::runtime_error(
                    "cannot load cryogenic transport table: "+transportPath);
            TransportProps::registerCryogenicTable(&transport);
        }

        double baseResidual=std::numeric_limits<double>::quiet_NaN();
        double candidateResidual=std::numeric_limits<double>::quiet_NaN();
        double bestCandidateResidual=
            std::numeric_limits<double>::infinity();
        double bestCandidateMaximumRelativeCorrection=0.0;
        double baseLowFrequencyResidual=
            std::numeric_limits<double>::quiet_NaN();
        double candidateLowFrequencyResidual=
            std::numeric_limits<double>::quiet_NaN();
        double bestCandidateLowFrequencyResidual=
            std::numeric_limits<double>::infinity();
        int residualAttempts=0;
        if(requireFineResidualDecrease)
        {
            SectorFrozenSteadyResidualContext baseContext(
                fineBase.domain,fineBase.mesh,registry,&gerg,
                fineBase.turbulence.get(),fineBase.metadata.time);
            SteadyResidualEvaluator::Result residual;
            if(!baseContext.evaluate(residualReferenceStep,residual))
                throw std::runtime_error(
                    "fine base residual evaluation failed at cell "
                    +std::to_string(fineBase.domain.failure().cell)+": "
                    +fineBase.domain.failure().detail);
            baseResidual=SectorFrozenSteadyResidualContext::normalizedRms(
                residual,fineBase.domain);
            const SteadyResidualEvaluator::Result restrictedResidual=
                SectorFASDefect::restrictResidual(
                    fineBase.domain,fineBase.mesh,residual,
                    coarseBase.domain,coarseBase.mesh);
            baseLowFrequencyResidual=
                SectorFrozenSteadyResidualContext::normalizedRms(
                    restrictedResidual,coarseBase.domain);
        }

        std::unique_ptr<Domain> corrected;
        SectorCoarseCorrection::Result result;
        double acceptedMaximumRelativeCorrection=maximumRelativeCorrection;
        bool accepted=false;
        const int attempts=requireFineResidualDecrease
            ?residualBacktracks+1:1;
        for(int attempt=0;attempt<attempts;++attempt)
        {
            const double trialMaximumRelativeCorrection=
                std::ldexp(maximumRelativeCorrection,-attempt);
            auto trial=std::make_unique<Domain>(
                fineBase.domain.Nx(),fineBase.domain.Ny(),
                fineBase.domain.Nz(),fineBase.mesh.length(),0.065,
                fineBase.mesh.sectorAngle(),fineBase.domain.nSpecies());
            trial->setSpecies(fineBase.domain.speciesNames());
            trial->thermoModel().setGERGTable(&gerg);
            if(!trial->replaceConservedField(
                    captureField(fineBase.domain),
                    captureTemperature(fineBase.domain)))
                throw std::runtime_error(
                    "fine corrected target clone failed");

            const SectorCoarseCorrection::Result trialResult=
                SectorCoarseCorrection::apply(
                    fineBase.domain,fineBase.mesh,
                    coarseBase.domain,coarseImproved.domain,coarseBase.mesh,
                    *trial,smoothingSweeps,24,
                    trialMaximumRelativeCorrection,
                    excludeBoundaryParentLayer,
                    reverseCorrection?-1.0:1.0,
                    freezeFineComposition,
                    conservativeParentComposition);
            if(requireFineResidualDecrease)
            {
                ++residualAttempts;
                SectorFrozenSteadyResidualContext trialContext(
                    *trial,fineBase.mesh,registry,&gerg,
                    fineBase.turbulence.get(),fineBase.metadata.time);
                SteadyResidualEvaluator::Result residual;
                if(!trialContext.evaluate(residualReferenceStep,residual))
                    continue;
                candidateResidual=
                    SectorFrozenSteadyResidualContext::normalizedRms(
                        residual,fineBase.domain);
                const SteadyResidualEvaluator::Result restrictedResidual=
                    SectorFASDefect::restrictResidual(
                        *trial,fineBase.mesh,residual,
                        coarseBase.domain,coarseBase.mesh);
                candidateLowFrequencyResidual=
                    SectorFrozenSteadyResidualContext::normalizedRms(
                        restrictedResidual,coarseBase.domain);
                if(candidateResidual<bestCandidateResidual)
                {
                    bestCandidateResidual=candidateResidual;
                    bestCandidateLowFrequencyResidual=
                        candidateLowFrequencyResidual;
                    bestCandidateMaximumRelativeCorrection=
                        trialMaximumRelativeCorrection;
                }
                if(!(candidateResidual<baseResidual*(1.0-1.0e-8)))
                    continue;
            }
            corrected=std::move(trial);
            result=trialResult;
            acceptedMaximumRelativeCorrection=
                trialMaximumRelativeCorrection;
            accepted=true;
            break;
        }
        if(!accepted)
            throw std::runtime_error(
                "fine residual safeguard rejected every coarse correction "
                "scale: base="+std::to_string(baseResidual)
                +" best="+std::to_string(bestCandidateResidual)
                +" bestMaxCorrection="
                +std::to_string(bestCandidateMaximumRelativeCorrection)
                +" lowFrequencyBase="
                +std::to_string(baseLowFrequencyResidual)
                +" lowFrequencyBest="
                +std::to_string(bestCandidateLowFrequencyResidual)
                +" attempts="+std::to_string(residualAttempts));

        const double massBefore=SectorDiagnostics(fineBase.mesh)
            .evaluate(fineBase.domain).massInventory;
        const double massAfter=SectorDiagnostics(fineBase.mesh)
            .evaluate(*corrected).massInventory;

        std::unique_ptr<SectorTurbulenceOperator> correctedTurbulence;
        if(fineBase.turbulence)
        {
            correctedTurbulence=std::make_unique<SectorTurbulenceOperator>(
                fineBase.mesh,registry,0.05,0.001,0.25,false,
                SectorTurbulenceOperator::TurbulenceMode::WALE);
            correctedTurbulence->restoreFields(
                fineBase.turbulence->kField(),
                fineBase.turbulence->omegaField());
        }
        int threads=1;
#ifdef _OPENMP
        threads=omp_get_max_threads();
#endif
        SectorCheckpointMetadata metadata=coarseImproved.metadata;
        metadata.time=fineBase.metadata.time;
        metadata.iteration=coarseImproved.metadata.iteration;
        metadata.threadCount=threads;
#ifdef RFAM_GIT_COMMIT
        metadata.gitCommit=RFAM_GIT_COMMIT;
#endif
        SectorCheckpoint::write(
            outputPath,*corrected,metadata,correctedTurbulence.get());

        const double relativeMassChange=std::abs(massAfter-massBefore)
            /std::max(std::abs(massBefore),1.0e-300);
        std::cout << "fine/coarse grids="
                  << fineBase.domain.Nx()-2 << "x"
                  << fineBase.domain.Ny()-2 << "x"
                  << fineBase.domain.Nz()-2 << "/"
                  << coarseBase.domain.Nx()-2 << "x"
                  << coarseBase.domain.Ny()-2 << "x"
                  << coarseBase.domain.Nz()-2 << "\n"
                  << "base_restriction_error=" << baseRestrictionError << "\n"
                  << "max_relative_correction/requested/accepted/unscaled/trust_scale="
                  << maximumRelativeCorrection << "/"
                  << acceptedMaximumRelativeCorrection << "/"
                  << result.maximumUnscaledRelativeCorrection << "/"
                  << result.trustRegionScale << "\n"
                  << "fine_residual_required/base/candidate/best/best_max_correction/attempts="
                  << (requireFineResidualDecrease?1:0) << "/"
                  << baseResidual << "/" << candidateResidual << "/"
                  << bestCandidateResidual << "/"
                  << bestCandidateMaximumRelativeCorrection << "/"
                  << residualAttempts << "\n"
                  << "low_frequency_residual_base/candidate/best="
                  << baseLowFrequencyResidual << "/"
                  << candidateLowFrequencyResidual << "/"
                  << bestCandidateLowFrequencyResidual << "\n"
                  << "correction_direction="
                  << result.correctionDirection << "\n"
                  << "exclude_boundary_parent_layer/masked_fine_cells="
                  << (excludeBoundaryParentLayer?1:0) << "/"
                  << result.boundaryMaskedFineCells << "\n"
                  << "global_scale/local_limited/recovery_backtracks/max_local_reduction="
                  << result.acceptedGlobalScale << "/"
                  << result.localLimitedCells << "/"
                  << result.recoveryBacktracks << "/"
                  << result.maximumLocalReduction << "\n"
                  << "composition_frozen/reconstructed/projected/max_projection="
                  << result.compositionFrozenCells << "/"
                  << result.compositionReconstructedCells << "/"
                  << result.compositionProjectedCells << "/"
                  << result.maximumCompositionProjection << "\n"
                  << "parent_composition_projected/max_projection="
                  << result.compositionProjectedParents << "/"
                  << result.maximumParentCompositionProjection << "\n"
                  << "inventory_defect/flow/species="
                  << result.maximumInventoryDefectRelativeToBase << "/"
                  << result.maximumFlowInventoryDefectRelativeToBase << "/"
                  << result.maximumSpeciesInventoryDefectRelativeToBase << "\n"
                  << "mass_before/after/relative_change="
                  << massBefore << "/" << massAfter << "/"
                  << relativeMassChange << "\n"
                  << "first_recovery_failure_cell/status/backend="
                  << result.firstRecoveryFailedCell << "/"
                  << result.firstRecoveryEquilibriumStatus << "/"
                  << result.firstRecoveryBackendStatus << "\n"
                  << "first_failure_rho_base/candidate/drho="
                  << result.firstFailureBaseDensity << "/"
                  << result.firstFailureCandidateDensity << "/"
                  << result.firstFailureDensityCorrection << "\n"
                  << "first_failure_e_base/candidate/drhoE="
                  << result.firstFailureBaseSpecificInternalEnergy << "/"
                  << result.firstFailureCandidateSpecificInternalEnergy << "/"
                  << result.firstFailureTotalEnergyCorrection << "\n"
                  << "output=" << outputPath << " mode=" << metadata.mode
                  << " git=" << metadata.gitCommit << "\n"
                  << "SECTOR_COARSE_CORRECTION_APPLY_PASS\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << "SECTOR_COARSE_CORRECTION_APPLY_FAIL: "
                  << error.what() << "\n";
        return 1;
    }
}
