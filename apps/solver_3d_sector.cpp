/*---------------------------------------------------------------------------*\
    solver_3d_sector - Vortex mapped 3D 120-degree sector driver.

    Supports isothermal or finite-rate reacting flow, optional k-omega SST
    RANS, mapped cylindrical metrics, rotational periodicity, and calibrated
    torch/LOX/C3H8 injection boundary conditions.

    Grid levels: smoke 14,400; pathfinder 72,000; intermediate 432,000 cells.
    Coordinates: x=axial, y=radial (0 to r_wall), z=azimuthal (0 to 2*pi/3).
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <memory>
#include <chrono>
#include <cstring>
#include <exception>
#include <fstream>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "src/variables/ConservedVars.H"
#include "src/framework/Domain.H"
#include "src/framework/FluxScheme.H"
#include "src/framework/Operator.H"
#include "src/framework/SolverBase.H"
#include "src/framework/NonReflectingOutlet.H"
#include "src/framework/IsothermalProjection.H"
#include "src/framework/ModelCapabilities.H"
#include "src/framework/RunTelemetry.H"
#include "src/transport/CryogenicTransportTable.H"
#include "src/transport/TransportProps.H"
#include "src/framework/SectorConvectionOperator.H"
#include "src/framework/SectorDiagnostics.H"
#include "src/framework/SectorInjectionBoundary.H"
#include "src/thermo/GERGResidualTable.H"
#include "src/thermo/GERGSaturationTable.H"
#include "src/hotstart/HotStartFieldLoader.H"
#include "src/framework/SectorNonReflectingOutlet.H"
#include "src/framework/SectorDiffusionOperator.H"
#include "src/framework/SectorPerformanceDiagnostics.H"
#include "src/framework/SectorWallBoundary.H"
#include "src/framework/ProfiledOperator.H"
#include "src/framework/DualTimeIntegrator.H"
#include "src/framework/SteadyLocalTimeIntegrator.H"
#include "src/framework/SteadyLUSGSIntegrator.H"
#include "src/framework/FixedConservedSourceOperator.H"
#include "src/framework/SectorCheckpoint.H"
#include "src/framework/SectorResidualField.H"
#include "src/framework/SectorPressureSchurPreconditioner.H"
#include "src/framework/SectorStructuredPressureSchurPreconditioner.H"
#include "src/framework/SectorOutputWriter.H"
#include "src/mesh/SectorMeshMetrics.H"
#include "src/turbulence/SectorTurbulenceOperator.H"
#include "src/thermo/SpeciesRegistry.H"
#ifdef RFAM_HAVE_CHEMISTRY
#include "src/combustion/SectorChemistryOperator.H"
#include "src/combustion/FiniteRate/TwoStepCO.H"
#include "src/combustion/SandiegoSectorChemistry.H"
#endif

using namespace Foam;
using namespace std;

int main(int argc, char** argv)
{
    cout << "=== RFAM Vortex mapped 3D 120-degree sector ===\n\n";

    // Telemetry covers the structured path too. It is wired at phase
    // granularity here: the sector driver's step bookkeeping lives inside the
    // integrator objects, so per-step metrics are a follow-up, but the
    // manifest, heartbeat, stall detection and degradation ledger are complete.
    Foam::telemetry::Options telemetryOptions;
    for (int i = 1; i + 1 < argc; ++i)
    {
        const std::string option = argv[i];
        if (option == "--no-telemetry") telemetryOptions.enabled = false;
        else if (option == "--telemetry-dir") telemetryOptions.rootDir = argv[++i];
        else if (option == "--telemetry-stall-seconds")
            telemetryOptions.stallSeconds = std::atoi(argv[++i]);
        else if (option == "--strict-no-downgrade")
            telemetryOptions.strictNoDowngrade = true;
    }
    Foam::telemetry::Session& telemetry = Foam::telemetry::Session::instance();
    telemetry.start(argc, argv, telemetryOptions);
    RFAM_PHASE(sector_startup);
    RFAM_PHASE(sector_mesh);
    RFAM_PHASE(sector_initialize);
    RFAM_PHASE(sector_advance);
    RFAM_PHASE(sector_output);
    RFAM_PHASE(sector_finish);
    RFAM_ENTER(sector_startup);

    // Mapped cylindrical grid: axial, radial, and azimuthal coordinates.
    enum class GridLevel { Coarse, Smoke, Pathfinder, Intermediate };
    GridLevel gridLevel=GridLevel::Smoke;
    bool useRans = false;
    bool useLes = false;
    bool reacting = false;
    bool freezeChemistry = false;
    bool cflAcousticScale = false;
    bool dualTime = false;
    bool steadyLocalTime = false;
    bool steadyLusgs = false;
    bool lusgsAnderson = false;
    bool lusgsFrozenEulerBlock = false;
    bool lusgsHydroDiagonalBlock = false;
    bool lusgsPressureSchurFallback = false;
    bool lusgsPressureSchurLeastSquaresGradient = false;
    bool lusgsPressureSchurNonOrthogonalCorrection = false;
    bool lusgsPressureSchurStructuredMomentum = false;
    bool lusgsSegregatedAcousticFallback = false;
    bool lusgsSegregatedFlowFallback = false;
    bool lusgsSegregatedEnergyFallback = false;
    bool lusgsResidualLineSearch = false;
    bool lusgsBidirectionalLineSearch = false;
    int lusgsAndersonDepth = 4;
    double dtMult = 20.0;
    double dualInnerCfl = 5.0;
    double dualResidualTolerance = 0.05;
    bool dualLocalPseudoTime = false;
    double dualLocalPseudoRatio = 100.0;
    bool dualAitken = false;
    bool dualAnderson = false;
    int dualAndersonDepth = 4;
    bool dualInnerInjectionRefresh = false;
    int dualMaxInnerIters = 80;
    int dualMaxRetries = 4;
    double steadyResidualTolerance = 0.05;
    double steadyLocalRatio = 100.0;
    double lusgsPseudoCfl = 5.0;
    double lusgsMaximumUpdate = 0.05;
    double lusgsEnergyScale = 1.0;
    double steadyOutletRelaxation = 0.0;
    double steadyRampOverride = -1.0;   // >=0: constant injection ramp (steady)
    int steadyRampSteps = 0;            // >0: ramp over N steady pseudo steps
    double icPressure = 4e6, icTemperature = 700.0;  // initial field state
    double torchTubeMm = 0.0;        // >0: submerged torch tube (torch at x=L)
    double injectorFaceMm = 0.0;    // >0: injected plane at x (CAD 9.9mm LOX/C3H8 face)
    double sectorDegrees = 120.0;   // 120 sector or 360 full annulus
    bool allowSinglePhaseReactingDiagnostic = false;
    bool blueprintColdFlow = false;
    bool lusgsNoSpeciesLimiter = false;
    bool allowFrozenToFiniteRate = false;
    bool chemistryExplicit = false;
    bool sandiegoMechanism = false;
    int maxSteps = 50;
    double cfl=0.3;
    int checkpointEvery=0;
    int monitorEvery=0;
    std::string checkpointPath="rfam_sector.chk",restartPath;
    std::string fasSourcePath;
    std::string hotstartFieldPath;
    std::string exportAfterImportPath;
    std::string fieldCsvPath,wallCsvPath;
    for (int a = 1; a < argc; ++a)
    {
        if (std::strcmp(argv[a], "--pathfinder") == 0)
            gridLevel=GridLevel::Pathfinder;
        else if (std::strcmp(argv[a], "--intermediate") == 0)
            gridLevel=GridLevel::Intermediate;
        else if (std::strcmp(argv[a], "--grid") == 0 && a + 1 < argc)
        {
            const std::string value=argv[++a];
            if(value=="coarse") gridLevel=GridLevel::Coarse;
            else if(value=="smoke") gridLevel=GridLevel::Smoke;
            else if(value=="pathfinder") gridLevel=GridLevel::Pathfinder;
            else if(value=="intermediate") gridLevel=GridLevel::Intermediate;
            else
            {
                cerr << "unknown grid level: " << value
                     << " (expected coarse, smoke, pathfinder, or intermediate)\n";
                return 2;
            }
        }
        else if (std::strcmp(argv[a], "--rans") == 0) useRans = true;
        else if (std::strcmp(argv[a], "--les") == 0) useLes = true;
        else if (std::strcmp(argv[a], "--reacting") == 0) reacting = true;
        else if (std::strcmp(argv[a], "--freeze-chemistry") == 0)
            freezeChemistry = true;
        else if (std::strcmp(
            argv[a], "--allow-single-phase-reacting-diagnostic") == 0)
            allowSinglePhaseReactingDiagnostic = true;
        else if (std::strcmp(argv[a], "--chemistry-explicit") == 0)
            chemistryExplicit = true;
        else if (std::strcmp(argv[a], "--sandiego") == 0)
            sandiegoMechanism = true;
        else if (std::strcmp(argv[a], "--steps") == 0 && a + 1 < argc)
            maxSteps = std::max(0, std::atoi(argv[++a]));
        else if(std::strcmp(argv[a],"--cfl")==0&&a+1<argc)
            cfl=std::clamp(std::atof(argv[++a]),0.01,0.9);
        else if(std::strcmp(argv[a],"--dual-time")==0) dualTime=true;
        else if(std::strcmp(argv[a],"--steady-local-time")==0)
            steadyLocalTime=true;
        else if(std::strcmp(argv[a],"--steady-lusgs")==0)
            steadyLusgs=true;
        else if(std::strcmp(argv[a],"--lusgs-anderson")==0)
            lusgsAnderson=true;
        else if(std::strcmp(argv[a],"--lusgs-frozen-euler-block")==0)
            lusgsFrozenEulerBlock=true;
        else if(std::strcmp(argv[a],"--lusgs-hydro-diagonal-block")==0)
            lusgsHydroDiagonalBlock=true;
        else if(std::strcmp(argv[a],"--lusgs-pressure-schur-fallback")==0)
            lusgsPressureSchurFallback=true;
        else if(std::strcmp(
            argv[a],"--lusgs-pressure-schur-least-squares-gradient")==0)
            lusgsPressureSchurLeastSquaresGradient=true;
        else if(std::strcmp(
            argv[a],"--lusgs-pressure-schur-nonorthogonal-correction")==0)
            lusgsPressureSchurNonOrthogonalCorrection=true;
        else if(std::strcmp(
            argv[a],"--lusgs-pressure-schur-structured-momentum")==0)
            lusgsPressureSchurStructuredMomentum=true;
        else if(std::strcmp(argv[a],"--lusgs-segregated-acoustic-fallback")==0)
            lusgsSegregatedAcousticFallback=true;
        else if(std::strcmp(argv[a],"--lusgs-segregated-flow-fallback")==0)
            lusgsSegregatedFlowFallback=true;
        else if(std::strcmp(argv[a],"--lusgs-segregated-energy-fallback")==0)
            lusgsSegregatedEnergyFallback=true;
        else if(std::strcmp(argv[a],"--lusgs-residual-line-search")==0)
            lusgsResidualLineSearch=true;
        else if(std::strcmp(
            argv[a],"--lusgs-bidirectional-line-search")==0)
            lusgsBidirectionalLineSearch=true;
        else if(std::strcmp(argv[a],"--lusgs-anderson-depth")==0&&a+1<argc)
            lusgsAndersonDepth=std::clamp(std::atoi(argv[++a]),1,8);
        else if(std::strcmp(argv[a],"--dt-mult")==0&&a+1<argc)
            dtMult=std::clamp(std::atof(argv[++a]),1.0,1000.0);
        else if(std::strcmp(argv[a],"--dual-inner-cfl")==0&&a+1<argc)
            dualInnerCfl=std::clamp(std::atof(argv[++a]),0.001,100.0);
        else if(std::strcmp(argv[a],"--dual-res-tol")==0&&a+1<argc)
            dualResidualTolerance=std::clamp(
                std::atof(argv[++a]),1.0e-8,0.5);
        else if(std::strcmp(argv[a],"--dual-local-time")==0)
            dualLocalPseudoTime=true;
        else if(std::strcmp(argv[a],"--dual-local-ratio")==0&&a+1<argc)
            dualLocalPseudoRatio=std::clamp(
                std::atof(argv[++a]),1.0,10000.0);
        else if(std::strcmp(argv[a],"--dual-aitken")==0)
            dualAitken=true;
        else if(std::strcmp(argv[a],"--dual-anderson")==0)
            dualAnderson=true;
        else if(std::strcmp(argv[a],"--dual-anderson-depth")==0&&a+1<argc)
            dualAndersonDepth=std::clamp(std::atoi(argv[++a]),1,8);
        else if(std::strcmp(argv[a],"--dual-inner-injection-refresh")==0)
            dualInnerInjectionRefresh=true;
        else if(std::strcmp(argv[a],"--dual-iters")==0&&a+1<argc)
            dualMaxInnerIters=std::max(1,std::atoi(argv[++a]));
        else if(std::strcmp(argv[a],"--dual-retries")==0&&a+1<argc)
            dualMaxRetries=std::clamp(std::atoi(argv[++a]),0,12);
        else if(std::strcmp(argv[a],"--steady-res-tol")==0&&a+1<argc)
            steadyResidualTolerance=std::clamp(
                std::atof(argv[++a]),1.0e-8,0.5);
        else if(std::strcmp(argv[a],"--steady-local-ratio")==0&&a+1<argc)
            steadyLocalRatio=std::clamp(
                std::atof(argv[++a]),1.0,10000.0);
        else if(std::strcmp(argv[a],"--lusgs-cfl")==0&&a+1<argc)
            lusgsPseudoCfl=std::clamp(
                std::atof(argv[++a]),0.01,1000.0);
        else if(std::strcmp(argv[a],"--lusgs-max-update")==0&&a+1<argc)
            lusgsMaximumUpdate=std::clamp(
                std::atof(argv[++a]),1.0e-4,0.5);
        else if(std::strcmp(argv[a],"--lusgs-energy-scale")==0&&a+1<argc)
            lusgsEnergyScale=std::clamp(
                std::atof(argv[++a]),0.01,1.0);
        else if(std::strcmp(argv[a],"--steady-outlet-relax")==0&&a+1<argc)
            steadyOutletRelaxation=std::clamp(
                std::atof(argv[++a]),0.0,1.0);
        else if(std::strcmp(argv[a],"--steady-ramp")==0&&a+1<argc)
            steadyRampOverride=std::atof(argv[++a]);
        else if(std::strcmp(argv[a],"--steady-ramp-steps")==0&&a+1<argc)
            steadyRampSteps=std::max(1,std::atoi(argv[++a]));
        else if(std::strcmp(argv[a],"--ic-pressure")==0&&a+1<argc)
            icPressure=std::atof(argv[++a]);
        else if(std::strcmp(argv[a],"--ic-temperature")==0&&a+1<argc)
            icTemperature=std::atof(argv[++a]);
        else if(std::strcmp(argv[a],"--blueprint-cold-flow")==0)
            blueprintColdFlow=true;
        else if(std::strcmp(argv[a],"--torch-tube-mm")==0&&a+1<argc)
            torchTubeMm=std::max(0.0,std::atof(argv[++a]));
        else if(std::strcmp(argv[a],"--injector-face-mm")==0&&a+1<argc)
            injectorFaceMm=std::max(0.0,std::atof(argv[++a]));
        else if(std::strcmp(argv[a],"--sector-degrees")==0&&a+1<argc)
            sectorDegrees=std::atof(argv[++a]);
        else if(std::strcmp(argv[a],"--lusgs-no-species-limiter")==0)
            lusgsNoSpeciesLimiter=true;
        else if(std::strcmp(argv[a],"--allow-frozen-to-finite-rate")==0)
            allowFrozenToFiniteRate=true;
        else if(std::strcmp(argv[a],"--checkpoint-every")==0&&a+1<argc)
            checkpointEvery=std::max(1,std::atoi(argv[++a]));
        else if(std::strcmp(argv[a],"--checkpoint")==0&&a+1<argc)
            checkpointPath=argv[++a];
        else if(std::strcmp(argv[a],"--monitor-every")==0&&a+1<argc)
            monitorEvery=std::max(1,std::atoi(argv[++a]));
        else if(std::strcmp(argv[a],"--restart")==0&&a+1<argc)
            restartPath=argv[++a];
        else if(std::strcmp(argv[a],"--fas-source")==0&&a+1<argc)
            fasSourcePath=argv[++a];
        else if(std::strcmp(argv[a],"--hotstart-field")==0&&a+1<argc)
            hotstartFieldPath=argv[++a];
        else if(std::strcmp(argv[a],"--export-after-import")==0&&a+1<argc)
            exportAfterImportPath=argv[++a];
        else if(std::strcmp(argv[a],"--field-csv")==0&&a+1<argc)
            fieldCsvPath=argv[++a];
        else if(std::strcmp(argv[a],"--wall-csv")==0&&a+1<argc)
            wallCsvPath=argv[++a];
    }

    const int selectedImplicitModes=(dualTime?1:0)
        +(steadyLocalTime?1:0)+(steadyLusgs?1:0);
    if(selectedImplicitModes>1)
    {
        cerr << "TIME_INTEGRATOR_CONFLICT: choose exactly one of --dual-time, "
             << "--steady-local-time, or --steady-lusgs\n";
        return 2;
    }
    if((steadyLocalTime||steadyLusgs)&&reacting&&!freezeChemistry)
    {
        cerr << "STEADY_REACTING_REQUIRES_FROZEN_CHEMISTRY: a steady solver "
             << "cannot advance split finite-rate chemistry\n";
        return 2;
    }
    if(lusgsAnderson&&!steadyLusgs)
    {
        cerr << "LUSGS_ANDERSON_REQUIRES_LUSGS: use --steady-lusgs with "
             << "--lusgs-anderson\n";
        return 2;
    }
    if(lusgsFrozenEulerBlock&&!steadyLusgs)
    {
        cerr << "LUSGS_BLOCK_REQUIRES_LUSGS: use --steady-lusgs with "
             << "--lusgs-frozen-euler-block\n";
        return 2;
    }
    if(lusgsHydroDiagonalBlock&&!steadyLusgs)
    {
        cerr << "LUSGS_HYDRO_BLOCK_REQUIRES_LUSGS: use --steady-lusgs with "
             << "--lusgs-hydro-diagonal-block\n";
        return 2;
    }
    if(lusgsPressureSchurLeastSquaresGradient
        &&!lusgsPressureSchurFallback)
    {
        cerr << "LUSGS_PRESSURE_SCHUR_LS_REQUIRES_PRESSURE_SCHUR: use "
             << "--lusgs-pressure-schur-fallback with the least-squares "
             << "gradient option\n";
        return 2;
    }
    if(lusgsPressureSchurNonOrthogonalCorrection
        &&(!lusgsPressureSchurFallback
            ||!lusgsPressureSchurLeastSquaresGradient))
    {
        cerr << "LUSGS_PRESSURE_SCHUR_NONORTH_REQUIRES_LS: use pressure "
             << "Schur fallback and least-squares gradient with the "
             << "nonorthogonal correction option\n";
        return 2;
    }
    if(lusgsPressureSchurStructuredMomentum
        &&!lusgsPressureSchurFallback)
    {
        cerr << "LUSGS_PRESSURE_SCHUR_STRUCTURED_REQUIRES_PRESSURE_SCHUR: "
             << "use --lusgs-pressure-schur-fallback with the structured "
             << "momentum option\n";
        return 2;
    }
    if(lusgsPressureSchurStructuredMomentum
        &&(lusgsPressureSchurLeastSquaresGradient
            ||lusgsPressureSchurNonOrthogonalCorrection))
    {
        cerr << "LUSGS_PRESSURE_SCHUR_VARIANT_CONFLICT: structured momentum "
             << "and local LS/nonorthogonal Schur variants are mutually "
             << "exclusive\n";
        return 2;
    }
    if((lusgsPressureSchurFallback||lusgsSegregatedAcousticFallback
        ||lusgsSegregatedFlowFallback||lusgsSegregatedEnergyFallback)
        &&(!steadyLusgs||!lusgsResidualLineSearch
           ||!lusgsBidirectionalLineSearch))
    {
        cerr << "LUSGS_SEGREGATED_FALLBACK_REQUIRES_BIDIRECTIONAL_LINE_SEARCH: "
             << "use --steady-lusgs --lusgs-residual-line-search "
             << "--lusgs-bidirectional-line-search with a segregated "
             << "fallback option\n";
        return 2;
    }
    if(lusgsResidualLineSearch&&!steadyLusgs)
    {
        cerr << "LUSGS_LINE_SEARCH_REQUIRES_LUSGS: use --steady-lusgs with "
             << "--lusgs-residual-line-search\n";
        return 2;
    }
    if(lusgsBidirectionalLineSearch
        &&(!steadyLusgs||!lusgsResidualLineSearch))
    {
        cerr << "LUSGS_BIDIRECTIONAL_REQUIRES_LINE_SEARCH: use "
             << "--steady-lusgs --lusgs-residual-line-search with "
             << "--lusgs-bidirectional-line-search\n";
        return 2;
    }
    if(!fasSourcePath.empty()
        &&(!steadyLusgs||gridLevel!=GridLevel::Coarse
           ||!reacting||!freezeChemistry||!useLes
           ||restartPath.empty()))
    {
        cerr << "FAS_SOURCE_MODE_MISMATCH: --fas-source requires --grid "
             << "coarse --steady-lusgs --reacting --freeze-chemistry --les "
             << "and --restart\n";
        return 2;
    }
    if (freezeChemistry && !reacting)
    {
        cerr << "FREEZE_CHEMISTRY_REQUIRES_REACTING: use --reacting "
             << "--freeze-chemistry to preserve the 9-species state\n";
        return 2;
    }
    if (freezeChemistry && (sandiegoMechanism || chemistryExplicit))
    {
        cerr << "FREEZE_CHEMISTRY_CONFLICT: frozen chemistry cannot select "
             << "a chemistry mechanism or explicit chemistry integrator\n";
        return 2;
    }
    if (reacting && maxSteps > 0 &&
        !ModelCapabilities::reactingSectorProductionReady() &&
        !allowSinglePhaseReactingDiagnostic)
    {
        cerr << "REACTING_ADVANCE_BLOCKED: phase equilibrium, latent heat, "
             << "phase-aware energy/sound speed/boundary states/characteristics/"
             << "transport/chemistry are not ready\n";
        return 3;
    }
    if (reacting && maxSteps > 0 && allowSinglePhaseReactingDiagnostic)
        cerr << (freezeChemistry
            ? "WARNING: frozen-chemistry sector is diagnostic only; chemistry "
              "sources are disabled and results are not valid for LOX design use\n"
            : "WARNING: single-phase reacting sector is diagnostic only; "
              "results are not valid for LOX design use\n");
    if (maxSteps > 0 &&
        !allowSinglePhaseReactingDiagnostic &&
        !ModelCapabilities::nonReactingSectorProductionReady())
    {
        cerr << "SECTOR_ADVANCE_BLOCKED: phase equilibrium, latent heat, "
             << "phase-aware energy/sound speed/boundary states/characteristics "
             << "and dense cryogenic transport "
             << "are not ready for the configured 130 K LOX boundary\n";
        return 4;
    }
    int activeNx=60,activeNy=20,activeNz=12;
    const char* gridName="smoke";
    if(gridLevel==GridLevel::Coarse)
    { activeNx=30;activeNy=10;activeNz=6;gridName="coarse"; }
    else if(gridLevel==GridLevel::Pathfinder)
    { activeNx=120;activeNy=30;activeNz=20;gridName="pathfinder"; }
    else if(gridLevel==GridLevel::Intermediate)
    { activeNx=240;activeNy=50;activeNz=36;gridName="intermediate"; }
    // Full-annulus runs scale the circumferential count so the angular
    // resolution matches the 120-degree convention (e.g. 36 -> 108).
    if(sectorDegrees>180.0)
        activeNz=static_cast<int>(std::lround(
            activeNz*sectorDegrees/120.0));
    const int Nx=activeNx+2,Ny=activeNy+2,Nz=activeNz+2;
    if (telemetry.active())
    {
        auto flagOf=[&](const char* key,bool value)
        {telemetry.defineField(key,value?"true":"false",nullptr,"resolved");};
        telemetry.defineField("grid",gridName,nullptr,"resolved");
        telemetry.defineField("nx",std::to_string(activeNx),nullptr,"resolved");
        telemetry.defineField("ny",std::to_string(activeNy),nullptr,"resolved");
        telemetry.defineField("nz",std::to_string(activeNz),nullptr,"resolved");
        telemetry.defineField("sector_degrees",std::to_string(sectorDegrees),nullptr,"resolved");
        telemetry.defineField("max_steps",std::to_string(maxSteps),nullptr,"resolved");
        flagOf("rans",useRans);
        flagOf("les",useLes);
        flagOf("reacting",reacting);
        flagOf("freeze_chemistry",freezeChemistry);
        flagOf("dual_time",dualTime);
        flagOf("steady_local_time",steadyLocalTime);
        flagOf("steady_lusgs",steadyLusgs);
        flagOf("cfl_acoustic_scale",cflAcousticScale);
        flagOf("single_phase_reacting_diagnostic",allowSinglePhaseReactingDiagnostic);
        telemetry.flushManifest();
        Foam::json::Value& manifest=telemetry.describe();
        manifest.set("thermo_model_id",
                     Foam::json::Value(std::string("sector-") +
                         (freezeChemistry ? "frozen-chemistry" : "live")));
        Foam::json::Value capabilities=Foam::json::Value::object();
        capabilities.set("nonReactingSectorProductionReady",
            Foam::json::Value(Foam::ModelCapabilities::nonReactingSectorProductionReady()));
        capabilities.set("reactingSectorProductionReady",
            Foam::json::Value(Foam::ModelCapabilities::reactingSectorProductionReady()));
        capabilities.set("ransTetraProductionReady",
            Foam::json::Value(Foam::ModelCapabilities::ransTetraProductionReady()));
        manifest.set("capabilities",capabilities);
        telemetry.flushManifest();
    }
    RFAM_LEAVE(sector_startup);
    RFAM_ENTER(sector_mesh);
    const std::string profilePath=vortexProfilePath();
    SectorMeshMetrics mesh(profilePath,
                           Nx,Ny,Nz,sectorDegrees);
    const double Lx=mesh.length(), Ly=0.065, Lz=mesh.sectorAngle();
    const double p0=4e6, T0=700.0;      // injection reference state
    auto sumY=[&](const std::vector<double>& Y)->double
    {
        double s=0.0;
        for(double y:Y) s+=y;
        return s;
    };

    SpeciesRegistry registry;
    const std::vector<std::string> species=reacting?registry.names()
        :std::vector<std::string>{"O2","C3H8"};
    std::vector<double> background(species.size(),0.0);
    if(reacting)
    {
        background[registry.index("CO2")]=0.6468;
        background[registry.index("H2O")]=0.3532;
    }
    else background[0]=1.0;
    Domain dom(Nx, Ny, Nz, Lx, Ly, Lz, static_cast<int>(species.size()));
    // Sector state conversion and transport always use one explicit species
    // layout.  Leaving the two-species case implicit would select the legacy
    // constant-Cp thermodynamic branch in part of the pipeline.
    dom.setSpecies(species);
    cout << "  Grid: " << gridName << " " << activeNx << "x" << activeNy << "x" << activeNz
         << " = " << activeNx*activeNy*activeNz << " active cells"
         << " (plus halos)\n";
    cout << "  Domain: L=" << Lx*1e3 << " mm, Rmax=" << Ly*1e3
         << " mm, sector=120 deg\n\n";

    // IC: either a uniform quiescent gas at chamber pressure, or (opt-in) a
    // quasi-1D isentropic nozzle blueprint for the frozen cold-flow state.
    // The blueprint makes the extrapolated supersonic outlet an honest
    // boundary (no 26 kg/s blowout, no closed-outlet stall): chamber at the
    // frozen-design stagnation state, M=1 at the throat, supersonic exit.
    std::vector<double> torchKernel(background.size(),0.0);
    if(reacting)
    {
        torchKernel[registry.index("O2")]=0.784;
        torchKernel[registry.index("C3H8")]=0.216;
    }
    else torchKernel[0]=1.0;
    std::vector<double> blueMach(Nx,0.0),bluePressure(Nx,0.0),blueTemp(Nx,0.0);
    if(blueprintColdFlow)
    {
        // quasi-1D isentropic blueprint; mixture constants for the CO2/H2O
        // background at 700 K (gamma 1.25, R 285.2 J/kgK) - IC only.
        const double gamma=1.25, gasR=285.2, tTot=700.0;
        const double sectorFrac=mesh.sectorAngle()
            /(2.0*3.14159265358979323846);   // sectorAngle() is in radians
        std::vector<double> area(Nx,0.0);
        double aStar=1e30; int throatI=1;
        for(int i=1;i<Nx-1;++i)
        {
            const double r=mesh.wallFace(i);
            area[i]=3.14159265358979323846*r*r*sectorFrac;
            if(area[i]<aStar){ aStar=area[i]; throatI=i; }
        }
        // frozen-design stagnation state consistent with the corrected
        // sector mass flow (3.85/3 kg/s) and the choked-throat relation.
        const double mdot=0.9007+0.3378+0.1348/3.0;
        const double cRatio=2.0/(gamma+1.0);
        const double expnt=(gamma+1.0)/(2.0*(gamma-1.0));
        const double pTot=mdot
            /(std::sqrt(gamma/(gasR*tTot))*std::pow(cRatio,expnt)*aStar);
        for(int i=1;i<Nx-1;++i)
        {
            const double ratio=area[i]/aStar;
            const double lo=(i<throatI)?1.0e-6:1.0+1.0e-9;
            const double hi=(i<throatI)?1.0-1.0e-9:8.0;
            const auto f=[&](double M)
            {
                return std::pow(cRatio*(1.0+0.5*(gamma-1.0)*M*M),expnt)/M
                    -ratio;
            };
            double loM=lo,hiM=hi;
            for(int it=0;it<80;++it)
            {
                const double mid=0.5*(loM+hiM);
                if(f(mid)>0.0) hiM=mid; else loM=mid;
            }
            const double M=0.5*(loM+hiM);
            const double T=tTot/std::pow(
                1.0+0.5*(gamma-1.0)*M*M,1.0);
            blueMach[i]=M;
            bluePressure[i]=pTot/std::pow(
                1.0+0.5*(gamma-1.0)*M*M,gamma/(gamma-1.0));
            blueTemp[i]=T;
        }
        const double exitM=blueMach[Nx-2], exitP=bluePressure[Nx-2];
        cout << "  Blueprint cold-flow: pTot=" << pTot/1e6
             << " MPa, Tt=" << tTot << " K, throat i=" << throatI
             << " A*=" << aStar*1e6 << " mm2, exit M=" << exitM
             << " exit p=" << exitP/1e3 << " kPa\n";
    }
    for (int i=0; i<Nx; ++i)
    for (int j=0; j<Ny; ++j)
    for (int k=0; k<Nz; ++k)
    {
        PrimitiveVars& w = dom.stagedW(i,j,k);
        if(blueprintColdFlow)
        {
            const int bi=std::min(std::max(i,1),Nx-2);
            w.p = bluePressure[bi]; w.T = blueTemp[bi];
            w.U[0] = blueMach[bi]
                *std::sqrt(1.25*285.2*blueTemp[bi]);
            w.U[1] = w.U[2] = 0.0;
            w.Y=background;
        }
        else
        {
            w.p = icPressure; w.T = icTemperature;
            w.U[0] = 0.0;
            w.U[1] = w.U[2] = 0.0;
            w.Y=background;
        }
        if(reacting&&mesh.xFace(i)<0.010&&mesh.radiusCenter(i,j)<0.003)
        {
            w.T=2400.0;
            w.Y=torchKernel;
        }
    }
    // Load and register the GERG-2008 residual table: the dense-cryogenic
    // thermodynamic backend for the 130 K LOX boundary (the PR isochoric
    // path has no single-phase liquid root there).
    static GERGResidualTable gergResidualTable;
    {
        const std::string gergPath =
            std::string(RFAM_SOURCE_DIR)
            + "/tests/reference/gerg_residual_table_v3.csv";
        if (gergResidualTable.load(gergPath))
        {
            dom.thermoModel().setGERGTable(&gergResidualTable);
            cout << "  GERG-2008 residual table loaded: dense cryogenic "
                 << "thermo backend registered [" << gergPath << "]\n";
        }
        else
            cerr << "WARNING: GERG residual table missing; dense LOX "
                 << "thermo stays fail-closed: " << gergPath << "\n";
    }
    // Saturation/HEM companion: the dome + minimal two-phase path rescues
    // LOX near-field states that land on the mechanically unstable
    // spinodal branch (e.g. rho~54.5/e~-231977 at 126 K: teqp single-phase
    // root has dp/drho<0, PR has no liquid root - without HEM the cell
    // fails with INVERSION_FAILED).
    static GERGSaturationTable gergSaturationTable;
    {
        const std::string satPath =
            std::string(RFAM_SOURCE_DIR)
            + "/tests/reference/gerg_saturation_table_v3.csv";
        if (gergSaturationTable.load(satPath))
        {
            dom.thermoModel().setSaturationTable(&gergSaturationTable);
            cout << "  GERG saturation/HEM table loaded: two-phase rescue "
                 << "registered [" << satPath << "]\n";
        }
        else
            cerr << "WARNING: GERG saturation table missing; spinodal "
                 << "cells stay fail-closed: " << satPath << "\n";
    }



    // 加载并注册稠密纯流体输运表（Phase 1c）。
    // 未找到表时稠密查询保持 fail-close（UNSUPPORTED）。
    static CryogenicTransportTable cryogenicTransportTable;
    {
        const std::string tablePath=
            std::string(RFAM_SOURCE_DIR)
            + "/tests/reference/cryogenic_transport_runtime_coolprop72_v1.csv";
        if(cryogenicTransportTable.load(tablePath))
        {
            TransportProps::registerCryogenicTable(&cryogenicTransportTable);
            cout << "  Cryogenic transport table loaded: O2/C3H8 dense-fluid "
                 << "backend registered\n";
        }
        else
            cerr << "WARNING: cryogenic transport table missing; dense LOX "
                 << "transport stays fail-closed: " << tablePath << "\n";
    }


    double hotstartSourceTime=std::numeric_limits<double>::quiet_NaN();
    int hotstartSourceIteration=-1;
    // External hot-start field import (PTY authority).  The field is
    // reconciled through the run's own GERG/PR backend and written into
    // the Domain before the first step; the JSON report is saved next to
    // the field CSV for audit.
    if(!hotstartFieldPath.empty())
    {
        HotStartFieldLoader loader(registry,&gergResidualTable,mesh);
        const HotStartFieldLoadResult field=loader.load(hotstartFieldPath);
        if(!field.valid)
            throw std::runtime_error(
                "hot-start field import failed: "+field.error);
        int imported=0;
        for(int n=0;n<static_cast<int>(field.reconciled.size());++n)
        {
            if(!field.reconciled[static_cast<size_t>(n)].valid) continue;
            const ReconciledCell& rc=field.reconciled[static_cast<size_t>(n)];
            const int i=field.ijk[static_cast<size_t>(n)]/10000;
            const int j=(field.ijk[static_cast<size_t>(n)]/100)%100;
            const int k=field.ijk[static_cast<size_t>(n)]%100;
            PrimitiveVars w(static_cast<int>(species.size()));
            w.p=rc.p; w.T=rc.T;
            w.U[0]=rc.rhoU[0]/rc.rho;
            w.U[1]=rc.rhoU[1]/rc.rho;
            w.U[2]=rc.rhoU[2]/rc.rho;
            w.Y=rc.rhoY;
            for(double& y:w.Y) y/=rc.rho;
            dom.stagedW(i,j,k)=w;
            ++imported;
        }
        const std::string reportPath=hotstartFieldPath+".reconcile.json";
        {
            std::ofstream out(reportPath);
            out<<field.report.toJson()<<"\n";
        }
        hotstartSourceTime=field.sourceTime;
        hotstartSourceIteration=field.sourceIteration;
        cout<<"  Hot-start field imported: "<<imported<<"/"
            <<field.loadedCells<<" cells -> "<<reportPath<<"\n";
        const double massChange=field.report.totalMassAfter>0.0
            ?std::abs(field.report.totalMassAfter
                      -field.report.totalMassBefore)
                /field.report.totalMassBefore:0.0;
        const double energyChange=field.report.totalEnergyBefore!=0.0
            ?std::abs(field.report.totalEnergyAfter
                      -field.report.totalEnergyBefore)
                /std::abs(field.report.totalEnergyBefore):0.0;
        cout<<"    report: massChange="<<massChange
            <<" energyChange="<<energyChange
            <<" singlePhase="<<field.report.singlePhaseCells
            <<" twoPhase="<<field.report.twoPhaseCells<<"\n";
    }
    dom.primitiveToConserved();
    if(!exportAfterImportPath.empty()&&!hotstartFieldPath.empty())
    {
        std::ofstream out(exportAfterImportPath);
        out << "# rfam.hotstart-field/v1\n";
        out << "x,y,z,p,T,u,v,w,rho,e";
        for (const std::string& name : species) out << "," << name;
        out << "\n";
        for (int i = 1; i < Nx - 1; ++i)
        for (int j = 1; j < Ny - 1; ++j)
        for (int k = 1; k < Nz - 1; ++k)
        {
            const auto& w = dom.W(i, j, k);
            const auto& q = dom.Q(i, j, k);
            const auto center = mesh.cellCenter(i, j, k);
            out << center.x << "," << center.y << "," << center.z << ","
                << w.p << "," << w.T << ","
                << w.U[0] << "," << w.U[1] << "," << w.U[2] << ","
                << q.rho << "," << q.e;
            for (double y : w.Y) out << "," << y;
            out << "\n";
        }
        cout << "  Exported imported field -> " << exportAfterImportPath << "\n";
    }
    SectorInjectionBoundary injection(mesh,p0,species,background,injectorFaceMm);
    injection.setGERGTable(&gergResidualTable);
    const SectorInjectionReport& injectionReport=injection.report();
    const double targetMassFlow=injectionReport.torchMassFlow
        +injectionReport.loxMassFlow+injectionReport.fuelMassFlow;
    SectorNonReflectingOutlet outlet(mesh,101325.0);
    SectorWallBoundary wall(700.0);
    const double rampInitial=0.1,rampDuration=2.0e-4;
    double currentBoundaryRamp=rampInitial;

    // Solver
    SolverBase solver(dom);
    DualTimeIntegrator* dualIntegrator=nullptr;
    SteadyLocalTimeIntegrator* steadyIntegrator=nullptr;
    SteadyLUSGSIntegrator* lusgsIntegrator=nullptr;
    if(dualTime)
    {
        auto innerInjectionBoundary=[&](Domain& innerDomain)
        {
            try
            {
                injection.applyInStepTransaction(
                    innerDomain,currentBoundaryRamp);
                return !innerDomain.tripped();
            }
            catch(const std::exception& error)
            {
                if(!innerDomain.tripped())
                    innerDomain.recordFailure(
                        Domain::FailureCode::OPERATOR_EXCEPTION,-1,
                        "dual-time injection boundary",error.what());
                return false;
            }
        };
        dualIntegrator=new DualTimeIntegrator(
            dtMult,dualInnerCfl,dualMaxInnerIters,
            dualResidualTolerance,dualMaxRetries,12,
            dualLocalPseudoTime,cfl,dualLocalPseudoRatio,
            dualAitken,0.05,2.0,dualAnderson,dualAndersonDepth,
            1.0e-10,2.0,
            dualInnerInjectionRefresh
                ?std::function<bool(Domain&)>(innerInjectionBoundary)
                :std::function<bool(Domain&)>(),
            dualInnerInjectionRefresh);
        solver.setIntegrator(unique_ptr<TimeIntegrator>(dualIntegrator));
        cout << "  Time integration: dual-time (physical dt x" << dtMult
             << ", pseudo-CFL " << dualInnerCfl
             << ", residual ratio <= " << dualResidualTolerance
             << ", local pseudo-time "
             << (dualLocalPseudoTime?"ON":"OFF")
             << ", Aitken " << (dualAitken?"ON":"OFF")
             << ", Anderson " << (dualAnderson?"ON":"OFF")
             << ", inner injection refresh "
             << (dualInnerInjectionRefresh?"ON":"OFF")
             << ", inner iters <= " << dualMaxInnerIters << ")\n";
    }
    else if(steadyLocalTime)
    {
        steadyIntegrator=new SteadyLocalTimeIntegrator(
            cfl,steadyResidualTolerance,steadyLocalRatio,24);
        solver.setIntegrator(unique_ptr<TimeIntegrator>(steadyIntegrator));
        cout << "  Time integration: STEADY local pseudo-time "
             << "(physical clock frozen, residual ratio <= "
             << steadyResidualTolerance << ", local ratio <= "
             << steadyLocalRatio << ", outlet continuation "
             << steadyOutletRelaxation << ")\n";
    }
    else if(steadyLusgs)
    {
        auto steadyBoundaryHook=[&](Domain& innerDomain)
        {
            try
            {
                injection.applyInStepTransaction(
                    innerDomain,currentBoundaryRamp);
                std::vector<Domain::CompleteCellReplacement> axis;
                axis.reserve(static_cast<std::size_t>(Nx*Nz));
                for(int i=0;i<Nx;++i)
                for(int k=0;k<Nz;++k)
                {
                    PrimitiveVars primitive=innerDomain.W(i,1,k);
                    primitive.U[1]=0.0;
                    axis.push_back(Domain::CompleteCellReplacement{
                        i,0,k,primitive,innerDomain.Q(i,1,k)});
                }
                if(!innerDomain.replaceCompleteCellsInStepTransaction(axis))
                    return false;
                wall.applyInStepTransaction(innerDomain);
                innerDomain.applyCyclicZInStepTransaction(120.0);
                if(innerDomain.tripped()) return false;
                outlet.applySteadyInStepTransaction(
                    innerDomain,steadyOutletRelaxation);
                innerDomain.applyCyclicZInStepTransaction(120.0);
                return !innerDomain.tripped();
            }
            catch(const std::exception& error)
            {
                if(!innerDomain.tripped())
                    innerDomain.recordFailure(
                        Domain::FailureCode::OPERATOR_EXCEPTION,-1,
                        "steady LU-SGS boundary",error.what());
                return false;
            }
        };
        PressureVelocityCorrectionBuilder pressureCorrectionBuilder;
        const std::vector<unsigned char> inletPressureDirichlet=
            injection.pressureCorrectionDirichletMask();
        if(lusgsPressureSchurFallback)
        {
            pressureCorrectionBuilder=[
                &mesh,inletPressureDirichlet,
                lusgsPressureSchurLeastSquaresGradient,
                lusgsPressureSchurNonOrthogonalCorrection,
                lusgsPressureSchurStructuredMomentum](
                const Domain& domain,
                const std::vector<ConservedVars>& residual,
                const std::vector<unsigned char>& active,
                const StructuredImplicitStencil& stencil,
                double pseudoCfl,
                PressureVelocityCorrection& correction)
            {
                try
                {
                    if(lusgsPressureSchurStructuredMomentum)
                    {
                        const auto result=
                            SectorStructuredPressureSchurPreconditioner::build(
                                domain,mesh,residual,active,stencil,
                                inletPressureDirichlet,pseudoCfl,
                                120,1.0e-8);
                        correction.pressure=result.pressure;
                        correction.velocity=result.velocity;
                        correction.linearSweeps=result.krylovIterations;
                        correction.linearResidualRatio=result.residualRatio;
                        correction.maximumRelativePressure=
                            result.maximumRelativePressure;
                        correction.maximumVelocityMach=
                            result.maximumVelocityMach;
                        correction.structuredKrylovIterations=
                            result.krylovIterations;
                        correction.structuredResidualRatio=
                            result.residualRatio;
                        correction.structuredMomentumSolves=
                            result.audit.momentumSolves;
                        correction.structuredMomentumSweeps=
                            result.audit.momentumSweeps;
                        correction.structuredMomentumMaximumResidualRatio=
                            result.audit.maximumMomentumResidualRatio;
                        correction.structuredAcousticPreconditionerSolves=
                            result.audit.acousticPreconditionerSolves;
                        correction.structuredAcousticPreconditionerSweeps=
                            result.audit.acousticPreconditionerSweeps;
                        correction.structuredAcousticPreconditionerMaximumResidualRatio=
                            result.audit.maximumAcousticPreconditionerResidualRatio;
                        correction.structuredPressureDiffusionFaces=
                            result.audit.pressureDiffusionFaces;
                        correction.structuredMinimumPressureDiffusionMobility=
                            std::isfinite(result.audit.minimumPressureDiffusionMobility)
                            ?result.audit.minimumPressureDiffusionMobility:0.0;
                        correction.structuredMaximumPressureDiffusionMobility=
                            result.audit.maximumPressureDiffusionMobility;
                        correction.structuredMaximumPressureDiffusionCoefficient=
                            result.audit.maximumPressureDiffusionCoefficient;
                        correction.structuredMaximumPressureDiffusionReciprocityError=
                            result.audit.maximumPressureDiffusionReciprocityError;
                        correction.leastSquaresGradientApplications=
                            result.audit.leastSquaresApplications;
                        correction.leastSquaresGradientFallbacks=
                            result.audit.leastSquaresFallbacks;
                        int pressureDirichletFaces=0;
                        for(unsigned char value:inletPressureDirichlet)
                            pressureDirichletFaces+=value!=0U?1:0;
                        correction.inletPressureDirichletFaces=
                            pressureDirichletFaces;
                        correction.inletDirichletFaces=0;
                        correction.inletNeumannFaces=
                            (domain.Ny()-2)*(domain.Nz()-2);
                        return true;
                    }
                    const auto result=SectorPressureSchurPreconditioner::build(
                        domain,mesh,residual,active,stencil,
                        inletPressureDirichlet,pseudoCfl,160,1.0e-8,
                        lusgsPressureSchurLeastSquaresGradient,
                        lusgsPressureSchurNonOrthogonalCorrection?12:0,
                        true);
                    correction.pressure=result.pressure;
                    correction.velocity=result.velocity;
                    correction.linearSweeps=result.sweeps;
                    correction.linearResidualRatio=result.residualRatio;
                    correction.maximumRelativePressure=
                        result.maximumRelativePressure;
                    correction.maximumVelocityMach=
                        result.maximumVelocityMach;
                    correction.minimumFacePseudoTime=
                        result.minimumFacePseudoTime;
                    correction.maximumFacePseudoTime=
                        result.maximumFacePseudoTime;
                    correction.maximumFaceReciprocityError=
                        result.maximumFaceReciprocityError;
                    correction.maximumFaceNonOrthogonalityDegrees=
                        result.maximumFaceNonOrthogonalityDegrees;
                    correction.leastSquaresGradientApplications=
                        result.leastSquaresGradientApplications;
                    correction.leastSquaresGradientFallbacks=
                        result.leastSquaresGradientFallbacks;
                    correction.nonOrthogonalCorrections=
                        result.nonOrthogonalCorrections;
                    correction.nonOrthogonalResidualRatio=
                        result.nonOrthogonalResidualRatio;
                    correction.inletPressureDirichletFaces=
                        result.inletPressureDirichletFaces;
                    correction.inletDirichletFaces=
                        result.inletDirichletFaces;
                    correction.inletNeumannFaces=
                        result.inletNeumannFaces;
                    return true;
                }
                catch(...)
                {
                    return false;
                }
            };
        }
        lusgsIntegrator=new SteadyLUSGSIntegrator(
            lusgsPseudoCfl,steadyResidualTolerance,24,true,
            lusgsMaximumUpdate,1.0e-4,1.0,
            lusgsAnderson,lusgsAndersonDepth,1.0e-10,2.0,
            lusgsResidualLineSearch,
            lusgsResidualLineSearch
                ?std::function<bool(Domain&)>(steadyBoundaryHook)
                :std::function<bool(Domain&)>(),
            lusgsEnergyScale,1.0,true,1.0e-5,20,
            lusgsBidirectionalLineSearch,lusgsHydroDiagonalBlock,
            lusgsSegregatedEnergyFallback,lusgsSegregatedFlowFallback,
            lusgsSegregatedAcousticFallback,
            std::move(pressureCorrectionBuilder));
        solver.setIntegrator(unique_ptr<TimeIntegrator>(lusgsIntegrator));
        if(lusgsNoSpeciesLimiter)
        {
            lusgsIntegrator->setUseSpeciesLimiter(false);
            cout << "  LU-SGS species limiter DISABLED (diagnostic)\n";
        }
        cout << "  Time integration: STEADY LU-SGS "
             << "(physical clock frozen, residual ratio <= "
             << steadyResidualTolerance << ", pseudo-CFL "
             << lusgsPseudoCfl << ", max update "
             << lusgsMaximumUpdate << ", energy scale "
             << lusgsEnergyScale << ", Anderson "
             << (lusgsAnderson?"ON":"OFF") << " depth "
             << lusgsAndersonDepth << ", frozen Euler block "
             << (lusgsFrozenEulerBlock?"ON":"OFF")
             << ", hydro diagonal block "
             << (lusgsHydroDiagonalBlock?"ON":"OFF")
             << ", pressure Schur fallback "
             << (lusgsPressureSchurFallback?"ON":"OFF")
             << ", pressure Schur least-squares gradient "
             << (lusgsPressureSchurLeastSquaresGradient?"ON":"OFF")
             << ", pressure Schur nonorthogonal correction "
             << (lusgsPressureSchurNonOrthogonalCorrection?"ON":"OFF")
             << ", pressure Schur structured momentum "
             << (lusgsPressureSchurStructuredMomentum?"ON":"OFF")
             << ", segregated acoustic fallback "
             << (lusgsSegregatedAcousticFallback?"ON":"OFF")
             << ", segregated flow fallback "
             << (lusgsSegregatedFlowFallback?"ON":"OFF")
             << ", segregated energy fallback "
             << (lusgsSegregatedEnergyFallback?"ON":"OFF")
             << ", residual line search "
             << (lusgsResidualLineSearch?"ON":"OFF")
             << ", bidirectional "
             << (lusgsBidirectionalLineSearch?"ON":"OFF")
             << ", outlet continuation "
             << steadyOutletRelaxation << ")\n";
    }
    if (std::isfinite(hotstartSourceTime))
    {
        solver.setRestartState(hotstartSourceTime,
            hotstartSourceIteration >= 0 ? hotstartSourceIteration : 0);
        cout << "  Hot-start clock synced: t=" << hotstartSourceTime
             << " iter=" << hotstartSourceIteration << "\n";
    }
    std::vector<ProfiledOperator*> profiles;
    auto addProfiled=[&](std::unique_ptr<Operator> op)
    {
        auto profiled=std::make_unique<ProfiledOperator>(std::move(op));
        profiles.push_back(profiled.get()); solver.addOperator(std::move(profiled));
    };
    addProfiled(unique_ptr<Operator>(new SectorConvectionOperator(
        unique_ptr<FluxScheme>(reacting
            ?static_cast<FluxScheme*>(new AUSMFluxAdapter(species))
            :static_cast<FluxScheme*>(new AUSMFluxAdapter())),mesh)));
    {
        auto* conv = dynamic_cast<SectorConvectionOperator*>(
            profiles.back()->inner());
        if (conv && cflAcousticScale)
            conv->setCflAcousticScaling(true);
        if (conv)
        {
            conv->setFrozenEulerBlockCoupling(lusgsFrozenEulerBlock);
            conv->setHydroDiagonalBlockCoupling(
                lusgsHydroDiagonalBlock);
        }
    }
    SectorTurbulenceOperator* turbulence=nullptr;
    if(useRans || useLes)
    {
        turbulence=new SectorTurbulenceOperator(
            mesh,registry,0.05,0.001,0.25,
            allowSinglePhaseReactingDiagnostic,
            useLes ? SectorTurbulenceOperator::TurbulenceMode::WALE
                   : SectorTurbulenceOperator::TurbulenceMode::SST);
        addProfiled(unique_ptr<Operator>(turbulence));
    }
    addProfiled(unique_ptr<Operator>(new SectorDiffusionOperator(
        mesh,registry,turbulence,0.3,
        allowSinglePhaseReactingDiagnostic)));
#ifdef RFAM_HAVE_CHEMISTRY
    SectorChemistryOperator* chemistryOperator=nullptr;
    SandiegoSectorChemistryOperator* sandiegoOperator=nullptr;
    std::unique_ptr<TwoStepCO> chemistry;
    // Dual-time mode: CVODE chemistry is a rigid split operator, excluded
    // from the inner iteration chain and applied once per physical step.
    std::unique_ptr<Operator> splitChemistry;
    if(reacting && !freezeChemistry)
    {
        if(sandiegoMechanism)
        {
            sandiegoOperator=new SandiegoSectorChemistryOperator(registry);
            if(dualTime)
                splitChemistry.reset(sandiegoOperator);
            else
                addProfiled(unique_ptr<Operator>(sandiegoOperator));
            cout << "  Chemistry: San Diego 59-species mechanism"
                 << " (CVODE cell reactor, 9<->59 bridge)";
            if(dualTime) cout << " [dual-time split operator]";
            cout << "\n";
        }
        else
        {
            chemistry=std::make_unique<TwoStepCO>();
            chemistryOperator=new SectorChemistryOperator(chemistry.get(),registry);
            if (chemistryExplicit)
            {
                chemistryOperator->setExplicitChemistry(true, 256);
                cout << "  Chemistry: EXPLICIT substeps (diagnostic, no CVODE)\n";
            }
            if(dualTime)
                splitChemistry.reset(chemistryOperator);
            else
                addProfiled(unique_ptr<Operator>(chemistryOperator));
        }
    }
#else
    if(reacting && !freezeChemistry)
    {
        cerr << "reacting mode unavailable: rebuild with SUNDIALS_ROOT\n";
        return 2;
    }
#endif
    if(freezeChemistry)
        cout << "  Chemistry: FROZEN (9-species state retained; source "
             << "operators disabled)\n";
    if(!reacting)
        addProfiled(unique_ptr<Operator>(new IsothermalProjectionOperator(700.0)));

    const std::string turbulenceSuffix=
        useLes?"-les":(useRans?"-rans":"-laminar");
    const std::string runMode=std::string(freezeChemistry
        ?"frozen-chemistry":(reacting?"finite-rate":"isothermal"))
        +turbulenceSuffix;
    const std::string finiteRateRunMode="finite-rate"+turbulenceSuffix;
    const std::uint64_t profileHash=SectorCheckpoint::hashFile(profilePath);
    if(!restartPath.empty())
    {
        const SectorCheckpointMetadata loaded=
            SectorCheckpoint::read(restartPath,dom,turbulence);
        const bool migratingFiniteRateToFrozen=
            freezeChemistry && loaded.mode==finiteRateRunMode;
        const bool migratingFrozenToFiniteRate=
            allowFrozenToFiniteRate && !freezeChemistry
            && loaded.mode=="frozen-chemistry"+turbulenceSuffix;
        if(loaded.mode!=runMode && !migratingFiniteRateToFrozen
            && !migratingFrozenToFiniteRate)
            throw std::runtime_error("checkpoint mode mismatch: loaded="
                +loaded.mode+" requested="+runMode);
        if(loaded.profileHash!=profileHash) throw std::runtime_error("checkpoint profile hash mismatch");
        solver.setRestartState(loaded.time,loaded.iteration);
        cout << "  Restarted " << restartPath << " at iter=" << loaded.iteration
             << " t=" << loaded.time << "s git=" << loaded.gitCommit;
        if(migratingFiniteRateToFrozen)
            cout << " mode-migration=" << loaded.mode << "->" << runMode;
        cout << "\n";

        if(!fasSourcePath.empty())
        {
            SectorResidualField::Data fasSource=
                SectorResidualField::read(fasSourcePath);
            const std::uint64_t restartHash=
                SectorCheckpoint::hashFile(restartPath);
            SectorResidualField::requireCompatible(
                fasSource,dom,profileHash,runMode,restartHash,
                SectorResidualField::Kind::FAS_SOURCE);
            if(fasSource.metadata.time!=loaded.time
                ||fasSource.metadata.iteration!=loaded.iteration)
                throw std::runtime_error(
                    "FAS source clock/iteration does not match restart");
            for(int cell=0;cell<dom.N();++cell)
            {
                const int i=cell%Nx;
                const int j=(cell/Nx)%Ny;
                const int k=cell/(Nx*Ny);
                const unsigned char expected=
                    i>0&&i<Nx-1&&j>0&&j<Ny-1&&k>0&&k<Nz-1
                    ?1U:0U;
                if(fasSource.active[static_cast<std::size_t>(cell)]
                    !=expected)
                    throw std::runtime_error(
                        "FAS source active mask does not match coarse interior");
            }
            const std::uint64_t fineHash=
                fasSource.metadata.fineStateCheckpointHash;
            const double referenceStep=fasSource.metadata.referenceStep;
            const std::string sourceGit=fasSource.metadata.gitCommit;
            solver.addOperator(
                std::make_unique<FixedConservedSourceOperator>(
                    std::move(fasSource.field),
                    std::move(fasSource.active)));
            cout << "  FAS source: " << fasSourcePath
                 << " coarseHash=" << restartHash
                 << " fineHash=" << fineHash
                 << " referenceStep=" << referenceStep
                 << " git=" << sourceGit << "\n";
        }
    }
    const int segmentStartIteration=solver.iteration();

    // Boundary conditions
    auto applyBC = [&](){
        if(steadyRampOverride>=0.0)
            currentBoundaryRamp=steadyRampOverride;
        else if(steadyRampSteps>0)
            currentBoundaryRamp=rampInitial+(1.0-rampInitial)
                *std::min(double(solver.iteration())/double(steadyRampSteps),1.0);
        else
            currentBoundaryRamp=rampInitial+(1.0-rampInitial)
                *std::min(solver.time()/rampDuration,1.0);
        injection.apply(dom,currentBoundaryRamp);
        if(torchTubeMm>0.0)
        {
            // Submerged torch tube: the torch core is placed at x=L so the
            // flame sits inside the chamber flow field (design change:
            // lengthened torch tube).  Replace the tube-end cylindrical
            // cells with the torch state (same state as the face torch).
            std::vector<Domain::CompleteCellReplacement> tube;
            for(int i=1;i<Nx-1;++i)
            for(int j=1;j<Ny-1;++j)
            for(int k=1;k<Nz-1;++k)
            {
                const double xc=mesh.cellCenter(i,j,k).x;
                // whole tube bore: all r<2mm cells from the board face to
                // the outlet plane hold the torch state (no dead cells in
                // the tube, which dual-time big-dt steps would blow up).
                if(xc>torchTubeMm*1e-3+mesh.length()/Nx) continue;
                if(mesh.radiusCenter(i,j)>=0.002) continue;
                const PrimitiveVars torchState=injection.torchTubePrimitive(
                    currentBoundaryRamp);
                const thermo::EquilibriumState torchQ=injection.torchTubeThermo(
                    dom,currentBoundaryRamp);
                tube.push_back(Domain::CompleteCellReplacement{
                    i,j,k,torchState,torchQ});
            }
            if(!tube.empty())
                dom.replaceCompleteCells(tube);
        }
        // Axis (j=0): symmetry
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i=0; i<Nx; ++i)
        for (int k=0; k<Nz; ++k)
        {
            PrimitiveVars w = dom.W(i,1,k);
            w.U[1] = 0.0;   // zero radial velocity at axis
            dom.setCell(i,0,k,w);
        }
        // Wall (j=Ny-1): no-slip isothermal
        wall.apply(dom);
        // Cyclic z-direction (120°)
        dom.applyCyclicZ(sectorDegrees);
    };

    // Time march. Print the exact checkpoint mode so audit logs and restart
    // compatibility use one vocabulary (including LES and frozen chemistry).
    cout << "  Mode: " << runMode << "\n";
    cout << "  Running " << maxSteps << " steps...\n";
    const auto wallStart = std::chrono::steady_clock::now();
    double previousDt=0.0;
    double boundarySeconds=0.0,outletSeconds=0.0,stepSeconds=0.0;
    bool chemistryFailed=false;
    bool stepRejected=false;
    SectorBalance lastAppliedBoundaryBalance{};
    SectorInjectionFluxReport lastAppliedInjectionFlux{};
    bool haveAppliedBoundaryBalance=false;
#ifdef RFAM_HAVE_CHEMISTRY
    auto reportChemistryFailure = [&]()
    {
        const bool failed = chemistryOperator
            ? chemistryOperator->lastFailures() > 0
            : (sandiegoOperator && sandiegoOperator->lastFailures() > 0);
        if(!failed) return;
        if(chemistryOperator)
        {
            const std::streamsize previousPrecision=cerr.precision();
            cerr << std::setprecision(17)
                 << "  CHEMISTRY_FAIL iter=" << solver.iteration()
                 << " cells=" << chemistryOperator->lastFailures()
                 << " cvode/state/temperature/mass/element="
                 << chemistryOperator->lastCvodeFailures()
                 << "/" << chemistryOperator->lastStateFailures()
                 << "/" << chemistryOperator->lastTemperatureFailures()
                 << "/" << chemistryOperator->lastMassFailures()
                 << "/" << chemistryOperator->lastElementFailures()
                 << " max_element_relerr="
                 << chemistryOperator->lastFailureElementError()
                 << " cell=" << chemistryOperator->lastFailureI() << ","
                 << chemistryOperator->lastFailureJ() << ","
                 << chemistryOperator->lastFailureK()
                 << " dt/rho/p/T=" << chemistryOperator->lastFailureDt() << "/"
                 << chemistryOperator->lastFailureRho() << "/"
                 << chemistryOperator->lastFailurePressure() << "/"
                 << chemistryOperator->lastFailureTemperature()
                 << " Ychem=";
            for(double y:chemistryOperator->lastFailureChemY()) cerr << y << ",";
            cerr << "\n" << std::setprecision(previousPrecision);
            return;
        }
        const std::streamsize previousPrecision=cerr.precision();
        cerr << std::setprecision(17)
             << "  CHEMISTRY_FAIL iter=" << solver.iteration()
             << " cells=" << sandiegoOperator->lastFailures()
             << " cell=" << sandiegoOperator->lastFailureI() << ","
             << sandiegoOperator->lastFailureJ() << ","
             << sandiegoOperator->lastFailureK()
             << " rho/T=" << sandiegoOperator->lastFailureRho() << "/"
             << sandiegoOperator->lastFailureTemperature() << " Y9=";
        for(double y:sandiegoOperator->lastFailureY9()) cerr << y << ",";
        cerr << "\n" << std::setprecision(previousPrecision);
    };
#endif
    for (int s=0; s<maxSteps && !dom.tripped(); ++s)
    {
        auto phaseStart=std::chrono::steady_clock::now();
        applyBC();
        boundarySeconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-phaseStart).count();
        if(dom.tripped()) break;
        if(previousDt<=0.0) previousDt=solver.computeTimestep(cfl);
        phaseStart=std::chrono::steady_clock::now();
        if(steadyIntegrator||lusgsIntegrator)
            outlet.applySteady(dom,steadyOutletRelaxation);
        else outlet.apply(dom,previousDt);
        outletSeconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-phaseStart).count();
        if(dom.tripped()) break;
        dom.applyCyclicZ(sectorDegrees);
        if(dom.tripped()) break;
        const bool captureBoundary=(monitorEvery>0
            &&(solver.iteration()+1)%monitorEvery==0)||s==maxSteps-1;
        if(captureBoundary)
        {
            lastAppliedBoundaryBalance=SectorDiagnostics(mesh).evaluate(dom);
            lastAppliedInjectionFlux=injection.numericalMassFlows(dom);
            haveAppliedBoundaryBalance=true;
        }
        const double dtAc = solver.computeTimestep(cfl);
        const double dt = dualTime ? dtMult * dtAc : dtAc;
        double acceptedDt=dt;
        phaseStart=std::chrono::steady_clock::now();
        bool stepAccepted=false;
        try
        {
            stepAccepted=solver.step(dt);
            if(stepAccepted)
                acceptedDt=solver.previousAcceptedDt().value_or(dt);
#ifdef RFAM_HAVE_CHEMISTRY
            // Dual-time: CVODE chemistry is a split operator applied once per
            // physical step after the flow inner iteration converges.
            if(stepAccepted && dualTime && splitChemistry)
            {
                dom.beginStep();
                splitChemistry->apply(dom,acceptedDt);
                if(dom.tripped() || !dom.commit()) stepAccepted=false;
            }
#endif
        }
        catch(const std::exception& error)
        {
#ifdef RFAM_HAVE_CHEMISTRY
            reportChemistryFailure();
#endif
            if(!dom.tripped())
                dom.recordFailure(Domain::FailureCode::OPERATOR_EXCEPTION,-1,
                    "solver_3d_sector",error.what());
        }
        catch(...)
        {
#ifdef RFAM_HAVE_CHEMISTRY
            reportChemistryFailure();
#endif
            if(!dom.tripped())
                dom.recordFailure(Domain::FailureCode::OPERATOR_EXCEPTION,-1,
                    "solver_3d_sector","unknown exception during time step");
        }
        stepSeconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-phaseStart).count();
        if(!stepAccepted)
        {
            stepRejected=true;
            if(!dom.tripped())
                dom.recordFailure(Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
                    "solver_3d_sector",
                    "time integrator rejected step without failure detail");
            break;
        }
        previousDt=acceptedDt;
#ifdef RFAM_HAVE_CHEMISTRY
        // Tolerate a small number of locally skipped chemistry cells in the
        // diagnostic San Diego path; a runaway failure count still trips.
        const long chemActive = sandiegoOperator
            ? sandiegoOperator->lastActive()
            : (chemistryOperator ? chemistryOperator->lastActiveCells() : 0);
        const long chemFailures = chemistryOperator
            ? chemistryOperator->lastFailures()
            : (sandiegoOperator ? sandiegoOperator->lastFailures() : 0);
        const bool chemistryFailures = chemistryOperator
            ? chemFailures > 0
            : chemFailures > std::max(24L, chemActive / 8L);
        if(chemistryFailures)
        {
            chemistryFailed=true;
            reportChemistryFailure();
            break;
        }
#endif
        if(monitorEvery>0&&solver.iteration()%monitorEvery==0)
        {
            double monitorTmin=1e30,monitorTmax=-1e30;
            for(int i=1;i<Nx-1;++i) for(int j=1;j<Ny-1;++j)
            for(int k=1;k<Nz-1;++k)
            { monitorTmin=std::min(monitorTmin,dom.W(i,j,k).T);
              monitorTmax=std::max(monitorTmax,dom.W(i,j,k).T); }
            const SectorBalance monitorBalance=SectorDiagnostics(mesh).evaluate(dom);
            cout << "  MONITOR iter=" << solver.iteration() << " t=" << solver.time()
                 << " mass=" << monitorBalance.massInventory << " mdot_in/out="
                 << lastAppliedBoundaryBalance.inletMassFlow << "/"
                 << lastAppliedBoundaryBalance.outletMassFlow
                 << " mdot_axis/wall=" << lastAppliedBoundaryBalance.axisOutwardMassFlow
                 << "/" << lastAppliedBoundaryBalance.wallOutwardMassFlow
                 << " mdot_patch=" << lastAppliedInjectionFlux.background << "/"
                 << lastAppliedInjectionFlux.torch << "/"
                 << lastAppliedInjectionFlux.lox << "/" << lastAppliedInjectionFlux.fuel
                 << " T=" << monitorTmin << "/" << monitorTmax;
            if(reacting)
            {
                const SectorSpeciesBalance speciesBalance=
                    SectorDiagnostics(mesh).evaluateSpecies(dom,registry);
                const SectorPerformance p=
                    SectorPerformanceDiagnostics(mesh,registry).evaluate(dom,turbulence);
                cout << " Pc=" << p.chamberPressure << " throat=" << p.throatMassFlow
                     << " cstar=" << p.cStar << " qmax/qavg="
                     << p.wallHeatFluxMaximum << "/" << p.wallHeatFluxAverage
                     << " mdot_target=" << targetMassFlow
                     << " throat/target=" << p.throatMassFlow/targetMassFlow
                     << " minY=";
                for(size_t s=0;s<species.size();++s)
                    cout << (s?",":"") << species[s] << ":"
                         << speciesBalance.minimumMassFraction[s];
                cout << " sumYerr=" << speciesBalance.maximumMassFractionSumError
                     << " elem=" << speciesBalance.elementInventory[0] << "/"
                     << speciesBalance.elementInventory[1] << "/"
                     << speciesBalance.elementInventory[2] << "/"
                     << speciesBalance.elementInventory[3];
#ifdef RFAM_HAVE_CHEMISTRY
                if(chemistryOperator) cout << " chem=" << chemistryOperator->lastActiveCells()
                    << "/" << chemistryOperator->lastRetriedCells()
                    << "/" << chemistryOperator->lastFailures()
                    << "/limited:" << chemistryOperator->lastPositivityLimitedCells()
                    << " qchem_min/max/avg=" << chemistryOperator->lastHeatReleaseMinimum()
                    << "/" << chemistryOperator->lastHeatReleaseMaximum()
                    << "/" << chemistryOperator->lastHeatReleaseAverage();
                else if(sandiegoOperator)
                    cout << " sandiego_chem=" << sandiegoOperator->lastActiveCells()
                         << "/" << sandiegoOperator->lastFailures()
                         << " qchem_min/max/avg="
                         << sandiegoOperator->lastHeatReleaseMinimum()
                         << "/" << sandiegoOperator->lastHeatReleaseMaximum()
                         << "/" << sandiegoOperator->lastHeatReleaseAverage();
#endif
            }
            if(dualIntegrator)
                cout << " dual_inner/res0/res/ratio/update/backtrack/retry/dt="
                     << dualIntegrator->lastInnerIters() << "/"
                     << dualIntegrator->lastInitialResidual() << "/"
                     << dualIntegrator->lastResidual() << "/"
                     << dualIntegrator->lastResidualRatio() << "/"
                     << dualIntegrator->lastRelativeUpdate() << "/"
                     << dualIntegrator->lastInnerBacktracks() << "/"
                     << dualIntegrator->lastStepRetries() << "/"
                     << dualIntegrator->lastAcceptedDt()
                     << " pseudo_dt_min/max="
                     << dualIntegrator->lastPseudoDtMinimum() << "/"
                     << dualIntegrator->lastPseudoDtMaximum()
                     << " aitken="
                     << dualIntegrator->lastAitkenRelaxation()
                     << " anderson_depth="
                     << dualIntegrator->lastAndersonDepthUsed();
            if(steadyIntegrator)
                cout << " steady_res0/res/ratio/update/backtrack="
                     << steadyIntegrator->initialResidual() << "/"
                     << steadyIntegrator->lastResidual() << "/"
                     << steadyIntegrator->lastResidualRatio() << "/"
                     << steadyIntegrator->lastRelativeUpdate() << "/"
                     << steadyIntegrator->lastBacktracks()
                     << " relaxation=" << steadyIntegrator->relaxation()
                     << " pseudo_dt_min/max="
                     << steadyIntegrator->lastPseudoDtMinimum() << "/"
                     << steadyIntegrator->lastPseudoDtMaximum();
            if(lusgsIntegrator)
            {
                cout << " lusgs_res0/res/ratio/update/backtrack="
                     << lusgsIntegrator->initialResidual() << "/"
                     << lusgsIntegrator->lastResidual() << "/"
                     << lusgsIntegrator->lastResidualRatio() << "/"
                     << lusgsIntegrator->lastRelativeUpdate() << "/"
                     << lusgsIntegrator->lastBacktracks()
                     << " relaxation=" << lusgsIntegrator->relaxation()
                     << " anderson_depth="
                     << lusgsIntegrator->lastAndersonDepthUsed()
                     << " block_faces/fallbacks="
                     << lusgsIntegrator->lastBlockFaceApplications() << "/"
                     << lusgsIntegrator->lastBlockFaceFallbacks()
                     << " hydro_blocks/fallbacks/commit/pivot/output/matrix/rhs/offdiag="
                     << lusgsIntegrator->lastHydroBlockApplications() << "/"
                     << lusgsIntegrator->lastHydroBlockFallbacks() << "/"
                     << lusgsIntegrator->lastHydroBlockCommitFallbacks() << "/"
                     << lusgsIntegrator->lastHydroBlockPivotFallbacks() << "/"
                     << lusgsIntegrator->lastHydroBlockOutputFallbacks() << "/"
                     << lusgsIntegrator->lastHydroBlockMatrixFallbacks() << "/"
                     << lusgsIntegrator->lastHydroBlockRhsFallbacks() << "/"
                     << lusgsIntegrator->lastHydroBlockMaximumNormalizedOffDiagonal()
                     << " subspace/schur_evals/schur_best/schur_sweeps/schur_linear_ratio/schur_dp_rel/schur_du_mach/schur_face_dtau_min/schur_face_dtau_max/schur_face_reciprocity/schur_nonorth_deg/schur_ls/schur_ls_fallback/schur_nonorth_corrections/schur_nonorth_linear_ratio/schur_inlet_pressure_D/schur_inlet_continuity_D/schur_inlet_continuity_N/acoustic_evals/acoustic_best/flow_evals/flow_best/energy_evals/energy_best="
                     << lusgsIntegrator->lastAcceptedResidualSubspace() << "/"
                     << lusgsIntegrator->lastPressureSchurResidualEvaluations() << "/"
                     << lusgsIntegrator->lastBestPressureSchurTrialResidual() << "/"
                     << lusgsIntegrator->lastPressureSchurLinearSweeps() << "/"
                     << lusgsIntegrator->lastPressureSchurLinearResidualRatio() << "/"
                     << lusgsIntegrator->lastPressureSchurMaximumRelativePressure() << "/"
                     << lusgsIntegrator->lastPressureSchurMaximumVelocityMach() << "/"
                     << lusgsIntegrator->lastPressureSchurMinimumFacePseudoTime() << "/"
                     << lusgsIntegrator->lastPressureSchurMaximumFacePseudoTime() << "/"
                     << lusgsIntegrator->lastPressureSchurMaximumFaceReciprocityError() << "/"
                     << lusgsIntegrator->lastPressureSchurMaximumFaceNonOrthogonalityDegrees() << "/"
                     << lusgsIntegrator->lastPressureSchurLeastSquaresGradientApplications() << "/"
                     << lusgsIntegrator->lastPressureSchurLeastSquaresGradientFallbacks() << "/"
                     << lusgsIntegrator->lastPressureSchurNonOrthogonalCorrections() << "/"
                     << lusgsIntegrator->lastPressureSchurNonOrthogonalResidualRatio() << "/"
                     << lusgsIntegrator->lastPressureSchurInletPressureDirichletFaces() << "/"
                     << lusgsIntegrator->lastPressureSchurInletDirichletFaces() << "/"
                     << lusgsIntegrator->lastPressureSchurInletNeumannFaces() << "/"
                     << lusgsIntegrator->lastAcousticResidualEvaluations() << "/"
                     << lusgsIntegrator->lastBestAcousticTrialResidual() << "/"
                     << lusgsIntegrator->lastFlowResidualEvaluations() << "/"
                     << lusgsIntegrator->lastBestFlowTrialResidual() << "/"
                     << lusgsIntegrator->lastEnergyResidualEvaluations() << "/"
                     << lusgsIntegrator->lastBestEnergyTrialResidual()
                     << " schur_structured=krylov/ratio/mm_solves/mm_sweeps/mm_ratio/acoustic_solves/acoustic_sweeps/acoustic_ratio:"
                     << lusgsIntegrator->lastPressureSchurStructuredKrylovIterations() << "/"
                     << lusgsIntegrator->lastPressureSchurStructuredResidualRatio() << "/"
                     << lusgsIntegrator->lastPressureSchurStructuredMomentumSolves() << "/"
                     << lusgsIntegrator->lastPressureSchurStructuredMomentumSweeps() << "/"
                     << lusgsIntegrator->lastPressureSchurStructuredMomentumMaximumResidualRatio() << "/"
                     << lusgsIntegrator->lastPressureSchurStructuredAcousticPreconditionerSolves() << "/"
                     << lusgsIntegrator->lastPressureSchurStructuredAcousticPreconditionerSweeps() << "/"
                     << lusgsIntegrator->lastPressureSchurStructuredAcousticPreconditionerMaximumResidualRatio()
                     << " schur_pressure_diffusion=faces/min_mobility/max_mobility/max_coefficient/reciprocity:"
                     << lusgsIntegrator->lastPressureSchurStructuredPressureDiffusionFaces() << "/"
                     << lusgsIntegrator->lastPressureSchurStructuredMinimumPressureDiffusionMobility() << "/"
                     << lusgsIntegrator->lastPressureSchurStructuredMaximumPressureDiffusionMobility() << "/"
                     << lusgsIntegrator->lastPressureSchurStructuredMaximumPressureDiffusionCoefficient() << "/"
                     << lusgsIntegrator->lastPressureSchurStructuredMaximumPressureDiffusionReciprocityError()
                     << " trial_res/full/best/evals/direction/scale="
                     << lusgsIntegrator->lastTrialResidual() << "/"
                     << lusgsIntegrator->lastFullTrialResidual() << "/"
                     << lusgsIntegrator->lastBestTrialResidual() << "/"
                     << lusgsIntegrator->lastResidualEvaluations() << "/"
                     << lusgsIntegrator->lastAcceptedResidualDirection() << "/"
                     << lusgsIntegrator->lastAcceptedResidualScale()
                     << " residual_limit="
                     << lusgsIntegrator->lastResidualAcceptanceLimit()
                     << " diagonal_min/max="
                     << lusgsIntegrator->lastDiagonalMinimum() << "/"
                     << lusgsIntegrator->lastDiagonalMaximum()
                     << " face_max="
                     << lusgsIntegrator->lastMaximumFaceCoupling()
                     << " species_limiter="
                     << lusgsIntegrator->lastProjectedCells() << "/"
                     << lusgsIntegrator->lastMaximumSpeciesProjection()
                     << " component_res=";
                for(double value:lusgsIntegrator->lastComponentResiduals())
                    cout << value << ",";
                cout << " full_trial_component_res=";
                for(double value:
                    lusgsIntegrator->lastFullTrialComponentResiduals())
                    cout << value << ",";
                cout << " best_schur_component_res=";
                for(double value:
                    lusgsIntegrator->lastBestPressureSchurComponentResiduals())
                    cout << value << ",";
                cout << " accepted_schur_component_res=";
                for(double value:
                    lusgsIntegrator->lastAcceptedPressureSchurComponentResiduals())
                    cout << value << ",";
                cout << " max_res_cell/component/value="
                     << lusgsIntegrator->lastMaximumResidualCell() << "/"
                     << lusgsIntegrator->lastMaximumResidualComponent() << "/"
                     << lusgsIntegrator->lastMaximumResidual();
            }
            cout << " dtlim=";
            for(size_t p=0;p<profiles.size();++p)
                cout << (p?",":"") << profiles[p]->name() << ":"
                     << profiles[p]->lastStableDt();
            cout << "\n" << std::flush;
        }
        if(checkpointEvery>0&&solver.iteration()%checkpointEvery==0&&!dom.tripped())
        {
            int threadCount=1;
#ifdef _OPENMP
            threadCount=omp_get_max_threads();
#endif
#ifdef RFAM_GIT_COMMIT
            const std::string gitCommit=RFAM_GIT_COMMIT;
#else
            const std::string gitCommit="unknown";
#endif
            SectorCheckpointMetadata metadata{solver.time(),solver.iteration(),
                threadCount,runMode,gitCommit,profileHash};
            SectorCheckpoint::write(checkpointPath,dom,metadata,turbulence);
            cout << "  checkpoint " << checkpointPath << " iter="
                 << solver.iteration() << " t=" << solver.time() << "s\n";
        }
        const bool steadyConverged=
            (steadyIntegrator&&steadyIntegrator->converged())
            ||(lusgsIntegrator&&lusgsIntegrator->converged());
        if(steadyConverged)
        {
            const double residualRatio=steadyIntegrator
                ?steadyIntegrator->lastResidualRatio()
                :lusgsIntegrator->lastResidualRatio();
            cout << "  STEADY_CONVERGED iter=" << solver.iteration()
                 << " physical_t=" << solver.time()
                 << " residual_ratio=" << residualRatio << "\n";
            break;
        }
    }
    const double wallSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wallStart).count();
    const int segmentIterations=std::max(0,solver.iteration()-segmentStartIteration);
    const double activeCellCount=static_cast<double>(activeNx)
        *static_cast<double>(activeNy)*static_cast<double>(activeNz);

    const bool selectedSteady=steadyIntegrator||lusgsIntegrator;
    const bool finalSteadyConverged=
        (steadyIntegrator&&steadyIntegrator->converged())
        ||(lusgsIntegrator&&lusgsIntegrator->converged());
    if(selectedSteady&&maxSteps>0&&!finalSteadyConverged&&!dom.tripped())
        dom.recordFailure(Domain::FailureCode::TRANSACTION_PROTOCOL,-1,
            lusgsIntegrator?"steady LU-SGS":"steady local time",
            "steady residual gate was not reached before the iteration limit");

    if(dom.tripped())
    {
        const Domain::FailureInfo& failure=dom.failure();
        cerr << "  DOMAIN_FAILURE code=" << Domain::failureCodeName(failure.code)
             << " cell=" << failure.cell;
        if(failure.cell>=0)
        {
            const int failureI=failure.cell%Nx;
            const int failureJ=(failure.cell/Nx)%Ny;
            const int failureK=failure.cell/(Nx*Ny);
            cerr << " (i,j,k=" << failureI << "," << failureJ << ","
                 << failureK << ")";
        }
        cerr << " component=" << failure.component
             << " detail=" << failure.detail
             << " backend=" << failure.thermoDiagnostics.backendStatus
             << " eqstatus=" << static_cast<int>(failure.equilibriumStatus)
             << " accepted/rejected=" << dom.acceptedSteps() << "/"
             << dom.rejectedSteps();
        if(dualIntegrator)
            cerr << " dual_inner/res0/res/ratio/update/backtrack/retry/dt="
                 << dualIntegrator->lastInnerIters() << "/"
                 << dualIntegrator->lastInitialResidual() << "/"
                 << dualIntegrator->lastResidual() << "/"
                 << dualIntegrator->lastResidualRatio() << "/"
                 << dualIntegrator->lastRelativeUpdate() << "/"
                 << dualIntegrator->lastInnerBacktracks() << "/"
                 << dualIntegrator->lastStepRetries() << "/"
                 << dualIntegrator->lastAcceptedDt()
                 << " pseudo_dt_min/max="
                 << dualIntegrator->lastPseudoDtMinimum() << "/"
                 << dualIntegrator->lastPseudoDtMaximum()
                 << " aitken="
                 << dualIntegrator->lastAitkenRelaxation()
                 << " anderson_depth="
                 << dualIntegrator->lastAndersonDepthUsed();
        if(steadyIntegrator)
            cerr << " steady_res0/res/ratio/update/backtrack="
                 << steadyIntegrator->initialResidual() << "/"
                 << steadyIntegrator->lastResidual() << "/"
                 << steadyIntegrator->lastResidualRatio() << "/"
                 << steadyIntegrator->lastRelativeUpdate() << "/"
                 << steadyIntegrator->lastBacktracks()
                 << " relaxation=" << steadyIntegrator->relaxation()
                 << " pseudo_dt_min/max="
                 << steadyIntegrator->lastPseudoDtMinimum() << "/"
                 << steadyIntegrator->lastPseudoDtMaximum();
        if(lusgsIntegrator)
        {
            cerr << " lusgs_res0/res/ratio/update/backtrack="
                 << lusgsIntegrator->initialResidual() << "/"
                 << lusgsIntegrator->lastResidual() << "/"
                 << lusgsIntegrator->lastResidualRatio() << "/"
                 << lusgsIntegrator->lastRelativeUpdate() << "/"
                 << lusgsIntegrator->lastBacktracks()
                 << " relaxation=" << lusgsIntegrator->relaxation()
                 << " anderson_depth="
                 << lusgsIntegrator->lastAndersonDepthUsed()
                 << " block_faces/fallbacks="
                 << lusgsIntegrator->lastBlockFaceApplications() << "/"
                 << lusgsIntegrator->lastBlockFaceFallbacks()
                 << " hydro_blocks/fallbacks/commit/pivot/output/matrix/rhs/offdiag="
                 << lusgsIntegrator->lastHydroBlockApplications() << "/"
                 << lusgsIntegrator->lastHydroBlockFallbacks() << "/"
                 << lusgsIntegrator->lastHydroBlockCommitFallbacks() << "/"
                 << lusgsIntegrator->lastHydroBlockPivotFallbacks() << "/"
                 << lusgsIntegrator->lastHydroBlockOutputFallbacks() << "/"
                 << lusgsIntegrator->lastHydroBlockMatrixFallbacks() << "/"
                 << lusgsIntegrator->lastHydroBlockRhsFallbacks() << "/"
                 << lusgsIntegrator->lastHydroBlockMaximumNormalizedOffDiagonal()
                 << " subspace/schur_evals/schur_best/schur_sweeps/schur_linear_ratio/schur_dp_rel/schur_du_mach/schur_face_dtau_min/schur_face_dtau_max/schur_face_reciprocity/schur_nonorth_deg/schur_ls/schur_ls_fallback/schur_nonorth_corrections/schur_nonorth_linear_ratio/schur_inlet_pressure_D/schur_inlet_continuity_D/schur_inlet_continuity_N/acoustic_evals/acoustic_best/flow_evals/flow_best/energy_evals/energy_best="
                 << lusgsIntegrator->lastAcceptedResidualSubspace() << "/"
                 << lusgsIntegrator->lastPressureSchurResidualEvaluations() << "/"
                 << lusgsIntegrator->lastBestPressureSchurTrialResidual() << "/"
                 << lusgsIntegrator->lastPressureSchurLinearSweeps() << "/"
                 << lusgsIntegrator->lastPressureSchurLinearResidualRatio() << "/"
                 << lusgsIntegrator->lastPressureSchurMaximumRelativePressure() << "/"
                 << lusgsIntegrator->lastPressureSchurMaximumVelocityMach() << "/"
                 << lusgsIntegrator->lastPressureSchurMinimumFacePseudoTime() << "/"
                 << lusgsIntegrator->lastPressureSchurMaximumFacePseudoTime() << "/"
                 << lusgsIntegrator->lastPressureSchurMaximumFaceReciprocityError() << "/"
                 << lusgsIntegrator->lastPressureSchurMaximumFaceNonOrthogonalityDegrees() << "/"
                 << lusgsIntegrator->lastPressureSchurLeastSquaresGradientApplications() << "/"
                 << lusgsIntegrator->lastPressureSchurLeastSquaresGradientFallbacks() << "/"
                 << lusgsIntegrator->lastPressureSchurNonOrthogonalCorrections() << "/"
                 << lusgsIntegrator->lastPressureSchurNonOrthogonalResidualRatio() << "/"
                 << lusgsIntegrator->lastPressureSchurInletPressureDirichletFaces() << "/"
                 << lusgsIntegrator->lastPressureSchurInletDirichletFaces() << "/"
                 << lusgsIntegrator->lastPressureSchurInletNeumannFaces() << "/"
                 << lusgsIntegrator->lastAcousticResidualEvaluations() << "/"
                 << lusgsIntegrator->lastBestAcousticTrialResidual() << "/"
                 << lusgsIntegrator->lastFlowResidualEvaluations() << "/"
                 << lusgsIntegrator->lastBestFlowTrialResidual() << "/"
                 << lusgsIntegrator->lastEnergyResidualEvaluations() << "/"
                 << lusgsIntegrator->lastBestEnergyTrialResidual()
                 << " schur_structured=krylov/ratio/mm_solves/mm_sweeps/mm_ratio/acoustic_solves/acoustic_sweeps/acoustic_ratio:"
                 << lusgsIntegrator->lastPressureSchurStructuredKrylovIterations() << "/"
                 << lusgsIntegrator->lastPressureSchurStructuredResidualRatio() << "/"
                 << lusgsIntegrator->lastPressureSchurStructuredMomentumSolves() << "/"
                 << lusgsIntegrator->lastPressureSchurStructuredMomentumSweeps() << "/"
                 << lusgsIntegrator->lastPressureSchurStructuredMomentumMaximumResidualRatio() << "/"
                 << lusgsIntegrator->lastPressureSchurStructuredAcousticPreconditionerSolves() << "/"
                 << lusgsIntegrator->lastPressureSchurStructuredAcousticPreconditionerSweeps() << "/"
                 << lusgsIntegrator->lastPressureSchurStructuredAcousticPreconditionerMaximumResidualRatio()
                 << " schur_pressure_diffusion=faces/min_mobility/max_mobility/max_coefficient/reciprocity:"
                 << lusgsIntegrator->lastPressureSchurStructuredPressureDiffusionFaces() << "/"
                 << lusgsIntegrator->lastPressureSchurStructuredMinimumPressureDiffusionMobility() << "/"
                 << lusgsIntegrator->lastPressureSchurStructuredMaximumPressureDiffusionMobility() << "/"
                 << lusgsIntegrator->lastPressureSchurStructuredMaximumPressureDiffusionCoefficient() << "/"
                 << lusgsIntegrator->lastPressureSchurStructuredMaximumPressureDiffusionReciprocityError()
                 << " trial_res/full/best/evals/direction/scale="
                 << lusgsIntegrator->lastTrialResidual() << "/"
                 << lusgsIntegrator->lastFullTrialResidual() << "/"
                 << lusgsIntegrator->lastBestTrialResidual() << "/"
                 << lusgsIntegrator->lastResidualEvaluations() << "/"
                 << lusgsIntegrator->lastAcceptedResidualDirection() << "/"
                 << lusgsIntegrator->lastAcceptedResidualScale()
                 << " residual_limit="
                 << lusgsIntegrator->lastResidualAcceptanceLimit()
                 << " diagonal_min/max="
                 << lusgsIntegrator->lastDiagonalMinimum() << "/"
                 << lusgsIntegrator->lastDiagonalMaximum()
                 << " face_max="
                 << lusgsIntegrator->lastMaximumFaceCoupling()
                 << " species_limiter="
                 << lusgsIntegrator->lastProjectedCells() << "/"
                 << lusgsIntegrator->lastMaximumSpeciesProjection()
                 << " component_res=";
            for(double value:lusgsIntegrator->lastComponentResiduals())
                cerr << value << ",";
            cerr << " full_trial_component_res=";
            for(double value:
                lusgsIntegrator->lastFullTrialComponentResiduals())
                cerr << value << ",";
            cerr << " best_schur_component_res=";
            for(double value:
                lusgsIntegrator->lastBestPressureSchurComponentResiduals())
                cerr << value << ",";
            cerr << " accepted_schur_component_res=";
            for(double value:
                lusgsIntegrator->lastAcceptedPressureSchurComponentResiduals())
                cerr << value << ",";
            cerr << " max_res_cell/component/value="
                 << lusgsIntegrator->lastMaximumResidualCell() << "/"
                 << lusgsIntegrator->lastMaximumResidualComponent() << "/"
                 << lusgsIntegrator->lastMaximumResidual();
        }
        cerr << "\n";
        if (failure.cell >= 0)
        {
            const int fi = failure.cell % Nx;
            const int fj = (failure.cell / Nx) % Ny;
            const int fk = failure.cell / (Nx * Ny);
            const ConservedVars& fu = dom.Unew(fi, fj, fk);
            cerr.precision(17);
            cerr << "  FAILURE STATE rho=" << fu.rho
                 << " e=" << (fu.rho > 0.0 ? fu.rhoE / fu.rho : 0.0)
                 << " rhoE=" << fu.rhoE << " Y[0]="
                 << (fu.rhoY.empty() ? -1.0 : fu.rhoY[0] / fu.rho)
                 << " oldT=" << dom.W(fi, fj, fk).T
                 << " rhoU=" << fu.rhoU[0] << "/" << fu.rhoU[1] << "/" << fu.rhoU[2]
                 << " rhoZ=" << fu.rhoZ << " rhoC=" << fu.rhoC
                 << " sumY=" << [&]() {
                        double s = 0.0;
                        for (double y : fu.rhoY) s += y;
                        return s / fu.rho;
                    }() << "\n";
            if (!fu.rhoY.empty())
            {
                cerr << "  FAILURE Y=";
                for (size_t si = 0; si < fu.rhoY.size(); ++si)
                    cerr << (si ? "," : "") << fu.rhoY[si] / fu.rho;
                cerr << "\n";
            }
        }
        bool reported=false;
        for(int i=0;i<Nx&&!reported;++i) for(int j=0;j<Ny&&!reported;++j)
        for(int k=0;k<Nz&&!reported;++k)
        {
            const ConservedVars& u=dom.Unew(i,j,k);
            double speciesMass=0.0,minY=1e30;
            int minimumSpecies=-1;
            bool invalidSpecies=false;
            for(size_t s=0;s<u.rhoY.size();++s)
            {
                const double rhoY=u.rhoY[s];
                speciesMass+=rhoY;
                if(std::isfinite(u.rho)&&u.rho>0.0&&rhoY/u.rho<minY)
                { minY=rhoY/u.rho; minimumSpecies=static_cast<int>(s); }
                if(!std::isfinite(rhoY)) invalidSpecies=true;
            }
            if(!std::isfinite(u.rho)||u.rho<=0.0||!std::isfinite(u.rhoE)
               ||invalidSpecies||minY<-1e-10
               ||(!species.empty()&&std::isfinite(u.rho)&&u.rho>0.0
                  &&std::abs(speciesMass-u.rho)>1e-6*u.rho))
            {
                const PrimitiveVars& w=dom.W(i,j,k);
                cerr << "  TRIP cell=" << i << "," << j << "," << k
                     << " Unew(rho,rhoE,rhoUx)=" << u.rho << "," << u.rhoE
                     << "," << u.rhoU[0] << " old(p,T,Ux)=" << w.p << ","
                     << w.T << "," << w.U[0] << " minY=" << minY
                     << " species=" << (minimumSpecies>=0?species[minimumSpecies]:"unknown")
                     << " sumYerr=" << (speciesMass/u.rho-1.0) << "\n";
                reported=true;
            }
        }
    }

    cout << "  iters=" << solver.iteration() << " t=" << solver.time()
         << "s tripped=" << (dom.tripped()?"YES":"NO") << "\n\n";
    cout << "  wall=" << wallSeconds << "s segment_iters=" << segmentIterations
         << " throughput="
         << (activeCellCount*static_cast<double>(segmentIterations)
             /std::max(wallSeconds,1e-12))
         << " cell-steps/s\n";
    cout << "  profile driver: BC=" << boundarySeconds << "s outlet="
         << outletSeconds << "s step=" << stepSeconds << "s\n";
    for(const ProfiledOperator* profile:profiles)
        cout << "  profile " << profile->name() << ": apply="
             << profile->applySeconds() << "s dt=" << profile->timestepSeconds()
             << "s calls=" << profile->applyCalls() << "/"
             << profile->timestepCalls() << "\n";

    // Refresh halos to the final interior state before evaluating boundary
    // fluxes; the time loop applies BCs at the beginning of each step.
    if(!dom.tripped()) applyBC();
    if(!fieldCsvPath.empty())
        SectorOutputWriter::writeFieldCsv(fieldCsvPath,dom,mesh);
    if(!wallCsvPath.empty())
        SectorOutputWriter::writeWallCsv(wallCsvPath,dom,mesh,registry,turbulence);
    SectorBalance balance=SectorDiagnostics(mesh).evaluate(dom);
    if(haveAppliedBoundaryBalance)
    {
        balance.inletMassFlow=lastAppliedBoundaryBalance.inletMassFlow;
        balance.outletMassFlow=lastAppliedBoundaryBalance.outletMassFlow;
        balance.axisOutwardMassFlow=lastAppliedBoundaryBalance.axisOutwardMassFlow;
        balance.wallOutwardMassFlow=lastAppliedBoundaryBalance.wallOutwardMassFlow;
        balance.periodicNetMassFlow=lastAppliedBoundaryBalance.periodicNetMassFlow;
        balance.periodicRelativeMismatch=
            lastAppliedBoundaryBalance.periodicRelativeMismatch;
    }
    const SectorInjectionFluxReport finalInjectionFlux=
        injection.numericalMassFlows(dom);
    cout << "  mass inventory=" << balance.massInventory << " kg\n"
         << "  mdot inlet/outlet=" << balance.inletMassFlow << " / "
         << balance.outletMassFlow << " kg/s\n"
         << "  mdot axis/wall outward=" << balance.axisOutwardMassFlow << " / "
         << balance.wallOutwardMassFlow << " kg/s\n"
         << "  mdot background/torch/LOX/C3H8="
         << finalInjectionFlux.background << " / "
         << finalInjectionFlux.torch << " / "
         << finalInjectionFlux.lox << " / "
         << finalInjectionFlux.fuel << " kg/s\n"
         << "  periodic mass mismatch=" << balance.periodicRelativeMismatch
         << " (net " << balance.periodicNetMassFlow << " kg/s)\n";
    const double finalRamp=rampInitial+(1.0-rampInitial)
        *std::min(solver.time()/rampDuration,1.0);
    cout << "  prescribed torch/LOX/C3H8=" << injectionReport.torchMassFlow << " / "
         << injectionReport.loxMassFlow << " / " << injectionReport.fuelMassFlow
         << " kg/s\n"
         << "  mass-flow ramp=" << finalRamp << " (0.1 -> 1 over "
         << rampDuration << " s)\n"
         << "  calibrated axial U=" << injectionReport.torchAxialVelocity << " / "
         << injectionReport.loxAxialVelocity << " / " << injectionReport.fuelAxialVelocity
         << " m/s on faces " << injectionReport.torchFaces << " / "
         << injectionReport.loxFaces << " / " << injectionReport.fuelFaces << "\n";
    bool ok = !dom.tripped()&&!chemistryFailed&&!stepRejected;
    if(reacting)
    {
        const SectorSpeciesBalance speciesBalance=
            SectorDiagnostics(mesh).evaluateSpecies(dom,registry);
        const SectorPerformance performance=
            SectorPerformanceDiagnostics(mesh,registry).evaluate(dom,turbulence);
        cout << "  Pc/throat mdot/c*=" << performance.chamberPressure << " Pa / "
             << performance.throatMassFlow << " kg/s / " << performance.cStar
             << " m/s\n"
             << "  qwall max/avg/power=" << performance.wallHeatFluxMaximum
             << " / " << performance.wallHeatFluxAverage << " W/m2 / "
             << performance.wallHeatPower << " W (sector)\n";
        cout << "  max |sum(Y)-1|=" << speciesBalance.maximumMassFractionSumError
             << " element inventory C/H/O/N=" << speciesBalance.elementInventory[0]
             << " / " << speciesBalance.elementInventory[1] << " / "
             << speciesBalance.elementInventory[2] << " / "
             << speciesBalance.elementInventory[3] << " mol-atoms\n";
        for(size_t s=0;s<species.size();++s)
            cout << "  min Y_" << species[s] << "="
                 << speciesBalance.minimumMassFraction[s] << "\n";
        if(speciesBalance.maximumMassFractionSumError>1e-6) ok=false;
        for(double value:speciesBalance.minimumMassFraction)
            if(value<-1e-10) ok=false;
    }

    // Checks
    // rho positive everywhere
    double Tmin=1e30, Tmax=-1e30;
    double maxNut=0.0;
    for (int i=1;i<Nx-1;++i)
    for (int j=1;j<Ny-1;++j)
    for (int k=1;k<Nz-1;++k)
    {
        if (dom.U(i,j,k).rho <= 0) ok=false;
        Tmin = std::min(Tmin, dom.W(i,j,k).T);
        Tmax = std::max(Tmax, dom.W(i,j,k).T);
        if(turbulence) maxNut=std::max(maxNut,turbulence->turbulentViscosity(i,j,k,dom));
    }
    cout << "  interior T range=[" << Tmin << ", " << Tmax << "] K\n";
    if(turbulence) cout << "  max nu_t=" << maxNut << " m2/s\n";
#ifdef RFAM_HAVE_CHEMISTRY
    if(chemistryOperator)
        cout << "  chemistry active/failures=" << chemistryOperator->lastActiveCells()
             << " / " << chemistryOperator->lastFailures()
             << " retried=" << chemistryOperator->lastRetriedCells()
             << " positivity-limited="
             << chemistryOperator->lastPositivityLimitedCells()
             << " qchem min/max/avg=" << chemistryOperator->lastHeatReleaseMinimum()
             << " / " << chemistryOperator->lastHeatReleaseMaximum()
             << " / " << chemistryOperator->lastHeatReleaseAverage() << " W/m3\n";
#endif
    if (!reacting&&(std::abs(Tmin-700.0) > 1e-4 || std::abs(Tmax-700.0) > 1e-4)) ok=false;

    if (ok)
        cout << "  [PASS] 3D sector: no trip, rho>0\n";
    else
        cout << "  [FAIL] 3D sector failed\n";

    RFAM_ENTER(sector_finish);
    if (telemetry.active())
    {
        Foam::json::Value summary=Foam::json::Value::object();
        summary.set("ok",Foam::json::Value(ok));
        summary.set("max_steps",Foam::json::Value(maxSteps));
        telemetry.finish(ok?"pass":"failed",ok?0:1,summary);
    }
    if (!ok)
        cout << "\nSECTOR_3D_FAIL\n";
    else if (maxSteps == 0)
        cout << "\n" << (restartPath.empty()
            ? "SECTOR_3D_LOAD_OK" : "SECTOR_3D_RELOAD_OK") << "\n";
    else
        cout << "\nSECTOR_3D_STEP_PASS\n";
    return ok ? 0 : 1;
}
