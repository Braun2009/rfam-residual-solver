// RFAM tetrahedral cold-flow driver (first unstructured path).
//
// Scope: audited nozzle-only fluid mesh, mass-flow inlets by design zone,
// slip adiabatic walls, optional equilibrium pressure outlet, LU-SGS pseudo-
// time relaxation with conservative conserved-state recovery per cell.
// No combustion, diffusion or turbulence model is enabled here.
#include "TetraMesh.H"
#include "TetraFiniteVolume.H"
#include "TetraLUSGS.H"
#include "CompleteCellState.H"
#include "EquilibriumThermoModel.H"
#include "CoolPropThermoModel.H"
#include "TableAcoustic.H"
#include "ThermoCounters.H"
#include "StateAcceptance.H"
#include "TetraStateRecovery.H"
#include "TetraPressureOutlet.H"
#include "TetraInletAllocation.H"
#include "TetraBlueprintWall.H"
#include "TetraMolecularAdmission.H"
#include "TetraTurbulenceOperator.H"
#include "TransportProps.H"
#include "SpeciesRegistry.H"
#include "TetraSteadyMonitor.H"
#include "TetraPhysicalStepRetry.H"
#include "TetraPhysicalCheckpoint.H"
#include "ThermoDataProvenance.H"
#include "RunTelemetry.H"
#include "ModelCapabilities.H"
#ifdef RFAM_HAVE_COOLPROP
#include "CoolPropTransport.H"
#include "AbstractState.h"
#endif
#include "GERGResidualTable.H"
#include "GERGSaturationTable.H"
#include "FluxScheme.H"
#ifdef RFAM_HAVE_CHEMISTRY
#include "TwoStepCO.H"
#include "SandiegoChemistry.H"
#include "SandiegoReducedChemistry.H"
#include "SandiegoNasa7Data.H"
#include "TetraChemistryCVODE.H"
#endif
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <array>
#include <iterator>
#include <map>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef RFAM_USE_OPENMP
#include <omp.h>
#endif

using namespace Foam;

namespace
{

constexpr int kWallZone=0;
constexpr int kOutletZone=1;

struct ZoneInfo
{
    std::vector<std::size_t> faces;
    double area=0.0;
    double targetMassFlow=0.0;   // kg/s, positive means into the domain
    double currentMassFlow=0.0;  // kg/s, actual boundary flow this step (Cd model)
    double entropyMargin=std::numeric_limits<double>::quiet_NaN(); // s_face - s_supply
    bool entropyWarned=false;
    double supplyPressure=4.0e6;
    double supplyEnthalpy=std::numeric_limits<double>::quiet_NaN();
    double supplyEntropy=std::numeric_limits<double>::quiet_NaN();
    double chokedMassFlow=std::numeric_limits<double>::quiet_NaN(); // isentropic maximum
    double temperature=300.0;
    std::vector<double> composition;
    double pressure=0.0;         // interior average, rebuilt each step
    double blockedArea=0.0;      // inlet faces above supply pressure
    double availableCapacity=0.0;
    std::string name;
};

struct Config
{
    std::string mesh;
    std::string dump;
    int dumpEvery=-1;
    bool dumpPhaseEvery=false; // pair periodic fields with EOS phase and step metadata
    std::string seedField;
    bool burntInlet=false;          // impose burnt gas at the injector face
    double burntTemperature=-1.0;   // burnt-gas total temperature [K]
    bool inletCompositionSet=false;
    std::string physicalResume,physicalCheckpoint;
    std::string thermoBackend="coolprop";
    std::string inletMode="mass";
    double outletPressure=0.0; // 0: explicit zero-gradient diagnostic mode
    double backflowTemperature=300.0;
    std::vector<double> backflowY{1.0,0.0}; // explicit O2 reservoir; not ambient air
    int steps=100;
    int monitorEvery=10;
    double cfl=0.5;
    double p0=4.0e6;
    double T0=400.0;
    double pinj=4.0e6;   // fallback torch/chamber-side inlet pressure
    double loxPinj=6.0e6;
    double fuelPinj=4.5e6;
    double torchPinj=4.0e6;
    std::vector<double> background{0.784,0.216};
    double loxMassFlow=2.8;
    double fuelMassFlow=1.05;
    // Cold start before torch ignition. A nonzero central premixed feed needs
    // an explicit boundary temperature; reacting hot torch is a different case.
    double torchMassFlow=0.0;
    double loxTemperature=130.0;
    double fuelTemperature=450.0;
    double torchTemperature=std::numeric_limits<double>::quiet_NaN();
    bool controlMassFlow=true;
    bool blueprint=false;   // 1D continuity fill instead of stagnant start
    bool blueprintPlateWall=false; // opt-in initial velocity tangential to injector plate
    bool blueprintInjector=false;  // pre-fill injector channels with propellant states
    int blueprintLayers=8;         // connectivity layers primed from each inlet zone
    int blueprintBlend=0;          // outermost filled layers blended into ambient
    bool lusgs=false;       // implicit LU-SGS steady acceleration
    bool physicalTime=false; // conservative global physical-time cold-flow advance
    bool lowMachPreconditioning=false;  // Turkel-Weiss low-Mach preconditioning
    double lowMachCutoff=0.05;          // minimum reference Mach for preconditioning
    double physicalDt=0.0;  // requested physical step; 0 selects the current CFL diagnostic
    int physicalInnerSteps=30; // LU-SGS iterations per accepted physical step
    int physicalMaxRetries=9;
    bool physicalProbeInner=false;
    bool outletRoundoffShortcut=false;
    double physicalInnerTolerance=1e-3;
    bool interfaceConsistentFlux=false;  // equilibrium interface density in the flux
    bool molecularViscous=false; // opt-in, single-phase CoolProp molecular transport only
    std::string turbulenceModel;    // "" or "sst" (Menter k-omega SST, steady pseudo-time only)
    double turbulenceIntensity=0.05;
    double turbulenceLength=1e-3;
    double inletCd=0.0;     // >0: zone mass flow = Cd * isentropic capacity at face p
    bool supplyEntropyStrict=false; // false: report the inlet entropy margin, do not veto
    double supplyEntropyTolerance=1e-7; // relative entropy margin tolerance (allow losses)
    bool inletSupplyCap=false;      // true: mdot = min(design, ideal supply capacity at face p)
    bool inletFaceSupplyCap=false;  // opt-in local face-pressure/supply-capacity diagnostic
    bool lusgsCompositionCoupling=false; // opt-in expensive binary derivative audit
    double lusgsDtScale=2.0; // >0 adds a local pseudo-time term to the diagonal
    double lusgsDtScaleStart=0.0; // >0 ramps from this (smaller step) to lusgsDtScale
    int lusgsRampSteps=200;  // ramp length in iterations
    double lusgsRelax=0.7;   // correction relaxation (stability vs convergence)
    double twoPhaseDtScale=1.0; // local pseudo-time factor inside the two-phase band
    double injectorDtScale=1.0;  // pseudo-time factor inside the injector block (z<face)
    double guardTMin=80.0,guardTMax=4000.0;   // accepted temperature range
    double guardPMin=1.0e4,guardPMax=2.0e7;   // accepted pressure range
    int massFlowIterations=24;
    double massFlowTolerance=1e-3;
    // Run telemetry: on by default. The reacting path is the one where a
    // forgotten flag changes the chemistry, so the manifest records the
    // mechanism that was actually constructed, not the one that was requested.
    bool telemetry=true;
    std::string telemetryDir;
    int telemetryStallSeconds=30;
    bool strictNoDowngrade=false;
    // Reacting (hot fire) path: six-species C3H8/O2 chemistry with the
    // multi-species legacy (hybrid Helmholtz/PR) backend. The cold path is
    // untouched when --reacting is absent.
    bool reacting=false;
    // Chemistry mechanism: "twostep" (6-species diagnostic closure) or
    // "sandiego" (59-species, 271-reaction production mechanism).
    std::string mechanism="twostep";
    double chemSafety=0.5;
    int chemMaxSubsteps=64;
    bool chemCVODE=false;
    double ignitionT=0.0;
    double ignitionX=0.010;
    double ignitionR=0.004;
};

void usage()
{
    std::cout<<
        "usage: rfam_v20_tetra --mesh <file.rfamtet> [options]\n"
        "  --steps N              steady iterations or accepted physical steps (default 100)\n"
        "  --cfl C                local CFL number (default 0.5)\n"
        "  --monitor-every N      monitor interval (default 10)\n"
        "  --p0 P --T0 T          initial state (default 4e6 Pa, 400 K)\n"
        "  --thermo coolprop|legacy   default CoolProp 7.2.0 HEOS; no silent fallback\n"
        "  --inlet-mode mass|ghost   default prescribed mass + supply enthalpy; ghost diagnostic\n"
        "  --lox-pinj P --fuel-pinj P --torch-pinj P  supply stagnation pressures [Pa]\n"
        "  --pinj P                override all supply pressures [Pa]\n"
        "  --lusgs                 enable implicit stepping\n"
        "  --physical-time         global physical cold-flow steps; strict CoolProp mass inlet + pressure outlet\n"
        "                          uniform initial gas, torch off, inviscid; no diagnostic seed or relaxed EOS\n"
        "  --physical-dt DT       requested physical step [s]; implicit convergence decides acceptance\n"
        "  --physical-inner-iters N  LU-SGS iterations per physical step (default 30)\n"
        "  --physical-max-retries N  physical-step attempts (default 9)\n"
        "  --physical-probe-inner    log inner iteration time and gate metrics; diagnostic only\n"
        "  --outlet-roundoff-shortcut  diagnostic one-ULP matched-pressure outlet path\n"
        "  --physical-inner-tol R  relative and absolute physical-step residual tolerance (default 1e-3)\n"
        "  --physical-checkpoint FILE  write each accepted state as FILE.sN\n"
        "  --physical-resume FILE      resume an exact accepted physical checkpoint\n"
        "  --lusgs-composition-coupling  diagnostic binary pressure/composition Jacobian (extra EOS calls)\n"
        "  --outlet-p P           equilibrium pressure outlet [Pa]; default 0 (zero gradient diagnostic)\n"
        "  --backflow-T T --backflow-yO2 Y  explicit O2/C3H8 reservoir (default 300 K, pure O2; not air)\n"
        "  --yO2 Y --yC3H8 Y      initial mass fractions\n"
        "  --lox-mdot M --fuel-mdot M --torch-mdot M   inlet mass flows [kg/s]\n"
        "  --lox-T T --fuel-T T --torch-T T            inlet temperatures [K]\n"
        "  Torch defaults off; nonzero torch flow requires --torch-T (cold mixture only).\n"
        "  --no-mass-flow-control     use the constant-area velocity estimate only\n"
        "  --inlet-cd C           injector discharge coefficient: zone mdot follows\n"
        "                         Cd*isentropic capacity at the local face pressure\n"
        "  --supply-entropy-strict  veto mass-flow inlets whose face entropy is below\n"
        "                         the supply entropy (default: report only)\n"
        "  --supply-entropy-tolerance R  relative entropy margin tolerance (default 1e-7);\n"
        "                         allows s_face - s_supply >= -R*max(|s_supply|, 1) for losses\n"
        "  --inlet-supply-cap     cap the prescribed mdot at the ideal supply capacity\n"
        "                         at the local face pressure (no discharge coefficient)\n"
        "  --inlet-face-supply-cap  opt-in: cap each inlet face using its own interior pressure\n"
        "  --dump FILE            write final/last-accepted field and status\n"
        "  --dump-every N         also dump every N accepted steps (FILE.sN)\n"
        "  --dump-phase-every    pair each periodic field with EOS phase CSV and status JSON\n"
        "  --seed-field FILE      initialize from rounded primitive field (not exact restart)\n"
        "  --blueprint-plate-wall  remove normal velocity at the ambient injector-plate wall (initial condition only)\n"
        "  --molecular-viscous  opt-in, single-phase CoolProp viscosity/conduction with no-slip adiabatic walls\n"
        "  --no-telemetry       do not write the run telemetry directory\n"
        "  --telemetry-dir DIR  telemetry root (default runs/telemetry)\n"
        "  --telemetry-stall-seconds N  no-progress threshold for a STALL event (default 30)\n"
        "  --strict-no-downgrade  abort the run on any recorded downgrade\n"
    "  --reacting            finite-rate chemistry on the multi-species legacy EOS\n"
    "  --mechanism NAME      twostep (6 species, diagnostic) | sandiego (59 species, 271 rxn)\n"
    "                        multi-species legacy EOS; requires --thermo legacy\n"
    "  --chem-safety F       chemistry sub-step safety factor (default 0.5)\n"
    "  --chem-max-substeps N explicit diagnostic limit or BDF internal-step budget (default 64)\n"
    "  --chem-cvode         full-time San Diego SUNDIALS BDF chemistry (opt-in)\n"
    "  --ignition-T T        initialize a hot ignition kernel at --ignition-x/r\n"
    "  --ignition-x X --ignition-r R  kernel centre/radius [m] (default 0.01/0.004)\n";
}

bool parse(int argc,char**argv,Config& c)
{
    for(int i=1;i<argc;++i)
    {
        const std::string a=argv[i];
        auto need=[&](const char* name)->const char*
        {
            if(i+1>=argc)throw std::runtime_error(std::string(name)+" requires a value");
            return argv[++i];
        };
        if(a=="--help"||a=="-h"){usage();return false;}
        else if(a=="--mesh")c.mesh=need("--mesh");
        else if(a=="--thermo")c.thermoBackend=need("--thermo");
        else if(a=="--inlet-mode")c.inletMode=need("--inlet-mode");
        else if(a=="--outlet-p")c.outletPressure=std::atof(need("--outlet-p"));
        else if(a=="--backflow-T")c.backflowTemperature=std::atof(need("--backflow-T"));
        else if(a=="--backflow-yO2")
        {c.backflowY[0]=std::atof(need("--backflow-yO2"));c.backflowY[1]=1.-c.backflowY[0];}
        else if(a=="--steps")c.steps=std::atoi(need("--steps"));
        else if(a=="--cfl")c.cfl=std::atof(need("--cfl"));
        else if(a=="--monitor-every")c.monitorEvery=std::atoi(need("--monitor-every"));
        else if(a=="--p0")c.p0=std::atof(need("--p0"));
        else if(a=="--T0")c.T0=std::atof(need("--T0"));
        else if(a=="--pinj")
        {
            c.pinj=std::atof(need("--pinj"));
            c.loxPinj=c.fuelPinj=c.torchPinj=c.pinj;
        }
        else if(a=="--lox-pinj")c.loxPinj=std::atof(need("--lox-pinj"));
        else if(a=="--fuel-pinj")c.fuelPinj=std::atof(need("--fuel-pinj"));
        else if(a=="--torch-pinj")c.torchPinj=std::atof(need("--torch-pinj"));
        else if(a=="--yO2")c.background[0]=std::atof(need("--yO2"));
        else if(a=="--yC3H8")c.background[1]=std::atof(need("--yC3H8"));
        else if(a=="--lox-mdot")c.loxMassFlow=std::atof(need("--lox-mdot"));
        else if(a=="--fuel-mdot")c.fuelMassFlow=std::atof(need("--fuel-mdot"));
        else if(a=="--torch-mdot")c.torchMassFlow=std::atof(need("--torch-mdot"));
        else if(a=="--lox-T")c.loxTemperature=std::atof(need("--lox-T"));
        else if(a=="--fuel-T")c.fuelTemperature=std::atof(need("--fuel-T"));
        else if(a=="--torch-T")c.torchTemperature=std::atof(need("--torch-T"));
        else if(a=="--no-mass-flow-control")c.controlMassFlow=false;
        else if(a=="--inlet-cd")c.inletCd=std::atof(need("--inlet-cd"));
        else if(a=="--supply-entropy-strict")c.supplyEntropyStrict=true;
        else if(a=="--supply-entropy-tolerance")c.supplyEntropyTolerance=std::atof(need("--supply-entropy-tolerance"));
        else if(a=="--inlet-supply-cap")c.inletSupplyCap=true;
        else if(a=="--inlet-face-supply-cap")c.inletFaceSupplyCap=true;
        else if(a=="--blueprint")c.blueprint=true;
        else if(a=="--blueprint-plate-wall")c.blueprintPlateWall=true;
        else if(a=="--blueprint-injector")c.blueprintInjector=true;
        else if(a=="--blueprint-layers")c.blueprintLayers=std::atoi(need("--blueprint-layers"));
        else if(a=="--blueprint-blend")c.blueprintBlend=std::atoi(need("--blueprint-blend"));
        else if(a=="--lusgs")c.lusgs=true;
        else if(a=="--physical-time")c.physicalTime=true;
        else if(a=="--low-mach-preconditioning")c.lowMachPreconditioning=true;
        else if(a=="--low-mach-cutoff")c.lowMachCutoff=std::atof(need("--low-mach-cutoff"));
        else if(a=="--physical-dt")c.physicalDt=std::atof(need("--physical-dt"));
        else if(a=="--physical-inner-iters")c.physicalInnerSteps=std::atoi(need("--physical-inner-iters"));
        else if(a=="--physical-max-retries")c.physicalMaxRetries=std::atoi(need("--physical-max-retries"));
        else if(a=="--physical-probe-inner")c.physicalProbeInner=true;
        else if(a=="--outlet-roundoff-shortcut")c.outletRoundoffShortcut=true;
        else if(a=="--physical-inner-tol")c.physicalInnerTolerance=std::atof(need("--physical-inner-tol"));
        else if(a=="--interface-consistent-flux")c.interfaceConsistentFlux=true;
        else if(a=="--reacting")c.reacting=true;
        else if(a=="--mechanism")c.mechanism=need("--mechanism");
        else if(a=="--chem-safety")c.chemSafety=std::atof(need("--chem-safety"));
        else if(a=="--chem-max-substeps")c.chemMaxSubsteps=std::atoi(need("--chem-max-substeps"));
        else if(a=="--chem-cvode")c.chemCVODE=true;
        else if(a=="--ignition-T")c.ignitionT=std::atof(need("--ignition-T"));
        else if(a=="--ignition-x")c.ignitionX=std::atof(need("--ignition-x"));
        else if(a=="--ignition-r")c.ignitionR=std::atof(need("--ignition-r"));
        else if(a=="--molecular-viscous")c.molecularViscous=true;
        else if(a=="--no-telemetry")c.telemetry=false;
        else if(a=="--telemetry")c.telemetry=true;
        else if(a=="--telemetry-dir")c.telemetryDir=need("--telemetry-dir");
        else if(a=="--telemetry-stall-seconds")c.telemetryStallSeconds=std::atoi(need("--telemetry-stall-seconds"));
        else if(a=="--strict-no-downgrade")c.strictNoDowngrade=true;
        else if(a=="--turbulence")c.turbulenceModel=need("--turbulence");
        else if(a=="--turb-intensity")c.turbulenceIntensity=std::atof(need("--turb-intensity"));
        else if(a=="--turb-length")c.turbulenceLength=std::atof(need("--turb-length"));
        else if(a=="--lusgs-composition-coupling")c.lusgsCompositionCoupling=true;
        else if(a=="--lusgs-dt-scale")c.lusgsDtScale=std::atof(need("--lusgs-dt-scale"));
        else if(a=="--lusgs-dt-scale-start")c.lusgsDtScaleStart=std::atof(need("--lusgs-dt-scale-start"));
        else if(a=="--lusgs-ramp-steps")c.lusgsRampSteps=std::atoi(need("--lusgs-ramp-steps"));
        else if(a=="--two-phase-dt-scale")c.twoPhaseDtScale=std::atof(need("--two-phase-dt-scale"));
        else if(a=="--injector-dt-scale")c.injectorDtScale=std::atof(need("--injector-dt-scale"));
        else if(a=="--lusgs-relax")c.lusgsRelax=std::atof(need("--lusgs-relax"));
        else if(a=="--guard-T")c.guardTMin=std::atof(need("--guard-T"));
        else if(a=="--guard-T-max")c.guardTMax=std::atof(need("--guard-T-max"));
        else if(a=="--guard-p")c.guardPMin=std::atof(need("--guard-p"));
        else if(a=="--guard-p-max")c.guardPMax=std::atof(need("--guard-p-max"));
        else if(a=="--dump")c.dump=need("--dump");
        else if(a=="--dump-every")c.dumpEvery=std::atoi(need("--dump-every"));
        else if(a=="--dump-phase-every")c.dumpPhaseEvery=true;
        else if(a=="--seed-field")c.seedField=need("--seed-field");
        else if(a=="--burnt-inlet")c.burntInlet=true;
        else if(a=="--burnt-T")c.burntTemperature=std::atof(need("--burnt-T"));
        else if(a=="--physical-checkpoint")c.physicalCheckpoint=need("--physical-checkpoint");
        else if(a=="--physical-resume")c.physicalResume=need("--physical-resume");
        else throw std::runtime_error("unknown argument: "+a);
    }
    if(c.mesh.empty())throw std::runtime_error("--mesh is required");
    if(c.thermoBackend!="coolprop"&&c.thermoBackend!="legacy")
        throw std::runtime_error("--thermo must be coolprop or legacy");
    if(c.inletMode!="mass"&&c.inletMode!="ghost")throw std::runtime_error("inlet mode must be mass or ghost");
    if(c.inletMode=="mass"&&!c.controlMassFlow)throw std::runtime_error("--no-mass-flow-control only applies to ghost diagnostic mode");
    if(!std::isfinite(c.inletCd)||c.inletCd<0.||c.inletCd>2.)
        throw std::runtime_error("--inlet-cd must be in [0,2]");
    if(c.inletCd>0.&&(c.thermoBackend!="coolprop"||c.inletMode!="mass"))
        throw std::runtime_error("--inlet-cd requires --thermo coolprop --inlet-mode mass");
    if(c.inletSupplyCap&&((c.thermoBackend!="coolprop"&&!c.reacting)||c.inletMode!="mass"))
        throw std::runtime_error("--inlet-supply-cap requires --thermo coolprop (or --reacting legacy) --inlet-mode mass");
    if(c.inletCd>0.&&c.inletSupplyCap)
        throw std::runtime_error("--inlet-cd and --inlet-supply-cap are mutually exclusive");
    if(c.inletFaceSupplyCap&&(!c.inletSupplyCap||c.inletCd>0.||c.inletMode!="mass"||(c.thermoBackend!="coolprop"&&!c.reacting)))
        throw std::runtime_error("--inlet-face-supply-cap requires CoolProp mass inlet and --inlet-supply-cap, without Cd");
    if(c.molecularViscous&&c.thermoBackend!="coolprop"&&!c.reacting)
        throw std::runtime_error("--molecular-viscous requires the CoolProp backend"
            " (or --reacting for the Wilke/Sutherland gas closure)");
    if(!std::isfinite(c.outletPressure)||c.outletPressure<0.
       ||!std::isfinite(c.backflowTemperature)||c.backflowTemperature<=0.
       ||!std::isfinite(c.backflowY[0])||c.backflowY[0]<0.||c.backflowY[0]>1.)
        throw std::runtime_error("invalid pressure outlet/reservoir inputs");
    if(c.outletPressure>0.&&c.thermoBackend!="coolprop"&&!c.reacting)
        throw std::runtime_error("equilibrium pressure outlet requires CoolProp entropy support");
    if(c.reacting&&c.mechanism!="twostep"&&c.mechanism!="sandiego"
       &&c.mechanism!="sandiego-reduced")
        throw std::runtime_error("--mechanism must be twostep, sandiego or sandiego-reduced");
    if(c.reacting&&c.thermoBackend=="coolprop")
        throw std::runtime_error("--reacting requires --thermo legacy (multi-species backend)");
    if(c.reacting&&c.physicalTime&&c.thermoBackend!="legacy")
        throw std::runtime_error("--reacting with --physical-time requires --thermo legacy");
    if(!c.reacting&&c.ignitionT>0.)
        throw std::runtime_error("--ignition-T requires --reacting");
    if(c.reacting&&(!std::isfinite(c.chemSafety)||c.chemSafety<=0.||c.chemSafety>1.))
        throw std::runtime_error("--chem-safety must be in (0,1]");
    if(c.reacting&&(c.chemMaxSubsteps<1||c.chemMaxSubsteps>100000))
        throw std::runtime_error("--chem-max-substeps out of range");
    if(c.chemCVODE&&(!c.reacting||c.mechanism!="sandiego"))
        throw std::runtime_error("--chem-cvode requires --reacting --mechanism sandiego");
    if(c.steps<0||c.monitorEvery<=0||!std::isfinite(c.cfl)||c.cfl<=0.)
        throw std::runtime_error("invalid iteration controls");
    if(c.physicalTime&&c.lusgs)
        throw std::runtime_error("--physical-time and --lusgs are mutually exclusive");
    if(c.physicalTime&&(c.inletMode!="mass"||c.outletPressure<=0.))
    {
        // The pressure outlet needs entropy; the legacy multi-species model does
        // not provide it. A reacting physical-time run therefore uses the
        // zero-gradient outlet, which is the correct boundary once the nozzle
        // exit is supersonic but is not a Riemann invariant during startup.
        if(!(c.reacting&&c.thermoBackend=="legacy"&&c.outletPressure<=0.))
            throw std::runtime_error("physical-time requires mass-flow inlets and either a"
                " CoolProp pressure outlet or, for reacting legacy runs, the"
                " zero-gradient outlet (omit --outlet-p)");
    }
    if(c.outletRoundoffShortcut&&(!c.physicalTime||c.outletPressure<=0.))
        throw std::runtime_error("outlet roundoff shortcut requires physical-time pressure outlet");
    if(c.physicalTime&&(c.cfl>0.5||c.guardPMax!=2e7))
        throw std::runtime_error("physical-time requires CFL <= 0.5 and the 20 MPa pressure guard");
    if(c.physicalTime&&(c.inletCd>0.||c.inletSupplyCap||c.inletFaceSupplyCap))
        throw std::runtime_error("physical-time design-flow startup does not admit Cd or supply-cap inlet diagnostics");
    if(c.physicalTime&&(c.blueprint||c.blueprintInjector||c.torchMassFlow!=0.))
        throw std::runtime_error("physical-time startup requires a uniform initial gas and torch off");
    if(c.physicalTime&&!c.seedField.empty())
        throw std::runtime_error("physical-time cannot start from a rounded diagnostic seed");
    if((!c.physicalResume.empty()||!c.physicalCheckpoint.empty())&&!c.physicalTime)
        throw std::runtime_error("physical checkpoint requires --physical-time");
    if(c.physicalTime&&(!std::isfinite(c.physicalDt)||c.physicalDt<0.
       ||c.physicalInnerSteps<=0||c.physicalMaxRetries<1||c.physicalMaxRetries>9
       ||!std::isfinite(c.physicalInnerTolerance)
       ||c.physicalInnerTolerance<=0.||c.physicalInnerTolerance>=1.))
        throw std::runtime_error("invalid physical-time controls");
    for(double v:{c.p0,c.T0,c.loxPinj,c.fuelPinj,c.torchPinj,c.loxTemperature,c.fuelTemperature})
        if(!std::isfinite(v)||v<=0.)throw std::runtime_error("invalid inlet/initial thermodynamic inputs");
    for(double v:{c.loxMassFlow,c.fuelMassFlow,c.torchMassFlow})
        if(!std::isfinite(v)||v<0.)throw std::runtime_error("invalid inlet mass flow");
    if(c.torchMassFlow>0.&&(!std::isfinite(c.torchTemperature)||c.torchTemperature<=0.))
        throw std::runtime_error("nonzero torch feed requires explicit --torch-T > 0 K");
    for(double y:c.background)if(!std::isfinite(y)||y<0.||y>1.)throw std::runtime_error("invalid initial composition");
    if(std::abs(c.background[0]+c.background[1]-1.)>1e-12)throw std::runtime_error("initial mass fractions must sum to one");
    if(c.lusgsCompositionCoupling&&!c.lusgs)
        throw std::runtime_error("--lusgs-composition-coupling requires --lusgs");
    if(!c.lusgs&&!c.physicalTime&&c.steps>0)
        throw std::runtime_error("use --lusgs or --physical-time");
    if(c.blueprintPlateWall&&(!c.blueprint||!c.seedField.empty()))
        throw std::runtime_error("--blueprint-plate-wall requires --blueprint without --seed-field");
    if(c.dumpPhaseEvery&&(c.dump.empty()||c.dumpEvery<=0))
        throw std::runtime_error("--dump-phase-every requires --dump and --dump-every N>0");
    return true;
}


#ifdef RFAM_HAVE_CHEMISTRY
// ----------------------------------------------------------------------- //
// Operator-split finite-rate chemistry on the tetra state (TwoStepCO).
//
//   d(rhoY_s)/dt = omega_s [kg/m3/s] ,  d(rhoE)/dt = qdot [W/m3]
//
// The update is explicit with a per-cell pseudo-time limited by the local
// reaction timescale (species consumption time and a temperature-rise cap).
// Cells whose re-flash fails keep their previous state: the chemistry never
// publishes a state the EOS has not committed (fail closed, counted).
// Two-phase cells are skipped (evaporation, not finite-rate chemistry,
// owns those states in this solver).
// ----------------------------------------------------------------------- //
struct ChemistryReport
{
    long long cells=0,active=0,limited=0,substeps=0,skippedTwoPhase=0,
              skippedCold=0,failed=0;
    long long failSource=0,failNonfinite=0,failDt=0,failRecover=0,failCommit=0;
    double maxHeatRelease=0.0,maxOmega=0.0;
    double minT=std::numeric_limits<double>::quiet_NaN(),
           maxT=std::numeric_limits<double>::quiet_NaN();
    double heatReleaseTotal=0.0;
    double heatReleaseEnergyJ=0.0;
    long long rhsCalls=0;
    double maxStateChange=0.0,maxElementDefect=0.0,maxEnergyDefect=0.0;
    double minTimeFraction=0.0;
};

std::string chemistryDiagnosticNumber(double value)
{
    std::ostringstream text;
    text<<std::scientific<<std::setprecision(10)<<value;
    return text.str();
}

void applyChemistryCVODE(
    TetraFiniteVolume& fv,const thermo::EquilibriumThermoModel& model,
    const Foam::ChemistryModel& chem,const std::vector<double>& localDt,
    double dtScale,const Config& cfg,ChemistryReport& report)
{
    const std::size_t cells=fv.cellCount();
    if(localDt.size()!=cells||!(dtScale>0.)||!std::isfinite(dtScale))
        throw std::invalid_argument("CVODE chemistry timestep dimension/value mismatch");
    int workers=1;
#ifdef RFAM_USE_OPENMP
    workers=omp_get_max_threads();
#endif
    std::vector<std::unique_ptr<TetraChemistryCVODE>> integrators;
    for(int i=0;i<workers;++i)
        integrators.emplace_back(std::make_unique<TetraChemistryCVODE>(chem,model));
    std::vector<std::optional<CompleteCellState>> staged(cells);
    long long active=0,skippedCold=0,skippedTwoPhase=0,failed=0,substeps=0,rhsCalls=0;
    double maxHeat=0.,heatTotal=0.,energyTotal=0.,minT=1e300,maxT=-1e300,
        maxChange=0.,maxElement=0.,maxEnergy=0.,minTime=1.;
    std::size_t bad=cells;std::string detail;
    std::atomic<bool> abortStage{false};
    #pragma omp parallel for schedule(dynamic,32) reduction(+:active,skippedCold,skippedTwoPhase,failed,substeps,rhsCalls,heatTotal,energyTotal) reduction(max:maxHeat,maxT,maxChange,maxElement,maxEnergy) reduction(min:minT,minTime)
    for(std::size_t c=0;c<cells;++c)
    {
        if(abortStage.load(std::memory_order_relaxed))continue;
        const auto& q=fv.Q(c);const auto& w=fv.W(c);const auto& u=fv.U(c);
        if(q.T<500.){++skippedCold;continue;}
        if(q.region==thermo::EquilibriumRegion::TWO_PHASE
           ||q.region==thermo::EquilibriumRegion::SINGLE_LIQUID
           ||q.region==thermo::EquilibriumRegion::CONDENSED_SOLID)
        {++skippedTwoPhase;continue;}
        try
        {
            int worker=0;
#ifdef RFAM_USE_OPENMP
            worker=omp_get_thread_num();
#endif
            const double requested=localDt[c]*dtScale;
            const auto candidate=integrators[worker]->advance(w,u,q,requested,cfg.chemMaxSubsteps);
            substeps+=candidate.steps;rhsCalls+=candidate.rhsCalls;
            if(!candidate.success())throw std::runtime_error(candidate.reason);
            const auto& end=*candidate.state;
            if(end.Q.T<cfg.guardTMin||end.Q.T>cfg.guardTMax
               ||end.Q.p<cfg.guardPMin||end.Q.p>cfg.guardPMax)
                throw std::runtime_error("CVODE final physical range guard");
            std::vector<double> omega;double heat=0.;
            chem.computeSourceTerms(end.Q.T,end.Q.p,end.U.rho,end.W.Y,omega,heat);
            if(!std::isfinite(heat))throw std::runtime_error("CVODE nonfinite final heat source");
            staged[c]=end;++active;
            heatTotal+=heat*fv.mesh().volume(c);
            energyTotal+=(end.U.rhoE-u.rhoE)*fv.mesh().volume(c);
            maxHeat=std::max(maxHeat,heat);minT=std::min(minT,end.Q.T);maxT=std::max(maxT,end.Q.T);
            minTime=std::min(minTime,candidate.reachedTime/requested);
            maxElement=std::max(maxElement,candidate.elementDefect);
            maxEnergy=std::max(maxEnergy,candidate.energyDefect);
            maxChange=std::max(maxChange,std::abs(end.Q.T-q.T)/q.T);
            maxChange=std::max(maxChange,std::abs(end.Q.p-q.p)/q.p);
            for(std::size_t s=0;s<w.Y.size();++s)
                maxChange=std::max(maxChange,std::abs(end.W.Y[s]-w.Y[s]));
        }
        catch(const std::exception& error)
        {
            ++failed;
            abortStage.store(true,std::memory_order_relaxed);
            #pragma omp critical(rfam_cvode_stage_failure)
            {
                if(c<bad){bad=c;detail=error.what();}
                std::cerr<<std::setprecision(17)<<"CVODE_CELL_REJECTED cell="<<c
                    <<" p="<<q.p<<" T="<<q.T<<" rho="<<u.rho<<" e="<<q.e
                    <<" requested_dt="<<localDt[c]*dtScale<<" reason="<<error.what()<<" Y=";
                for(double y:w.Y)std::cerr<<y<<",";
                std::cerr<<std::endl;
            }
        }
    }
    report=ChemistryReport{};report.cells=cells;report.active=active;
    report.failed=failed;report.failRecover=failed;report.substeps=substeps;report.rhsCalls=rhsCalls;
    report.skippedCold=skippedCold;report.skippedTwoPhase=skippedTwoPhase;
    report.maxHeatRelease=maxHeat;report.heatReleaseTotal=heatTotal;report.heatReleaseEnergyJ=energyTotal;
    report.minT=active?minT:std::numeric_limits<double>::quiet_NaN();
    report.maxT=active?maxT:std::numeric_limits<double>::quiet_NaN();
    report.minTimeFraction=active?minTime:0.;report.maxStateChange=maxChange;
    report.maxElementDefect=maxElement;report.maxEnergyDefect=maxEnergy;
    if(bad<cells)throw std::runtime_error("CVODE chemistry stage rejected at cell "+std::to_string(bad)+": "+detail);
    for(std::size_t c=0;c<cells;++c)
        if(staged[c])fv.exchangeState(c,staged[c]->W,staged[c]->U,staged[c]->Q);
}

void applyChemistrySources(
    TetraFiniteVolume& fv,
    const thermo::EquilibriumThermoModel& model,
    const Foam::ChemistryModel& chem,
    const std::vector<double>& localDt,
    double dtScale,
    const Config& cfg,
    ChemistryReport& report)
{
    const std::size_t cells=fv.cellCount();
    const std::size_t speciesCount=fv.species().size();
    long long active=0,limited=0,substeps=0,skippedTwoPhase=0,skippedCold=0,failed=0;
    long long failSource=0,failNonfinite=0,failDt=0,failRecover=0,failCommit=0;
    double maxHeat=0.0,maxOmega=0.0,heatTotal=0.0,energyTotal=0.0,
           minT=1e300,maxT=-1e300;

    #pragma omp parallel for schedule(static) reduction(+:active,limited,substeps,skippedTwoPhase,skippedCold,failed,failSource,failNonfinite,failDt,failRecover,failCommit,heatTotal,energyTotal) reduction(max:maxHeat,maxOmega,maxT) reduction(min:minT)
    for(std::size_t c=0;c<cells;++c)
    {
        const PrimitiveVars w=fv.W(c);
        const ConservedVars u=fv.U(c);
        const thermo::EquilibriumState q=fv.Q(c);
        if(q.region==thermo::EquilibriumRegion::TWO_PHASE){++skippedTwoPhase;continue;}
        // Cold screen: below 500 K the Arrhenius rate is numerically zero.
        if(q.T<500.0){++skippedCold;continue;}

        std::vector<double> omega;
        double qdot=0.0;
        try
        {
            chem.computeSourceTerms(q.T,q.p,u.rho,w.Y,omega,qdot);
        }
        catch(const std::exception&)
        {
            ++failed;++failSource;continue;
        }
        bool finite=true;
        for(double value:omega)if(!std::isfinite(value)){finite=false;break;}
        if(!finite||!std::isfinite(qdot)){++failed;++failNonfinite;continue;}

        const double totalDt=localDt[c]*dtScale;
        if(!(totalDt>0.0)){++failed;++failDt;continue;}

        double dtLimit=totalDt;
        for(std::size_t sp=0;sp<omega.size()&&sp<w.Y.size();++sp)
        {
            const double source=std::abs(omega[sp]);
            if(source<=1e-12)continue;
            const double available=u.rho*std::max(w.Y[sp],1e-12);
            dtLimit=std::min(dtLimit,cfg.chemSafety*available/source);
        }
        if(qdot>1e-6&&q.thermodynamics)
        {
            const double cv=q.thermodynamics->Cv;
            if(std::isfinite(cv)&&cv>0.0)
                dtLimit=std::min(dtLimit,cfg.chemSafety*250.0*u.rho*cv/qdot);
        }
        const bool wasLimited=dtLimit<totalDt;
        double effectiveDt=std::min(totalDt,dtLimit);
        if(!(effectiveDt>0.0))effectiveDt=totalDt;
        const double ratio=totalDt/std::max(dtLimit,1e-30);
        const long long sub=static_cast<long long>(
            std::ceil(std::min(ratio,double(cfg.chemMaxSubsteps))));

        std::vector<double> rhoY=u.rhoY;
        for(std::size_t sp=0;sp<rhoY.size()&&sp<omega.size();++sp)
            rhoY[sp]+=effectiveDt*omega[sp];
        double rhoE=u.rhoE+effectiveDt*qdot;

        double sumSpecies=0.0;
        for(double& value:rhoY){value=std::max(value,0.0);sumSpecies+=value;}
        if(!(sumSpecies>0.0)||u.rho<=0.0){++failed;continue;}
        if(std::abs(sumSpecies-u.rho)>1e-14*u.rho)
        {
            const double factor=u.rho/sumSpecies;
            for(double& value:rhoY)value*=factor;
        }

        const double kinetic=0.5*u.rho*(w.U[0]*w.U[0]+w.U[1]*w.U[1]+w.U[2]*w.U[2]);
        const double internalEnergy=(rhoE-kinetic)/u.rho;
        if(!std::isfinite(internalEnergy)){++failed;continue;}

        // The EOS takes NORMALISED mass fractions; rhoY carries species
        // densities (sum = rho). Convert before the flash (2026-10-06: this
        // was the INVALID_INPUT failure that rejected every cell).
        std::vector<double> Ynorm(rhoY.size(),0.0);
        for(std::size_t sp=0;sp<rhoY.size();++sp)
            Ynorm[sp]=rhoY[sp]/u.rho;

        thermo::EquilibriumSeed seed;
        seed.pressure=q.p;
        seed.temperature=q.T;
        seed.previousTemperature=q.T;
        seed.previousRegion=q.region;
        const auto redone=model.recoverRhoE(u.rho,internalEnergy,Ynorm,&seed);
        if(!redone.success())
        {
            ++failed;++failRecover;
            if(failRecover<=3)
                #pragma omp critical(rfam_chem_diag)
                std::cerr<<"CHEM_RECOVER_FAIL cell="<<c
                         <<" p="<<q.p<<" T="<<q.T<<" rho="<<u.rho
                         <<" e="<<internalEnergy
                         <<" Y0="<<rhoY[0]<<" Y1="<<rhoY[1]<<" Y3="<<rhoY[3]
                         <<" status="<<redone.diagnostics.backendStatus<<std::endl;
            continue;
        }

        PrimitiveVars wNew(speciesCount);
        wNew.p=redone.state->p;
        wNew.T=redone.state->T;
        wNew.U=w.U;
        wNew.Y=Ynorm;
        wNew.Z=w.Z;
        wNew.C=w.C;
        const auto committed=CompleteStateConversion::fromPrimitiveAndEquilibrium(wNew,*redone.state);
        if(!committed.success())
        {
            ++failed;++failCommit;
            if(failCommit<=3)
                #pragma omp critical(rfam_chem_diag)
                std::cerr<<"CHEM_COMMIT_FAIL cell="<<c
                         <<" p="<<wNew.p<<" T="<<wNew.T
                         <<" status="<<committed.diagnostics.backendStatus<<std::endl;
            continue;
        }

        fv.setState(c,committed.state->W,committed.state->U,committed.state->Q);
        ++active;
        if(wasLimited)++limited;
        substeps+=sub;
        heatTotal+=qdot*fv.mesh().volume(c);
        energyTotal+=qdot*effectiveDt*fv.mesh().volume(c);
        maxHeat=std::max(maxHeat,qdot);
        for(double value:omega)maxOmega=std::max(maxOmega,std::abs(value));
        minT=std::min(minT,redone.state->T);
        maxT=std::max(maxT,redone.state->T);
    }

    report.cells=static_cast<long long>(cells);
    report.active=active;
    report.limited=limited;
    report.substeps=substeps;
    report.skippedTwoPhase=skippedTwoPhase;
    report.skippedCold=skippedCold;
    report.failed=failed;
    report.failSource=failSource;
    report.failNonfinite=failNonfinite;
    report.failDt=failDt;
    report.failRecover=failRecover;
    report.failCommit=failCommit;
    report.maxHeatRelease=maxHeat;
    report.maxOmega=maxOmega;
    report.minT=active>0?minT:std::numeric_limits<double>::quiet_NaN();
    report.maxT=active>0?maxT:std::numeric_limits<double>::quiet_NaN();
    report.heatReleaseTotal=heatTotal;
    report.heatReleaseEnergyJ=energyTotal;
}
#endif

std::string tetraMeshFingerprint(const std::string& path)
{
    std::ifstream input(path,std::ios::binary);
    if(!input)throw std::runtime_error("cannot hash tetra mesh");
    const std::string bytes(std::istreambuf_iterator<char>{input},{});
    if(input.bad())throw std::runtime_error("tetra mesh hash read failed");
    return thermoProvenance::sha256Hex(bytes);
}

std::string physicalBoundaryIdentity(const Config& cfg)
{
    std::string material;
    const auto number=[&](const char* key,double value)
    {thermoProvenance::appendDouble(material,key,value);};
    const auto flag=[&](const char* key,bool value)
    {thermoProvenance::appendSigned(material,key,value?1:0);};
    number("initial-p",cfg.p0);number("initial-T",cfg.T0);
    number("initial-Y-O2",cfg.background[0]);number("initial-Y-C3H8",cfg.background[1]);
    number("lox-supply-p",cfg.loxPinj);number("lox-supply-T",cfg.loxTemperature);
    number("fuel-supply-p",cfg.fuelPinj);number("fuel-supply-T",cfg.fuelTemperature);
    number("lox-mdot",cfg.loxMassFlow);number("fuel-mdot",cfg.fuelMassFlow);
    number("torch-mdot",cfg.torchMassFlow);number("outlet-p",cfg.outletPressure);
    number("backflow-T",cfg.backflowTemperature);number("backflow-Y-O2",cfg.backflowY[0]);
    number("physical-dt",cfg.physicalDt);number("physical-inner-tol",cfg.physicalInnerTolerance);
    number("cfl",cfg.cfl);number("lusgs-relax",cfg.lusgsRelax);
    number("two-phase-scale",cfg.twoPhaseDtScale);number("injector-scale",cfg.injectorDtScale);
    number("guard-T-min",cfg.guardTMin);number("guard-T-max",cfg.guardTMax);
    number("guard-p-min",cfg.guardPMin);number("guard-p-max",cfg.guardPMax);
    flag("outlet-roundoff-shortcut",cfg.outletRoundoffShortcut);
    flag("interface-consistent-flux",cfg.interfaceConsistentFlux);
    thermoProvenance::appendString(material,"inlet-mode",cfg.inletMode);
    thermoProvenance::appendString(material,"thermo-backend",cfg.thermoBackend);
    return thermoProvenance::sha256Hex(material);
}

} // namespace

int main(int argc,char**argv)
{
    Config cfg;
    try{if(!parse(argc,argv,cfg))return 2;}
    catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 2;}

    Foam::telemetry::Options telemetryOptions;
    telemetryOptions.enabled=cfg.telemetry;
    telemetryOptions.rootDir=cfg.telemetryDir;
    telemetryOptions.stallSeconds=cfg.telemetryStallSeconds;
    telemetryOptions.strictNoDowngrade=cfg.strictNoDowngrade;
    Foam::telemetry::Session& telemetry=Foam::telemetry::Session::instance();
    telemetry.start(argc,argv,telemetryOptions);
    {
        auto field=[&](const char* key,const std::string& value,const char* flag=nullptr)
        {telemetry.defineField(key,value,flag);};
        auto number=[&](const char* key,double value,const char* flag=nullptr)
        {field(key,std::to_string(value),flag);};
        auto integer=[&](const char* key,long long value,const char* flag=nullptr)
        {field(key,std::to_string(value),flag);};
        auto flagOf=[&](const char* key,bool value)
        {field(key,value?"true":"false");};
        field("mesh",cfg.mesh);
        field("thermo",cfg.thermoBackend);
        field("inlet_mode",cfg.inletMode);
        // The mechanism is the field whose absence from the old status file made
        // "which chemistry did this run use" unanswerable after the fact.
        field("mechanism",cfg.mechanism);
        flagOf("reacting",cfg.reacting);
        flagOf("chem_cvode",cfg.chemCVODE);
        integer("steps",cfg.steps);
        number("cfl",cfg.cfl);
        integer("monitor_every",cfg.monitorEvery);
        number("p0",cfg.p0);
        number("T0",cfg.T0);
        number("lox_pinj",cfg.loxPinj);
        number("fuel_pinj",cfg.fuelPinj);
        number("torch_pinj",cfg.torchPinj);
        number("lox_mdot",cfg.loxMassFlow);
        number("fuel_mdot",cfg.fuelMassFlow);
        number("torch_mdot",cfg.torchMassFlow);
        number("lox_T",cfg.loxTemperature);
        number("fuel_T",cfg.fuelTemperature);
        number("torch_T",cfg.torchTemperature);
        number("outlet_p",cfg.outletPressure);
        number("inlet_cd",cfg.inletCd);
        flagOf("molecular_viscous",cfg.molecularViscous);
        field("turbulence",cfg.turbulenceModel);
        number("turb_intensity",cfg.turbulenceIntensity);
        number("turb_length",cfg.turbulenceLength);
        flagOf("lusgs",cfg.lusgs);
        number("lusgs_dt_scale",cfg.lusgsDtScale);
        number("lusgs_relax",cfg.lusgsRelax);
        number("two_phase_dt_scale",cfg.twoPhaseDtScale);
        number("injector_dt_scale",cfg.injectorDtScale);
        flagOf("physical_time",cfg.physicalTime);
        number("physical_dt",cfg.physicalDt);
        integer("physical_inner_iters",cfg.physicalInnerSteps);
        integer("physical_max_retries",cfg.physicalMaxRetries);
        flagOf("low_mach_preconditioning",cfg.lowMachPreconditioning);
        number("low_mach_cutoff",cfg.lowMachCutoff);
        flagOf("outlet_roundoff_shortcut",cfg.outletRoundoffShortcut);
        flagOf("inlet_supply_cap",cfg.inletSupplyCap);
        flagOf("inlet_face_supply_cap",cfg.inletFaceSupplyCap);
        flagOf("supply_entropy_strict",cfg.supplyEntropyStrict);
        flagOf("burnt_inlet",cfg.burntInlet);
        number("burnt_T",cfg.burntTemperature);
        flagOf("lusgs_composition_coupling",cfg.lusgsCompositionCoupling);
        field("dump",cfg.dump);
        field("physical_resume",cfg.physicalResume);
        telemetry.flushManifest();
    }
    RFAM_PHASE(react_setup);
    RFAM_PHASE(react_mesh_load);
    RFAM_PHASE(react_thermo_setup);
    RFAM_PHASE(react_chemistry_setup);
    RFAM_PHASE(react_wall_distance);
    RFAM_PHASE(react_initial);
    RFAM_PHASE(react_step);
    RFAM_PHASE(react_chemistry_stage);
    RFAM_PHASE(react_monitor);
    RFAM_PHASE(react_dump);
    RFAM_PHASE(react_finish);
    RFAM_ENTER(react_setup);

#ifdef RFAM_USE_OPENMP
    std::cout<<"openmp_threads="<<omp_get_max_threads()<<"\n";
#endif
    try
    {
        if(cfg.physicalTime)
        {
            if(Foam::projectionAcceptEnabled()
               ||Foam::stateAcceptanceRel()!=Foam::STATE_ACCEPTANCE_REL_DEFAULT
               ||Foam::speciesRepairRel()!=Foam::SPECIES_REPAIR_REL_DEFAULT)
                throw std::runtime_error("physical-time requires strict default state recovery (no projection or relaxed acceptance)");
            cfg.supplyEntropyStrict=true;
        }
        RFAM_ENTER(react_mesh_load);
        TetraMesh mesh;mesh.read(cfg.mesh);
        RFAM_LEAVE(react_mesh_load);
        TetraPhysicalCheckpoint::Metadata checkpointIdentity;
        if(!cfg.physicalResume.empty()||!cfg.physicalCheckpoint.empty())
        {
            checkpointIdentity.meshSha256=tetraMeshFingerprint(cfg.mesh);
            std::string identity=thermoProvenance::sha256Hex(
                physicalBoundaryIdentity(cfg)+tetraMeshFingerprint(argv[0]));
            const char* tablePath=std::getenv("RFAM_TETRA_TABLE");
            if(tablePath&&*tablePath)
                identity=thermoProvenance::sha256Hex(identity+tetraMeshFingerprint(tablePath));
            else
                identity=thermoProvenance::sha256Hex(identity+"no-property-table");
            checkpointIdentity.boundaryIdentity=identity;
        }
        std::cout<<"mesh cells="<<mesh.cellCount()<<" faces="<<mesh.faceCount()
                 <<" points="<<mesh.pointCount()<<"\n";

        // Mesh scale diagnostics. Small cells are never frozen to conceal a
        // failed thermodynamic recovery or make a residual look converged.
        std::vector<double> volumes;
        volumes.reserve(mesh.cellCount());
        for(std::size_t c=0;c<mesh.cellCount();++c)volumes.push_back(mesh.volume(c));
        std::sort(volumes.begin(),volumes.end());
        const double medianVolume=volumes[volumes.size()/2];
        const double sliverVolume=1e-3*medianVolume;
        std::size_t sliverCount=0;
        double sliverVolumeSum=0.0;
        for(double v:volumes)if(v<sliverVolume){++sliverCount;sliverVolumeSum+=v;}
        std::cout<<"volume statistics: min="<<volumes.front()<<" median="<<medianVolume
                 <<" max="<<volumes.back()
                 <<" sliver_threshold="<<sliverVolume
                 <<" sliver_cells="<<sliverCount
                 <<" sliver_volume_fraction="<<sliverVolumeSum/3.0e-3
                 <<std::endl;

        const std::vector<std::string> coldSpecies{"O2","C3H8"};
        // Reacting order MUST match the selected mechanism exactly:
        //   twostep  -> Foam::TwoStepCO order (C3H8,O2,CO,CO2,H2O,N2)
        //   sandiego -> Foam::SandiegoChemistry order (SandiegoNasa7Data)
        const std::vector<std::string> twoStepSpecies{"C3H8","O2","CO","CO2","H2O","N2"};
        std::vector<std::string> sandiegoSpecies;
        for(int i=0;i<SandiegoNasa7::nSpecies;++i)
            sandiegoSpecies.push_back(SandiegoNasa7::species[i].name);
        const std::vector<std::string> reducedSpecies=
            Foam::SandiegoReducedChemistry::defaultSubset();
        const std::vector<std::string>& reactingSpecies=
            (cfg.mechanism=="sandiego")?sandiegoSpecies
            :(cfg.mechanism=="sandiego-reduced")?reducedSpecies
            :twoStepSpecies;
        const std::vector<std::string>& species=
            cfg.reacting?reactingSpecies:coldSpecies;
#ifndef RFAM_HAVE_CHEMISTRY
        if(cfg.reacting)throw std::runtime_error(
            "--reacting requires a build with RFAM_HAVE_CHEMISTRY");
#endif
        RFAM_ENTER(react_thermo_setup);
        auto model=cfg.thermoBackend=="coolprop"?thermo::makeCoolPropThermoModel():
            thermo::makeLegacySinglePhaseThermoModel(species);
        // Species-order composition mapper: the CLI keeps the two-component
        // (O2,C3H8) vocabulary; reacting runs expand it into the chemistry order.
        const auto speciesIndexOf=[&](const std::string& name)->int
        {
            for(std::size_t i=0;i<species.size();++i)
                if(species[i]==name)return static_cast<int>(i);
            return -1;
        };
        const int iO2=speciesIndexOf("O2");
        const int iC3H8=speciesIndexOf("C3H8");
        if(cfg.reacting&&(iO2<0||iC3H8<0))
            throw std::runtime_error("reacting species basis lacks O2/C3H8");
        auto mixtureY=[&](double yO2,double yC3H8)->std::vector<double>
        {
            if(!cfg.reacting)return std::vector<double>{yO2,yC3H8};
            std::vector<double> y(species.size(),0.0);
            y[iC3H8]=yC3H8;y[iO2]=yO2;
            return y;
        };
        // Chamber background: reacting runs start from inert combustion
        // products (Fluent pilot precedent) so only the injected propellant
        // can burn; the cold path keeps the torch-blend background.
        std::vector<double> backgroundY=mixtureY(cfg.background[0],cfg.background[1]);
        if(cfg.reacting)
        {
            backgroundY.assign(species.size(),0.0);
            backgroundY[speciesIndexOf("CO2")]=0.6468;
            backgroundY[speciesIndexOf("H2O")]=0.3532;
        }
        const std::vector<double> backflowSpeciesY=
            mixtureY(cfg.backflowY[0],cfg.backflowY[1]);
        std::cout<<"thermo_model="<<model->modelId()<<"\n";
        Foam::telemetry::Session::instance().describe().set(
            "thermo_model_id",Foam::json::Value(model->modelId()));
        RFAM_LEAVE(react_thermo_setup);
#ifdef RFAM_HAVE_CHEMISTRY
        std::unique_ptr<Foam::ChemistryModel> chemOperator;
        RFAM_ENTER(react_chemistry_setup);
        if(cfg.reacting)
        {
            if(cfg.mechanism=="sandiego")
                chemOperator.reset(new Foam::SandiegoChemistry());
            else if(cfg.mechanism=="sandiego-reduced")
                chemOperator.reset(new Foam::SandiegoReducedChemistry(
                    reactingSpecies));
            else
                chemOperator.reset(new Foam::TwoStepCO());
            chemOperator->print();
            std::cout<<"chemistry_mechanism="<<cfg.mechanism
                     <<" species="<<chemOperator->nSpecies()<<std::endl;
            if(telemetry.active())
            {
                // Recorded from the constructed operator, not from the flag:
                // the two can differ if a later change makes the flag advisory.
                Foam::json::Value chemistry=Foam::json::Value::object();
                chemistry.set("mechanism",Foam::json::Value(cfg.mechanism));
                chemistry.set("integrator",Foam::json::Value(cfg.chemCVODE?"cvode":"explicit"));
                chemistry.set("n_species",Foam::json::Value(static_cast<long long>(chemOperator->nSpecies())));
                Foam::json::Value speciesList=Foam::json::Value::array();
                for(const auto& name:species)speciesList.push(Foam::json::Value(name));
                chemistry.set("species",speciesList);
                Foam::json::Value& manifest=telemetry.describe();
                manifest.set("chemistry",chemistry);
                manifest.set("thermo_model_id",Foam::json::Value(model->modelId()));
                Foam::json::Value capabilities=Foam::json::Value::object();
                capabilities.set("twoPhaseFlash",Foam::json::Value(Foam::ModelCapabilities::twoPhaseFlash));
                capabilities.set("phaseAwareChemistry",Foam::json::Value(Foam::ModelCapabilities::phaseAwareChemistry));
                capabilities.set("phaseAwareSoundSpeed",Foam::json::Value(Foam::ModelCapabilities::phaseAwareSoundSpeed));
                capabilities.set("ransTetraProductionReady",Foam::json::Value(Foam::ModelCapabilities::ransTetraProductionReady()));
                manifest.set("capabilities",capabilities);
                Foam::json::Value meshInfo=Foam::json::Value::object();
                meshInfo.set("path",Foam::json::Value(cfg.mesh));
                meshInfo.set("cells",Foam::json::Value(static_cast<long long>(mesh.cellCount())));
                meshInfo.set("faces",Foam::json::Value(static_cast<long long>(mesh.faceCount())));
                manifest["files"].set("mesh",meshInfo);
                telemetry.flushManifest();
            }
        }
        RFAM_LEAVE(react_chemistry_setup);
#endif
        ChemistryReport chemistryReport;
        checkpointIdentity.modelId=model->modelId();

        // Dense-cryogenic backend for the 130 K LOX boundary: the PR isochoric
        // path has no single-phase liquid root there, so the same GERG-2008
        // residual table that the structured solver registers is used here.
        static GERGResidualTable gergTable;
        static GERGSaturationTable saturationTable;
        if(cfg.thermoBackend=="legacy")
        {
            const std::string path=std::string(RFAM_SOURCE_DIR)
                +"/tests/reference/gerg_residual_table_v3.csv";
            if(gergTable.load(path))
            {
                model->setGERGTable(&gergTable);
                std::cout<<"gerg_table=loaded "<<path<<std::endl;
            }
            else
                std::cout<<"gerg_table=missing (dense cryogenic states will fail closed)"
                         <<std::endl;
            // Dome + minimal two-phase companion: cold mixed cells that land
            // on the unstable spinodal branch have no single-phase root.
            const std::string satPath=std::string(RFAM_SOURCE_DIR)
                +"/tests/reference/gerg_saturation_table_v3.csv";
            if(saturationTable.load(satPath))
            {
                model->setSaturationTable(&saturationTable);
                std::cout<<"gerg_saturation=loaded "<<satPath<<std::endl;
            }
            else
                std::cout<<"gerg_saturation=missing (spinodal cells fail closed)"
                         <<std::endl;
        }

        // The conserved field must be generated by the SAME model that will
        // invert it. Building U with StateConversion::primitiveToConserved
        // uses a different (legacy EOS) energy convention, which makes every
        // later revertRhoE recovery systematically inconsistent.
        auto consistentState=[&](const PrimitiveVars& w,const thermo::EquilibriumState& q)
        {
            const auto cell=CompleteStateConversion::fromPrimitiveAndEquilibrium(w,q);
            if(!cell.success())
                throw std::runtime_error("model-consistent state construction failed: "
                    +std::string(thermo::equilibriumStatusName(cell.status)));
            return *cell.state;
        };

        std::map<int,ZoneInfo> zones;
        for(std::size_t f=0;f<mesh.faceCount();++f)
        {
            const auto& face=mesh.face(f);
            if(face.neighbour>=0)continue;
            auto& zone=zones[face.boundaryZone];
            zone.faces.push_back(f);
            const auto a=mesh.faceAreaVector(f);
            zone.area+=std::sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
        }
        if(!zones.count(kOutletZone))throw std::runtime_error("mesh has no outlet zone");
        auto configure=[&](const std::vector<int>& ids,const std::string& name,
                           double mdot,double supplyPressure,double temperature,
                           const std::vector<double>& composition)
        {
            // Each zone in the family carries an equal share of the design flow.
            const double share=mdot/static_cast<double>(ids.size());
            for(int id:ids)
            {
                auto it=zones.find(id);
                if(it==zones.end())throw std::runtime_error("missing inlet zone "+std::to_string(id));
                it->second.targetMassFlow=share;
                it->second.currentMassFlow=share;
                it->second.supplyPressure=supplyPressure;
                it->second.temperature=temperature;
                it->second.composition=composition;
                it->second.name=name+std::to_string(id);
            }
        };
        configure({10,11,12},"lox",cfg.loxMassFlow,cfg.loxPinj,cfg.loxTemperature,
                  mixtureY(1.0,0.0));
        configure({20,21,22},"fuel",cfg.fuelMassFlow,cfg.fuelPinj,cfg.fuelTemperature,
                  mixtureY(0.0,1.0));
        // Route B has one central equivalent patch. Its cold-flow composition
        // is the mass-weighted core-torch blend, not ambient background gas.
        const double torchLox=0.0980,torchFuel=0.0367;
        configure({30},"torch",cfg.torchMassFlow,cfg.torchPinj,cfg.torchTemperature,
                  mixtureY(torchLox/(torchLox+torchFuel),
                           torchFuel/(torchLox+torchFuel)));
        // Burnt-gas injector face: skip the injector entirely and impose the
        // equilibrium combustion products of the design O/F on every propellant
        // zone. This removes the two things this case cannot model faithfully --
        // the liquid-to-gas injector transition and the (2.4x oversized) LOX
        // annulus area -- and leaves a single-phase steady chamber, which is what
        // the dome boundary-layer question actually needs. The inflow is
        // single-phase gas by construction, so no injector momentum is claimed.
        if(cfg.burntInlet)
        {
            std::vector<double> burnt(species.size(),0.0);
            const auto setByName=[&](const char* name,double y)
            {
                const int g=speciesIndexOf(name);
                if(g<0)throw std::runtime_error(std::string("--burnt-inlet needs ")+name
                    +" in the mechanism");
                burnt[static_cast<std::size_t>(g)]=y;
            };
            // C3H8 + 5 O2 -> 3 CO2 + 4 H2O at the design O/F; the deficit or
            // excess oxygen is carried as O2 (lean) or as CO/H2 (rich is not
            // represented by the two-step mechanism, so rich runs stay CO2/H2O).
            const double of=cfg.loxMassFlow/std::max(cfg.fuelMassFlow,1e-30);
            const double ofStoich=159.99/44.097;
            const double f=std::min(of/ofStoich,1.0);
            const double mCO2=3.0*44.0095*f, mH2O=4.0*18.01528*f;
            const double mO2ex=(1.0-f)*cfg.loxMassFlow;
            const double total=mCO2+mH2O+mO2ex;
            setByName("CO2",mCO2/total);
            setByName("H2O",mH2O/total);
            setByName("O2",mO2ex/total);
            const double T=cfg.burntTemperature;
            if(!std::isfinite(T)||T<=0.)
                throw std::runtime_error("--burnt-inlet requires --burnt-T (burnt-gas total temperature)");
            for(int id:{10,11,12,20,21,22})
            {
                auto it=zones.find(id);
                if(it==zones.end())continue;
                it->second.composition=burnt;
                it->second.temperature=T;
                it->second.name="burnt"+std::to_string(id);
            }
            std::cout<<"burnt_inlet=1 O/F="<<of<<" T_total="<<T<<" Y_CO2="
                     <<burnt[static_cast<std::size_t>(speciesIndexOf("CO2"))]
                     <<" Y_H2O="<<burnt[static_cast<std::size_t>(speciesIndexOf("H2O"))]
                     <<" Y_O2="<<burnt[static_cast<std::size_t>(speciesIndexOf("O2"))]
                     <<std::endl;
        }

        std::vector<int> inletZones;
        for(auto& kv:zones)
        {
            if(kv.first==kWallZone||kv.first==kOutletZone)continue;
            inletZones.push_back(kv.first);
            auto& zone=kv.second;
            if(zone.targetMassFlow==0.)continue;
            const auto supply=model->initializePTY(zone.supplyPressure,zone.temperature,zone.composition);
            if(!supply.success()||!tetraAdmitsAcoustics(*supply.state))
                throw std::runtime_error("supply thermodynamics failed for "+zone.name+": "+supply.diagnostics.backendStatus);
            zone.supplyEnthalpy=supply.state->h;
            if(supply.state->entropyAvailable)zone.supplyEntropy=supply.state->s;
            if(cfg.inletCd>0.&&std::isfinite(zone.supplyEntropy))
            {
                // Isentropic sonic-point scan: the Cd injector model must cap
                // the flow at the choked value when the local face pressure
                // falls below the critical pressure.
                double best=0.;bool usable=false;
                for(int sample=1;sample<=60;++sample)
                {
                    const double p=zone.supplyPressure*(1.0-0.985*sample/60.0);
                    const auto isen=model->initializePSY(p,zone.supplyEntropy,zone.composition);
                    if(!isen.success()||!isen.state)continue;
                    const double drop=zone.supplyEnthalpy-isen.state->h;
                    if(!(drop>0.))continue;
                    usable=true;
                    best=std::max(best,isen.state->rho*std::sqrt(2.0*drop)*zone.area);
                }
                if(!usable)
                    throw std::runtime_error("choked capacity scan failed for "+zone.name);
                zone.chokedMassFlow=best;
            }
            std::cout<<"supply "<<zone.name<<" p_Pa="<<zone.supplyPressure<<" T_K="<<zone.temperature
                     <<" rho_kg_m3="<<supply.state->rho<<" h0_J_kg="<<zone.supplyEnthalpy
                     <<" sound_m_s="<<supply.state->acoustic.soundSpeed()<<std::endl;
        }
        std::cout<<"zones:\n";
        for(auto& kv:zones)
            std::cout<<"  zone "<<kv.first<<(kv.second.name.empty()?"":" ("+kv.second.name+")")
                     <<" faces="<<kv.second.faces.size()
                     <<" area_mm2="<<kv.second.area*1e6
                     <<" mdot="<<kv.second.targetMassFlow
                     <<" choked_mdot="<<(cfg.inletCd>0.?kv.second.chokedMassFlow:0.)
                     <<" T="<<kv.second.temperature<<"\n";
        if(cfg.inletCd>0.)
            std::cout<<"inlet_cd="<<cfg.inletCd
                     <<" : each main inlet delivers Cd*isentropic capacity at its local"
                     <<" face pressure (design targets remain the zone mdot values)\n";
        if(cfg.inletSupplyCap)
            std::cout<<"inlet_supply_cap=1 : each main inlet delivers min(design mdot,"
                     <<" ideal isentropic capacity at its local face pressure)\n";
        if(cfg.inletFaceSupplyCap)
            std::cout<<"inlet_face_supply_cap=1 : each face is checked at its own pressure;"
                     <<" capacity-limited faces redistribute zone flow\n";

        TetraFiniteVolume fv(mesh,species,
            std::unique_ptr<FluxScheme>(new AUSMFluxAdapter(species)));
        fv.setProgressCallback([](std::size_t done,std::size_t total)
        {
            Foam::telemetry::progress(done,total,"faces");
        });
        // Gas-phase molecular transport for the multi-species legacy backend.
        // Declared outside the provider block because the provider lambda, which
        // captures it by address, outlives that block.
        SpeciesRegistry transportRegistry;
        if(cfg.molecularViscous)
        {
            if(cfg.thermoBackend=="coolprop")
            {
#ifdef RFAM_HAVE_COOLPROP
                fv.setMolecularTransportProvider([](const PrimitiveVars& w,
                                                     const thermo::EquilibriumState& q)
                {
                    if(!tetraAdmitsMolecularTransport(q.region,-1.))
                        throw std::runtime_error("molecular-viscous two-phase closure not admitted");
                    // Transport flashes use their own thread-local backends, not
                    // the state cache used by conserved-state recovery.
                    static thread_local std::unique_ptr<CoolProp::AbstractState> oxygen(
                        CoolProp::AbstractState::factory("HEOS","Oxygen"));
                    static thread_local std::unique_ptr<CoolProp::AbstractState> propane(
                        CoolProp::AbstractState::factory("HEOS","Propane"));
                    static thread_local std::unique_ptr<CoolProp::AbstractState> mixture(
                        CoolProp::AbstractState::factory("HEOS","Oxygen&Propane"));
                    const bool pureO2=w.Y[0]==1.,pureC3=w.Y[1]==1.;
                    auto& state=pureO2?*oxygen:(pureC3?*propane:*mixture);
                    state.unspecify_phase();
                    if(!pureO2&&!pureC3)state.set_mass_fractions(w.Y);
                    state.update(CoolProp::PT_INPUTS,w.p,w.T);
                    if(!tetraAdmitsMolecularTransport(q.region,state.Q()))
                        throw std::runtime_error("molecular-viscous phase/transport closure not admitted");
                    const auto props=transport::fromState(state);
                    return std::make_pair(props.mu,props.lambda);
                });
                std::cout<<"molecular_viscous=on single_phase_only=1 wall=no-slip-adiabatic\n";
#else
                throw std::runtime_error("molecular-viscous requires native CoolProp build");
#endif
            }
            else
            {
                // Multi-species legacy gas closure (Wilke mixing of per-species
                // Sutherland viscosities, Eucken conductivity). This is what
                // lets the reacting path carry a boundary layer at all: CoolProp
                // holds only the O2/C3H8 binary and cannot follow burnt gas.
                SpeciesRegistry& registry=transportRegistry;
                for(const auto& name:species)
                    if(!registry.has(name))
                        throw std::runtime_error("--molecular-viscous legacy closure has no"
                            " registry entry for "+name);
                // Report how much of a wide mechanism borrows Sutherland constants
                // from its nearest-molecular-weight homologue, so the approximation
                // is visible in every run log instead of being inherited silently.
                {
                    std::vector<double> probe(species.size(),
                        1.0/static_cast<double>(species.size()));
                    double muProbe=0.0;
                    std::size_t substituted=0;
                    if(TransportProps::tryMuMixGeneral(species,probe,1500.0,registry,
                                                       muProbe,substituted))
                        std::cout<<"molecular_viscous=wilke species="<<species.size()
                                 <<" sutherland_substituted="<<substituted
                                 <<" mu_probe_Pas="<<muProbe<<std::endl;
                }
                fv.setMolecularTransportProvider(
                    [&species,&registry,&model](const PrimitiveVars& w,
                                                const thermo::EquilibriumState& q)
                    {
                        // Wilke needs a normalised composition; the reacting
                        // recovery closes Y only to its own repair tolerance, so
                        // normalise here instead of fail-closing on a state the
                        // chemistry step itself just accepted.
                        double sumY=0.;
                        for(double yi:w.Y) sumY+=yi;
                        if(!(std::abs(sumY-1.)<=1e-3))
                            throw std::runtime_error("molecular-viscous legacy gas closure:"
                                " composition does not sum to one (sum="
                                +std::to_string(sumY)+")");
                        PrimitiveVars wn=w;
                        for(double& yi:wn.Y) yi/=sumY;
                        const auto r=TransportProps::evaluateCommittedLegacyGas(
                            species,wn,q,*model,registry);
                        if(!r.success())
                        {
                            static std::atomic<int> warned{0};
                            if(warned.fetch_add(1)<3)
                            {
                                std::string comp;
                                double sy=0.;
                                for(std::size_t i=0;i<species.size()&&i<w.Y.size();++i)
                                { comp+=species[i]+"="+std::to_string(w.Y[i])+"";
                                  sy+=w.Y[i]; }
                                std::cerr<<"LEGACY_TRANSPORT_DIAG p="<<w.p<<" T="<<w.T
                                    <<" rho="<<q.rho<<" sumY="<<sy<<" ["
                                    <<comp<<"]\n";
                            }
                            throw std::runtime_error(std::string(
                                "molecular-viscous legacy gas closure not admitted: ")
                                +TransportProps::statusName(r.status));
                        }
                        return std::make_pair(r.point.mu,r.point.kappa);
                    });
                std::cout<<"molecular_viscous=on closure=wilke-sutherland species="
                         <<species.size()<<" wall=no-slip-adiabatic\n";
            }
        }
        // Menter k-omega SST on the unstructured tetra path. The eddy viscosity
        // enters through the same mu/lambda channel as the molecular transport,
        // so it also reaches the LU-SGS diagonal. Steady pseudo-time only: the
        // physical-time branch rejects --molecular-viscous, which SST requires.
        std::unique_ptr<TetraTurbulenceOperator> turbulence;
        if(!cfg.turbulenceModel.empty())
        {
            if(cfg.turbulenceModel!="sst")
                throw std::invalid_argument("--turbulence supports only 'sst'");
            if(!cfg.molecularViscous)
                throw std::invalid_argument("--turbulence sst requires --molecular-viscous");
            TetraTurbulenceOperator::Settings ts;
            ts.intensity=cfg.turbulenceIntensity;
            ts.lengthScale=cfg.turbulenceLength;
            ts.wallZone=kWallZone;
            turbulence.reset(new TetraTurbulenceOperator(fv,ts));
            TetraTurbulenceOperator* raw=turbulence.get();
            fv.setTurbulentTransportProvider(
                [raw](std::size_t cell,const PrimitiveVars& w,
                      const thermo::EquilibriumState& q)
                { return raw->transportAddition(cell,w,q); });
            std::cout<<"turbulence=sst closure=menter-1994 transports=k,omega"
                     <<" wall=no-slip intensity="<<cfg.turbulenceIntensity
                     <<" length_scale_m="<<cfg.turbulenceLength<<std::endl;
        }
        if(cfg.interfaceConsistentFlux)
        {
            // Interface-consistent flux: the mass flux uses the equilibrium
            // density of the interface state (p_half, h_upwind, Y_upwind) taken
            // from the property table, so the mass/energy/species fluxes all
            // belong to one equilibrium state. Cell states remain exact through
            // the recovery; only the flux's projection error drops.
            const thermo::TableAcoustic* table=thermo::tableAcoustic();
            if(!table)
                throw std::runtime_error("--interface-consistent-flux requires RFAM_TETRA_TABLE");
            fv.setConsistentDensityHook([table](double yO2,double p,double h,double& rho)
                {return table->stateFromPHY(yO2,p,h,rho);});
            std::cout<<"interface_consistent_flux=on table="<<table->path()<<std::endl;
        }

        PrimitiveVars initial(species.size());
        initial.p=cfg.p0;initial.T=cfg.T0;initial.U={0.0,0.0,0.0};initial.Y=backgroundY;
        const auto initialQ=model->initializePTY(initial.p,initial.T,initial.Y);
        if(!initialQ.success())
            throw std::runtime_error("initial state rejected: "
                +std::string(thermo::equilibriumStatusName(initialQ.status)));
        // Optional 1D continuity blueprint: fill the domain with the design
        // mass flow instead of a stagnant state. Without this the explicit
        // local-time-stepping run needs O(1e4) pseudo steps just to fill the
        // chamber; with it the near-injector region is the only part that has
        // to develop. Slab areas come from the mesh itself.
        if(cfg.blueprint)
        {
            const double totalMassFlow=cfg.loxMassFlow+cfg.fuelMassFlow+cfg.torchMassFlow;
            const auto plateNormals=cfg.blueprintPlateWall
                ?tetraPlateWallNormals(mesh,0.2709527893066406)
                :std::vector<std::array<double,3>>{};
            std::size_t plateCells=0;
            constexpr int kBins=400;
            std::vector<double> slabVolume(kBins,0.0);
            double zMin=1e30,zMax=-1e30;
            for(std::size_t c=0;c<mesh.cellCount();++c)
            {
                const double z=mesh.cellCentroid(c).z;
                zMin=std::min(zMin,z);zMax=std::max(zMax,z);
            }
            const double binWidth=(zMax-zMin)/kBins;
            for(std::size_t c=0;c<mesh.cellCount();++c)
            {
                const double z=mesh.cellCentroid(c).z;
                const int bin=std::min(kBins-1,std::max(0,int((z-zMin)/binWidth)));
                slabVolume[bin]+=mesh.volume(c);
            }
            std::size_t filled=0;
            for(std::size_t c=0;c<mesh.cellCount();++c)
            {
                const double z=mesh.cellCentroid(c).z;
                const int bin=std::min(kBins-1,std::max(0,int((z-zMin)/binWidth)));
                const double area=std::max(slabVolume[bin]/binWidth,1e-8);
                const double speed=totalMassFlow/(initialQ.state->rho*area);
                PrimitiveVars w=initial;
                w.U={0.0,0.0,std::min(speed,300.0)};
                if(cfg.blueprintPlateWall&&removeBlueprintPlateImpact(w,initial.Y,
                    initial.p,initial.T,plateNormals[c]))++plateCells;
                const auto cell=consistentState(w,*initialQ.state);
                fv.setState(c,cell.W,cell.U,cell.Q);
                ++filled;
            }
            std::cout<<"blueprint_initialized cells="<<filled
                     <<" plate_wall_projected="<<plateCells
                     <<" total_mass_flow="<<totalMassFlow
                     <<" mean_axial_speed="<<totalMassFlow
                        /(initialQ.state->rho*(2.958603075278e-3/0.4065))
                     <<std::endl;
        }
        else
        {
            const auto cell=consistentState(initial,*initialQ.state);
            for(std::size_t c=0;c<mesh.cellCount();++c)
                fv.setState(c,cell.W,cell.U,cell.Q);
        }
        if(cfg.blueprintInjector)
        {
            // Zone-driven pre-fill: the cold-flow injector channels are primed
            // with their propellant before the chamber transient starts. Each
            // inlet zone's upstream cell cluster is found by a breadth-first
            // walk over the mesh connectivity (the fluid passages), so the
            // assignment follows the actual channels instead of a geometric
            // guess, and the state comes from the zone's own supply (3 % below
            // the supply pressure, marginally colder: the inlet availability
            // check rejects a face backpressure at or above the supply, and the
            // isentropic capacity flash needs the primed state just inside the
            // supply isentrope).
            std::size_t filled=0;
            std::size_t zoneCount=0;
            std::vector<char> visited(mesh.cellCount(),0);
            for(int id:{10,11,12,20,21,22,30})
            {
                auto it=zones.find(id);
                if(it==zones.end())continue;
                ZoneInfo& zone=it->second;
                if(!(zone.targetMassFlow>0.))continue;   // torch off stays ambient
                // Pre-fill pressure sits just above the chamber, at 5 % of the
                // available supply-to-chamber drop. Priming close to the supply
                // (the first attempt used 0.97*p_supply) leaves only a few per
                // cent of margin, and the jet dynamic head (rho*u^2/2 ~ 0.2 MPa
                // for the fuel) then pushes the channel static pressure past
                // the supply during the first steps, which the inlet
                // availability check correctly rejects.
                const double pressure=cfg.p0+0.05*(zone.supplyPressure-cfg.p0);
                const double temperature=0.996*zone.temperature;
                std::vector<double> composition=zone.composition;
                if(composition.size()!=species.size())continue;
                std::vector<std::size_t> frontier;
                for(std::size_t f:zone.faces)
                {
                    const std::size_t owner=static_cast<std::size_t>(mesh.face(f).owner);
                    if(!visited[owner]){visited[owner]=1;frontier.push_back(owner);}
                }
                std::size_t cluster=0;
                for(int layer=0;layer<cfg.blueprintLayers;++layer)
                {
                    std::vector<std::size_t> next;
                    // Linear blend over the outermost filled layers: a 22:1
                    // liquid/gas density jump confined to one cell layer makes
                    // the first pseudo-time steps very stiff (the recovery
                    // thrashes and the step cost explodes). Blending the
                    // outermost layers into the ambient state spreads the
                    // interface over a few cells -- physically an injector
                    // primed with a residual gas blanket -- and keeps the flux
                    // inside the smooth regime.
                    double blendWeight=0.;
                    if(cfg.blueprintBlend>0&&layer>=cfg.blueprintLayers-cfg.blueprintBlend)
                        blendWeight=double(layer-(cfg.blueprintLayers-cfg.blueprintBlend)+1)
                                    /double(cfg.blueprintBlend);
                    // Per-cell flashes are independent (the thermo model keeps
                    // thread-local states), so the pre-fill is parallel; only
                    // the breadth-first walk itself stays serial.
                    std::atomic<int> layerFilled{0};
                    #pragma omp parallel for schedule(dynamic,32)
                    for(std::size_t i=0;i<frontier.size();++i)
                    {
                        const std::size_t c=frontier[i];
                        PrimitiveVars ambient=fv.W(c);
                        PrimitiveVars w(species.size());
                        w.p=(1.-blendWeight)*pressure+blendWeight*ambient.p;
                        w.T=(1.-blendWeight)*temperature+blendWeight*ambient.T;
                        w.Y=composition;
                        for(std::size_t s=0;s<species.size();++s)
                            w.Y[s]=(1.-blendWeight)*composition[s]+blendWeight*ambient.Y[s];
                        w.U={0.0,0.0,0.0};
                        const auto q=model->initializePTY(w.p,w.T,w.Y);
                        if(q.success())
                        {
                            // Primed channel velocity = the steady injection
                            // velocity mdot/(rho*A) along the channel axis.
                            // A stagnant liquid pre-fill is not a physical
                            // starting state: LOX is nearly incompressible, so
                            // a few pseudo-time steps of mass inflow with no
                            // outflow compress the channel by a per cent and
                            // the static pressure spikes past the supply (A9
                            // tripped the inlet availability check at step 10
                            // exactly this way).
                            const double uInjection=zone.targetMassFlow>0.
                                ?zone.targetMassFlow/(q.state->rho*zone.area):0.;
                            w.U={0.0,0.0,(1.-blendWeight)*uInjection};
                            const auto cell=consistentState(w,*q.state);
                            fv.setState(c,cell.W,cell.U,cell.Q);
                            layerFilled.fetch_add(1,std::memory_order_relaxed);
                        }
                    }
                    filled+=static_cast<std::size_t>(layerFilled.load());
                    cluster+=static_cast<std::size_t>(layerFilled.load());
                    for(std::size_t c:frontier)
                    {
                        for(std::size_t f:mesh.cellFaces(c))
                        {
                            const auto& face=mesh.face(f);
                            for(int other:{face.owner,face.neighbour})
                            {
                                if(other<0)continue;
                                const std::size_t o=static_cast<std::size_t>(other);
                                if(!visited[o]){visited[o]=1;next.push_back(o);}
                            }
                        }
                    }
                    frontier.swap(next);
                    if(frontier.empty())break;
                }
                ++zoneCount;
                std::cout<<"injector_prefill_zone "<<zone.name<<" cells="<<cluster
                         <<" p="<<pressure<<" T="<<temperature<<std::endl;
            }
            std::cout<<"injector_prefill cells="<<filled
                     <<" zones="<<zoneCount<<" layers="<<cfg.blueprintLayers<<std::endl;
        }
#ifdef RFAM_HAVE_CHEMISTRY
        if(cfg.reacting&&cfg.ignitionT>0.0)
        {
            // Hot ignition kernel (steady-start pilot): seed a hot premixed
            // pocket on the torch axis. Without it the cold chamber would
            // need thermal runaway from the 700 K initialization and the
            // pseudo-time run would stall before ignition.
            std::size_t kernelCells=0;
            // Torch axis from the zone's own face owners: the mesh axis is NOT
            // at (y,z)=(0,0) in this geometry, so the kernel is placed relative
            // to the torch zone centroid, not the global origin.
            double axisX=0.0,axisY=0.0,axisZ=0.0;
            double outletZ=0.0;
            {
                const auto it=zones.find(30);
                if(it!=zones.end()&&!it->second.faces.empty())
                {
                    double sx=0.,sy=0.,sz=0.;
                    for(std::size_t f:it->second.faces)
                    {
                        const auto centre=mesh.cellCentroid(
                            static_cast<std::size_t>(mesh.face(f).owner));
                        sx+=centre.x;sy+=centre.y;sz+=centre.z;
                    }
                    const double n=static_cast<double>(it->second.faces.size());
                    axisX=sx/n;axisY=sy/n;axisZ=sz/n;
                }
                const auto out=zones.find(kOutletZone);
                if(out!=zones.end()&&!out->second.faces.empty())
                {
                    double sz=0.;
                    for(std::size_t f:out->second.faces)
                        sz+=mesh.cellCentroid(
                            static_cast<std::size_t>(mesh.face(f).owner)).z;
                    outletZ=sz/static_cast<double>(out->second.faces.size());
                }
                std::cout<<"ignition_axis x="<<axisX<<" y="<<axisY
                         <<" z="<<axisZ<<" outlet_z="<<outletZ<<std::endl;
            }
            // Chamber direction: the injector plate sits at the torch face and
            // the chamber runs toward the outlet zone.
            const double axialSign=(outletZ<axisZ)?-1.0:1.0;
            const double kernelZ=axisZ+axialSign*std::abs(cfg.ignitionX);
            const double radialTol=cfg.ignitionR;
            for(std::size_t c=0;c<mesh.cellCount();++c)
            {
                const auto centre=mesh.cellCentroid(c);
                const double radius=std::hypot(centre.x-axisX,centre.y-axisY);
                if(radius>radialTol)continue;
                if(std::abs(centre.z-kernelZ)>0.02)continue;
                const PrimitiveVars w=fv.W(c);
                PrimitiveVars wHot=w;
                wHot.T=cfg.ignitionT;
                // The kernel must carry a REACTIVE mixture: the chamber
                // background is inert products, so a hot pocket of products
                // cannot ignite. Stoichiometric C3H8/O2 (0.216/0.784 by mass)
                // matches the torch blend and sits inside the flammable range.
                wHot.Y.assign(fv.species().size(),0.0);
                if(cfg.reacting)
                {
                    wHot.Y[iC3H8]=0.216;
                    wHot.Y[iO2]=0.784;
                }
                else
                {
                    wHot.Y=w.Y;
                }
                const auto hot=model->initializePTY(wHot.p,cfg.ignitionT,wHot.Y);
                if(!hot.success())continue;
                // The constructor requires the primitive to carry the SAME
                // (p,T,Y) as the equilibrium state it is paired with.
                wHot.p=hot.state->p;
                wHot.T=hot.state->T;
                const auto cell=consistentState(wHot,*hot.state);
                fv.setState(c,cell.W,cell.U,cell.Q);
                ++kernelCells;
            }
            std::cout<<"ignition_kernel cells="<<kernelCells
                     <<" T="<<cfg.ignitionT<<" x="<<cfg.ignitionX
                     <<" r="<<cfg.ignitionR<<std::endl;
        }
#endif
        double seedMaxRhoRel=0.,seedMaxHRel=0.;
        if(!cfg.seedField.empty())
        {
            std::ifstream input(cfg.seedField,std::ios::binary);
            char magic[9]{};std::uint64_t count=0,ns=0;
            input.read(magic,9);input.read(reinterpret_cast<char*>(&count),8);input.read(reinterpret_cast<char*>(&ns),8);
            if(!input||std::string(magic,9)!="RFAMTETF1"||count!=mesh.cellCount()||ns!=species.size())
                throw std::runtime_error("seed field header/mesh/species mismatch");
            std::vector<float> centres(count*3),fields(count*8),ys(count*ns);
            input.read(reinterpret_cast<char*>(centres.data()),centres.size()*sizeof(float));
            input.read(reinterpret_cast<char*>(fields.data()),fields.size()*sizeof(float));
            input.read(reinterpret_cast<char*>(ys.data()),ys.size()*sizeof(float));
            if(!input)throw std::runtime_error("truncated seed field");
            char extra;if(input.read(&extra,1))throw std::runtime_error("seed field has trailing data");
            std::size_t bad=count;std::string reason;
            std::vector<CompleteCellState> seeded(count,CompleteCellState(species.size()));
            #pragma omp parallel for schedule(static)
            for(std::size_t c=0;c<count;++c)
            {
                try
                {
                    const auto centre=mesh.cellCentroid(c);
                    if(std::abs(centre.x-centres[3*c])>1e-7||std::abs(centre.y-centres[3*c+1])>1e-7
                       ||std::abs(centre.z-centres[3*c+2])>1e-7)throw std::runtime_error("seed mesh coordinates differ");
                    PrimitiveVars w(species.size());w.p=fields[8*c];w.T=fields[8*c+1];
                    for(int d=0;d<3;++d)w.U[d]=fields[8*c+2+d];
                    double sum=0.;
                    for(std::size_t s=0;s<ns;++s){w.Y[s]=ys[c*ns+s];sum+=w.Y[s];}
                    if(!std::isfinite(sum)||std::abs(sum-1.)>2e-6)throw std::runtime_error("invalid rounded seed composition");
                    for(double& value:w.Y)value/=sum;
                    auto result=model->initializePTY(w.p,w.T,w.Y);
                    if(!result.success())
                    {
                        // The rounded seed can land exactly on the saturation
                        // line, where CoolProp refuses a PT flash. Rebuild from
                        // the conserved pair instead: the (rho,e) recovery owns
                        // the lever rule for that case. e = h - p/rho.
                        const double rhoSeed=fields[8*c+5],hSeed=fields[8*c+7];
                        const double eSeed=rhoSeed>0.?hSeed-w.p/rhoSeed:0.;
                        result=model->recoverRhoE(rhoSeed,eSeed,w.Y);
                    }
                    if(!result.success())throw std::runtime_error(result.diagnostics.backendStatus);
                    if(!tetraAdmitsAcoustics(*result.state))throw std::runtime_error("unadmitted seed state");
                    // The recovery is the authority: its (rho,e) closure was just
                    // verified. The dumped h is float32-rounded and, in the
                    // projection-accepted development runs, carries the last
                    // step's projection residual, so the old 2e-5 check
                    // rejected valid seeds (A46, cell 83956). Keep a gross 5 %
                    // sanity bound against unit/reference errors and record the
                    // actual mismatch for the status.
                    const double rhoRel=std::abs(result.state->rho-fields[8*c+5])
                        /std::max(result.state->rho,1.);
                    const double hRel=std::abs(result.state->h-fields[8*c+7])
                        /std::max(std::abs(result.state->h),1e5);
                    // In the record-based projection mode the dumped h carries
                    // the last step's projection residual (measured 2-11 % in
                    // the near-critical band), so it cannot be a rejection
                    // criterion: the recovered equilibrium state is the
                    // authority and the mismatch is recorded. Outside that mode
                    // the strict 2e-5 seed contract stays (unit/reference
                    // errors must fail closed).
                    const double seedGate=projectionAcceptEnabled()?1.0:2e-5;
                    if(rhoRel>seedGate||hRel>seedGate)
                        throw std::runtime_error("seed model/energy reference mismatch");
                    seedMaxRhoRel=std::max(seedMaxRhoRel,rhoRel);
                    seedMaxHRel=std::max(seedMaxHRel,hRel);
                    seeded[c]=consistentState(w,*result.state);
                }
                catch(const std::exception& e)
                {
                    #pragma omp critical(rfam_seed_failure)
                    {if(c<bad){bad=c;reason=e.what();}}
                }
            }
            if(bad<count)throw std::runtime_error("seed field rejected at cell "+std::to_string(bad)+": "+reason);
            for(std::size_t c=0;c<count;++c)fv.setState(c,seeded[c].W,seeded[c].U,seeded[c].Q);
            std::cout<<"FIELD_SEED_INITIALIZED cells="<<count<<" exact_restart=0 source="<<cfg.seedField
                     <<" max_rho_rel="<<seedMaxRhoRel<<" max_h_rel="<<seedMaxHRel<<std::endl;
        }
        double initialMass=0.0;
        for(std::size_t c=0;c<mesh.cellCount();++c)
            initialMass+=fv.Q(c).rho*mesh.volume(c);
        std::cout<<"initial: p="<<cfg.p0<<"Pa T="<<cfg.T0<<"K rho="<<initialQ.state->rho
                 <<"kg/m3 sound="<<initialQ.state->acoustic.soundSpeed()
                 <<"m/s mass="<<initialMass<<"kg\n";

        std::vector<double> inletFaceMassFlow(mesh.faceCount(),0.0);
        std::vector<PrimitiveVars> boundaryW(mesh.faceCount(),initial);
        std::vector<thermo::EquilibriumState> boundaryQ(mesh.faceCount(),*initialQ.state);
        const TetraBoundaryProvider provider=[&](std::size_t,std::size_t face,
                                                 PrimitiveVars& w,thermo::EquilibriumState& q)
        {w=boundaryW[face];q=boundaryQ[face];};

        fv.setBoundaryFluxProvider(
            [&](std::size_t,std::size_t f,ConservedVars& flux)
            {
                const int id=mesh.face(f).boundaryZone;
                if(id==kOutletZone&&cfg.outletPressure>0.)
                {
                    // Characteristic construction returns the boundary trace,
                    // not a ghost extrapolation to feed a second Riemann solve.
                    flux=TetraPressureOutlet::physicalFlux(boundaryW[f],boundaryQ[f],fv.faceNormal(f));
                    return true;
                }
                if(id==kWallZone||id==kOutletZone||cfg.inletMode!="mass")return false;
                const auto& zone=zones.at(id);if(!(zone.targetMassFlow>0.))return false;
                const auto& w=boundaryW[f];const auto& n=fv.faceNormal(f);
                if(cfg.inletFaceSupplyCap&&inletFaceMassFlow[f]==0.)
                {
                    flux=ConservedVars(species.size());
                    for(int d=0;d<3;++d)flux.rhoU[d]=w.p*n[d];
                    return true;
                }
                const double mass=cfg.inletFaceSupplyCap
                    ?-inletFaceMassFlow[f]/fv.faceArea(f)
                    :-zone.currentMassFlow/zone.area;
                flux=ConservedVars(species.size());flux.rho=mass;
                flux.rhoE=mass*zone.supplyEnthalpy;
                for(int d=0;d<3;++d)flux.rhoU[d]=mass*w.U[d]+w.p*n[d];
                for(std::size_t s=0;s<species.size();++s)flux.rhoY[s]=mass*zone.composition[s];
                flux.rhoZ=mass*w.Z;flux.rhoC=mass*w.C;return true;
            });

        // Zone face normal (unit) cached once: inlets inject along -n.
        std::map<int,std::vector<std::array<double,3>>> zoneNormals;
        for(auto& kv:zones)
        {
            auto& normals=zoneNormals[kv.first];
            normals.reserve(kv.second.faces.size());
            for(std::size_t f:kv.second.faces)
            {
                auto a=mesh.faceAreaVector(f);
                const double magnitude=std::sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
                if(!(magnitude>0))throw std::runtime_error("zero-area boundary face");
                normals.push_back({a.x/magnitude,a.y/magnitude,a.z/magnitude});
            }
        }

        auto fillInlet=[&](ZoneInfo& zone,int zoneId,double speed)
        {
            const auto& normals=zoneNormals[zoneId];
            for(std::size_t i=0;i<zone.faces.size();++i)
            {
                const std::size_t f=zone.faces[i];
                const std::size_t owner=static_cast<std::size_t>(mesh.face(f).owner);
                const auto& n=normals[i];
                boundaryW[f]=fv.W(owner);
                boundaryW[f].p=zone.supplyPressure;
                boundaryW[f].T=zone.temperature;
                boundaryW[f].Y=zone.composition;
                boundaryW[f].U={-speed*n[0],-speed*n[1],-speed*n[2]};
            }
        };
        auto zoneFlow=[&](ZoneInfo& zone,int zoneId)->double
        {
            if(cfg.inletFaceSupplyCap&&zoneId!=kWallZone&&zoneId!=kOutletZone)
                return -zone.currentMassFlow;
            double total=0.0;   // positive = outflow of the domain
            for(std::size_t f:zone.faces)
            {
                const std::size_t owner=static_cast<std::size_t>(mesh.face(f).owner);
                total+=fv.faceMassFlow(f,fv.W(owner),boundaryW[f],fv.Q(owner),boundaryQ[f]);
            }
            return total;
        };

        std::unique_ptr<TetraPressureOutlet> pressureOutlet;
        if(cfg.outletPressure>0.)
        {
            pressureOutlet=std::make_unique<TetraPressureOutlet>(*model,cfg.outletPressure,cfg.backflowTemperature,backflowSpeciesY,
                cfg.outletRoundoffShortcut);
            std::cout<<"outlet=equilibrium-pressure p_Pa="<<cfg.outletPressure
                     <<" backflow_T_K="<<cfg.backflowTemperature<<" backflow_Y_O2="<<cfg.backflowY[0]
                     <<" backflow_is_air=0"<<std::endl;
        }
        else std::cout<<"outlet=zero-gradient diagnostic=1"<<std::endl;
        double boundaryInletSeconds=0.,boundaryWallSeconds=0.,boundaryOutletSeconds=0.;
        std::size_t outletExactPressureOwners=0,outletNearOutflow=0,outletNearBackflow=0;
        std::size_t outletNearBackflowAbove=0,outletNearBackflowBelow=0;
        double outletMaxRelativePressureDifference=0.,outletMaxNearBackflowSpeed=0.;
        auto buildBoundary=[&]()
        {
            const auto inletStart=std::chrono::steady_clock::now();
            for(auto& kv:zones)
            {
                ZoneInfo& zone=kv.second;
                if(kv.first==kWallZone||kv.first==kOutletZone)continue;
                if(zone.targetMassFlow==0.)
                {
                    for(std::size_t f:zone.faces)
                    {
                        const auto owner=mesh.face(f).owner;auto w=fv.W(owner);
                        const auto& n=fv.faceNormal(f);
                        const double un=w.U[0]*n[0]+w.U[1]*n[1]+w.U[2]*n[2];
                        for(int d=0;d<3;++d)w.U[d]-=2.*un*n[d];
                        boundaryW[f]=w;boundaryQ[f]=fv.Q(owner);
                    }
                    continue;
                }
                double sum=0.0;
                for(std::size_t f:zone.faces)
                    sum+=fv.W(static_cast<std::size_t>(mesh.face(f).owner)).p*fv.faceArea(f);
                zone.pressure=zone.faces.empty()?cfg.p0:sum/zone.area;
                if(cfg.inletFaceSupplyCap)
                {
                    std::vector<double> areas(zone.faces.size()),capacities(zone.faces.size());
                    zone.blockedArea=0.;zone.availableCapacity=0.;
                    zone.entropyMargin=std::numeric_limits<double>::infinity();
                    for(std::size_t i=0;i<zone.faces.size();++i)
                    {
                        const auto f=zone.faces[i];
                        const double pFace=fv.W(mesh.face(f).owner).p;
                        areas[i]=fv.faceArea(f);
                        if(pFace>=zone.supplyPressure)
                        {
                            zone.blockedArea+=areas[i];
                            continue;
                        }
                        thermo::EquilibriumSeed faceSeed;
                        faceSeed.temperature=zone.temperature;
                        const auto isen=model->initializePSY(pFace,zone.supplyEntropy,
                                                             zone.composition,&faceSeed);
                        if(!isen.success()||!isen.state)
                            throw std::runtime_error("inlet face isentropic capacity failed for "+zone.name
                                +" face "+std::to_string(f)+": "+isen.diagnostics.backendStatus);
                        const double drop=zone.supplyEnthalpy-isen.state->h;
                        if(drop>0.)
                            capacities[i]=isen.state->rho*std::sqrt(2.*drop)*areas[i];
                        zone.availableCapacity+=capacities[i];
                    }
                    const auto allocation=allocateTetraInletFlow(areas,capacities,zone.targetMassFlow);
                    zone.currentMassFlow=0.;
                    for(std::size_t i=0;i<zone.faces.size();++i)
                    {
                        const auto f=zone.faces[i];
                        const std::size_t owner=static_cast<std::size_t>(mesh.face(f).owner);
                        const double faceFlow=allocation[i];
                        inletFaceMassFlow[f]=faceFlow;
                        zone.currentMassFlow+=faceFlow;
                        if(faceFlow==0.)
                        {
                            // Diagnostic no-return valve: pressure traction only.
                            // This is not a model of the upstream feed system.
                            boundaryW[f]=fv.W(owner);
                            boundaryQ[f]=fv.Q(owner);
                            continue;
                        }
                        const double pFace=fv.W(owner).p;
                        const double mass=faceFlow/areas[i];
                        double t=zone.temperature;bool closed=false;
                        thermo::EquilibriumState q;
                        for(int iter=0;iter<25;++iter)
                        {
                            const auto result=model->initializePTY(pFace,t,zone.composition);
                            if(!result.success()||!tetraAdmitsAcoustics(*result.state))
                                throw std::runtime_error("inlet face static flash failed for "+zone.name
                                    +" face "+std::to_string(f)+": "+result.diagnostics.backendStatus);
                            q=*result.state;
                            const double u=mass/q.rho;
                            const double residual=q.h+.5*u*u-zone.supplyEnthalpy;
                            if(std::abs(residual)<=1e-8*std::max(std::abs(zone.supplyEnthalpy),1e5))
                            {closed=true;break;}
                            if(!q.thermodynamics)
                                throw std::runtime_error("two-phase inlet face needs a dedicated PH flash");
                            const auto& td=*q.thermodynamics;
                            const double derivative=td.Cp+u*u*td.dpdT_rho/(q.rho*td.dpdrho_T);
                            if(!std::isfinite(derivative)||derivative<=0.)
                                throw std::runtime_error("invalid inlet face PH derivative");
                            t+=std::clamp(-residual/derivative,-.05*t,.05*t);
                        }
                        if(!closed)
                            throw std::runtime_error("inlet face total enthalpy failed to close for "+zone.name);
                        if(q.entropyAvailable&&std::isfinite(zone.supplyEntropy))
                        {
                            const double margin=q.s-zone.supplyEntropy;
                            zone.entropyMargin=std::min(zone.entropyMargin,margin);
                            if(margin<-cfg.supplyEntropyTolerance*std::max(std::abs(zone.supplyEntropy),1.)
                                &&cfg.supplyEntropyStrict)
                                throw std::runtime_error("inlet face below supply entropy for "+zone.name);
                        }
                        const auto& n=zoneNormals[kv.first][i];
                        const double speed=mass/q.rho;
                        PrimitiveVars w(species.size());w.p=q.p;w.T=q.T;w.Y=zone.composition;
                        w.U={-speed*n[0],-speed*n[1],-speed*n[2]};
                        boundaryW[f]=w;boundaryQ[f]=q;
                    }
                    continue;
                }
                if(cfg.inletMode=="mass")
                {
                    // Prescribed mdot + stagnation h; pressure is extrapolated
                    // from the domain. Solve static h(T,p)+u^2/2=h_supply with
                    // u=mdot/(rho*A). Supply pressure is an availability check.
                    // With --inlet-cd the zone mass flow is not prescribed: it
                    // is Cd times the isentropic capacity at the local face
                    // pressure (capped at the choked flow), so the injector
                    // adapts to the chamber pressure and the supply entropy
                    // remains realizable by construction.
                    if(zone.pressure>=zone.supplyPressure)
                        throw std::runtime_error("inlet backpressure exceeds supply for "+zone.name);
                    double targetMassFlow=zone.targetMassFlow;
                    if(cfg.inletCd>0.||cfg.inletSupplyCap)
                    {
                        const auto isen=model->initializePSY(zone.pressure,zone.supplyEntropy,zone.composition);
                        if(!isen.success()||!isen.state)
                            throw std::runtime_error("isentropic inlet capacity failed for "+zone.name
                                +": "+isen.diagnostics.backendStatus);
                        const double drop=zone.supplyEnthalpy-isen.state->h;
                        double capacity=0.;
                        if(drop>0.)capacity=isen.state->rho*std::sqrt(2.0*drop)*zone.area;
                        if(std::isfinite(zone.chokedMassFlow))
                            capacity=std::min(capacity,zone.chokedMassFlow);
                        if(cfg.inletCd>0.)
                        {
                            targetMassFlow=cfg.inletCd*capacity;
                            if(!(targetMassFlow>0.))
                                throw std::runtime_error("Cd inlet capacity vanished for "+zone.name);
                        }
                        else targetMassFlow=std::min(zone.targetMassFlow,capacity);
                    }
                    zone.currentMassFlow=targetMassFlow;
                    const double mass=targetMassFlow/zone.area;
                    double t=zone.temperature;bool closed=false;
                    thermo::EquilibriumState q;
                    for(int iter=0;iter<25;++iter)
                    {
                        const auto result=model->initializePTY(zone.pressure,t,zone.composition);
                        if(!result.success()||!tetraAdmitsAcoustics(*result.state))
                            throw std::runtime_error("static inlet flash failed for "+zone.name+": "+result.diagnostics.backendStatus);
                        q=*result.state;const double u=mass/q.rho;
                        const double residual=q.h+.5*u*u-zone.supplyEnthalpy;
                        if(std::abs(residual)<=1e-8*std::max(std::abs(zone.supplyEnthalpy),1e5)){closed=true;break;}
                        if(!q.thermodynamics)throw std::runtime_error("two-phase supply PH inversion requires dedicated flash");
                        const auto& td=*q.thermodynamics;
                        const double derivative=td.Cp+u*u*td.dpdT_rho/(q.rho*td.dpdrho_T);
                        if(!std::isfinite(derivative)||derivative<=0.)throw std::runtime_error("invalid inlet PH derivative");
                        t+=std::clamp(-residual/derivative,-.05*t,.05*t);
                    }
                    if(!closed)throw std::runtime_error("inlet total enthalpy failed to close for "+zone.name);
                    // The supply-entropy margin of the imposed mass-flow inlet.
                    // The mass-flow inlet (Fluent-equivalent) does not claim the
                    // injector is lossless, so by default a negative margin is a
                    // reported diagnostic, not a veto; --supply-entropy-strict
                    // restores the fail-closed behaviour.
                    if(std::isfinite(zone.supplyEntropy)&&q.entropyAvailable)
                        zone.entropyMargin=q.s-zone.supplyEntropy;
                    else zone.entropyMargin=std::numeric_limits<double>::quiet_NaN();
                    if(std::isfinite(zone.entropyMargin)&&zone.entropyMargin<-cfg.supplyEntropyTolerance*std::max(std::abs(zone.supplyEntropy),1.))
                    {
                        if(cfg.supplyEntropyStrict)
                        {
                            double capacity=std::numeric_limits<double>::quiet_NaN();
                            try
                            {
                                const auto isen=model->initializePSY(zone.pressure,zone.supplyEntropy,zone.composition);
                                if(isen.success()&&isen.state)
                                {
                                    const double drop=zone.supplyEnthalpy-isen.state->h;
                                    if(drop>0.)capacity=isen.state->rho*std::sqrt(2.*drop)*zone.area;
                                }
                            }
                            catch(const std::exception&){ }
                            std::ostringstream detail;
                            detail<<"inlet flow/area/backpressure violates supply entropy for "<<zone.name
                                  <<"; equilibrium inlet model cannot realize requested mass flux"
                                  <<" target_mdot_kg_s="<<targetMassFlow
                                  <<" capacity_estimate_kg_s="<<capacity
                                  <<" capacity_ratio="<<(std::isfinite(capacity)?capacity/targetMassFlow:std::numeric_limits<double>::quiet_NaN())
                                  <<" area_m2="<<zone.area
                                  <<" interior_p_Pa="<<zone.pressure
                                  <<" supply_p_Pa="<<zone.supplyPressure
                                  <<" static_p_Pa="<<q.p
                                  <<" static_T_K="<<q.T
                                  <<" static_rho_kg_m3="<<q.rho
                                  <<" static_speed_m_s="<<(mass/q.rho)
                                  <<" supply_h_J_kg="<<zone.supplyEnthalpy
                                  <<" static_h_J_kg="<<q.h
                                  <<" supply_s_J_kgK="<<zone.supplyEntropy
                                  <<" static_s_J_kgK="<<q.s
                                  <<" entropy_margin_J_kgK="<<zone.entropyMargin;
                            throw std::runtime_error(detail.str());
                        }
                        if(!zone.entropyWarned)
                        {
                            zone.entropyWarned=true;
                            std::cerr<<"WARNING inlet "<<zone.name<<" face state is below the supply entropy by "
                                     <<-zone.entropyMargin<<" J/(kg*K); mass-flow inlet imposed, margin reported"<<std::endl;
                        }
                    }
                    const double speed=mass/q.rho;
                    for(std::size_t i=0;i<zone.faces.size();++i)
                    {
                        const auto f=zone.faces[i];const auto& n=zoneNormals[kv.first][i];
                        PrimitiveVars w(species.size());w.p=q.p;w.T=q.T;w.Y=zone.composition;
                        w.U={-speed*n[0],-speed*n[1],-speed*n[2]};
                        boundaryW[f]=w;boundaryQ[f]=q;
                    }
                    continue;
                }
                PrimitiveVars ghost(species.size());
                ghost.p=zone.supplyPressure > 0.0 ? zone.supplyPressure : cfg.pinj;
                ghost.T=zone.temperature;ghost.Y=zone.composition;
                const auto ghostQ=model->initializePTY(ghost.p,ghost.T,ghost.Y);
                if(!ghostQ.success())
                    throw std::runtime_error("inlet ghost state rejected for zone "+zone.name
                        +": "+thermo::equilibriumStatusName(ghostQ.status));
                for(std::size_t f:zone.faces)boundaryQ[f]=*ghostQ.state;

                const double rhoMean=ghostQ.state->rho;
                double speed=zone.targetMassFlow/(rhoMean*std::max(zone.area,1e-30));
                if(cfg.controlMassFlow)
                {
                    // Bisection on the numerical patch mass flow. The kernel
                    // flux is upwind, so the achieved flow is not exactly
                    // rho*u*A and must be measured rather than assumed.
                    double lower=0.0;
                    fillInlet(zone,kv.first,lower);
                    double lowerFlow=-zoneFlow(zone,kv.first);
                    if(lowerFlow>zone.targetMassFlow*(1.+cfg.massFlowTolerance))
                        throw std::runtime_error("static inlet pressure and mass-flow target are incompatible for "+zone.name
                            +"; zero-speed inflow="+std::to_string(lowerFlow)
                            +" target="+std::to_string(zone.targetMassFlow));
                    double upper=std::max(speed,1.0);
                    double upperFlow=0.0;
                    bool bracketed=false;
                    for(int expand=0;expand<40;++expand)
                    {
                        fillInlet(zone,kv.first,upper);
                        upperFlow=-zoneFlow(zone,kv.first);
                        if(upperFlow>=zone.targetMassFlow){bracketed=true;break;}
                        lower=upper;lowerFlow=upperFlow;upper*=2.0;
                    }
                    if(!bracketed)
                        throw std::runtime_error("inlet mass-flow target not bracketed for zone "+zone.name);
                    const double tolerance=cfg.massFlowTolerance*std::abs(zone.targetMassFlow);
                    for(int it=0;it<cfg.massFlowIterations;++it)
                    {
                        const double middle=0.5*(lower+upper);
                        fillInlet(zone,kv.first,middle);
                        const double flow=-zoneFlow(zone,kv.first);
                        if(std::abs(flow-zone.targetMassFlow)<=tolerance){speed=middle;break;}
                        if(flow<zone.targetMassFlow){lower=middle;lowerFlow=flow;}
                        else {upper=middle;upperFlow=flow;}
                        speed=middle;
                    }
                    (void)lowerFlow; (void)upperFlow;
                }
                fillInlet(zone,kv.first,speed);
                if(cfg.controlMassFlow&&std::abs(-zoneFlow(zone,kv.first)-zone.targetMassFlow)>
                   cfg.massFlowTolerance*zone.targetMassFlow)
                    throw std::runtime_error("inlet mass-flow control did not converge for "+zone.name);
                {
                    std::size_t f=zone.faces.front();
                    const std::size_t owner=static_cast<std::size_t>(mesh.face(f).owner);
                    const double rho=fv.Q(owner).rho;
                    const auto& gw=boundaryW[f];
                    const double gspeed=std::sqrt(gw.U[0]*gw.U[0]+gw.U[1]*gw.U[1]+gw.U[2]*gw.U[2]);
                    const auto& gq=boundaryQ[f];
                    static std::atomic<bool> reported{false};
                    if(!reported.exchange(true))
                        std::cout<<"inlet ghost "<<zone.name
                                 <<": p_int="<<zone.pressure
                                 <<" T_in="<<zone.temperature
                                 <<" rho_in="<<gq.rho
                                 <<" sound_in="<<gq.acoustic.soundSpeed()
                                 <<" speed="<<gspeed
                                 <<" Mach_in="<<gspeed/gq.acoustic.soundSpeed()
                                 <<" mdot_achieved="<<-zoneFlow(zone,kv.first)
                                 <<std::endl;
                }
            }
            const auto inletDone=std::chrono::steady_clock::now();
            boundaryInletSeconds=std::chrono::duration<double>(inletDone-inletStart).count();
            for(int id:{kWallZone,kOutletZone})
            {
                const auto groupStart=std::chrono::steady_clock::now();
                auto it=zones.find(id);
                if(it==zones.end())continue;
                const auto& normals=zoneNormals[id];
                std::size_t bad=it->second.faces.size();std::string detail;
                std::size_t exactPressureOwners=0,nearOutflow=0,nearBackflow=0;
                std::size_t nearBackflowAbove=0,nearBackflowBelow=0;
                double maxRelativePressureDifference=0.,maxNearBackflowSpeed=0.;
                #pragma omp parallel for schedule(static) reduction(+:exactPressureOwners,nearOutflow,nearBackflow,nearBackflowAbove,nearBackflowBelow) reduction(max:maxRelativePressureDifference,maxNearBackflowSpeed)
                for(std::size_t i=0;i<it->second.faces.size();++i)
                {
                    try
                    {
                    const std::size_t f=it->second.faces[i];
                    const std::size_t owner=static_cast<std::size_t>(mesh.face(f).owner);
                    PrimitiveVars w=fv.W(owner);
                    thermo::EquilibriumState q=fv.Q(owner);
                    if(id==kOutletZone&&pressureOutlet)
                    {
                        exactPressureOwners+=w.p==cfg.outletPressure;
                        maxRelativePressureDifference=std::max(maxRelativePressureDifference,
                            std::abs(w.p-cfg.outletPressure)/cfg.outletPressure);
                        const auto& normal=fv.faceNormal(f);
                        const double un=w.U[0]*normal[0]+w.U[1]*normal[1]+w.U[2]*normal[2];
                        const double ulp=std::nextafter(cfg.outletPressure,
                            std::numeric_limits<double>::infinity())-cfg.outletPressure;
                        if(std::abs(w.p-cfg.outletPressure)<=ulp)
                        {
                            nearOutflow+=un>=0.;
                            nearBackflow+=un<0.;
                            if(un<0.)
                            {
                                nearBackflowAbove+=w.p>cfg.outletPressure;
                                nearBackflowBelow+=w.p<cfg.outletPressure;
                                maxNearBackflowSpeed=std::max(maxNearBackflowSpeed,-un);
                            }
                        }
                    }
                    if(id==kWallZone)
                    {
                        const auto& n=normals[i];
                        const double un=w.U[0]*n[0]+w.U[1]*n[1]+w.U[2]*n[2];
                        w.U[0]-=2.0*un*n[0];
                        w.U[1]-=2.0*un*n[1];
                        w.U[2]-=2.0*un*n[2];
                    }
                    else if(pressureOutlet)
                    {
                        pressureOutlet->apply(fv.W(owner),fv.Q(owner),fv.faceNormal(f),w,q);
                    }
                    else
                    {
                        // Outflow: zero gradient. Backflow: hold a quiescent
                        // ambient state instead of ingesting chamber gas
                        // through the nozzle exit.
                        const auto& n=normals[i];
                        const double un=w.U[0]*n[0]+w.U[1]*n[1]+w.U[2]*n[2];
                        if(un<0.0)
                        {
                            w=initial;
                            w.U={0.0,0.0,0.0};
                            q=*initialQ.state;
                        }
                    }
                    boundaryW[f]=w;
                    boundaryQ[f]=q;
                    }
                    catch(const std::exception& e)
                    {
                        #pragma omp critical(rfam_outlet_failure)
                        {if(i<bad){bad=i;detail=e.what();}}
                    }
                }
                if(bad<it->second.faces.size())
                    throw std::runtime_error("boundary face "+std::to_string(it->second.faces[bad])+": "+detail);
                const double groupSeconds=std::chrono::duration<double>(
                    std::chrono::steady_clock::now()-groupStart).count();
                if(id==kWallZone)boundaryWallSeconds=groupSeconds;
                else
                {
                    boundaryOutletSeconds=groupSeconds;
                    outletExactPressureOwners=exactPressureOwners;
                    outletNearOutflow=nearOutflow;
                    outletNearBackflow=nearBackflow;
                    outletNearBackflowAbove=nearBackflowAbove;
                    outletNearBackflowBelow=nearBackflowBelow;
                    outletMaxNearBackflowSpeed=maxNearBackflowSpeed;
                    outletMaxRelativePressureDifference=maxRelativePressureDifference;
                }
            }
        };

        // The implicit step stages cells; this callback stages the boundary
        // side effects. A rejected inlet/outlet restores both parts of the
        // accepted state, including per-zone pressure and ghost composition.
        bool boundaryRejected=false;
        bool physicalBoundaryTrial=false;
        bool boundaryTrialRejected=false;
        double boundarySnapshotSeconds=0.,boundaryBuildSeconds=0.;
        const auto validateBoundaryTransaction=[&]()
        {
            const auto boundaryStart=std::chrono::steady_clock::now();
            struct SavedBoundary
            {
                std::size_t face;
                PrimitiveVars w;
                thermo::EquilibriumState q;
            };
            std::vector<SavedBoundary> saved;
            std::vector<std::pair<std::size_t,double>> savedFaceFlow;
            struct SavedZone
            {
                double pressure,flow,blockedArea,availableCapacity,entropyMargin;
            };
            std::map<int,SavedZone> savedPressure;
            for(const auto& item:zones)
            {
                const auto& z=item.second;
                savedPressure[item.first]={z.pressure,z.currentMassFlow,z.blockedArea,
                                           z.availableCapacity,z.entropyMargin};
                for(auto f:item.second.faces)
                {
                    saved.push_back({f,boundaryW[f],boundaryQ[f]});
                    if(cfg.inletFaceSupplyCap)savedFaceFlow.push_back({f,inletFaceMassFlow[f]});
                }
            }
            const auto snapshotDone=std::chrono::steady_clock::now();
            try{buildBoundary();}
            catch(...)
            {
                if(physicalBoundaryTrial)boundaryTrialRejected=true;
                else boundaryRejected=true;
                for(auto& value:saved)
                {
                    boundaryW[value.face]=std::move(value.w);
                    boundaryQ[value.face]=std::move(value.q);
                }
                for(const auto& item:savedFaceFlow)inletFaceMassFlow[item.first]=item.second;
                for(const auto& item:savedPressure)
                {
                    auto& z=zones.at(item.first);
                    z.pressure=item.second.pressure;
                    z.currentMassFlow=item.second.flow;
                    z.blockedArea=item.second.blockedArea;
                    z.availableCapacity=item.second.availableCapacity;
                    z.entropyMargin=item.second.entropyMargin;
                }
                throw;
            }
            const auto buildDone=std::chrono::steady_clock::now();
            boundarySnapshotSeconds=std::chrono::duration<double>(snapshotDone-boundaryStart).count();
            boundaryBuildSeconds=std::chrono::duration<double>(buildDone-snapshotDone).count();
        };

        std::vector<ConservedVars> residual;
        std::vector<double> dt;
        double physicalTime=0.0;
        double physicalDtUsed=0.0;
        int initialAcceptedStep=0;
        if(!cfg.physicalResume.empty())
        {
            const auto loaded=TetraPhysicalCheckpoint::read(
                cfg.physicalResume,fv,checkpointIdentity);
            if(loaded.steps>static_cast<std::uint64_t>(std::numeric_limits<int>::max()-cfg.steps))
                throw std::runtime_error("physical checkpoint step count exceeds driver range");
            initialAcceptedStep=static_cast<int>(loaded.steps);
            physicalTime=loaded.physicalTime;
            physicalDtUsed=loaded.lastDt;
            std::cout<<"PHYSICAL_CHECKPOINT_RESTORED steps="<<initialAcceptedStep
                     <<" physical_time_s="<<physicalTime
                     <<" last_dt_s="<<physicalDtUsed<<std::endl;
        }
        double physicalInnerResidualRatio=1.0;
        double physicalInnerMaxResidualRatio=1.0;
        double physicalInnerResidualAbsolute=0.0;
        double physicalInnerMaxResidualAbsolute=0.0;
        double physicalBudgetError=0.0;
        double physicalCflLimit=0.0;
        bool physicalBoundaryInputRejected=false;
        PhysicalStepRetryAudit physicalRetryAudit;
        int physicalInnerAcceptedIterations=0;
        bool physicalInnerConverged=false;
        const auto addScaledConserved=[](ConservedVars& target,
                                          const ConservedVars& source,double scale)
        {
            target.rho+=scale*source.rho;
            target.rhoE+=scale*source.rhoE;
            target.rhoZ+=scale*source.rhoZ;
            target.rhoC+=scale*source.rhoC;
            for(int d=0;d<3;++d)target.rhoU[d]+=scale*source.rhoU[d];
            for(std::size_t s=0;s<target.rhoY.size();++s)
                target.rhoY[s]+=scale*source.rhoY[s];
        };
        const auto physicalBoundaryBudgetError=[&](const std::vector<ConservedVars>& equation)
        {
            ConservedVars sum(species.size()),scale(species.size());
            for(const auto& r:equation)addScaledConserved(sum,r,1.);
            for(std::size_t f=0;f<mesh.faceCount();++f)
            {
                if(mesh.face(f).neighbour>=0)continue;
                const auto& flux=fv.lastFaceFlux(f);
                scale.rho+=std::abs(flux.rho*fv.faceArea(f));
                scale.rhoE+=std::abs(flux.rhoE*fv.faceArea(f));
                scale.rhoZ+=std::abs(flux.rhoZ*fv.faceArea(f));
                scale.rhoC+=std::abs(flux.rhoC*fv.faceArea(f));
                for(int d=0;d<3;++d)
                    scale.rhoU[d]+=std::abs(flux.rhoU[d]*fv.faceArea(f));
                for(std::size_t s=0;s<species.size();++s)
                    scale.rhoY[s]+=std::abs(flux.rhoY[s]*fv.faceArea(f));
            }
            double error=std::max({std::abs(sum.rho)/std::max(scale.rho,1e-12),
                                   std::abs(sum.rhoE)/std::max(scale.rhoE,1e-12),
                                   std::abs(sum.rhoZ)/std::max(scale.rhoZ,1e-12),
                                   std::abs(sum.rhoC)/std::max(scale.rhoC,1e-12)});
            for(int d=0;d<3;++d)
                error=std::max(error,std::abs(sum.rhoU[d])/std::max(scale.rhoU[d],1e-12));
            for(std::size_t s=0;s<species.size();++s)
                error=std::max(error,std::abs(sum.rhoY[s])/std::max(scale.rhoY[s],1e-12));
            return error;
        };
        std::atomic<std::size_t> frozenCells{0};
        std::atomic<std::size_t> limitedCells{0};
        std::atomic<std::size_t> stiffenedCells{0};
        double frozenVolume=0.0;
        const double frozenVolumeLimit=1.0e-3*2.958603075278e-3;
        double referenceMass=-1.0;
        double referenceResidual=-1.0;
        int acceptedSteps=initialAcceptedStep;
        bool advanceFailed=false,numericallyConverged=false;
        double windowChange=0.;std::size_t windowLimited=0;
        const double totalInflow=cfg.loxMassFlow+cfg.fuelMassFlow+cfg.torchMassFlow;
        double oxygenInflow=0.,fuelInflow=0.;
        double energyReference=std::max({std::abs(initialQ.state->e),
            initialQ.state->acoustic.pressureDensityDerivative,1e5});
        for(int id:inletZones)
        {
            const auto& zone=zones.at(id);
            oxygenInflow+=zone.targetMassFlow*zone.composition[0];
            fuelInflow+=zone.targetMassFlow*zone.composition[1];
            if(zone.targetMassFlow>0.)energyReference=std::max(energyReference,std::abs(zone.supplyEnthalpy));
        }
        std::unique_ptr<TetraSteadyMonitor> steadyMonitor;
        if(totalInflow>0.)steadyMonitor=std::make_unique<TetraSteadyMonitor>(
            totalInflow,initialQ.state->acoustic.soundSpeed(),energyReference,oxygenInflow,fuelInflow);
        TetraSteadyMonitor::Report steadyReport;
        // Binary field writer shared by the periodic (--dump-every) and final
        // dumps. Cell-centred: p, T, U, rho, sound speed, h and the species
        // mass fractions, matching the RFAMTETF1 dump the report tools read.
        const auto writeFieldDump=[&](const std::string& path)
        {
            std::ofstream out(path,std::ios::binary);
            if(!out)throw std::runtime_error("cannot open dump file "+path);
            const std::uint64_t cells=mesh.cellCount();
            const std::uint64_t speciesCount=species.size();
            out.write("RFAMTETF1",9);
            out.write(reinterpret_cast<const char*>(&cells),8);
            out.write(reinterpret_cast<const char*>(&speciesCount),8);
            std::vector<float> buffer(static_cast<std::size_t>(cells)*3);
            for(std::size_t c=0;c<cells;++c)
            {
                const auto centre=mesh.cellCentroid(c);
                buffer[3*c+0]=static_cast<float>(centre.x);
                buffer[3*c+1]=static_cast<float>(centre.y);
                buffer[3*c+2]=static_cast<float>(centre.z);
            }
            out.write(reinterpret_cast<const char*>(buffer.data()),
                      std::streamsize(buffer.size()*sizeof(float)));
            std::vector<float> fields(static_cast<std::size_t>(cells)*8,0.0f);
            for(std::size_t c=0;c<cells;++c)
            {
                const auto& w=fv.W(c);
                fields[8*c+0]=static_cast<float>(w.p);
                fields[8*c+1]=static_cast<float>(w.T);
                fields[8*c+2]=static_cast<float>(w.U[0]);
                fields[8*c+3]=static_cast<float>(w.U[1]);
                fields[8*c+4]=static_cast<float>(w.U[2]);
                fields[8*c+5]=static_cast<float>(fv.Q(c).rho);
                fields[8*c+6]=static_cast<float>(fv.Q(c).acoustic.soundSpeed());
                fields[8*c+7]=static_cast<float>(fv.Q(c).h);
            }
            out.write(reinterpret_cast<const char*>(fields.data()),
                      std::streamsize(fields.size()*sizeof(float)));
            std::vector<float> composition(static_cast<std::size_t>(cells)*speciesCount);
            for(std::size_t c=0;c<cells;++c)
                for(std::size_t s=0;s<speciesCount;++s)
                    composition[c*speciesCount+s]=static_cast<float>(fv.W(c).Y[s]);
            out.write(reinterpret_cast<const char*>(composition.data()),
                      std::streamsize(composition.size()*sizeof(float)));
            if(!out)throw std::runtime_error("field dump write failed "+path);
            std::cout<<"FIELD_DUMP_WRITTEN "<<path<<" cells="<<cells
                     <<" species="<<speciesCount<<std::endl;
        };
        const auto writePhaseSidecar=[&](const std::string& fieldPath)
        {
            std::ofstream phases(fieldPath+".phase.csv");
            if(!phases)throw std::runtime_error("cannot write phase metadata");
            phases<<std::setprecision(17)
                <<"cell,region,vapor_mass_fraction,vapor_volume_fraction,acoustic_evidence,projected\n";
            for(std::size_t c=0;c<mesh.cellCount();++c)
            {
                const auto& q=fv.Q(c);
                phases<<c<<","<<static_cast<int>(q.region)<<","<<q.vaporMassFraction
                      <<","<<q.vaporVolumeFraction<<","<<static_cast<int>(q.acoustic.evidence)
                      <<","<<q.conservedStateProjected<<"\n";
            }
            if(!phases)throw std::runtime_error("phase metadata write failed");
        };
        const auto writeBoundaryBudget=[&](const std::string& fieldPath)
        {
            std::ofstream budget(fieldPath+".boundary.csv");
            if(!budget)throw std::runtime_error("cannot write boundary budget");
            budget<<std::setprecision(17)
                <<"zone,area_m2,target_inflow_kg_s,outward_mass_kg_s,outward_O2_kg_s,"
                  "outward_C3H8_kg_s,outward_energy_W\n";
            for(const auto& item:zones)
            {
                const auto& zone=item.second;ConservedVars total(species.size());
                for(std::size_t f:zone.faces)
                {
                    const auto& flux=fv.lastFaceFlux(f);const double area=fv.faceArea(f);
                    total.rho+=area*flux.rho;total.rhoE+=area*flux.rhoE;
                    for(std::size_t s=0;s<species.size();++s)total.rhoY[s]+=area*flux.rhoY[s];
                }
                budget<<item.first<<","<<zone.area<<","<<zone.targetMassFlow<<","<<total.rho
                      <<","<<total.rhoY[0]<<","<<total.rhoY[1]<<","<<total.rhoE<<"\n";
            }
            if(!budget)throw std::runtime_error("boundary budget write failed");
        };
        RFAM_ENTER(react_step);
        for(int localStep=0;localStep<=cfg.steps;++localStep)
        {
            const int step=initialAcceptedStep+localStep;
            telemetry.subPhase("step",__FILE__,__LINE__);
            if(localStep>0&&cfg.physicalTime)
            {
                physicalRetryAudit.beginStep();
                bool acceptedPhysicalStep=false;
                std::string failureDetail;
                bool nonRetryableBoundaryFailure=false;
                double trialDt=cfg.physicalDt;
                for(int retry=0;retry<cfg.physicalMaxRetries&&!acceptedPhysicalStep;++retry)
                {
                    std::vector<ConservedVars> reference(mesh.cellCount());
                    std::vector<CompleteCellState> old(mesh.cellCount());
                    for(std::size_t c=0;c<mesh.cellCount();++c)
                    {
                        reference[c]=fv.U(c);
                        old[c].W=fv.W(c);old[c].U=fv.U(c);old[c].Q=fv.Q(c);
                    }
                    // The turbulence fields are advanced inside the trial step and
                    // are not part of the cell-state transaction, so snapshot them
                    // for the reject path.
                    std::vector<double> oldK,oldOmega;
                    if(turbulence){oldK=turbulence->kField();oldOmega=turbulence->omegaField();}
                    try
                    {
                        validateBoundaryTransaction();
                        fv.computeResidual(provider,residual);
                        fv.localTimeStep(cfg.cfl,dt,provider);
                        const double stableDt=*std::min_element(dt.begin(),dt.end());
                        physicalCflLimit=stableDt;
                        if(!(stableDt>0.)||!std::isfinite(stableDt))
                            throw std::runtime_error("physical-time global CFL step is invalid");
                        const double requestedDt=cfg.physicalDt>0.?cfg.physicalDt:stableDt;
                        trialDt=trialDt>0.?trialDt:requestedDt;
                        if(retry==0&&cfg.physicalDt>stableDt)
                            std::cerr<<"PHYSICAL_DT_ABOVE_CFL requested="<<cfg.physicalDt
                                     <<" global_cfl_diagnostic="<<stableDt<<std::endl;
                        if(!(trialDt>0.)||!std::isfinite(trialDt))
                            throw std::runtime_error("physical-time step is invalid");
                        TetraLUSGS physicalSolver(fv);
                        physicalSolver.setGuard(cfg.guardTMin,cfg.guardTMax,cfg.guardPMin,cfg.guardPMax);
                        TetraLUSGS::Report physicalReport;
                        double initialResidual=0.,initialMaxResidual=0.;bool innerConverged=false;
                        int usedInner=0;double finalRatio=1.,finalMaxRatio=1.;
                        double finalAbsolute=0.,finalMaxAbsolute=0.,finalBudget=0.;
                        for(int inner=0;inner<cfg.physicalInnerSteps;++inner)
                        {
                            const auto innerStart=std::chrono::steady_clock::now();
                            physicalBoundaryTrial=true;
                            physicalReport=physicalSolver.step(provider,cfg.cfl,1.0,cfg.lusgsRelax,*model,
                                validateBoundaryTransaction,cfg.twoPhaseDtScale,cfg.injectorDtScale,
                                &reference,trialDt);
                            physicalBoundaryTrial=false;
                            // Chemistry and turbulence advance inside the trial so a
                            // rejected step rolls both back with the cell states.
#ifdef RFAM_HAVE_CHEMISTRY
                            if(cfg.reacting&&chemOperator&&!cfg.chemCVODE)
                                applyChemistrySources(fv,*model,*chemOperator,
                                                      std::vector<double>(mesh.cellCount(),trialDt),
                                                      1.0,cfg,chemistryReport);
                            if(cfg.chemCVODE)
                                physicalReport.maxDelta=std::max(physicalReport.maxDelta,
                                                                 chemistryReport.maxStateChange);
#endif
                            if(turbulence)
                            {
                                if(!turbulence->wallAdjacentCells()) turbulence->initialize();
                                turbulence->advance(std::vector<double>(mesh.cellCount(),trialDt));
                            }
                            const double current=physicalReport.residualNorm;
                            const double currentMax=physicalReport.residualMaxNorm;
                            if(inner==0)
                            {
                                initialResidual=std::max(current,1e-30);
                                initialMaxResidual=std::max(currentMax,1e-30);
                            }
                            finalRatio=current/initialResidual;
                            finalMaxRatio=currentMax/initialMaxResidual;
                            usedInner=inner+1;
                            if(cfg.physicalProbeInner)
                            {
                                const double seconds=std::chrono::duration<double>(
                                    std::chrono::steady_clock::now()-innerStart).count();
                                const auto centre=mesh.cellCentroid(physicalReport.worstResidualCell);
                                std::cerr<<"PHYSICAL_INNER_DIAGNOSTIC step="<<step
                                         <<" attempt="<<retry+1<<" inner="<<usedInner
                                         <<" trial_dt_s="<<trialDt<<" elapsed_s="<<seconds
                                         <<" residual_s="<<physicalReport.residualSeconds
                                         <<" linear_s="<<physicalReport.linearSeconds
                                         <<" recovery_s="<<physicalReport.recoverySeconds
                                         <<" boundary_s="<<physicalReport.boundarySeconds
                                         <<" boundary_snapshot_s="<<boundarySnapshotSeconds
                                         <<" boundary_build_s="<<boundaryBuildSeconds
                                         <<" inlet_s="<<boundaryInletSeconds
                                         <<" wall_s="<<boundaryWallSeconds
                                         <<" outlet_s="<<boundaryOutletSeconds
                                         <<" outlet_exact_pressure_owners="<<outletExactPressureOwners
                                         <<" outlet_near_outflow="<<outletNearOutflow
                                         <<" outlet_near_backflow="<<outletNearBackflow
                                         <<" outlet_near_backflow_above="<<outletNearBackflowAbove
                                         <<" outlet_near_backflow_below="<<outletNearBackflowBelow
                                         <<" outlet_max_near_backflow_speed="<<outletMaxNearBackflowSpeed
                                         <<" outlet_roundoff_hits="<<(pressureOutlet?pressureOutlet->roundoffHits():0)
                                         <<" outlet_max_relative_dp="<<outletMaxRelativePressureDifference
                                         <<" residual_rms="<<current
                                         <<" residual_max="<<currentMax
                                         <<" ratio_rms="<<finalRatio
                                         <<" ratio_max="<<finalMaxRatio
                                         <<" max_delta="<<physicalReport.maxDelta
                                         <<" worst_cell="<<physicalReport.worstResidualCell
                                         <<" worst_component="<<physicalReport.worstResidualComponent
                                         <<" worst_component_value="<<physicalReport.worstComponentValue
                                         <<" worst_cell_inlet_spectral_fraction="
                                         <<physicalReport.worstCellInletSpectralFraction
                                         <<" worst_cell_volume_m3="<<mesh.volume(physicalReport.worstResidualCell)
                                         <<" worst_cell_centroid_m=("<<centre.x<<","<<centre.y<<","<<centre.z<<")"
                                         <<" tolerance="<<cfg.physicalInnerTolerance
                                         <<std::endl;
                            }
                            if(finalRatio<=cfg.physicalInnerTolerance
                               &&finalMaxRatio<=cfg.physicalInnerTolerance
                               &&physicalReport.maxDelta<=cfg.physicalInnerTolerance)
                            {
                                const auto post=physicalSolver.measurePhysicalResidual(
                                    provider,reference,trialDt);
                                finalRatio=post.first/initialResidual;
                                finalMaxRatio=post.second/initialMaxResidual;
                                finalAbsolute=post.first;
                                finalMaxAbsolute=post.second;
                                finalBudget=physicalBoundaryBudgetError(physicalSolver.residual());
                                if(cfg.physicalProbeInner)
                                    std::cerr<<"PHYSICAL_GATE_DIAGNOSTIC step="<<step
                                             <<" attempt="<<retry+1<<" inner="<<usedInner
                                             <<" post_ratio_rms="<<finalRatio
                                             <<" post_ratio_max="<<finalMaxRatio
                                             <<" absolute_rms="<<finalAbsolute
                                             <<" absolute_max="<<finalMaxAbsolute
                                             <<" budget_error="<<finalBudget
                                             <<" max_delta="<<physicalReport.maxDelta
                                             <<" tolerance="<<cfg.physicalInnerTolerance
                                             <<std::endl;
                                if(finalRatio<=cfg.physicalInnerTolerance
                                   &&finalMaxRatio<=cfg.physicalInnerTolerance
                                   &&finalAbsolute<=cfg.physicalInnerTolerance
                                   &&finalMaxAbsolute<=cfg.physicalInnerTolerance
                                   &&finalBudget<=cfg.physicalInnerTolerance)
                                {innerConverged=true;break;}
                            }
                        }
                        if(!innerConverged)
                            throw std::runtime_error("physical-time inner LU-SGS did not converge");
                        physicalTime+=trialDt;
                        physicalDtUsed=trialDt;
                        physicalInnerResidualRatio=finalRatio;
                        physicalInnerMaxResidualRatio=finalMaxRatio;
                        physicalInnerResidualAbsolute=finalAbsolute;
                        physicalInnerMaxResidualAbsolute=finalMaxAbsolute;
                        physicalBudgetError=finalBudget;
                        physicalInnerAcceptedIterations=usedInner;
                        physicalInnerConverged=true;
                        acceptedSteps=step;
                        acceptedPhysicalStep=true;
                    }
                    catch(const std::exception& e)
                    {
                        failureDetail=e.what();
                        const bool retryableTrialBoundary=boundaryTrialRejected;
                        physicalBoundaryTrial=false;
                        boundaryTrialRejected=false;
                        const auto decision=classifyPhysicalStepFailure(
                            failureDetail,retryableTrialBoundary);
                        nonRetryableBoundaryFailure=decision.boundaryInputRejected();
                        for(std::size_t c=0;c<mesh.cellCount();++c)
                            fv.exchangeState(c,old[c].W,old[c].U,old[c].Q);
                        if(turbulence&&!oldK.empty())turbulence->restoreFields(oldK,oldOmega);
                        if(nonRetryableBoundaryFailure)
                        {
                            physicalBoundaryInputRejected=true;
                            std::cerr<<"PHYSICAL_BOUNDARY_INPUT_REJECTED step="<<step
                                     <<" reason="<<failureDetail<<std::endl;
                            break;
                        }
                        try{validateBoundaryTransaction();}
                        catch(const std::exception& rollbackError)
                        {
                            failureDetail+="; accepted-state boundary rebuild failed: "+std::string(rollbackError.what());
                            if(isNonRetryablePhysicalBoundaryFailure(rollbackError.what()))
                            {
                                physicalBoundaryInputRejected=true;
                                boundaryRejected=true;
                                std::cerr<<"PHYSICAL_BOUNDARY_INPUT_REJECTED step="<<step
                                         <<" reason="<<failureDetail<<std::endl;
                            }
                            advanceFailed=true;
                            break;
                        }
                        if(retry+1<cfg.physicalMaxRetries)
                        {
                            trialDt*=0.5;
                            physicalRetryAudit.recordRetry(retryableTrialBoundary);
                            std::cerr<<"PHYSICAL_STEP_RETRY step="<<step
                                     <<" retry="<<retry+1<<" dt="<<trialDt
                                     <<" reason="<<failureDetail<<std::endl;
                        }
                        else
                            std::cerr<<"PHYSICAL_STEP_ATTEMPTS_EXHAUSTED step="<<step
                                     <<" attempts="<<cfg.physicalMaxRetries
                                     <<" reason="<<failureDetail<<std::endl;
                    }
                }
                if(acceptedPhysicalStep&&!cfg.physicalCheckpoint.empty())
                {
                    auto metadata=checkpointIdentity;
                    metadata.steps=static_cast<std::uint64_t>(step);
                    metadata.physicalTime=physicalTime;
                    metadata.lastDt=physicalDtUsed;
                    const auto checkpoint=cfg.physicalCheckpoint+".s"+std::to_string(step);
                    TetraPhysicalCheckpoint::write(checkpoint,fv,metadata);
                    std::cout<<"PHYSICAL_CHECKPOINT_WRITTEN "<<checkpoint<<std::endl;
                }
                if(!acceptedPhysicalStep)
                {
                    physicalInnerConverged=false;
                    advanceFailed=true;
                    std::cerr<<"TETRA_UPDATE_REJECTED physical_step="<<step
                             <<" "<<failureDetail<<std::endl;
                    if(!cfg.dump.empty())
                    {
                        std::ofstream failure(cfg.dump+".failure.txt");
                        failure<<"failed_step="<<step<<" accepted_steps="<<acceptedSteps
                               <<" physical_time_s="<<std::setprecision(17)<<physicalTime
                               <<"\n"<<failureDetail<<"\n";
                    }
                    break;
                }
            }
            else if(localStep>0&&cfg.lusgs)
            {
                static std::unique_ptr<TetraLUSGS> lusgs;
                if(!lusgs)
                {
                    lusgs.reset(new TetraLUSGS(fv));
                    lusgs->setCompositionCoupling(cfg.lusgsCompositionCoupling);
                }
                lusgs->setGuard(cfg.guardTMin,cfg.guardTMax,cfg.guardPMin,cfg.guardPMax);
                // Pseudo-time ramp: start with a small, robust step and grow it
                // as the field settles (standard steady-acceleration practice).
                double dtScale=cfg.lusgsDtScale;
                if(cfg.lusgsDtScaleStart>0.0)
                {
                    const double fraction=std::min(1.0,
                        static_cast<double>(step)/std::max(cfg.lusgsRampSteps,1));
                    dtScale=cfg.lusgsDtScaleStart
                        +(cfg.lusgsDtScale-cfg.lusgsDtScaleStart)*fraction;
                }
                TetraLUSGS::Report report;
                try
                {
                    const auto validateChemistryTransaction=[&]()
                    {
#ifdef RFAM_HAVE_CHEMISTRY
                        if(cfg.chemCVODE)
                        {
                            std::vector<double> chemistryDt;fv.localTimeStep(cfg.cfl,chemistryDt);
                            applyChemistryCVODE(fv,*model,*chemOperator,chemistryDt,dtScale,cfg,chemistryReport);
                        }
#endif
                        validateBoundaryTransaction();
                    };
                    telemetry.subPhase("lusgs",__FILE__,__LINE__);
                    Foam::telemetry::setStep(step,cfg.steps,acceptedSteps);
                    report=lusgs->step(provider,cfg.cfl,dtScale,cfg.lusgsRelax,*model,
                                       validateChemistryTransaction,cfg.twoPhaseDtScale,cfg.injectorDtScale);
                    stiffenedCells.store(report.stiffened);
                    acceptedSteps=step;
#ifdef RFAM_HAVE_CHEMISTRY
                    if(cfg.reacting&&chemOperator&&!cfg.chemCVODE)
                        applyChemistrySources(fv,*model,*chemOperator,dt,dtScale,cfg,
                                              chemistryReport);
                    if(cfg.chemCVODE)report.maxDelta=std::max(report.maxDelta,chemistryReport.maxStateChange);
#endif
                }
                catch(const std::exception& e)
                {
                    advanceFailed=true;
                    std::cerr<<"TETRA_UPDATE_REJECTED step="<<step<<" "<<e.what()<<std::endl;
                    if(!cfg.dump.empty())
                    {
                        std::ofstream failure(cfg.dump+".failure.txt");
                        failure<<"failed_step="<<step<<" accepted_steps="<<acceptedSteps<<"\n"<<e.what()<<"\n";
                    }
                    break;
                }
                residual=lusgs->residual();
                fv.localTimeStep(cfg.cfl,dt);
                if(turbulence)
                {
                    telemetry.subPhase("turbulence",__FILE__,__LINE__);
                    if(!turbulence->wallAdjacentCells())
                    {
                        RFAM_ENTER(react_wall_distance);
                        turbulence->initialize();
                        RFAM_LEAVE(react_wall_distance);
                    }
                    turbulence->advance(dt);
                }
                frozenCells.store(report.frozen);
                limitedCells.store(report.limited);
                windowChange=std::max(windowChange,report.maxDelta);
                windowLimited=std::max(windowLimited,report.limited);
                if(report.frozenVolume>0.0)frozenVolume=report.frozenVolume;
                if(frozenVolume>frozenVolumeLimit)
                    throw std::runtime_error("frozen cell volume "
                        +std::to_string(frozenVolume)
                        +" m3 exceeds the negligible limit");
            }
            else
            {
                try{buildBoundary();}
                catch(const std::exception& e)
                {
                    boundaryRejected=true;
                    if(cfg.physicalTime)
                    {
                        const std::string detail=e.what();
                        advanceFailed=true;
                        if(isNonRetryablePhysicalBoundaryFailure(detail))
                        {
                            physicalBoundaryInputRejected=true;
                            std::cerr<<"PHYSICAL_BOUNDARY_INPUT_REJECTED initial reason="
                                     <<detail<<std::endl;
                        }
                        if(!cfg.dump.empty())
                        {
                            std::ofstream failure(cfg.dump+".failure.txt");
                            failure<<"failed_step=0 accepted_steps=0 physical_time_s=0\n"
                                   <<detail<<"\n";
                        }
                        break;
                    }
                    throw;
                }
                telemetry.subPhase("residual",__FILE__,__LINE__);
                fv.computeResidual(provider,residual,nullptr);
                telemetry.subPhase("local_dt",__FILE__,__LINE__);
                if(cfg.physicalTime)fv.localTimeStep(cfg.cfl,dt,provider);
                else fv.localTimeStep(cfg.cfl,dt);
                if(turbulence)
                {
                    if(!turbulence->wallAdjacentCells()) turbulence->initialize();
                    turbulence->advance(dt);
                }
                std::vector<double> sorted=dt;
                std::sort(sorted.begin(),sorted.end());
                auto below=[&](double limit)
                {return std::count_if(sorted.begin(),sorted.end(),[&](double v){return v<limit;});};
                std::cout<<"dt statistics: min="<<sorted.front()
                         <<" median="<<sorted[sorted.size()/2]
                         <<" max="<<sorted.back()
                         <<" cells<1e-9="<<below(1e-9)
                         <<" cells<1e-8="<<below(1e-8)
                         <<" cells<1e-7="<<below(1e-7)
                         <<" total="<<sorted.size()<<std::endl;
                if(cfg.physicalTime)
                {
                    const auto c=static_cast<std::size_t>(std::min_element(dt.begin(),dt.end())-dt.begin());
                    const auto centre=mesh.cellCentroid(c);
                    std::cout<<"PHYSICAL_CFL_LIMIT cell="<<c<<" dt_s="<<dt[c]
                             <<" volume_m3="<<mesh.volume(c)
                             <<" centroid_m=("<<centre.x<<","<<centre.y<<","<<centre.z<<")"
                             <<" boundary_zones=";
                    for(auto f:mesh.cellFaces(c))
                        if(mesh.face(f).neighbour<0)std::cout<<mesh.face(f).boundaryZone<<",";
                    std::cout<<std::endl;
                }
            }
            if(step%cfg.monitorEvery!=0&&step!=cfg.steps)continue;
            // Monitor the accepted state with freshly rebuilt boundary values.
            fv.computeResidual(provider,residual);
            if(steadyMonitor)
            {
                steadyReport=steadyMonitor->evaluate(fv,residual,windowChange,windowLimited,step);
                // An unspecified external pressure is not a physical steady
                // boundary condition for this subsonic startup problem.
                numericallyConverged=!cfg.physicalTime&&steadyReport.converged&&pressureOutlet&&cfg.inletMode=="mass"
                    &&(!cfg.inletFaceSupplyCap||std::all_of(inletZones.begin(),inletZones.end(),[&](int id)
                    {
                        const auto& z=zones.at(id);
                        return z.currentMassFlow>=z.targetMassFlow*(1.-1e-8);
                    }));
            }
            windowChange=0.;windowLimited=0;

            double massResidual=0.0,energyScale=0.0,massSum=0.0;
            std::size_t twoPhaseCells=0;double twoPhaseVolume=0.;
            double weightedResidual=0.0,totalVolume=0.0;
            double pChamber=0.0,tMin=1e30,tMax=-1e30,pMin=1e30,pMax=-1e30,machMax=0.0;
            std::size_t pmaxCell=0;
            std::size_t chamberCells=0;
            for(std::size_t c=0;c<mesh.cellCount();++c)
            {
                massResidual+=residual[c].rho*residual[c].rho;
                energyScale+=std::abs(residual[c].rhoE);
                weightedResidual+=residual[c].rho*residual[c].rho/mesh.volume(c);
                totalVolume+=mesh.volume(c);
                const auto& w=fv.W(c);
                if(fv.Q(c).region==thermo::EquilibriumRegion::TWO_PHASE)
                {++twoPhaseCells;twoPhaseVolume+=mesh.volume(c);}
                massSum+=fv.Q(c).rho*mesh.volume(c);
                const auto centre=mesh.cellCentroid(c);
                const double axialMm=(centre.z-0.271)*1e3;
                if(axialMm>30.0&&axialMm<250.0&&std::hypot(centre.x-0.049215348,centre.y-0.048350049)<0.040)
                {pChamber+=w.p;++chamberCells;}
                tMin=std::min(tMin,w.T);tMax=std::max(tMax,w.T);
                if(w.p>pMax){pMax=w.p;pmaxCell=c;}
                pMin=std::min(pMin,w.p);
                const double speed=std::sqrt(w.U[0]*w.U[0]+w.U[1]*w.U[1]+w.U[2]*w.U[2]);
                machMax=std::max(machMax,speed/fv.Q(c).acoustic.soundSpeed());
            }
            const double norm=std::sqrt(massResidual/std::max(energyScale,1e-30));
            if(referenceMass<0)referenceMass=massSum;
            if(referenceResidual<0&&localStep>0)referenceResidual=norm;
            double lox=0.0,fuel=0.0,torch=0.0;
            for(int id:inletZones)
            {
                const double f=-zoneFlow(zones[id],id);
                if(id>=10&&id<=12)lox+=f;
                else if(id>=20&&id<=22)fuel+=f;
                else torch+=f;
            }
            std::cout<<"STEP "<<step
                     <<" res="<<norm
                     <<" res_vol="<<std::sqrt(weightedResidual/std::max(totalVolume,1e-30))
                     <<(referenceResidual>0?(" ratio="+std::to_string(norm/referenceResidual)):"")
                     <<" mass_kg="<<massSum
                     <<" dMass="<<(massSum-referenceMass)
                     <<" mdot_in="<<lox<<"/"<<fuel<<"/"<<torch
                     <<" mdot_out="<<zoneFlow(zones[kOutletZone],kOutletZone)
                     <<" pChamber_Pa="<<(chamberCells?pChamber/chamberCells:0.0)
                     <<" p=["<<pMin<<","<<pMax<<"]"
                     <<" pmax_at=("<<mesh.cellCentroid(pmaxCell).z-0.271<<","
                     <<std::hypot(mesh.cellCentroid(pmaxCell).x-0.0492153,mesh.cellCentroid(pmaxCell).y-0.04835)<<")"
                     <<" pmax_vol="<<mesh.volume(pmaxCell)
                     <<" T=["<<tMin<<","<<tMax<<"]"
                     <<" Mach_max="<<machMax
                     <<" two_phase_cells="<<twoPhaseCells
                     <<" two_phase_volume_fraction="<<twoPhaseVolume/totalVolume
                     <<(turbulence
                        ?(" nut_max="+std::to_string(turbulence->maxNut())
                          +" turb_wall_cells="+std::to_string(turbulence->wallAdjacentCells()))
                        :"")
#ifdef RFAM_HAVE_CHEMISTRY
                     <<(cfg.reacting
                        ?(" chem_active="+std::to_string(chemistryReport.active)
                          +" chem_limited="+std::to_string(chemistryReport.limited)
                          +" chem_failed="+std::to_string(chemistryReport.failed)
                          +" chem_fail_src="+std::to_string(chemistryReport.failSource)
                          +" chem_fail_nf="+std::to_string(chemistryReport.failNonfinite)
                          +" chem_fail_dt="+std::to_string(chemistryReport.failDt)
                          +" chem_fail_rec="+std::to_string(chemistryReport.failRecover)
                          +" chem_fail_commit="+std::to_string(chemistryReport.failCommit)
                          +" chem_skip_two="+std::to_string(chemistryReport.skippedTwoPhase)
                          +" chem_skip_cold="+std::to_string(chemistryReport.skippedCold)
                          +" chem_substeps="+std::to_string(chemistryReport.substeps)
                          +" chem_rhs_calls="+std::to_string(chemistryReport.rhsCalls)
                          +" chem_time_fraction_min="+chemistryDiagnosticNumber(chemistryReport.minTimeFraction)
                          +" chem_element_defect="+chemistryDiagnosticNumber(chemistryReport.maxElementDefect)
                          +" chem_energy_defect_J_kg="+chemistryDiagnosticNumber(chemistryReport.maxEnergyDefect)
                          +" qdot_max="+chemistryDiagnosticNumber(chemistryReport.maxHeatRelease)
                          +" qdot_total_W="+chemistryDiagnosticNumber(chemistryReport.heatReleaseTotal)
                          +" chem_energy_J="+chemistryDiagnosticNumber(chemistryReport.heatReleaseEnergyJ)
                          +" chem_T=["+chemistryDiagnosticNumber(chemistryReport.minT)+","
                          +chemistryDiagnosticNumber(chemistryReport.maxT)+"]")
                        :std::string())
#endif
                     <<" bflux_fallback="<<fv.boundaryFluxFallbacks.load()
                     <<" dtScale="<<(cfg.lusgsDtScaleStart>0.0?std::min(1.0,double(step)/std::max(cfg.lusgsRampSteps,1)):-1.0)
                     <<" frozen="<<frozenCells.load()
                     <<" frozen_vol_frac="<<frozenVolume/2.958603075278e-3
                     <<" limited="<<limitedCells.load()
                     <<" stiffened="<<stiffenedCells.load()
                     <<" iflux_fb="<<fv.consistentDensityFallbacks()
                     <<" dt_min="<<(dt.empty()?0.0:*std::min_element(dt.begin(),dt.end()))
                     <<" steady_L1_max="<<steadyReport.maxLocal
                     <<" steady_budget_max="<<steadyReport.maxBudget
                     <<" state_change="<<steadyReport.change
                     <<" steady_samples="<<steadyReport.stableSamples
                     <<" physical_time_s="<<physicalTime
                     <<" physical_dt_s="<<(cfg.physicalTime?physicalDtUsed:0.0)
                     <<" physical_step_retries="<<physicalRetryAudit.stepRetries
                     <<" physical_boundary_trial_retries="<<physicalRetryAudit.boundaryTrialRetries
                     <<" physical_inner_iters="<<physicalInnerAcceptedIterations
                     <<" physical_inner_residual_ratio="<<physicalInnerResidualRatio
                     <<" physical_inner_max_residual_ratio="<<physicalInnerMaxResidualRatio
                     <<" physical_inner_residual_absolute="<<physicalInnerResidualAbsolute
                     <<" physical_inner_max_residual_absolute="<<physicalInnerMaxResidualAbsolute
                     <<" physical_budget_error="<<physicalBudgetError
                     <<" physical_cfl_limit_s="<<physicalCflLimit
                     <<" physical_dt_over_cfl="<<(physicalCflLimit>0.?physicalDtUsed/physicalCflLimit:0.)
                     <<(Foam::thermo::thermoCountersEnabled()
                        ?(" thermo_pt="+std::to_string(Foam::thermo::thermoCounters().pt.load())
                          +" isent="+std::to_string(Foam::thermo::thermoCounters().isentrope.load())
                          +" ident="+std::to_string(Foam::thermo::thermoCounters().identity.load())
                          +" recov="+std::to_string(Foam::thermo::thermoCounters().recover.load())
                          +" recovFail="+std::to_string(Foam::thermo::thermoCounters().recoverNewtonFail.load())
                          +" spRepair="+std::to_string(Foam::thermo::thermoCounters().speciesRepairEvents.load())
                          +" pureFallback="+std::to_string(Foam::thermo::thermoCounters().pureFallbackEvents.load())
                          +" pureSnapMax="+std::to_string(Foam::thermo::thermoCounters().maxPureSnapFraction.load())
                          +" spRepairMax="+std::to_string(Foam::thermo::thermoCounters().maxSpeciesRepairRel.load())
                          +" maxRelRho="+std::to_string(Foam::thermo::thermoCounters().maxRelRho.load())
                          +" maxRelE="+std::to_string(Foam::thermo::thermoCounters().maxRelEnergy.load()))
                        :std::string())
                     <<std::endl;
            if(telemetry.active())
            {
                const double outletFlow=zoneFlow(zones[kOutletZone],kOutletZone);
                Foam::telemetry::setStep(step,cfg.steps,acceptedSteps);
                Foam::telemetry::metric(Foam::telemetry::Metric::Residual,norm);
                if(referenceResidual>0.)
                    Foam::telemetry::metric(Foam::telemetry::Metric::ResidualRatio,norm/referenceResidual);
                Foam::telemetry::metric(Foam::telemetry::Metric::TemperatureMin,tMin);
                Foam::telemetry::metric(Foam::telemetry::Metric::TemperatureMax,tMax);
                Foam::telemetry::metric(Foam::telemetry::Metric::PressureMin,pMin);
                Foam::telemetry::metric(Foam::telemetry::Metric::PressureMax,pMax);
                Foam::telemetry::metric(Foam::telemetry::Metric::ChamberPressure,
                                        chamberCells?pChamber/chamberCells:0.0);
                Foam::telemetry::metric(Foam::telemetry::Metric::MachMax,machMax);
                Foam::telemetry::metric(Foam::telemetry::Metric::MdotLox,lox);
                Foam::telemetry::metric(Foam::telemetry::Metric::MdotFuel,fuel);
                Foam::telemetry::metric(Foam::telemetry::Metric::MdotTorch,torch);
                Foam::telemetry::metric(Foam::telemetry::Metric::MdotOut,outletFlow);
                Foam::telemetry::metric(Foam::telemetry::Metric::MassKg,massSum);
                Foam::telemetry::metric(Foam::telemetry::Metric::MassError,massSum-referenceMass);
                Foam::telemetry::metric(Foam::telemetry::Metric::PhysicalTime,physicalTime);
                Foam::telemetry::count(Foam::telemetry::Count::Frozen,frozenCells.load());
                Foam::telemetry::count(Foam::telemetry::Count::Limited,limitedCells.load());
                Foam::telemetry::count(Foam::telemetry::Count::Stiffened,stiffenedCells.load());
                Foam::telemetry::count(Foam::telemetry::Count::TwoPhaseCells,twoPhaseCells);
                Foam::telemetry::count(Foam::telemetry::Count::BoundaryFallback,
                                       fv.boundaryFluxFallbacks.load());
                Foam::telemetry::count(Foam::telemetry::Count::InterfaceFallback,
                                       fv.consistentDensityFallbacks());
                Foam::telemetry::count(Foam::telemetry::Count::ChemAccepted,
                                       static_cast<unsigned long long>(chemistryReport.active));
                Foam::telemetry::count(Foam::telemetry::Count::ChemRejected,
                                       static_cast<unsigned long long>(chemistryReport.failed));
                Foam::telemetry::count(Foam::telemetry::Count::ChemSkippedCold,
                                       static_cast<unsigned long long>(chemistryReport.skippedCold));
                Foam::telemetry::count(Foam::telemetry::Count::ChemSkippedTwoPhase,
                                       static_cast<unsigned long long>(chemistryReport.skippedTwoPhase));
                if(chemistryReport.failed>0)
                    Foam::telemetry::degrade("chemistry.cell_rejected",
                        Foam::telemetry::Severity::Approximation,
                        std::to_string(chemistryReport.failed)
                        +" cells rejected by the chemistry integrator");
                if(frozenCells.load()>0)
                    Foam::telemetry::degrade("solve.frozen_cells",
                        Foam::telemetry::Severity::Approximation,
                        std::to_string(frozenCells.load())+" cells held by the physical guard");
                if(!cfg.physicalTime&&cfg.outletPressure<=0.0)
                    Foam::telemetry::degrade("outlet.zero_gradient",
                        Foam::telemetry::Severity::Approximation,
                        "outlet is a zero-gradient diagnostic boundary, not a physical one");
                Foam::json::Value record=Foam::json::Value::object();
                record.set("step",Foam::json::Value(static_cast<long long>(step)));
                record.set("residual",Foam::json::Value(norm));
                if(referenceResidual>0.)
                    record.set("ratio",Foam::json::Value(norm/referenceResidual));
                record.set("mass_kg",Foam::json::Value(massSum));
                record.set("mdot_out",Foam::json::Value(outletFlow));
                record.set("chem_cells",Foam::json::Value(static_cast<long long>(chemistryReport.active)));
                record.set("chem_rejected",Foam::json::Value(static_cast<long long>(chemistryReport.failed)));
                telemetry.event("STEP",record);
            }

            if(cfg.inletFaceSupplyCap)
                for(int id:inletZones)
                {
                    const auto& zone=zones.at(id);
                    if(zone.targetMassFlow>0.)
                        std::cout<<"INLET_FACE_CAP step="<<step<<" zone="<<id
                                 <<" blocked_area_frac="<<zone.blockedArea/zone.area
                                 <<" available_kg_s="<<zone.availableCapacity
                                 <<" delivered_kg_s="<<zone.currentMassFlow<<std::endl;
                }
            if(cfg.dumpEvery>0&&!cfg.dump.empty()&&step>0&&step%cfg.dumpEvery==0)
            {
                const auto path=cfg.dump+".s"+std::to_string(step);
                writeFieldDump(path);
                if(cfg.dumpPhaseEvery)
                {
                    writePhaseSidecar(path);
                    writeBoundaryBudget(path);
                    std::ofstream metadata(path+".status.json");
                    if(!metadata)throw std::runtime_error("cannot write periodic status metadata");
                    metadata<<"{\"steps\":"<<step<<",\"converged\":"
                            <<(numericallyConverged?"true":"false")<<","
                            <<"\"production_validated\":false,\"phase_sidecar_available\":true,"
                            <<"\"boundary_budget_available\":true,"
                            <<"\"physical_time_mode\":"<<(cfg.physicalTime?"true":"false")<<","
                            <<"\"physical_time_s\":"<<std::setprecision(17)<<physicalTime<<","
                            <<"\"physical_dt_s\":"<<physicalDtUsed<<"," 
                            <<"\"physical_max_retries\":"<<cfg.physicalMaxRetries<<"," 
                            <<"\"physical_probe_inner\":"<<(cfg.physicalProbeInner?"true":"false")<<"," 
                            <<"\"outlet_roundoff_shortcut\":"<<(cfg.outletRoundoffShortcut?"true":"false")<<"," 
                            <<"\"physical_step_retries\":"<<physicalRetryAudit.stepRetries<<"," 
                            <<"\"physical_boundary_trial_retries\":"<<physicalRetryAudit.boundaryTrialRetries<<"," 
                            <<"\"physical_retry_count_total\":"<<physicalRetryAudit.totalRetries<<"," 
                            <<"\"physical_boundary_trial_retry_count_total\":"<<physicalRetryAudit.totalBoundaryTrialRetries<<"," 
                            <<"\"physical_inner_iters\":"<<physicalInnerAcceptedIterations<<"," 
                            <<"\"physical_inner_tolerance\":"<<cfg.physicalInnerTolerance<<"," 
                            <<"\"physical_inner_residual_ratio\":"<<std::setprecision(17)<<physicalInnerResidualRatio<<","
                            <<"\"physical_inner_max_residual_ratio\":"<<physicalInnerMaxResidualRatio<<","
                            <<"\"physical_inner_residual_absolute\":"<<physicalInnerResidualAbsolute<<","
                            <<"\"physical_inner_max_residual_absolute\":"<<physicalInnerMaxResidualAbsolute<<","
                            <<"\"physical_budget_error\":"<<physicalBudgetError<<","
                            <<"\"physical_cfl_limit_s\":"<<physicalCflLimit<<","
                            <<"\"physical_dt_over_cfl\":"<<(physicalCflLimit>0.?physicalDtUsed/physicalCflLimit:0.)<<","
                            <<"\"physical_inner_converged\":"<<(physicalInnerConverged?"true":"false")<<","
                            <<"\"physical_boundary_input_rejected\":"<<(physicalBoundaryInputRejected?"true":"false")<<","
                            <<"\"projection_accept\":"<<(Foam::projectionAcceptEnabled()?"true":"false")
                            <<",\"molecular_viscous\":"<<(cfg.molecularViscous?"true":"false")
                            <<",\"interface_consistent_flux\":"<<(cfg.interfaceConsistentFlux?"true":"false")
                            <<",\"two_phase_dt_scale\":"<<cfg.twoPhaseDtScale
                            <<",\"injector_dt_scale\":"<<cfg.injectorDtScale
                            <<"}\n";
                }
            }
            if(numericallyConverged)break;
        }
        // Recompute accepted-state fluxes after a rejected step as well.
        fv.computeResidual(provider,residual);

        if(!cfg.dump.empty())
        {
            writeFieldDump(cfg.dump);
            const std::uint64_t cells=mesh.cellCount();
            const std::uint64_t speciesCount=species.size();
            writePhaseSidecar(cfg.dump);
            writeBoundaryBudget(cfg.dump);
            std::cout<<"FIELD_DUMP_WRITTEN "<<cfg.dump<<" cells="<<cells
                     <<" species="<<speciesCount<<std::endl;
        }
        std::cout<<(advanceFailed?"TETRA_COLD_FLOW_FAILED":numericallyConverged?
                    "TETRA_NUMERICAL_STEADY_PASS":"TETRA_COLD_FLOW_ITERATION_LIMIT")
                 <<" accepted_steps="<<acceptedSteps
                 <<" converged="<<numericallyConverged<<" production_validated=0"<<std::endl;
        double supplyEntropyMin=1e300;bool supplyEntropyCompatible=true;
        for(int id:inletZones)
        {
            const auto& zone=zones.at(id);
            if(!(zone.targetMassFlow>0.)||!std::isfinite(zone.entropyMargin))continue;
            supplyEntropyMin=std::min(supplyEntropyMin,zone.entropyMargin);
            if(zone.entropyMargin<-cfg.supplyEntropyTolerance*std::max(std::abs(zone.supplyEntropy),1.))
                supplyEntropyCompatible=false;
        }
        const bool designFlowSatisfied=std::all_of(inletZones.begin(),inletZones.end(),[&](int id)
        {
            const auto& z=zones.at(id);
            return z.currentMassFlow>=z.targetMassFlow*(1.-1e-8);
        });
        const bool boundaryStateValid=!boundaryRejected&&supplyEntropyCompatible
            &&(!cfg.inletFaceSupplyCap||designFlowSatisfied);
        if(!cfg.dump.empty())
        {
            std::ofstream status(cfg.dump+".status.json");
            status<<std::setprecision(17)<<"{\"converged\":"<<(numericallyConverged?"true":"false")
                  <<",\"production_validated\":false,\"steps\":"<<acceptedSteps
                  <<",\"chemistry_integrator\":\""<<(cfg.chemCVODE?"cvode-bdf-full-time":"limited-explicit-diagnostic")<<"\""
                  <<",\"physical_time_mode\":"<<(cfg.physicalTime?"true":"false")
                  <<",\"physical_time_s\":"<<physicalTime
                  <<",\"physical_dt_s\":"<<physicalDtUsed
                  <<",\"physical_max_retries\":"<<cfg.physicalMaxRetries
                  <<",\"physical_probe_inner\":"<<(cfg.physicalProbeInner?"true":"false")
                  <<",\"outlet_roundoff_shortcut\":"<<(cfg.outletRoundoffShortcut?"true":"false")
                  <<",\"outlet_roundoff_hits\":"<<(pressureOutlet?pressureOutlet->roundoffHits():0)
                  <<",\"physical_step_retries\":"<<physicalRetryAudit.stepRetries
                  <<",\"physical_boundary_trial_retries\":"<<physicalRetryAudit.boundaryTrialRetries
                  <<",\"physical_retry_count_total\":"<<physicalRetryAudit.totalRetries
                  <<",\"physical_boundary_trial_retry_count_total\":"<<physicalRetryAudit.totalBoundaryTrialRetries
                  <<",\"physical_inner_iters\":"<<physicalInnerAcceptedIterations
                  <<",\"physical_inner_tolerance\":"<<cfg.physicalInnerTolerance
                  <<",\"physical_inner_residual_ratio\":"<<physicalInnerResidualRatio
                  <<",\"physical_inner_max_residual_ratio\":"<<physicalInnerMaxResidualRatio
                  <<",\"physical_inner_residual_absolute\":"<<physicalInnerResidualAbsolute
                  <<",\"physical_inner_max_residual_absolute\":"<<physicalInnerMaxResidualAbsolute
                  <<",\"physical_budget_error\":"<<physicalBudgetError
                  <<",\"physical_cfl_limit_s\":"<<physicalCflLimit
                  <<",\"physical_dt_over_cfl\":"<<(physicalCflLimit>0.?physicalDtUsed/physicalCflLimit:0.)
                  <<",\"physical_inner_converged\":"<<(physicalInnerConverged?"true":"false")
                  <<",\"physical_boundary_input_rejected\":"<<(physicalBoundaryInputRejected?"true":"false")
                  <<",\"lusgs_composition_coupling\":"<<(cfg.lusgsCompositionCoupling?"true":"false")
                  <<",\"inlet_face_supply_cap\":"<<(cfg.inletFaceSupplyCap?"true":"false")
                  <<",\"blueprint_plate_wall\":"<<(cfg.blueprintPlateWall?"true":"false")
                  <<",\"molecular_viscous\":"<<(cfg.molecularViscous?"true":"false")
                  <<",\"inlet_cd\":"<<cfg.inletCd
                  <<",\"state_acceptance_rel\":"<<Foam::stateAcceptanceRel()
                  <<",\"projection_accept\":"<<(Foam::projectionAcceptEnabled()?"true":"false")
                  <<",\"species_repair_rel\":"<<Foam::speciesRepairRel()
                  <<",\"inlet_supply_cap\":"<<(cfg.inletSupplyCap?"true":"false")
                  <<",\"inlet_design_flow_satisfied\":"<<(designFlowSatisfied?"true":"false")
                  <<",\"interface_consistent_flux\":"<<(cfg.interfaceConsistentFlux?"true":"false")
                  <<",\"cfl\":"<<cfg.cfl
                  <<",\"lusgs_dt_scale\":"<<cfg.lusgsDtScale
                  <<",\"lusgs_dt_scale_start\":"<<cfg.lusgsDtScaleStart
                  <<",\"lusgs_ramp_steps\":"<<cfg.lusgsRampSteps
                  <<",\"lusgs_relax\":"<<cfg.lusgsRelax
                  <<",\"two_phase_dt_scale\":"<<cfg.twoPhaseDtScale
                  <<",\"injector_dt_scale\":"<<cfg.injectorDtScale
                  <<",\"blueprint_injector\":"<<(cfg.blueprintInjector?"true":"false")
                  <<",\"blueprint_layers\":"<<cfg.blueprintLayers
                  <<",\"blueprint_blend\":"<<cfg.blueprintBlend
                  <<",\"steady_monitor_available\":"<<(steadyMonitor?"true":"false")
                  <<",\"steady_L1_max\":"<<steadyReport.maxLocal
                  <<",\"steady_budget_max\":"<<steadyReport.maxBudget
                  <<",\"steady_stable_samples\":"<<steadyReport.stableSamples
                  <<",\"supply_entropy_checked\":"<<(cfg.thermoBackend=="coolprop"&&cfg.inletMode=="mass"?"true":"false")
                  <<",\"boundary_state_valid\":"<<(boundaryStateValid?"true":"false")
                  <<",\"boundary_transaction_rejected\":"<<(boundaryRejected?"true":"false")
                  <<",\"supply_entropy_enforced\":"<<(cfg.supplyEntropyStrict?"true":"false")
                  <<",\"supply_entropy_compatible\":"<<(supplyEntropyCompatible?"true":"false")
                  <<",\"supply_entropy_min_J_kgK\":"<<(supplyEntropyMin<1e299?std::to_string(supplyEntropyMin):"null")
                  <<",\"update_rejected\":"<<(advanceFailed?"true":"false")
                  <<",\"thermo_backend\":\""<<cfg.thermoBackend
                  <<"\",\"thermo_model_id\":\""<<model->modelId()
                  <<"\",\"model\":\""<<(cfg.molecularViscous
                        ?"single-phase molecular viscous, no-slip adiabatic tetra diagnostic"
                        :"inviscid slip-wall tetra prototype")
                  <<"\",\"torch_mdot_kg_s\":"
                  <<cfg.torchMassFlow<<",\"inlet_mode\":\""<<cfg.inletMode
                  <<"\",\"seeded_from_rounded_field\":"<<(!cfg.seedField.empty()?"true":"false")
                  <<",\"exact_restart\":"<<(!cfg.physicalResume.empty()?"true":"false")
                  <<",\"outlet_model\":\""
                  <<(pressureOutlet?"equilibrium-pressure":"zero-gradient diagnostic")
                  <<"\",\"outlet_pressure_Pa\":"<<cfg.outletPressure
                  <<",\"backflow_temperature_K\":"<<cfg.backflowTemperature
                  <<",\"backflow_Y_O2\":"<<cfg.backflowY[0]<<",\"backflow_is_air\":false}\n";
        }
        RFAM_ENTER(react_finish);
        if(telemetry.active())
        {
            Foam::json::Value summary=Foam::json::Value::object();
            summary.set("advance_failed",Foam::json::Value(advanceFailed));
            summary.set("numerically_converged",Foam::json::Value(numericallyConverged));
            summary.set("accepted_steps",Foam::json::Value(static_cast<long long>(acceptedSteps)));
            summary.set("boundary_state_valid",Foam::json::Value(boundaryStateValid));
            summary.set("mechanism",Foam::json::Value(cfg.mechanism));
            summary.set("reacting",Foam::json::Value(cfg.reacting));
            summary.set("physical_time_s",Foam::json::Value(physicalTime));
            telemetry.finish(advanceFailed?"failed"
                             :(numericallyConverged?"pass":"iteration_limit"),
                             advanceFailed?1:0,summary);
        }
        return advanceFailed?1:0;
    }
    catch(const std::exception& e)
    {
        std::cerr<<"TETRA_COLD_FLOW_FAILED: "<<e.what()<<"\n";
        if(Foam::telemetry::Session::instance().active())
        {
            Foam::json::Value summary=Foam::json::Value::object();
            summary.set("error",Foam::json::Value(std::string(e.what())));
            summary.set("mechanism",Foam::json::Value(cfg.mechanism));
            Foam::telemetry::Session::instance().finish("failed",1,summary);
        }
        return 1;
    }
}
