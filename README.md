# realFluidAllMach

**Research implementation toward an auditable real-fluid reacting-flow solver**

The repository name describes the target architecture. The current code does
not yet contain a Turkel/Weiss residual/Jacobian preconditioner, a validated
two-phase model, or a monolithic production CHT solve. Status claims below are
limited to implemented and tested capabilities.

## Philosophy

This is NOT "another PIMPLE variant" or "density-based with a patch."

The target is a **unified all-Mach conservative solver**. The current
standalone implementation provides tested building blocks toward that target:

- **Conserved variables**: ρ, ρu, ρE, ρY (foundation)
- **Primitive variables**: p, T, u, Y (coupled, not sequential)
- **Real-fluid EOS**: PR/SRK with thermodynamic consistency (not scattered inline functions)
- **All-speed numerical flux**: AUSM+-up with low-Mach correction is
  implemented; HLLC and KT remain architecture options
- **Combustion**: FPV/FGM as **source terms** (not T lookup override)
- **CHT**: Monolithic fluid-solid-coolant coupling (not partitioned exchange)

## Design Principles

### 1. Thermodynamic Consistency
One state object holds all consistent properties:
```cpp
RealFluidState state(p, T, Y);
// Computes once, caches: rho, h, e, Cp, Cv, Z, a², derivatives
// All from unified PR EOS + departure functions
```

No more:
- `rho()` in one place
- `Cp()` somewhere else  
- `Cv()` negative somewhere
- `H()` inconsistent

### 2. Energy Equation is Primary
**NOT**: `h,p,Y → T` inversion as main path (病态)

**YES**: `ρE` + EOS + energy consistency → all states simultaneously

Inversion `T(h,p)` exists only as **auxiliary tool**, not main loop.

### 3. Combustion as Source, Not Override
**NOT**: `T = fpvTable.lookup(f, C)`

**YES**: 
```cpp
qdot = fpvTable.getHeatRelease(f, C)
omega_i = fpvTable.getSpeciesSource(f, C)
// T determined by energy equation + EOS
```

### 4. All-Speed Flux Today; Matrix Preconditioning Gated

The implemented AUSM+-up flux includes its standard low-Mach correction. An
optional `ExperimentalCflAcousticScaler` changes only the explicit CFL
estimate and is disabled by default. It does not transform the flux residual,
time derivative, eigenstructure, or Jacobian and is not Turkel/Weiss all-Mach
preconditioning.

A real-fluid preconditioning matrix, dual-time integration, and low-Mach
accuracy/convergence validation remain future work. Production and validation
reports must not infer those capabilities from the CFL scaler. See
`docs/all_mach_preconditioning_status.md`.

### 5. CHT as Monolithic System
One global nonlinear residual:
```
F(fluid, solid, coolant) = 0
q_gas = q_solid = q_coolant
T_wall continuous
```

Not: "fluid solver → interface BC → solid solver → exchange"

## Project Status

**Current baseline**: standalone C++17 core with real-fluid thermo, an AUSM+-up
all-speed flux, reacting-flow operators, CHT components, and Vortex bridge
tests. The experimental CFL scaler has a scope-contract test; full all-Mach
matrix preconditioning is not implemented or validated.

The production OpenFOAM executable described by the long-term architecture is
not yet present. The repository currently builds and validates the standalone
core. Registered CTest targets plus `scripts/validation_ladder.sh` are the
authoritative Core contract suite; `scripts/local_test.sh` is a complementary
standalone/legacy integration sweep and must not replace the validation ladder.

Freeze reports should use `tools/artifact_composite.py` instead of listing
hand-computed hashes for every source and generated artifact. Generation and
verification are silent on success; a mismatch returns one structured error.
See `docs/artifact_composite.md` for the closed path and manifest contract.

**Stage 0 foundation**: RealFluidState unit tests
- PR cubic root selection
- Gibbs energy root selection  
- Cv/Cp positivity
- Sound speed positivity
- T(h,p) bracketed Newton

**Stage 1** (Next): Cold flow real-fluid solver
- Conserved variable framework
- AUSM+-up flux
- injector boundary/mapping diagnostics; transcritical advancement remains
  blocked until phase-aware thermo and transport are validated

**Future Stages**: residual/Jacobian all-Mach preconditioning -> phase-aware
combustion -> validated CHT coupling -> full thrust chamber

## Repository Structure

```
realFluidAllMach/
├── src/              # Core solver source
│   ├── thermo/       # Real-fluid thermodynamic kernel (地基)
│   ├── variables/    # Conserved/primitive fields
│   ├── numerics/     # Flux, time integration, experimental CFL scaling
│   ├── combustion/   # FPV/FGM source-term models
│   ├── cht/          # Fluid-solid-coolant strong coupling
│   └── solver/       # Main solver loop
├── apps/             # Executable drivers (solver_3d_tetra, solver_3d_sector, …)
│   └── _unbuilt/     # drivers referenced by neither build system
├── tests/            # C++ tests: CTest suite + the local regression suite
│   ├── thermo/  data/  reference/  *.cmake
│   └── _unbuilt/     # tests referenced by neither build system
├── cases/            # Test/validation/production cases
├── config/           # Species data, chemistry tables
├── scripts/          # build / test / validation-ladder drivers
├── tools/            # Python pre/post-processing and audits
├── docs/             # Architecture & design docs (docs/history = dated records)
├── work/  runs/  rfam_work/   # run output and scratch (mostly untracked)
└── _archive/         # superseded material (read-only)
```

## Key Design Files

- `docs/architecture.md`: Overall system design
- `docs/thermo_pr_design.md`: PR EOS kernel philosophy
- `docs/governing_equations.md`: Conserved form + all-Mach preconditioning
- `src/thermo/RealFluidState.H`: **地基** - unified state object

## Building

```bash
cd realFluidAllMach
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

## Running Tests

```bash
cd build
ctest --output-on-failure

# Full standalone regression suite (WSL/Linux)
cd ..
bash scripts/local_test.sh
```

OpenMP CPU parallelism is enabled by default when the compiler supports it.
Disable it with `-DRFAM_ENABLE_OPENMP=OFF`; select the runtime thread count with
`OMP_NUM_THREADS`.

The current V20 isothermal driver smoke/pathfinder run is:

```bash
OMP_NUM_THREADS=8 ./build/rfam_v20_sector_smoke --grid pathfinder --steps 50
```

The sector driver supports `smoke` (60x20x12), `pathfinder` (120x30x20),
and `intermediate` (240x50x36) grid levels. The legacy `--pathfinder` flag
remains available; `--intermediate` is an alias for `--grid intermediate`.

## Running Cases

The unstructured tetrahedral solvers (cold flow and reacting, both with
optional Menter k-omega SST) are documented in a single place:

**`docs/TETRA_SOLVER_GUIDE.md`** -- build, full flag reference, turbulence model
limits, transport closures, and the configuration constraints that change the
answer (liquid LOX, cold-chamber design flow, pressure-outlet cost).

```bash
# cold flow + SST
./build_wsl/rfam_v20_tetra --mesh <grid.rfamtet> --thermo coolprop \
  --inlet-mode mass --molecular-viscous --turbulence sst --lusgs

# hot fire + SST (needs SUNDIALS; A406 has it)
./build-a406/rfam_v20_tetra_react --mesh <grid.rfamtet> --thermo legacy \
  --reacting --mechanism twostep --inlet-mode mass \
  --molecular-viscous --turbulence sst --lusgs
```

Structured cases still run through the case driver:

```bash
cd cases/verification/transcritical_injection_2D
python ../../../tools/run/run_case.py
```

## Dependencies

- C++17
- A C++17 compiler and CMake 3.16+
- SUNDIALS for the CVODE combustion regression tests
- OpenFOAM-11, Eigen3 and HDF5 are planned production-solver dependencies;
  they are not required by the current standalone core build
- Python 3.9+ (tools)

## License

MIT (to be determined)

## Contact

Built for real rocket thrust chamber simulation.
Not a toy, not a patch, a platform.

## Acknowledgments

Born from lessons learned in:
- OF11 PengRobinsonGas Cv崩溃 (2026-06-20)
- CHT多区网格6个bug叠加
- FPV温度跳变导致密度塌陷
- 临界点附近数值奇异性

每一个坑都是设计输入。
