# RFAM + preCICE thermal coupling smoke test

This is a two-process integration test using real preCICE 3 (socket transport).
`Fluid` advances a small, low-Mach 1D single-phase oxygen strip with the RFAM
AUSM+-up flux. `Solid` advances a 1D conducting slab. The interface uses a film
model at one coupling vertex. The program exchanges integrated heat power (W),
with a conservative mapping, and solid surface temperature (K), with a consistent
mapping. The fluid accounts for the energy crossing the two ends separately.

The fluid assembles `R = div(F) - S` and iterates the backward-Euler equation
`U - U_start + dt*R(U) = 0`. The Mach-dependent pressure-mode scaling acts on the
inner corrections only. The stop criterion uses the original physical-step
residual. The pressure-mode multiplier is `(1 + sigma)/(1 + sigma/beta^2)`, with
`sigma = dt * acoustic_rate` and `beta = max(Mach, mach_floor)` for subsonic
flow. This retains the physical-time term and tends to unity for short physical
steps. With a Mach floor of `1.0` or a local Mach >= 1 the scaling is the identity.
The implementation uses single-phase EOS derivatives at frozen composition; it
is not the Turkel–Weiss implementation.

Each preCICE time window takes one physical step in this example. On an implicit
coupling retry, each participant restores its own window-start checkpoint. Only
accepted windows add to the reported time and energy. The Python test launches
both participants in a private working directory, checks energy balance and real
checkpoint restores, and compares results with the inner scaling on and off.

```sh
cd /home/coffee/rfam-cur
# Replace the prefix if using a different preCICE installation.
cmake -S . -B build -DRFAM_ENABLE_PRECICE=ON \
    -DCMAKE_PREFIX_PATH="$HOME/.local/opt/precice-rfam"
cmake --build build -j4 --target rfam_precice_thermal \
    test_residual_integrator test_mach_residual_integrator test_cht2d
OMP_NUM_THREADS=1 ctest --test-dir build \
    -R '^(test_cht2d|test_(mach_)?residual_integrator|test_precice_residual_coupling)$' \
    --output-on-failure
```

The default tests include a 15-case single-phase sweep (Mach 0.005, 0.05, 0.3,
0.9, 1.2; acoustic CFL 0.025, 0.18, 0.4). For larger steps, run the diagnostic
below. A rejected step is reported as such and checked for rollback; the
stress diagnostic does not require all steps to be accepted.

```sh
OMP_NUM_THREADS=1 build/test_mach_residual_integrator --stress
```

On the 3-window preCICE case the worst inner count fell from 616 with the original
`beta^2` pressure multiplier to 108 with physical-time weighting, versus 20 with
scaling disabled. At Mach 0.05 and acoustic CFL 0.8 the unscaled case converges
but the scaled case does not reach tolerance within 3000 inner iterations. At
Mach 0.9 and CFL 1.6 the counts are 69 unscaled and 48 scaled; at Mach 1.2 both
are 77. These results are not evidence of a general speedup.

This tests integration and single-phase iterative behavior, not a coolant
pressure drop, high-Mach accuracy, finite-rate chemistry, or the tetrahedral
production path.
