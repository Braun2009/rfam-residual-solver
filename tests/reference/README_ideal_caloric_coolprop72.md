# CoolProp 7.2.0 ideal-caloric reference

`ideal_caloric_coolprop72.csv` is a frozen mass-specific ideal-gas caloric
reference for the low-temperature range not covered by the current NASA-7
tables:

- O2: 80-300 K
- C3H8: 80-350 K

The O2 grid uses 5 K spacing. The C3H8 grid uses 2.5 K spacing from 80 to
150 K, where the ideal-gas heat-capacity curve has stronger curvature, and
5 K spacing from 150 to 350 K. The interpolation accuracy criteria are not
relaxed for the extended range.

Generate or verify it with exactly CoolProp 7.2.0:

```text
python tools/generate_ideal_caloric_coolprop72.py
python tools/generate_ideal_caloric_coolprop72.py --check
```

The generator queries only `HMASS_IDEALGAS` and `CP0MASS` through `PropsSI`.
Each enthalpy is stored as
`HMASS_IDEALGAS(T) - HMASS_IDEALGAS(298.15 K)` for the same species. The
298.15 K point is an explicit table node. Floating-point fields use 17-digit
scientific notation, ASCII, and LF line endings; `--check` requires a
byte-for-byte match.

The checked-in C++ backend uses cubic Hermite interpolation with the frozen
enthalpy and heat capacity as endpoint value and derivative. Consequently the
returned `cp` is the analytic derivative of the returned sensible enthalpy,
not a separately interpolated curve. The generator audits source derivatives,
pressure independence, interpolation error, positive heat capacity, and
monotone enthalpy before writing the file. The C++ test independently checks
all frozen nodes, every knot, and dense points across both ranges.

This data supplements only the ideal caloric term used with Peng-Robinson. It
does not supply liquid or two-phase transport, phase equilibrium, latent heat,
or HEOS residual properties. It is not connected to a production capability
gate by this change.
