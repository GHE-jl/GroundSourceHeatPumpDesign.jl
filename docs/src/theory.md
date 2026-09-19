```@meta
CurrentModule = GroundSourceHeatPumpDesign
```

# Modeling theory

## Heat pump performance interpolation

Manufacturer capacity (CAP) and coefficient of performance (COP) data are given at discrete
source (and, for the 2D overload, load) entering fluid temperatures (EFTs).
[`heat_pump_performance`](@ref) recovers a value at an arbitrary EFT by interpolation, then clamps
outside the provided grid (no extrapolation).

**1D** — linear interpolation between the two grid points bracketing `Tsrc`:

```math
w = \frac{T_{src} - T_i}{T_{i+1} - T_i}, \qquad
\mathrm{CAP}(T_{src}) = (1-w)\,\mathrm{CAP}_i + w\,\mathrm{CAP}_{i+1}
```

with the same weight `w` applied to `COP`.

**2D** — bilinear interpolation on the rectangle bracketing `(Tsrc, Tload)`, with `f` standing for
either `CAP` or `COP` and `f₁₁, f₁₂, f₂₁, f₂₂` its values at the rectangle's four corners:

```math
w_x = \frac{T_{src} - T_{src,i}}{T_{src,i+1} - T_{src,i}}, \qquad
w_y = \frac{T_{load} - T_{load,j}}{T_{load,j+1} - T_{load,j}}
```

```math
f(T_{src}, T_{load}) = (1-w_x)(1-w_y)f_{11} + (1-w_x)w_y f_{12} + w_x(1-w_y)f_{21} + w_x w_y f_{22}
```

## Building load to ground load conversion

[`ground_load`](@ref)'s sign convention (see the README's "Conventions" section) follows directly
from an energy balance across the heat pump. In heating mode, the building load `Qb > 0` is met by
compressor work `W` plus heat drawn from the ground loop (`-Qg`, since `Qg < 0` when the ground is
being cooled): `Qb = W - Qg`. With `COP_heating = Qb / W`:

```math
Q_g = -Q_b\left(1 - \frac{1}{\mathrm{COP}_{heating}}\right)
```

In cooling mode, the building load `Qb < 0` (magnitude `|Qb|` of heat removed from the building) is
rejected to the ground loop along with the compressor work: `Qg = |Qb| + W`. With
`COP_cooling = |Qb| / W` and `|Qb| = -Qb`:

```math
Q_g = -Q_b\left(1 + \frac{1}{\mathrm{COP}_{cooling}}\right)
```

`ground_load` has three methods addressing the same underlying difficulty in increasing order of
completeness: the true COP and the installed capacity both depend on the source EFT, which in turn
depends on the ground load being computed — a circular dependency (see
[`gshp_simulation`](@ref)'s docstring for how the simulation resolves it without a full outer
iteration).

1. **Fixed COP** applies the equations above directly, with a single nominal COP and no capacity
   limit — appropriate once the operating COP is already known (e.g. read off the curve below at a
   previously-computed EFT), or as a first design estimate.
2. **Fixed COP, capacity-capped** clamps `Qb` to `[-n_HP * CAP_cooling_unit, n_HP *
   CAP_heating_unit]` before applying method 1, reporting the clamped-off remainder as `Qb_unmet`.
   This is the method used for the non-iterative simulation passes (see [`gshp_simulation`](@ref)),
   since the true EFT (and thus true COP) isn't known until the ground response has been simulated.
3. **Performance-curve-based** replaces the fixed COP with a lookup from
   [`heat_pump_performance`](@ref) at the current EFT(s), capping the heat pump's contribution at
   `n_HP` times the resulting CAP. This is only exact once the EFT is already known (e.g. as a
   diagnostic after a simulation pass, or a nominal case with the source EFT already given).

Both capacity-aware methods take the single heat pump unit's CAP (a fixed value, or the
temperature-dependent curve) and an `n_HP` keyword (default `1`) that scales it to the total
installed capacity -- the same `n_HP` used by [`n_active_heat_pumps`](@ref) below.

## Heat pump staging

[`n_active_heat_pumps`](@ref) reports how many of the `n_HP` installed units are needed to deliver
a given load `Qb_hp` at a source EFT `Tsrc`, following Pasquier et al. (2013):

```math
n(t) = \mathrm{clamp}\!\left(\frac{|Q_{b,hp}(t)|}{\mathrm{CAP}(T_{src}(t))},\ 0,\ n_{HP}\right)
```

This is a **fractional**, time-averaged count rather than an integer: within a time step, a heat
pump typically cycles on and off rather than switching exactly on the step boundary, so the
fraction of the step it spends running is a better estimate of its duty than rounding `n(t)` up or
down would be.

## GSHP simulation, without a full outer iteration

Simulating the ground temperature needs the ground load, which itself depends on the source EFT
through the heat pump's COP/CAP — the same circular dependency noted above, now at the scale of a
full time series rather than one time step. The fully correct, iterative approach (e.g. Pasquier et
al. 2013) resolves it by iterating over the whole pipeline (load conversion, staging, ground
response) until everything is self-consistent.

[`gshp_simulation`](@ref) keeps each case's resolution of that circularity as small and local as
it can, rather than iterating the whole pipeline:
- With **fixed** COP/CAP (case 1), there is no circularity to resolve — COP/CAP don't depend on the
  source EFT, so a single pass suffices.
- With **temperature-dependent** COP/CAP curves (case 2), a first pass converts the building load
  with a fixed nominal COP (method 2 above) to get an initial `Qg`/`T_out`; a small fixed-point loop
  then alternates the curve-based conversion (method 3 above, at the previous iteration's `T_out`)
  with the ground response until `T_out` stops moving — typically well under the 5-iteration safety
  cap (mirroring `fNSCiter.m`'s identically-capped loop, which recomputes COP/CAP from the current
  source EFT every pass). [`n_active_heat_pumps`](@ref) then reports the number of units operating
  and the COP actually achieved at that converged EFT, self-consistently.
- With a **non-stationary** flow rate (case 3), the active flow rate at each time step is picked
  from the ground load alone, never from the simulated temperature, so it stays a single pass.

See the README's "Simulation cases" section for which cases are implemented and which are planned.
