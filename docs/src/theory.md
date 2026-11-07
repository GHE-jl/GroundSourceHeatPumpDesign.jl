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

## Part-load degradation

The fractional unit of `n(t)` cycles on and off within the time step, and its COP is lower than at
steady state. With the part-load ratio `PLR = n - ⌊n⌋` and a degradation coefficient `C_d`,
[`part_load_COP`](@ref) gives the effective COP of the `n` running units:

```math
\mathrm{PLF} = 1 - C_d\,(1 - \mathrm{PLR}), \qquad
\mathrm{COP}_{eff} = \mathrm{COP}\,\frac{n}{\lfloor n \rfloor + \mathrm{PLR}/\mathrm{PLF}}
```

The `⌊n⌋` units at full load keep their COP, and only the cycling one is degraded. With `C_d = 0`
(the default of every [`gshp_simulation`](@ref) case), `COP_eff = COP`.

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
  from the ground load alone, never from the simulated temperature, so it stays a single pass. With
  the COP/CAP curves (case 3b), the same fixed-point loop as case 2 is used; after the first
  iteration, a flow rate can only increase, so that the switching between the discrete flow rates
  does not prevent the convergence.

The fluid temperatures come from the ground temperature at the borehole wall and the effective
borehole resistance (or, for the non-stationary cases, from the outlet response of each flow rate):

```math
T_{out} = T_0 + Q_g * h_{out}, \qquad T_{in} = T_{out} + \frac{Q_g}{n_b \dot V \rho_f c_f}
```

with `*` the temporal superposition and, for `outlet_model = :resistance`,
`h_out = (Rb* + g) / (n_b H) - 1 / (2 n_b V̇ ρ_f c_f)`.

## Source-side mitigation

With operating limits `Tlim` (cases 4 and 5), the temperature at `Tlim_at` (seen from the ground:
`:inlet` = `T_in`, the default; `:outlet` = `T_out`; `:mean`) is computed, and at the time steps
where it leaves the limits the ground load is corrected following Pasquier et al. (2013):

```math
\Delta Q_g(t) = \frac{T_{lim} - T(t)}{a}, \qquad
a = \frac{R_b^* + g(\Delta t)}{n_b H} + \frac{w - 1/2}{n_b \dot V \rho_f c_f}
```

where `a` is the temperature response to a unit load over one time step, and `w = 1, 0, 1/2` for
`:inlet`, `:outlet`, `:mean`. The correction only shaves the load toward zero (never reverses it).
Since the later time steps feel the corrected load, the temperature is recomputed and the
correction repeated until the temperature sits on the limit, which takes a few passes. The shaved
part `Qs_unmet = Qg_hp - Qg` must be covered on the source side (e.g. a boiler or a dry cooler).

In case 5, before the mitigation, the largest flow rate is forced at the time steps where the
temperature is within `T_margin` of the limits (as in `fNSCiter.m`): the bigger flow rate reduces
the `T_in`/`T_out` split, so that less load has to be shaved.

## Effective borehole resistance at the fluid temperature

With the `fluid` keyword, the fluid properties and `Rb*` are recalculated at the mean fluid
temperature `T_f = (T_in + T_out)/2` of the previous pass. Since `Rb*` only enters the fluid
temperature instantaneously, its change is added as a resistance at each time step:

```math
T_f(t) = T_0 + Q_g * \frac{R_b^*(T_0) + g}{n_b H} + Q_g(t)\,\frac{R_b^*(T_f(t)) - R_b^*(T_0)}{n_b H}
```

and the simulation is repeated until `T_f` changes by less than `tol` (a few passes). For the
non-stationary flow, the same holds with `Rb*` of the active flow rate.

## Pumping power and performance indicators

[`pumping_power`](@ref) gives the pressure drop of one U-loop by the Darcy-Weisbach equation, with
the minor losses of one 180° bend and two 90° elbows (`K = 1.5 + 2 × 0.9`):

```math
\Delta P = \left(f\,\frac{2H}{2 r_i} + K\right)\frac{\rho_f v^2}{2}, \qquad
P_{pump} = \frac{\Delta P\, n_b \dot V}{\eta}
```

[`report_gshp_simulation`](@ref) then gives the electricity of the heat pumps
`W_hp = |Q_{b,hp}| / COP`, the seasonal performance factors `SPF_hp = Σ|Q_{b,hp}| / ΣW_hp` and
`SPF_sys = Σ|Q_{b,hp}| / Σ(W_hp + P_pump)`, and the annual ground load imbalance `ΣQ_g / n_years`
(positive when more heat is rejected than extracted).
