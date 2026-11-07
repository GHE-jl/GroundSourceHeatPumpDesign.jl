```@meta
CurrentModule = GroundSourceHeatPumpDesign
```

# Tutorial

This walks through the building blocks the package provides today: reading a heat pump's
performance off manufacturer data, converting a building load into the corresponding ground load,
and combining both into a minimal GSHP simulation.

## 1. Heat pump performance from manufacturer data

Manufacturer catalogs give capacity (CAP) and coefficient of performance (COP) as a function of
the source (and sometimes load) entering fluid temperature (EFT). [`heat_pump_performance`](@ref)
interpolates that table at an arbitrary operating temperature. See `script/script_COP.jl` for a
runnable version of this section, including the 2D surface and a plotted paired-vector sweep.

1D form — CAP/COP as a function of source EFT only:

```julia
using GroundSourceHeatPumpDesign

Tsrc_grid = [-5.0, 0.0, 5.0, 10.0, 15.0]     # source EFT grid [°C]
CAP_vec   = [8.0, 9.0, 10.0, 11.0, 12.0]      # heating capacity [kW]
COP_vec   = [3.0, 3.4, 3.8, 4.2, 4.6]

CAP, COP = heat_pump_performance(7.5, Tsrc_grid, CAP_vec, COP_vec)
```

2D form — CAP/COP as a function of both source and load EFT (bilinear interpolation):

```julia
Tload_grid = [30.0, 40.0, 50.0]               # load EFT grid [°C]
CAP_table  = [8.0 7.0 6.0; 9.0 8.0 7.0; 10.0 9.0 8.0; 11.0 10.0 9.0; 12.0 11.0 10.0]
COP_table  = [3.0 2.6 2.2; 3.4 3.0 2.6; 3.8 3.4 3.0; 4.2 3.8 3.4; 4.6 4.2 3.8]

CAP, COP = heat_pump_performance(7.5, 35.0, Tsrc_grid, Tload_grid, CAP_table, COP_table)
```

Both forms also accept vectors of source/load temperatures and return vectors of CAP/COP.

## 2. Building load to ground load

[`ground_load`](@ref) converts a building thermal load to the corresponding ground load using the
sign convention: a positive building load (heating) extracts heat from the ground loop (negative
ground load); a negative building load (cooling) rejects heat to the ground loop (positive ground
load). It is one function name with methods of increasing complexity.

The simplest method takes a fixed COP, with no capacity limit:

```julia
Qb = 8.0                          # 8 kW building heating load
Qg = ground_load(Qb, 4.0, 5.0)    # COP_heating = 4.0, COP_cooling = 5.0 → Qg = -6.0 kW
```

Adding `CAP_heating_unit`/`CAP_cooling_unit` clamps `Qb` to the installed capacity first, and
reports the unmet remainder alongside `Qg`. These are the capacity of *one* heat pump unit; the
optional `n_HP` keyword (default `1`) scales them to the total installed capacity:

```julia
Qg, Qb_unmet, Qb_hp = ground_load(Qb, 4.0, 5.0, 10.0, 10.0)          # 1 unit, CAP = 10.0 kW
Qg, Qb_unmet, Qb_hp = ground_load(Qb, 4.0, 5.0, 10.0, 10.0; n_HP=5)  # 5 units, CAP = 50.0 kW total
```

Replacing the fixed COP with `Tsrc` (and, for the 2D case, `Tload`) plus performance-curve tuples
combines both steps: it looks up CAP/COP from the performance tables above at the current EFT(s),
caps the heat pump's contribution at that capacity (again scaled by `n_HP`), and reports any load
left unmet:

```julia
perf_heating = (Tsrc_grid, Tload_grid, CAP_table, COP_table)
perf_cooling = (Tsrc_grid, Tload_grid, CAP_table, COP_table)

Qg, Qb_unmet, Qb_hp = ground_load(8.0, 7.5, 35.0, perf_heating, perf_cooling; n_HP=5)
```

All of the above also accept vectors (`Qb`, and `Tsrc`/`Tload` where applicable) for a time series.
See `script/script_thermal_loads.jl` for a runnable version of this section, using a synthetic
year of building load (from `GroundHeatExchanger.jl`'s `ground_load_profile`) to compare all three
methods.

Once a source EFT time series `Tsrc` is known (e.g. from a ground temperature simulation),
[`n_active_heat_pumps`](@ref) reports how many of the installed units are running at each time
step and the COP actually achieved, given the load the heat pumps actually deliver (`Qb .- Qb_unmet`,
for `Qb`/`Qb_unmet`/`Tsrc` vectors matching in length):

```julia
n_HP = 5
n_active, COP_actual = n_active_heat_pumps(Qb .- Qb_unmet, Tsrc, n_HP, perf_heating, perf_cooling)
```

This is a diagnostic on top of the conversion (used for reporting), not a further step in it — it
does not feed back into `Qg`.

## 3. A minimal GSHP simulation, without a full outer iteration

The source EFT isn't known until the ground temperature has been simulated, and simulating the
ground temperature needs the ground load, which itself depends on the source EFT through the
heat pump's COP — a circular dependency that the full, iterative approach (e.g. Pasquier et al.
2013) resolves by iterating over the whole pipeline. This package instead keeps each case's
resolution as small and local as it can:
- **Case 1** (fixed COP): genuinely single-pass — COP/CAP are fixed numbers, so there's no
  circularity to resolve in the first place.
- **Case 2** (temperature-dependent COP/CAP curves): starts from a fixed nominal COP (the curve's
  mean) to get a first `T_out` estimate, then a small fixed-point loop alternates the curve-based
  COP (at the previous iteration's `T_out`) with the ground response until `T_out` stops moving —
  typically well under `max_iter=5` iterations (see `fNSCiter.m`'s identically-capped loop). This
  is *not* a full outer iteration over building-load-conversion, staging, etc.: only the
  Qg-EFT-COP link is resolved.
- **Case 3** (non-stationary flow): the active flow rate at each time step is picked from the
  ground load alone (never from the simulated temperature), so it stays a genuine single pass.
  Case 3b adds case 2's COP/CAP curves.
- **Case 4** (source-side mitigation): any case with `Tlim`, the ground load is shaved where the
  fluid would leave the limits (a few correction passes, not a full outer iteration).
- **Case 5**: case 3b with `Tlim`, the largest flow rate near the limits, then the mitigation.

See the README's "Simulation cases" section for the details.

[`gshp_simulation`](@ref) runs the whole pipeline — building load to ground load (capped at the
installed capacity), then the ground response (g-function, effective borehole resistance, and
stationary convolution, all from
[`GroundHeatExchanger.jl`](https://github.com/GHE-jl/GroundHeatExchanger.jl)) — and reports the
ground temperature entering (`T_in`) and leaving (`T_out`) the field, see the README's
"Conventions" section for the nomenclature. It is one function name with multiple methods, one per
case of increasing complexity (see the README's "Simulation cases" section for what's implemented
and what's planned); within a case, a lower-complexity method takes an already-computed g-function
and effective borehole resistance, and a higher one takes the field's raw geometry/soil/pipe inputs
and computes them before delegating to the first — reuse the first directly when calling
`gshp_simulation` repeatedly for the same field.

Case 1, fixed (temperature-independent) COP/CAP — see `script/script_gshp_1_fixedCOP.jl`:

```julia
n_HP = 5
res = gshp_simulation(t, Qb, 3.8, 4.7, 12_000.0, 13_000.0, n_HP,
    xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0)
```

Case 2, temperature-dependent COP/CAP curves — see `script/script_gshp_2_variableCOP.jl`. The
converged `res.T_out` is consistent with the curve-based COP actually achieved there
(`res.COP_actual`, also reported by `res.n_active`'s underlying interpolation), and `res.n_iter`
reports how many fixed-point iterations that took:

```julia
perf_heating = (Tsrc_grid, CAP_vec, COP_vec)
perf_cooling = (Tsrc_grid, CAP_vec, COP_vec)
res = gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP,
    xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0)
```

Case 4, operating limits with source-side mitigation — see `script/script_gshp_4_mitigation.jl`.
The limits are checked at the ground inlet `T_in` by default (`Tlim_at = :inlet`, the heat-pump
LFT), and the shaved ground load is returned as `res.Qs_unmet`:

```julia
res = gshp_simulation(t, Qb, 3.8, 4.7, 12_000.0, 13_000.0, n_HP,
    xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0; Tlim = [1.0, 30.0])
```

Case 5, everything together — see `script/script_gshp_5_all.jl`. The flow rate per borehole is
picked among `V_candidates` to target `ΔT_target` between the inlet and outlet, and `Cd` adds the
part-load degradation:

```julia
V_candidates = [25.0, 30.0, 35.0] ./ 60_000 # m3/s per borehole
res = gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP,
    xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V_candidates, 2.5, T0;
    Cd = 0.15, Tlim = [1.0, 30.0], T_margin = 1.0)
```

Adding `fluid = (:water, 100)` (or e.g. `(:MPG, 20)`) recalculates the effective borehole
resistance at the mean fluid temperature, as in `script/script_gshp_5_all.jl`.

`res` is a named tuple (`Qg`, `Qb_hp`, `Qb_unmet`, `Qs_unmet`, `T_in`, `T_out`, `n_active`,
`COP_actual`, `Rbe`, `g`, `P_pump`, `n_iter_fluid`, plus `n_iter`/`converged` for the curve-based
cases and `V_t`/`h_out` for the non-stationary ones). [`report_gshp_simulation`](@ref) then
summarizes the operating-temperature-limit exceedances, long-term ground temperature drift,
energy/peak-power use, electricity, ground load imbalance and seasonal performance factors, as in
the scripts above:

```julia
report_gshp_simulation(res, Qb, [1.0, 30.0], 8760)
```
