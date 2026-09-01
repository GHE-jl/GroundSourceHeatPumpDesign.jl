```@meta
CurrentModule = GroundSourceHeatPumpDesign
```

# Tutorial

This walks through the two pieces the package provides today: reading a heat pump's performance
off manufacturer data, and converting a building load into the corresponding ground load.

## 1. Heat pump performance from manufacturer data

Manufacturer catalogs give capacity (CAP) and coefficient of performance (COP) as a function of
the source (and sometimes load) entering fluid temperature (EFT). [`heat_pump_performance`](@ref)
interpolates that table at an arbitrary operating temperature.

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

[`Q_COP`](@ref) converts a building thermal load to the corresponding ground load using the sign
convention: a positive building load (heating) extracts heat from the ground loop (negative ground
load); a negative building load (cooling) rejects heat to the ground loop (positive ground load).

```julia
Qb = 8.0                      # 8 kW building heating load
Qg = Q_COP(Qb, 4.0, 5.0)      # COP_heating = 4.0, COP_cooling = 5.0 → Qg = -6.0 kW
```

[`ground_load_from_heat_pump`](@ref) combines both steps: given a building load and the current
source/load EFTs, it looks up CAP/COP from the performance tables above, caps the heat pump's
contribution at its capacity, and reports any load left unmet:

```julia
perf_heating = (Tsrc_grid, Tload_grid, CAP_table, COP_table)
perf_cooling = (Tsrc_grid, Tload_grid, CAP_table, COP_table)

Q_ground, Q_unmet = ground_load_from_heat_pump(8.0, 7.5, 35.0, perf_heating, perf_cooling)
```

## 3. A minimal, non-iterative GSHP simulation

The source EFT isn't known until the ground temperature has been simulated, and simulating the
ground temperature needs the ground load, which itself depends on the source EFT through the
heat pump's COP — a circular dependency that the full, iterative approach (e.g. Pasquier et al.
2013) resolves by iterating. This package instead offers the simple, single-pass alternative used
for a first design pass (in the spirit of GHE Analysis): convert the building load with a *fixed*
nominal COP, cap it at the installed capacity, simulate the ground temperature once (no
iteration), and only then read off the heat pump performance actually seen at the resulting source
EFT — as a diagnostic, not fed back into the simulation.

[`gshp_simulation`](@ref) runs the whole pipeline — building load to ground load (capped at the
installed capacity), then the ground response (g-function, effective borehole resistance, and
stationary convolution, all from
[`GroundHeatExchanger.jl`](https://github.com/GHE-jl/GroundHeatExchanger.jl)) — and reports the
ground temperature entering (`T_in`) and leaving (`T_out`) the field, see the README's
"Conventions" section for the nomenclature. It is one function name with multiple methods, one per
case of increasing complexity (see the README's "Simulation cases" section for what's implemented
and what's planned); the two implemented so far differ only in how the heat pump's COP/CAP are
described.

Case 1, fixed (temperature-independent) COP/CAP — see `script/script_gshp_1_fixedCOP.jl`:

```julia
n_HP = 5
res = gshp_simulation(t, Qb, 3.8, 4.7, 12_000.0, 13_000.0, n_HP,
    xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0)
```

Case 2, temperature-dependent COP/CAP curves — see `script/script_gshp_2_variableCOP.jl`. Here the
one-pass conversion still uses a fixed nominal COP (the mean of the curve), but
`res.n_active`/`res.COP_actual` are reported at the resulting, actual `T_out`:

```julia
perf_heating = (Tsrc_grid, CAP_vec, COP_vec)
perf_cooling = (Tsrc_grid, CAP_vec, COP_vec)
res = gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP,
    xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0)
```

`res` is a named tuple (`Qg`, `Qb_hp`, `Qb_unmet`, `T_in`, `T_out`, `n_active`, `Rbe`, `g`, plus
`COP_actual` for case 2). [`unmet_hours`](@ref), [`annual_drift`](@ref) and [`thermal_energy`](@ref)
then summarize the operating-temperature-limit exceedances, long-term ground temperature drift, and
energy/peak-power use for a final report, as in the two scripts above.
