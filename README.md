# GroundSourceHeatPumpDesign.jl

[![CI](https://github.com/GHE-jl/GroundSourceHeatPumpDesign.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GHE-jl/GroundSourceHeatPumpDesign.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://github.com/GHE-jl/GroundSourceHeatPumpDesign.jl/actions/workflows/Documenter.yml/badge.svg?branch=main)](https://GHE-jl.github.io/GroundSourceHeatPumpDesign.jl)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPLv3-blue.svg)](LICENSE)

Heat pump modelling for ground-source heat pump (GSHP) systems, part of the
[GHE.jl](https://GHE-jl.github.io/GHE.jl-docs) ecosystem. It sits at the same level as
[GroundHeatExchangerSizing.jl](https://github.com/GHE-jl/GroundHeatExchangerSizing.jl) and
[ThermalResponseTest.jl](https://github.com/GHE-jl/ThermalResponseTest.jl), and like them depends
on [`GroundHeatExchanger.jl`](https://github.com/GHE-jl/GroundHeatExchanger.jl) for the ground-side
physics (g-function, effective borehole resistance, temporal superposition).

## Quick start

```julia
using GroundSourceHeatPumpDesign

# Manufacturer performance data (1D: source EFT only) for a single heat pump unit
Tsrc_grid_heat = [-5.0, 0.0, 5.0, 10.0, 15.0]        # source EFT grid [°C]
CAP_heat_vec   = [8.0, 9.0, 10.0, 11.0, 12.0] .* 1000 # heating capacity [W]
COP_heat_vec   = [3.0, 3.4, 3.8, 4.2, 4.6]
Tsrc_grid_cool = [10.0, 20.0, 30.0, 40.0]
CAP_cool_vec   = [13.0, 12.0, 11.0, 10.0] .* 1000
COP_cool_vec   = [5.8, 5.0, 4.3, 3.7]
perf_heating = (Tsrc_grid_heat, CAP_heat_vec, COP_heat_vec)
perf_cooling = (Tsrc_grid_cool, CAP_cool_vec, COP_cool_vec)

Qb, Tsrc = 8000.0, 7.5   # 8 kW building heating load, 7.5°C source EFT

# Ground load, temperature-dependent COP/CAP, capacity-capped
Qg, Qb_unmet, Qb_hp = ground_load(Qb, Tsrc, perf_heating, perf_cooling)
```

## What it does

- **Heat pump performance interpolation**: Capacity (CAP) and coefficient of performance (COP)
  from manufacturer data, via 1D or 2D (source/load entering fluid temperature) interpolation.
- **Building-to-ground load conversion**: Converts a building heating/cooling load into the
  corresponding ground load using the heat pump's COP. One function name (`ground_load`) with
  methods of increasing complexity — fixed COP, fixed COP with capacity capping, and
  performance-curve-based (temperature-dependent CAP/COP, capacity-capped) — plus
  `n_active_heat_pumps` for the resulting heat pump staging/COP as a diagnostic. See
  `script/script_thermal_loads.jl` for a worked example of all of the above, and
  `script/script_COP.jl` for the underlying temperature -> CAP/COP lookup with the heat pump
  curves plotted.
- **`gshp_simulation`**: A minimal ground source heat pump (GSHP) simulation combining the load-conversion
  tools above with `GroundHeatExchanger.jl`'s g-function, effective borehole resistance and
  temporal superposition. Each case resolves the load/temperature circularity as locally as it can
  rather than iterating the whole pipeline — cases 1 and 3 are a single pass, cases 2 and 5 are a
  small, quickly-converging fixed-point loop over the source EFT-dependent COP (see "Simulation
  cases" below). Any load beyond the installed heat pump capacity is reported as unmet
  (`Qb_unmet`), to be covered by an auxiliary system. With operating limits `Tlim`, the ground load
  is also shaved so that the fluid stays within them, and the shaved part is reported as
  `Qs_unmet`.

See the [documentation](https://GHE-jl.github.io/GroundSourceHeatPumpDesign.jl) for a tutorial and
the full API reference.

## Scripts

Run from the package root with `julia --project=script/ script/<name>.jl`. First-time setup:
```
julia --project=script/ -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

| Script | What it shows |
|---|---|
| `script_COP.jl` | Heat pump CAP/COP interpolation from manufacturer data (`heat_pump_performance`) — 1D (source EFT) and 2D (source + load EFT) overloads, with the curves/surface plotted |
| `script_thermal_loads.jl` | Building-to-ground load conversion (`ground_load`), from fixed COP to temperature-dependent, capacity-capped, plus heat pump staging (`n_active_heat_pumps`) — increasing complexity, on a synthetic building load |
| `script_gshp_1_fixedCOP.jl` | Full GSHP simulation, case 1: fixed COP/CAP, stationary ground response |
| `script_gshp_2_variableCOP.jl` | Full GSHP simulation, case 2: temperature-dependent COP/CAP curves, stationary ground response |
| `script_gshp_3_nonstationaryFlow.jl` | Full GSHP simulation, case 3: fixed COP/CAP, non-stationary (load-based, discrete) flow rate, compared against the stationary case 1 |
| `script_gshp_4_mitigation.jl` | Full GSHP simulation, case 4: source-side mitigation to respect the operating limits, compared against case 1 |
| `script_gshp_5_all.jl` | Full GSHP simulation, case 5: COP/CAP curves with part-load degradation, non-stationary flow rate, flow override near the limits and mitigation, compared against case 2 |

## Simulation cases

`gshp_simulation` is one function name with multiple methods, one per case, dispatched on the
arguments (a fixed COP/CAP is passed as plain numbers; a temperature-dependent one as
`heat_pump_performance`-style curves). Within a case, there are two methods in increasing order of
complexity: one takes an already-computed g-function and effective borehole resistance (the base
method — reuse it directly when calling `gshp_simulation` repeatedly for the same field, since
building those is the expensive part); the other takes the field's raw geometry/soil/pipe inputs,
computes them once, and delegates to the first. Implemented so far:

1. **Fixed COP** (`script_gshp_1_fixedCOP.jl`) — temperature-independent COP/CAP, stationary
   (constant flow rate) ground response, no source-side load shaving.
2. **Variable COP** (`script_gshp_2_variableCOP.jl`) — temperature-dependent COP/CAP
   (`heat_pump_performance` curves), stationary ground response, no source-side load shaving. The
   building-to-ground conversion starts from a fixed nominal COP (the source EFT isn't known until
   the ground response is simulated), then a small fixed-point loop alternates the curve-based
   conversion (at the previous iteration's `T_out`) with the ground response until `T_out` stops
   moving — typically well under the 5-iteration safety cap, mirroring `fNSCiter.m`'s
   identically-capped loop. Heat pump staging/COP are then reported at the converged, actual
   source EFT.
3. **Non-stationary flow, fixed COP** (`script_gshp_3_nonstationaryFlow.jl`) — same building load
   conversion as case 1, but the ground response switches between a small set of controllable
   per-borehole flow rates instead of a single constant `V`, via `GroundHeatExchanger.jl`'s
   `convolution_ns` state-vector machinery. The active flow rate at each time step is picked from
   the ground load alone (never from the simulated temperature — see `_flow_rate_states`), so this
   stays a single, non-iterative pass. A bigger flow rate always shrinks the `T_in`/`T_out` split
   for the same load; whether that also lowers a given `T_out` extreme depends on the load's sign
   there too (see the script for a worked comparison against the stationary case). Case 3b combines
   it with case 2's COP/CAP curves (same fixed-point loop).
4. **Source-side mitigation** (`script_gshp_4_mitigation.jl`) — any of the cases above with the
   `Tlim = [Tlim_low, Tlim_high]` keyword. At the time steps where the fluid would leave the limits,
   the ground load is shaved toward zero with the correction scheme of Pasquier et al. (2013),
   `ΔQg = (Tlim - T) / a`, with `a` the temperature response to a unit load over one time step. It
   is repeated until the temperature sits on the limit (no full outer iteration). The shaved part
   `Qs_unmet = Qg_hp - Qg` must be covered by a source-side auxiliary system. The limits are checked
   at `Tlim_at = :inlet` (default, `T_in`), `:outlet` (`T_out`) or `:mean`.
5. **All together** (`script_gshp_5_all.jl`) — case 3b with `Tlim`. Before the mitigation, the
   largest flow rate is forced at the time steps where the fluid is within `T_margin` of the limits
   (as in `fNSCiter.m`), so that the flow is used before the load is shaved.

Options common to every case:

- `Cd` — part-load degradation coefficient (default `0`, none). The fractional unit of the staging
  count cycles on and off, which lowers its COP: `PLF = 1 - Cd (1 - PLR)`, see `part_load_COP`.
  With `Cd > 0`, the degraded COP is used in the building-to-ground conversion.
- `outlet_model` (non-stationary flow) — `:resistance` (default, outlet response from the effective
  borehole resistance of each flow rate) or `:gfunction` (outlet g-function of each
  flow rate, Beaudry et al. 2024; needs the grout heat capacity `Cg`).
- `fluid` (raw-input methods) — e.g. `(:MPG, 20)`. The fluid properties and the effective borehole
  resistance `Rb*` are recalculated at the mean fluid temperature of the previous pass, and the
  simulation is repeated until that temperature stops moving (typically 3 passes). Near 0 °C with
  a glycol mix, the viscosity increases markedly and the flow can become laminar, raising `Rb*`. The
  fluid heat capacity of the inlet/outlet split stays at the given value. With
  `outlet_model = :resistance` only.
- `η_pump` (raw-input methods) — pump efficiency of the pumping power `P_pump` (default `0.5`), see
  `pumping_power`: Darcy-Weisbach pressure drop over the U-loop (`2H`), plus one 180° bend and two
  90° elbows per U-loop, for all the boreholes in parallel. It is zero when the heat pumps are off.

`report_gshp_simulation(res, Qb, Tlim, steps_per_year)` reports every case from the output of
`gshp_simulation`: temperatures and hours outside the limits, energy and peak power of the
building, heat pumps, auxiliary systems and ground, the electricity of the heat pumps
(`|Qb_hp| / COP_actual`) and of the pump, the annual ground load imbalance, and the seasonal
performance factors of the heat pumps and of the heat pumps with the pumping.

## Conventions

- **Sides**: The building is the **load** side, and the ground is the **source** side.
- **Sign convention**: A building load `Qb` is positive when the building is heating (winter) and
  negative when it is cooling (summer). The corresponding ground load `Qg` is negative when the
  ground is cooling (heat extracted from it, building heating) and positive when the ground is
  heating (heat rejected to it, building cooling). See `ground_load`.
- **Nomenclature**: The different fluid circuit of a ground source heat pump systems are:
  - Load side (Building -> condensor of heat pump): `EWT_b = LWT_l` and `LWT_b = EWT_l` (`b` for the building, `l` for load at the condensor). The
    heat pump's leaving load-side temperature is the the building entry temperature, and vice versa.
  - Source side (Evaporator of heat pump -> ground): `LWT_s = EWT_g = T_in` and `EWT_s = LWT_g = T_out`
    (`s` for source at the evaporator, `g` for ground). The fluid leaving the heat pump's source side is the entering GHE (`T_in`), and vice-versa. This package considers a 100% efficiency if a secondary heat exchanger is present in the ground loop.
  - The COP/CAP curves are indexed by `EWT_s` (`T_out`). The operating temperature limits are
    checked, seen from the ground, at `Tlim_at = :inlet` (default, `T_in = LWT_s`), `:outlet`
    (`T_out = EWT_s`) or `:mean`.
- **Auxiliary vs. mitigation**: "Auxiliary" is unmet **building**-side demand covered by a separate
  system (boiler, electric resistance, cooling tower, air-source heat pump, etc.). The function `ground_load` gives the `Qb_unmet` vector. "Mitigation" is a **source**-side correction applied to keep
  the source fluid within the heat pump's operating envelope. It gives the `Qs_unmet` vector that must be covered by a separate system on the source side (typically a boiler).

## Installation

This package is not yet registered in the Julia General registry. Install directly from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/GHE-jl/GroundSourceHeatPumpDesign.jl")
```

## Dependencies

### Library

| Package | Used for |
|---|---|
| [GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl) | Ground-side physics used by `gshp_simulation`: g-function (`gfunction`), effective borehole resistance, and stationary fluid/inlet/outlet temperature |

### Scripts only

| Package | Used in |
|---|---|
| [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl) | Plotting in all `script_*.jl` examples |

## References

- Pasquier, P., Marcotte, D., Bernier, M., & Kummert, M. (2013). Simulation of ground-coupled heat
  pump systems using a spectral approach. *Proceedings of BS2013*, 2691–2698.
  Source of the fractional, time-averaged heat pump staging count used by `n_active_heat_pumps`,
  and the correction scheme used for source-side mitigation (case 4 above).
- Beaudry, G., Pasquier, P., & Nguyen, A. (2024). New formulations and experimental validation of
  non-stationary convolutions for the fast simulation of time-variant flowrates in ground heat
  exchangers. *Science and Technology for the Built Environment*, 30(3), 208–219.
  <https://doi.org/10.1080/23744731.2023.2279468>
  Non-stationary (state-vector) convolution and outlet g-function used by the non-stationary
  flow cases (3, 3b and 5).
- Bernier, M. A., Pinel, P., Labib, R., & Paillot, R. (2004). A Multiple Load Aggregation Algorithm
  for Annual Hourly Simulations of GCHP Systems. *HVAC&R Research*, 10(4), 471–487.
  <https://doi.org/10.1080/10789669.2004.10391115>
  Synthetic annual building load profile used in the example scripts
  (`ground_load_profile` from `GroundHeatExchanger.jl`).

## License

AGPL-3.0, see [LICENSE](LICENSE).
