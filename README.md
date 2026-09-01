# GroundSourceHeatPumpDesign.jl

[![CI](https://github.com/GHE-jl/GroundSourceHeatPumpDesign.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GHE-jl/GroundSourceHeatPumpDesign.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://github.com/GHE-jl/GroundSourceHeatPumpDesign.jl/actions/workflows/Documenter.yml/badge.svg?branch=main)](https://GHE-jl.github.io/GroundSourceHeatPumpDesign.jl)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

Heat pump modelling for ground-source heat pump (GSHP) systems, part of the
[GHE.jl](https://GHE-jl.github.io/GHE.jl-docs) ecosystem. It sits at the same level as
[GroundHeatExchangerSizing.jl](https://github.com/GHE-jl/GroundHeatExchangerSizing.jl) and
[ThermalResponseTest.jl](https://github.com/GHE-jl/ThermalResponseTest.jl), and like them depends
on [`GroundHeatExchanger.jl`](https://github.com/GHE-jl/GroundHeatExchanger.jl) for the ground-side
physics (g-function, effective borehole resistance, temporal superposition).

## What it does

- **Heat pump performance interpolation**: Capacity (CAP) and coefficient of performance (COP)
  from manufacturer data, via 1D or 2D (source/load entering fluid temperature) interpolation.
- **Building-to-ground load conversion**: Converts a building heating/cooling load into the
  corresponding ground load using the heat pump's COP, capped at the installed heat pumps'
  combined capacity. See `script/script_heat_pump.jl` for a worked example going from a fixed
  nominal COP to a temperature-dependent one, with the heat pump curves plotted.
- **`gshp_simulation`**: A minimal, non-iterative ground source heat pump (GSHP) simulation combining the load-conversion
  tools above with `GroundHeatExchanger.jl`'s g-function, effective borehole resistance and
  temporal superposition. There is no load shaving or optimization yet — any load beyond the
  installed heat pump capacity is simply reported as unmet, to be covered by an auxiliary system.
  Implemented as multiple methods, in increasing order of complexity — see the "Simulation cases"
  section below.

See the [documentation](https://GHE-jl.github.io/GroundSourceHeatPumpDesign.jl) for a tutorial and
the full API reference, and the `script/` folder for runnable, self-contained examples:
`script_heat_pump.jl` (heat pump performance and load conversion only), `script_gshp_1_fixedCOP.jl`
and `script_gshp_2_variableCOP.jl` (full simulation, numbered by increasing complexity — see below).

## Simulation cases

`gshp_simulation` is one function name with multiple methods, one per case, dispatched on the
arguments (a fixed COP/CAP is passed as plain numbers; a temperature-dependent one as
`heat_pump_performance`-style curves). Implemented so far:

1. **Fixed COP** (`script_gshp_1_fixedCOP.jl`) — temperature-independent COP/CAP, stationary
   (constant flow rate) ground response, no source-side load shaving.
2. **Variable COP** (`script_gshp_2_variableCOP.jl`) — temperature-dependent COP/CAP
   (`heat_pump_performance` curves), stationary ground response, no source-side load shaving. The
   building-to-ground conversion still uses a fixed nominal COP for its single pass (the source EFT
   isn't known until the ground response is simulated), but heat pump staging/COP are then reported
   at the resulting, actual source EFT.

Planned, not yet implemented:

3. **Source-side load shaving** (either COP style) — cap `Qg` further so the source EFT stays
   within the heat pumps' operating limits, along the lines of the correction scheme in Pasquier et
   al. (2013); needs a couple of extra `fluid_temperature` passes but still no full outer iteration.
4. **Non-stationary (time-varying) flow rate** (either COP style, with or without shaving) — swap
   the stationary `ground_response`/`fluid_temperature` call for `GroundHeatExchanger.jl`'s
   `convolution_ns`/`impulse_func_ns`/state-vector machinery, with a flow-rate control law (e.g. a
   `V(EWT)` map, as used in real hybrid/standing-column-well projects).

## Conventions

- **Sides**: The building is the **load** side, and the ground is the **source** side.
- **Sign convention**: A building load `Qb` is positive when the building is heating (winter) and
  negative when it is cooling (summer). The corresponding ground load `Qg` is negative when the
  ground is cooling (heat extracted from it, building heating) and positive when the ground is
  heating (heat rejected to it, building cooling). See `Q_COP`.
- **Nomenclature**: The different fluid circuit of a ground source heat pump systems are:
  - Load side (Building -> condensor of heat pump): `EWT_b = LWT_l` and `LWT_b = EWT_l` (`b` for the building, `l` for load at the condensor). The
    heat pump's leaving load-side temperature is the the building entry temperature, and vice versa.
  - Source side (Evaporator of heat pump -> ground): `LWT_s = EWT_g = T_in` and `EWT_s = LWT_g = T_out`
    (`s` for source at the evaporator, `g` for ground). The fluid leaving the heat pump's source side is the entering GHE (`T_in`), and vice-versa. This package considers a 100% efficiency if a secondary heat exchanger is present in the ground loop.
  - The heat pump's operating temperature limits are checked against `EWT_s` (`T_out`), matching the `Tsrc` grid that indexes `heat_pump_performance`.
- **Auxiliary vs. mitigation**: "Auxiliary" is unmet **building**-side demand covered by a separate
  system (boiler, electric resistance, cooling tower, air-source heat pump, etc.). The function `ground_load_capped`' gives the `Qb_unmet` vector. "Mitigation" is a **source**-side correction applied to keep
  the source fluid within the heat pump's operating envelope. It gives the `Qs_unmet` vector that must be covered by a separate system on the source side (typically a boiler).

## Installation

This package is not yet registered in the Julia General registry. Install directly from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/GHE-jl/GroundSourceHeatPumpDesign.jl")
```

## License

GPL-3.0, see [LICENSE](LICENSE).
