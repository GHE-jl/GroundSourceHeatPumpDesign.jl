```@meta
CurrentModule = GroundSourceHeatPumpDesign
```

# GroundSourceHeatPumpDesign.jl

*Heat pump modelling for ground-source heat pump (GSHP) systems, part of the GHE.jl ecosystem.*

This package sits at the same level as
[`GroundHeatExchangerSizing.jl`](https://GHE-jl.github.io/GroundHeatExchangerSizing.jl) and
[`ThermalResponseTest.jl`](https://GHE-jl.github.io/ThermalResponseTest.jl), and like them depends
on [`GroundHeatExchanger.jl`](https://github.com/GHE-jl/GroundHeatExchanger.jl) for the ground-side
physics. Its own heat pump/load layer has no ground-physics dependency and can be used on its own.

!!! note "Early / in-development package"
    This package is at an early stage: it currently provides heat pump performance interpolation,
    building-to-ground load conversion and a GSHP simulation workflow without a full outer
    iteration.
    Expect the public API to change without notice until the package reaches a tagged release.

## What it does today

- **Heat pump performance interpolation** ([`heat_pump_performance`](@ref)) — capacity (CAP) and
  coefficient of performance (COP) from manufacturer performance data, via 1D linear interpolation
  on the source entering fluid temperature (EFT) or 2D bilinear interpolation on both source and
  load EFTs. Values are clamped outside the provided grid.
- **Building-to-ground load conversion** ([`ground_load`](@ref)) — converts a building
  heating/cooling load into the corresponding ground load using the heat pump's COP, respecting the
  sign convention that a building heating load extracts heat from the ground loop and a building
  cooling load rejects heat to it. One function name with methods of increasing complexity: fixed
  COP, fixed COP with capacity capping, and performance-curve-based (temperature-dependent CAP/COP,
  capacity-capped).
- **[`gshp_simulation`](@ref)** — a GSHP simulation combining the above with
  [`GroundHeatExchanger.jl`](https://github.com/GHE-jl/GroundHeatExchanger.jl)'s g-function and
  effective borehole resistance: fixed or curve-based COP/CAP, part-load degradation, constant or
  non-stationary flow rate, source-side mitigation to respect the operating limits, and effective
  borehole resistance at the fluid temperature. Plus [`n_active_heat_pumps`](@ref),
  [`part_load_COP`](@ref), [`pumping_power`](@ref), [`unmet_hours`](@ref), [`annual_drift`](@ref),
  [`thermal_energy`](@ref) and [`report_gshp_simulation`](@ref) (electricity, ground load
  imbalance and seasonal performance factors) for reporting. Implemented as multiple methods of
  increasing complexity — see the `script/script_gshp_*.jl` examples, and the [Conventions](https://github.com/GHE-jl/GroundSourceHeatPumpDesign.jl#conventions)
  section of the README for the sign convention and nomenclature used throughout.

## Where to next

- **[Tutorial](@ref)** — a worked example combining the pieces above.
- **[Modeling theory](@ref)** — the interpolation, load-conversion and staging equations, and the
  non-iterative simulation's rationale.
- **[API reference](@ref)** — the complete docstring reference for every exported function.
- **[References](@ref)** — the bibliography underpinning the implementation.
