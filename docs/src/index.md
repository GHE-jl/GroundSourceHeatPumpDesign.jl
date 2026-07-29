```@meta
CurrentModule = GroundSourceHeatPumpDesign
```

# GroundSourceHeatPumpDesign.jl

*Heat pump modelling for ground-source heat pump (GSHP) systems, part of the GHE.jl ecosystem.*

This package is a standalone sibling of
[`GroundHeatExchangerSizing.jl`](https://GHE-jl.github.io/GroundHeatExchangerSizing.jl) and
[`ThermalResponseTest.jl`](https://GHE-jl.github.io/ThermalResponseTest.jl) — it does not depend on
either — focused on the heat pump side of the system (performance and load conversion) rather than
the ground loop side.

!!! note "Early / in-development package"
    This package is at an early stage: it currently provides heat pump performance interpolation
    and building-to-ground load conversion. Expect the public API to change without notice until
    the package reaches a tagged release.

## What it does today

- **Heat pump performance interpolation** ([`heat_pump_performance`](@ref)) — capacity (CAP) and
  coefficient of performance (COP) from manufacturer performance data, via 1D linear interpolation
  on the source entering fluid temperature (EFT) or 2D bilinear interpolation on both source and
  load EFTs. Values are clamped outside the provided grid.
- **Building-to-ground load conversion** ([`Q_COP`](@ref), [`ground_load_from_heat_pump`](@ref)) —
  converts a building heating/cooling load into the corresponding ground load using the heat
  pump's COP, respecting the sign convention that a building heating load extracts heat from the
  ground loop and a building cooling load rejects heat to it.

## Where to next

- **[Tutorial](@ref)** — a worked example combining both pieces above.
- **[API reference](@ref)** — the complete docstring reference for every exported function.
