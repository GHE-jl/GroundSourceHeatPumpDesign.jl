# GroundSourceHeatPumpDesign.jl

[![CI](https://github.com/GHE-jl/GroundSourceHeatPumpDesign.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GHE-jl/GroundSourceHeatPumpDesign.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://github.com/GHE-jl/GroundSourceHeatPumpDesign.jl/actions/workflows/Documenter.yml/badge.svg?branch=main)](https://GHE-jl.github.io/GroundSourceHeatPumpDesign.jl)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPLv3-blue.svg)](LICENSE)

Heat pump modelling for ground-source heat pump (GSHP) systems, part of the
[GHE.jl](https://GHE-jl.github.io/GHE.jl-docs) ecosystem. It is a standalone sibling of
[GroundHeatExchangerSizing.jl](https://github.com/GHE-jl/GroundHeatExchangerSizing.jl) and
[ThermalResponseTest.jl](https://github.com/GHE-jl/ThermalResponseTest.jl) — it does not depend on
either — focused on the heat pump side of the system (performance and load conversion) rather than
the ground loop side.

> **Early / in development.** This package currently provides heat pump performance interpolation
> and building-to-ground load conversion. The public API may change without notice until the
> package reaches a tagged release.

## What it does today

- **Heat pump performance interpolation** — capacity (CAP) and coefficient of performance (COP)
  from manufacturer data, via 1D or 2D (source/load entering fluid temperature) interpolation.
- **Building-to-ground load conversion** — converts a building heating/cooling load into the
  corresponding ground load using the heat pump's COP.

See the [documentation](https://GHE-jl.github.io/GroundSourceHeatPumpDesign.jl) for a tutorial and
the full API reference.

## Installation

Not yet registered in the Julia General registry. Install directly from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/GHE-jl/GroundSourceHeatPumpDesign.jl")
```

## License

AGPL-3.0, see [LICENSE](LICENSE).
