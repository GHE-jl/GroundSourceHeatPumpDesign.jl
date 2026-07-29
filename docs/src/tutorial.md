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
