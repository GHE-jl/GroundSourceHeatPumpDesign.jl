# References

The methods implemented in this package draw on the following sources.

## Heat pump staging and simulation

- **Pasquier, P., Marcotte, D., Bernier, M., & Kummert, M.** (2013). Simulation of ground-coupled
  heat pump systems using a spectral approach. *Proceedings of BS2013*, 2691–2698.
  Source of the fractional, time-averaged heat pump staging count used by
  [`n_active_heat_pumps`](@ref), and the correction scheme planned for source-side load shaving
  (see the README's "Simulation cases" section, case 3).

## Building load profile

- **Bernier, M. A., Pinel, P., Labib, R., & Paillot, R.** (2004). A Multiple Load Aggregation
  Algorithm for Annual Hourly Simulations of GCHP Systems. *HVAC&R Research*, 10(4), 471–487.
  <https://doi.org/10.1080/10789669.2004.10391115>
  Synthetic annual building load profile used in the example scripts
  (`ground_load_profile` from [`GroundHeatExchanger.jl`](https://github.com/GHE-jl/GroundHeatExchanger.jl)).
