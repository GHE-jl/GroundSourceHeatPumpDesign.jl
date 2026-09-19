module GroundSourceHeatPumpDesign

include("heat_pump.jl")
include("thermal_loads.jl")
include("reporting.jl")
include("gshp_simulation.jl")

# Functions from heat_pump.jl
export heat_pump_performance

# Functions from thermal_loads.jl
export ground_load, n_active_heat_pumps

# Functions from reporting.jl
export unmet_hours, annual_drift, thermal_energy, report_gshp_simulation

# Functions from gshp_simulation.jl
export gshp_simulation

end
