module GroundSourceHeatPumpDesign

include("heat_pump.jl")
include("thermal_loads.jl")

# Functions from heat_pump.jl
export heat_pump_performance

# Functions from thermal_loads.jl
export Q_COP, ground_load_from_heat_pump

end
