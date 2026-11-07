# GSHP simulation, case 4: source-side mitigation. Same inputs as script_gshp_1_fixedCOP.jl (fixed
# COP/CAP, constant flow rate), simulated twice:
#   - "case 1": no operating limits, the ground temperatures are free to leave [Tlim_low, Tlim_high]
#   - "case 4": `Tlim` given, the ground load is shaved at the time steps where the fluid would
#     leave the limits, and the shaved part (`Qs_unmet`) is covered by a source-side auxiliary system
# The limits are checked at the ground inlet T_in (`Tlim_at = :inlet` by default, the heat-pump LFT;
# seen from the ground). The same `Tlim` keyword also works with the COP/CAP curves of case 2.

using Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate() # Once per machine, to install dependencies in the local environment

using GroundSourceHeatPumpDesign
using GroundHeatExchanger: borefield_rectangle, fluid_property, ground_load_profile
using CairoMakie

# Simulation horizon
n_years = 10
steps_per_year = 8760
nt = steps_per_year * n_years
t = collect(3600.0:3600.0:(3600.0 * nt)) # s, for the ground response and convolution
years = (t ./ 3600 .- 1) ./ steps_per_year # years, aligned with t

# Building load: one year of hourly loads (synthetic profile, Bernier 2004), repeated over the
# design period. Positive = building heating, negative = building cooling.
Qb_year = -ground_load_profile(collect(0.0:1.0:(steps_per_year - 1)), 80_000.0) # W
Qb = repeat(Qb_year, n_years)

# Fixed nominal heat pump performance, single unit, and count
COP_heating, COP_cooling = 3.8, 4.7
CAP_heating_unit, CAP_cooling_unit = 12_000.0, 13_000.0 # W
n_HP = 5

# Borehole field and ground properties
xy = borefield_rectangle(3, 3, 6.0) # 3x3 field, 6 m spacing [m]
H = 120.0 # borehole depth [m]
D = 2.0 # buried depth [m]
rb = 0.08 # borehole radius [m]
ks, Cs = 3.0, 2.11e6 # ground thermal conductivity [W/mK], volumetric heat capacity [J/m3K]
T0 = 10.0 # undisturbed ground temperature [degC]

# Pipe, grout and fluid properties (single U-loop, 20% glycol mix)
s = 0.05 # shank spacing [m]
ro, ri = 0.022, 0.017 # pipe outer/inner radius [m]
kg, kp = 1.6, 0.4 # grout, pipe thermal conductivity [W/mK]
kf, cf, ρf, μf = fluid_property(T0, :MPG, percentage=20)
V = 5e-4 # flow rate per borehole [m3/s]

# Operating temperature limits
Tlim = [1.0, 30.0] # degC

# Case 1 (no limits) and case 4 (source-side mitigation)
res1 = gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating_unit, CAP_cooling_unit, n_HP,
    xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0)
res4 = gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating_unit, CAP_cooling_unit, n_HP,
    xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0; Tlim)

# Tables: temperatures, hours outside the limits and energy/power
println("Case 1: no operating limits")
report_gshp_simulation(res1, Qb, Tlim, steps_per_year)
println("\nCase 4: source-side mitigation")
report_gshp_simulation(res4, Qb, Tlim, steps_per_year)

# Fig 1. Ground inlet temperature, with and without mitigation
fig1 = Figure(size=(900, 450))
ax = Axis(fig1[1, 1], xlabel="Time [years]", ylabel="T_in [degC]",
    title="Ground inlet temperature (heat-pump LFT)")
lines!(ax, years, res1.T_in, color=(:gray, 0.6), label="Case 1, no limits")
lines!(ax, years, res4.T_in, color=:steelblue, label="Case 4, mitigation")
hlines!(ax, Tlim, color=:firebrick, linestyle=:dash, label="Tlim")
axislegend(ax, position=:rb)
display(fig1)

# Fig 2. Ground load and its part shaved by the mitigation
fig2 = Figure(size=(900, 450))
ax = Axis(fig2[1, 1], xlabel="Time [years]", ylabel="Load [kW]", title="Ground load")
lines!(ax, years, res1.Qg ./ 1000, color=(:gray, 0.6), label="Heat pump ground load")
lines!(ax, years, res4.Qs_unmet ./ 1000, color=:firebrick,
    label="Shaved, covered on the source side")
axislegend(ax, position=:rb)
display(fig2)
