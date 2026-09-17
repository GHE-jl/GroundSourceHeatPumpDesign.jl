# thermal_loads: building to ground load conversion, showcased in increasing order of
# complexity (see ground_load's docstring for the three method tiers this follows):
#   1. Fixed nominal COP, no capacity limit
#   2. Fixed nominal COP, capacity-capped
#   3. Performance-curve-based (temperature-dependent COP/CAP, capacity-capped)
# plus n_active_heat_pumps() as a post diagnostic on top of method 3.
# The source EFT is just an illustrative sweep, and the building load is
# a synthetic year of hourly demand (GroundHeatExchanger.jl's ground_load_profile, Bernier 2004).

using Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate() # Once per machine, to install dependencies in the local environment

using GroundSourceHeatPumpDesign
using GroundHeatExchanger: ground_load_profile
using CairoMakie

# Manufacturer performance data (1D: source EFT only) for a single heat pump unit
Tsrc_grid_heat = [-5.0, 0.0, 5.0, 10.0, 15.0] # [°C] source EFT grid for heating
CAP_heat_vec = [8.0, 9.0, 10.0, 11.0, 12.0] .* 1000 # [W] Capacity corresponding to Tsrc_grid_heat
COP_heat_vec = [3.0, 3.4, 3.8, 4.2, 4.6] # [-] COP corresponding to Tsrc_grid_heat
Tsrc_grid_cool = [10.0, 20.0, 30.0, 40.0] # [°C] source EFT grid for cooling
CAP_cool_vec = [13.0, 12.0, 11.0, 10.0] .* 1000 # [W] Capacity corresponding to Tsrc_grid_cool
COP_cool_vec = [5.8, 5.0, 4.3, 3.7] # [-] COP corresponding to Tsrc_grid_cool

perf_heating = (Tsrc_grid_heat, CAP_heat_vec, COP_heat_vec) # tuple for heat_pump_performance
perf_cooling = (Tsrc_grid_cool, CAP_cool_vec, COP_cool_vec) # tuple for heat_pump_performance

n_HP = 5 # [-] number of heat pump units installed

# Synthetic building load: one hourly year (Bernier 2004 profile, as in script_gshp_1_fixedCOP.jl)
t_h = collect(0.0:1.0:(8760.0 - 1)) # [h]
Qb = -ground_load_profile(t_h, 80_000.0) # [W] >0 heating, <0 cooling

# Illustrative source EFT sweep (cold source when heating, warm source when cooling).
Tsrc = 12.0 .- 10.0 .* cos.(2π .* t_h ./ 8760) # [°C]

# Plot of inputs: the raw building load and source EFT before any conversion, i.e. the annual
# heating/cooling demand shape (top) and the sinusoidal source temperature sweep (cold in winter,
# warm in summer) that drives every variable-COP/CAP result further down (bottom).
fig0 = Figure()
ax1 = Axis(fig0[1, 1], xlabel="Time [h]", ylabel="Building load [kW]",
    title="Synthetic building load (Bernier 2004 profile)")
lines!(ax1, t_h, Qb ./ 1000)
ax2 = Axis(fig0[2, 1], xlabel="Time [h]", ylabel="Temperature [°C]",
    title="Illustrative source EFT")
lines!(ax2, t_h, Tsrc)
display(fig0)

# Method 1: fixed COP, no capacity limit

# Fixed nominal COP: a single design COP, independent of Tsrc
COP_heating_nominal = sum(COP_heat_vec) / length(COP_heat_vec) # [-] mean COP over the grid = 3.8
COP_cooling_nominal = sum(COP_cool_vec) / length(COP_cool_vec) # [-] mean COP over the grid = 4.7

# Conversion with fixed COP
Qg_fixed = ground_load(Qb, COP_heating_nominal, COP_cooling_nominal) # [W]

# Variable COP: COP read off the curve at the source EFT at each time step with
# heat_pump_performance, then ground_load's vector-COP method (still method 1) applies it point
# by point.
_, COP_heating_series = heat_pump_performance(Tsrc, Tsrc_grid_heat, CAP_heat_vec, COP_heat_vec)
_, COP_cooling_series = heat_pump_performance(Tsrc, Tsrc_grid_cool, CAP_cool_vec, COP_cool_vec)

# Conversion with variable COP
Qg_variable = ground_load(Qb, COP_heating_series, COP_cooling_series) # [W]

# Fixed vs. variable COP, both uncapped (method 1, so no unmet load in either case): shows how
# letting COP track Tsrc shifts the ground load's magnitude relative to a single nominal COP.
# Larger heat extraction where Tsrc is coldest (worse heating COP) and larger heat rejection where
# Tsrc is warmest (worse cooling COP), tracking the Tsrc sweep from fig0.
fig1 = Figure()
ax1 = Axis(fig1[1, 1], xlabel="Time [h]", ylabel="Ground load [kW]",
    title="Method 1: fixed vs. variable COP (no capacity limit)")
lines!(ax1, t_h, Qg_fixed ./ 1000, label="Fixed nominal COP")
lines!(ax1, t_h, Qg_variable ./ 1000, label="Variable COP (illustrative Tsrc)")
axislegend(ax1, position=:cb)
ax2 = Axis(fig1[2, 1], xlabel="Time [h]", ylabel="Ground load [kW]",
    title="Method 1: fixed vs. variable COP (no capacity limit)")
lines!(ax2, t_h, (Qg_fixed .- Qg_variable) ./ 1000, label="Difference (fixed - variable)")
axislegend(ax2, position=:cb)
display(fig1)

# Method 2: fixed nominal COP, capacity-capped

Qg_capped, Qb_unmet = ground_load(Qb, COP_heating_nominal, COP_cooling_nominal,
    maximum(CAP_heat_vec), maximum(CAP_cool_vec); n_HP)
Qb_hp = Qb .- Qb_unmet

println("Method 2 (capacity-capped, fixed COP):")
println("  Unmet hours: $(count(!=(0.0), Qb_unmet)) / $(length(Qb))")
println("  Unmet loads: $(sum(Qb_unmet) / 1000) kW")

# Effect of the installed capacity limit on top of a fixed nominal COP (method 2): top panel
# compares raw building demand against what the installed heat pumps can actually cover; middle
# panel is the resulting unmet load, non-zero only where demand exceeds the installed
# n_HP * CAP_heating_unit/CAP_cooling_unit (the peaks of the winter/summer demand); bottom panel is
# the corresponding (capped) ground load.
fig2 = Figure()
ax1 = Axis(fig2[1, 1], ylabel="Building load [kW]",
    title="Method 2: capacity-capped (fixed nominal COP)")
lines!(ax1, t_h, Qb ./ 1000, color=(:gray, 0.6), label="Building demand")
lines!(ax1, t_h, Qb_hp ./ 1000, color=:steelblue, label="Covered by heat pump(s)")
axislegend(ax1, position=:rt)
ax2 = Axis(fig2[2, 1], ylabel="Unmet building load [kW]")
lines!(ax2, t_h, Qb_unmet ./ 1000, color=:firebrick)
ax3 = Axis(fig2[3, 1], xlabel="Time [h]", ylabel = "Ground load [kW]")
lines!(ax3, t_h, Qg_capped ./ 1000, color=:goldenrod)
display(fig2)

# Method 3: performance-curve-based (variable COP + capacity together)

Qg_curve, Qb_unmet_curve = ground_load(Qb, Tsrc, perf_heating, perf_cooling; n_HP)

println("Method 3 (performance-curve-based):")
println("  Unmet hours: $(count(!=(0.0), Qb_unmet_curve)) / $(length(Qb))")
println("  Unmet loads: $(sum(Qb_unmet_curve) / 1000) kW")

# Isolates the effect of variable, temperature-dependent COP/CAP (method 3) against a fixed
# nominal COP/CAP (method 2), both capacity-capped.
fig3 = Figure()
ax1 = Axis(fig3[1, 1], xlabel="Time [h]", ylabel="Ground load [kW]",
    title="Capacity-capped: fixed nominal COP (method 2) vs. performance curve (method 3)")
lines!(ax1, t_h, Qg_capped ./ 1000, label="Method 2 (fixed nominal COP)")
lines!(ax1, t_h, Qg_curve ./ 1000, label="Method 3 (performance curve)")
axislegend(ax1, position=:cb)
ax2 = Axis(fig3[2, 1], xlabel="Time [h]", ylabel="Ground load [kW]")
lines!(ax2, t_h, (Qg_capped .- Qg_curve) ./ 1000, label="Difference (method 2 - method 3)")
axislegend(ax2, position=:cb)
display(fig3)

# Note: the difference between Qg_capped (method 2) and Qg_curve (method 3) has two causes, both
# from replacing a single averaged nominal COP/CAP with the real, temperature-dependent curve:
#   - COP varies with Tsrc: a below-average source EFT gives a below-average (worse) COP, so the
#     same building load pulls more heat from the ground than the fixed-COP estimate, and vice
#     versa above average.
#   - CAP varies with Tsrc: method 2's cap (n_HP * maximum(CAP_heat_vec)) is only reached at the
#     single best grid point. The real, temperature-dependent capacity is lower everywhere else, so
#     method 3 caps (and reports unmet load) earlier (compare the "Unmet hours"/"Unmet loads"
#     printed above for methods 2 and 3).

# n_active_heat_pumps: staging/COP diagnostic on top of method 3

Qb_hp_curve = Qb .- Qb_unmet_curve # [W] building load actually delivered by the heat pumps
n_active, COP_actual = n_active_heat_pumps(Qb_hp_curve, Tsrc, n_HP, perf_heating, perf_cooling)

println("Mean active heat pumps: $(round(sum(n_active) / length(n_active); digits=2)) / $n_HP")

# Heat pump staging diagnostic on top of method 3: top panel is the fractional number of units
# (0 to n_HP) needed to deliver Qb_hp_curve at each time step -- it tracks the building demand
# shape, saturating at n_HP whenever load is capped; bottom panel is the actual COP achieved at
# that time step's Tsrc, from the same interpolation.
fig4 = Figure(size=(900, 450))
ax1 = Axis(fig4[1, 1], ylabel="Active heat pumps [-]",
    title="Heat pump staging (illustrative source EFT)")
lines!(ax1, t_h, n_active)
ax2 = Axis(fig4[2, 1], xlabel="Time [h]", ylabel="COP [-]")
lines!(ax2, t_h, COP_actual)
display(fig4)
