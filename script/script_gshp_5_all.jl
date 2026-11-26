# GSHP simulation, case 5: all the considerations together, compared with case 2.
#   - "case 2": COP/CAP curves, constant flow rate, no operating limits
#   - "case 5": COP/CAP curves with a part-load degradation (Cd), non-stationary flow rate picked
#     from the ground load among three set points, largest flow rate within `T_margin` of the
#     limits, source-side mitigation to respect them, and the effective borehole resistance
#     recalculated at the mean fluid temperature (`fluid` keyword)
# The limits are checked at the ground inlet T_in (`Tlim_at = :inlet` by default, the heat-pump LFT;
# seen from the ground). The outlet model can also be the outlet g-function of each flow
# rate (`outlet_model = :gfunction`, needs the grout heat capacity `Cg`).

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

# Heat pump performance (single unit, source EFT only), count and part-load degradation
Tsrc_grid_heat = [-5.0, 0.0, 5.0, 10.0, 15.0]
CAP_heat_vec = [8.0, 9.0, 10.0, 11.0, 12.0] .* 1000 # W
COP_heat_vec = [3.0, 3.4, 3.8, 4.2, 4.6]
Tsrc_grid_cool = [10.0, 20.0, 30.0, 40.0]
CAP_cool_vec = [13.0, 12.0, 11.0, 10.0] .* 1000 # W
COP_cool_vec = [5.8, 5.0, 4.3, 3.7]
perf_heating = (Tsrc_grid_heat, CAP_heat_vec, COP_heat_vec)
perf_cooling = (Tsrc_grid_cool, CAP_cool_vec, COP_cool_vec)
n_HP = 8
Cd = 0.15 # part-load degradation coefficient [-]

# Borehole field and ground properties
xy = borefield_rectangle(3, 3, 6.0) # 3x3 field, 6 m spacing [m]
H = 120.0 # borehole depth [m]
D = 2.0 # buried depth [m]
rb = 0.08 # borehole radius [m]
ks, Cs = 3.0, 2.11e6 # ground thermal conductivity [W/mK], volumetric heat capacity [J/m3K]
T0 = 10.0 # undisturbed ground temperature [degC]

# Pipe, grout and fluid properties (single U-loop, water)
s = 0.05 # shank spacing [m]
ro, ri = 0.022, 0.017 # pipe outer/inner radius [m]
kg, kp = 1.6, 0.4 # grout, pipe thermal conductivity [W/mK]
kf, cf, ρf, μf = fluid_property(T0, :water)
fluid = (:water, 100) # fluid of `fluid_property`, for the Rb* recalculation of case 5

# Flow rates: constant for case 2, three set points for case 5
V = 30.0 / 60_000 # 30 L/min per borehole [m3/s]
V_candidates = [25.0, 30.0, 35.0] ./ 60_000 # 25/30/35 L/min per borehole [m3/s]
ΔT_target = 2.5 # K, field inlet/outlet split targeted when picking the flow rate

# Operating temperature limits, and margin below which the largest flow rate is forced
Tlim = [1.0, 30.0] # degC
T_margin = 1.0 # K

# Case 2 and case 5
res2 = gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP, xy, rb, H, D, ks, Cs, s, ro, ri,
    kg, kp, kf, cf, ρf, μf, V, T0)
res5 = gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP, xy, rb, H, D, ks, Cs, s, ro, ri,
    kg, kp, kf, cf, ρf, μf, V_candidates, ΔT_target, T0; Cd, Tlim, T_margin, fluid)

# Tables: temperatures, hours outside the limits and energy/power
println("Case 2: COP/CAP curves, constant flow rate, no limits ($(res2.n_iter) iterations)")
report_gshp_simulation(res2, Qb, Tlim, steps_per_year)
println("\nCase 5: Cd = $Cd, staged flow rate, override and mitigation ($(res5.n_iter) iterations)")
report_gshp_simulation(res5, Qb, Tlim, steps_per_year)
println("\nTime share of each flow rate in case 5:")
for V_i in V_candidates
    println("  $(round(V_i * 60_000; digits=1)) L/min: " *
        "$(round(100 * count(==(V_i), res5.V_t) / nt; digits=1)) %")
end

# Fig 1. Ground inlet temperature
fig1 = Figure(size=(900, 450))
ax = Axis(fig1[1, 1], xlabel="Time [years]", ylabel="T_in [degC]",
    title="Ground inlet temperature (heat-pump LFT)")
lines!(ax, years, res2.T_in, color=(:gray, 0.6), label="Case 2")
lines!(ax, years, res5.T_in, color=:steelblue, label="Case 5")
hlines!(ax, Tlim, color=:firebrick, linestyle=:dash, label="Tlim")
axislegend(ax, position=:rb)
display(fig1)

# Fig 2. Flow rate and effective COP over the last year
last_year = (nt - steps_per_year + 1):nt
fig2 = Figure(size=(900, 500))
ax1 = Axis(fig2[1, 1], ylabel="Flow rate [L/min]", title="Case 5, last simulated year")
lines!(ax1, years[last_year], res5.V_t[last_year] .* 60_000)
ax2 = Axis(fig2[2, 1], xlabel="Time [years]", ylabel="COP [-]")
lines!(ax2, years[last_year], res2.COP_actual[last_year], color=(:gray, 0.6), label="Case 2")
lines!(ax2, years[last_year], res5.COP_actual[last_year], color=:steelblue,
    label="Case 5 (with Cd)")
axislegend(ax2, position=:rb)
display(fig2)
