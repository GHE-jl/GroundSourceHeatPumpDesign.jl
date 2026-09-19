# GSHP simulation, case 1 (simplest): fixed (temperature-independent) COP/CAP, stationary
# (constant flow rate) ground response, no source-side load shaving. This is `gshp_simulation`'s
# fixed-COP method (see its docstring).

using Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate() # Once per machine, to install dependencies in the local environment

using GroundSourceHeatPumpDesign
using GroundHeatExchanger: borefield_rectangle, fluid_property,
    ground_load_profile, resistance_ULoop_effective
using CairoMakie

# Simulation horizon
n_years = 10
steps_per_year = 8760
nt = steps_per_year * n_years
t = collect(3600.0:3600.0:(3600.0 * nt)) # s, for the ground response and convolution
t_h = t ./ 3600 .- 1 # h, aligned with t (t_h[1] = 0)

# Building load: one year of hourly loads (synthetic profile, Bernier 2004), repeated over the
# design period. Positive = building heating, negative = building cooling.
Qb_year = -ground_load_profile(collect(0.0:1.0:(steps_per_year - 1)), 80_000.0) # W
Qb = repeat(Qb_year, n_years)

# Fixed nominal heat pump performance (no source EFT dependence), single unit, and count
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

# Operating temperature limits, checked against T_out (EWT_s, entering the heat pumps)
Tlim = [1.0, 30.0] # degC

res = gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating_unit, CAP_cooling_unit, n_HP,
    xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0)

# Table: temperatures, unmet hours and energy/power for the last simulated year
println("Fixed COP: heating $COP_heating, cooling $COP_cooling")
report_gshp_simulation(res.T_in, res.T_out, Qb, res.Qb_hp, res.Qb_unmet, res.Qg, Tlim,
    steps_per_year)

# Figures
years = t_h ./ steps_per_year

# Fig 1. Ground temperature (T_in entering, T_out leaving the field)
block_mean(x, n) = [sum(@view x[((i - 1) * n + 1):(i * n)]) / n for i in 1:(length(x) ÷ n)]
t_yearly = 1:n_years
T_in_yearly = block_mean(res.T_in, steps_per_year)
T_out_yearly = block_mean(res.T_out, steps_per_year)

fig1 = Figure()
ax = Axis(fig1[1, 1], xlabel="Time [years]", ylabel="Temperature [degC]",
    title="Ground temperature (T_in entering, T_out leaving the field)")
lines!(ax, years, res.T_in, color=:steelblue, label="T_in, hourly")
lines!(ax, years, res.T_out, color=:firebrick, label="T_out, hourly")
lines!(ax, t_yearly, T_in_yearly, color=:steelblue, linestyle=:dash, linewidth=3,
    label="T_in, yearly mean")
lines!(ax, t_yearly, T_out_yearly, color=:firebrick, linestyle=:dash, linewidth=3,
    label="T_out, yearly mean")
axislegend(ax, position=:rb, nbanks=2)
display(fig1)

# Fig 2. Building thermal load and unmet load (covered by the auxiliary system)
fig2 = Figure()
ax1 = Axis(fig2[1, 1], ylabel="Load [kW]", title="Building thermal load")
lines!(ax1, years, Qb ./ 1000, color=(:gray, 0.6), label="Building demand")
lines!(ax1, years, res.Qb_hp ./ 1000, color=:steelblue, label="Covered by heat pump(s)")
axislegend(ax1, position=:rt)
ax2 = Axis(fig2[2, 1], xlabel="Time [years]", ylabel="Auxiliary load [kW]",
    title="Unmet load (covered by the auxiliary system)")
lines!(ax2, years, res.Qb_unmet ./ 1000, color=:firebrick)
display(fig2)

# Fig 3. Ground thermal load (extraction/rejection)
fig3 = Figure()
ax = Axis(fig3[1, 1], xlabel="Time [years]", ylabel="Ground load [kW]", title="Ground thermal load")
lines!(ax, years, res.Qg ./ 1000)
display(fig3)

# Fig 4. Heat pump operation
# No COP_actual here: with fixed COP/CAP it never departs from the nominal COP_heating/COP_cooling
# above, so only the number of active units is worth plotting.
fig4 = Figure()
ax = Axis(fig4[1, 1], xlabel="Time [years]", ylabel="Active heat pumps [-]",
    title="Heat pump operation")
lines!(ax, years, res.n_active)
display(fig4)

# Fig 5. Effective borehole thermal resistance
V_range = range(0.2, 1.0, length=20) ./ 1000
Rbe_range = [resistance_ULoop_effective(Vi, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
             for Vi in V_range]

fig5 = Figure()
ax = Axis(fig5[1, 1], xlabel="Flow rate per borehole [L/s]", ylabel="Rb,eff [m.K/W]",
    title="Effective borehole thermal resistance")
lines!(ax, V_range .* 1000, Rbe_range)
scatter!(ax, [V * 1000], [res.Rbe], color=:firebrick, label="Design point")
axislegend(ax)
display(fig5)

# Fig 6. g-function (ground response)
fig6 = Figure()
ax = Axis(fig6[1, 1], xlabel="Time [s]", ylabel="g [degC.m/W]", title="g-function", xscale=log10)
lines!(ax, t, res.g)
display(fig6)