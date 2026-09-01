# GSHP simulation, case 2: temperature-dependent COP/CAP (heat_pump_performance curves),
# stationary (constant flow rate) ground response, no source-side load shaving. This is
# `gshp_simulation`'s curve-based method — see its docstring for how the chicken-and-egg problem
# (the source EFT isn't known until the ground response is simulated) is resolved with a single,
# non-iterative pass. See script_gshp_1_fixedCOP.jl for the simpler, fixed-COP counterpart, and the
# README's "Conventions" section for the sign convention and T_in/T_out nomenclature.

using Pkg; Pkg.activate(@__DIR__)
Pkg.instantiate()

using GroundSourceHeatPumpDesign
using GroundHeatExchanger: borefield_rectangle, water_k, water_cp, water_ρ, water_μ,
    ground_load_profile, resistance_ULoop_effective
using CairoMakie

# Simulation horizon
n_years = 10
steps_per_year = 8760
nt = steps_per_year * n_years
t = collect(3600.0:3600.0:(3600.0 * nt)) # s, for the ground response and convolution
t_h = t ./ 3600 .- 1 # h, aligned with t (t_h[1] = 0)

# Building load: one year of hourly loads (synthetic profile, Bernier 2004), repeated over the
# design period. Positive = heating, negative = cooling.
Qb_year = -ground_load_profile(collect(0.0:1.0:(steps_per_year - 1)), 80_000.0) # W
Qb = repeat(Qb_year, n_years)

# Heat pump performance (single unit, source EFT only) and count
Tsrc_grid_heat = [-5.0, 0.0, 5.0, 10.0, 15.0]
CAP_heat_vec = [8.0, 9.0, 10.0, 11.0, 12.0] .* 1000 # W
COP_heat_vec = [3.0, 3.4, 3.8, 4.2, 4.6]
Tsrc_grid_cool = [10.0, 20.0, 30.0, 40.0]
CAP_cool_vec = [13.0, 12.0, 11.0, 10.0] .* 1000 # W
COP_cool_vec = [5.8, 5.0, 4.3, 3.7]
perf_heating = (Tsrc_grid_heat, CAP_heat_vec, COP_heat_vec)
perf_cooling = (Tsrc_grid_cool, CAP_cool_vec, COP_cool_vec)
n_HP = 5

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
kf, cf, ρf, μf = water_k(T0), water_cp(T0), water_ρ(T0), water_μ(T0)
V = 0.5 / 1000 # flow rate per borehole [m3/s]

# Operating temperature limits, checked against T_out (EWT_s, entering the heat pumps)
Tlim = [1.0, 30.0] # degC

res = gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP, xy, rb, H, D, ks, Cs, s, ro, ri,
    kg, kp, kf, cf, ρf, μf, V, T0)

# Table: temperatures, unmet hours and energy/power for the last simulated year
last_year = (nt - steps_per_year + 1):nt
println("Ground temperature, year $n_years (last simulated year):")
println("  T_in  min/max = $(round(minimum(res.T_in[last_year]); digits=1)) / " *
        "$(round(maximum(res.T_in[last_year]); digits=1)) degC")
println("  T_out min/max = $(round(minimum(res.T_out[last_year]); digits=1)) / " *
        "$(round(maximum(res.T_out[last_year]); digits=1)) degC")
println("  Hours with T_out outside [$(Tlim[1]), $(Tlim[2])] degC: " *
        "$(unmet_hours(res.T_out[last_year], Tlim))")
println("  Annual drift (T_out, last year - first year): " *
        "$(round(annual_drift(res.T_out, steps_per_year); digits=2)) degC")

MWh(Q) = thermal_energy(Q, 1.0) / 1e6
println("\nEnergy (MWh) and peak power (kW) over $n_years years:")
println("  Building heating : $(round(MWh(max.(Qb, 0.0)); digits=1)) MWh, " *
        "peak $(round(maximum(Qb) / 1000; digits=1)) kW")
println("  Building cooling : $(round(-MWh(min.(Qb, 0.0)); digits=1)) MWh, " *
        "peak $(round(-minimum(Qb) / 1000; digits=1)) kW")
println("  Heat pump heating: $(round(MWh(max.(res.Qb_hp, 0.0)); digits=1)) MWh, " *
        "peak $(round(maximum(res.Qb_hp) / 1000; digits=1)) kW")
println("  Heat pump cooling: $(round(-MWh(min.(res.Qb_hp, 0.0)); digits=1)) MWh, " *
        "peak $(round(-minimum(res.Qb_hp) / 1000; digits=1)) kW")
println("  Auxiliary heating: $(round(MWh(max.(res.Qb_unmet, 0.0)); digits=1)) MWh, " *
        "peak $(round(maximum(res.Qb_unmet) / 1000; digits=1)) kW")
println("  Auxiliary cooling: $(round(-MWh(min.(res.Qb_unmet, 0.0)); digits=1)) MWh, " *
        "peak $(round(-minimum(res.Qb_unmet) / 1000; digits=1)) kW")
println("  Ground extraction: $(round(-MWh(min.(res.Qg, 0.0)); digits=1)) MWh, " *
        "peak $(round(-minimum(res.Qg) / 1000; digits=1)) kW")
println("  Ground rejection : $(round(MWh(max.(res.Qg, 0.0)); digits=1)) MWh, " *
        "peak $(round(maximum(res.Qg) / 1000; digits=1)) kW")

# Figures
years = t_h ./ steps_per_year

block_mean(x, n) = [sum(@view x[((i - 1) * n + 1):(i * n)]) / n for i in 1:(length(x) ÷ n)]
t_monthly = range(0, n_years, length=nt ÷ 730) # 8760 / 12
t_yearly = 1:n_years
T_in_monthly, T_in_yearly = block_mean(res.T_in, 730), block_mean(res.T_in, steps_per_year)
T_out_monthly, T_out_yearly = block_mean(res.T_out, 730), block_mean(res.T_out, steps_per_year)

fig1 = Figure(size=(900, 450))
ax = Axis(fig1[1, 1], xlabel="Time [years]", ylabel="Temperature [degC]",
    title="Ground temperature (T_in entering, T_out leaving the field)")
lines!(ax, years, res.T_in, color=(:steelblue, 0.25), label="T_in, hourly")
lines!(ax, years, res.T_out, color=(:firebrick, 0.25), label="T_out, hourly")
lines!(ax, t_monthly, T_in_monthly, color=:steelblue, linewidth=2, label="T_in, monthly mean")
lines!(ax, t_monthly, T_out_monthly, color=:firebrick, linewidth=2, label="T_out, monthly mean")
lines!(ax, t_yearly, T_in_yearly, color=:steelblue, linestyle=:dash, linewidth=3,
    label="T_in, yearly mean")
lines!(ax, t_yearly, T_out_yearly, color=:firebrick, linestyle=:dash, linewidth=3,
    label="T_out, yearly mean")
axislegend(ax, position=:rb, nbanks=2)
display(fig1)

fig2 = Figure(size=(900, 500))
ax1 = Axis(fig2[1, 1], ylabel="Load [kW]", title="Building thermal load")
lines!(ax1, years, Qb ./ 1000, color=(:gray, 0.6), label="Building demand")
lines!(ax1, years, res.Qb_hp ./ 1000, color=:steelblue, label="Covered by heat pump(s)")
axislegend(ax1, position=:rt)
ax2 = Axis(fig2[2, 1], xlabel="Time [years]", ylabel="Auxiliary load [kW]",
    title="Unmet load (covered by the auxiliary system)")
lines!(ax2, years, res.Qb_unmet ./ 1000, color=:firebrick)
display(fig2)

fig3 = Figure(size=(900, 350))
ax = Axis(fig3[1, 1], xlabel="Time [years]", ylabel="Ground load [kW]", title="Ground thermal load")
lines!(ax, years, res.Qg ./ 1000)
display(fig3)

fig4 = Figure(size=(900, 450))
ax1 = Axis(fig4[1, 1], ylabel="Active heat pumps [-]", title="Heat pump operation")
lines!(ax1, years, res.n_active)
ax2 = Axis(fig4[2, 1], xlabel="Time [years]", ylabel="COP [-]")
lines!(ax2, years, res.COP_actual)
display(fig4)

V_range = range(0.2, 1.0, length=20) ./ 1000
Rbe_range = [resistance_ULoop_effective(Vi, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
             for Vi in V_range]
fig5 = Figure(size=(600, 400))
ax = Axis(fig5[1, 1], xlabel="Flow rate per borehole [L/s]", ylabel="Rb,eff [m.K/W]",
    title="Effective borehole thermal resistance")
lines!(ax, V_range .* 1000, Rbe_range)
scatter!(ax, [V * 1000], [res.Rbe], color=:firebrick, label="Design point")
axislegend(ax)
display(fig5)

fig6 = Figure(size=(600, 400))
ax = Axis(fig6[1, 1], xlabel="Time [s]", ylabel="g [degC.m/W]", title="g-function", xscale=log10)
lines!(ax, t, res.g)
display(fig6)
