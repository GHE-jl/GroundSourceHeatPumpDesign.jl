# GSHP simulation, case 3: fixed COP/CAP, non-stationary (time-varying) flow rate, no source-side
# shaving. This is `gshp_simulation`'s non-stationary method — see its docstring, and
# `_flow_rate_states`'s, for the load-based flow-rate control law that keeps this a single,
# non-iterative pass. See script_gshp_1_fixedCOP.jl for the stationary (constant-flow) counterpart
# this is compared against below, and the README's "Conventions" section for the sign convention and
# T_in/T_out nomenclature.
#
# To isolate the effect of switching flow rates, both runs below share the exact same building load,
# heat pump sizing and borehole field — only the flow-rate strategy differs:
#   - "stationary": constant flow rate, held at the middle candidate (30 L/min per borehole)
#   - "non-stationary": switches between 25/30/35 L/min per borehole, picked from the ground load
#
# What's guaranteed: the inlet/outlet split |T_in - T_out| shrinks monotonically with a bigger flow
# rate for the same ground load (it's proportional to Qg/V), so the non-stationary run's worst-case
# split is always <= the stationary run's held at the same flow. What's NOT guaranteed: whether that
# also lowers the peak T_out specifically. T_out is set by two competing effects of a bigger V — a
# smaller |T_in-T_out| split (pulls T_out toward the mean fluid temperature Tf) and a smaller
# effective borehole resistance Rbₑ (pulls Tf itself down for a heat-rejecting/cooling load, up for a
# heat-extracting/heating one) — and which one wins depends on the sign and size of the load. Run the
# numbers below rather than assuming a single direction.

using Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate() # Once per machine, to install dependencies in the local environment

using GroundSourceHeatPumpDesign
using GroundHeatExchanger: borefield_rectangle, fluid_property,
    ground_load_profile
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

# Pipe, grout and fluid properties (single U-loop, water)
s = 0.05 # shank spacing [m]
ro, ri = 0.022, 0.017 # pipe outer/inner radius [m]
kg, kp = 1.6, 0.4 # grout, pipe thermal conductivity [W/mK]
kf, cf, ρf, μf = fluid_property(T0, :water)

# Controllable flow-rate set points and target inlet/outlet split used to pick between them
V_candidates = [25.0, 30.0, 35.0] ./ 60_000 # 25/30/35 L/min per borehole -> m3/s
V_nominal = V_candidates[2] # middle candidate, used for the stationary comparison run
ΔT_target = 2.5 # K, field inlet/outlet split targeted when picking the active flow rate

# Operating temperature limits, checked against T_out (EWT_s, entering the heat pumps) -- reporting
# only: neither run below enforces this, `unmet_hours` is just applied to the raw T_out afterwards.
Tlim = [1.0, 30.0] # degC

res_stationary = gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating_unit,
    CAP_cooling_unit, n_HP, xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V_nominal, T0)

res_ns = gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating_unit, CAP_cooling_unit, n_HP,
    xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V_candidates, ΔT_target, T0)

# Table: temperatures and unmet hours for the last simulated year, stationary vs. non-stationary
last_year = (nt - steps_per_year + 1):nt
println("Flow-rate candidates [L/min]: ", round.(V_candidates .* 60_000, digits=1),
    " (stationary run held at $(round(V_nominal * 60_000, digits=1)) L/min)")
for (name, res) in ("stationary" => res_stationary, "non-stationary" => res_ns)
    ΔT = abs.(res.T_in[last_year] .- res.T_out[last_year])
    println("\n$name:")
    println("  T_out min/max = $(round(minimum(res.T_out[last_year]); digits=2)) / " *
        "$(round(maximum(res.T_out[last_year]); digits=2)) degC")
    println("  |T_in - T_out| max = $(round(maximum(ΔT); digits=2)) degC")
    println("  Hours with T_out outside [$(Tlim[1]), $(Tlim[2])] degC: " *
        "$(unmet_hours(res.T_out[last_year], Tlim))")
end
println("\nThe inlet/outlet split shrinks with the non-stationary run (guaranteed, it's the flow-rate " *
    "control law's whole point). Whether the heating-season T_out minimum rises and/or the " *
    "cooling-season T_out maximum falls too depends on the load sign/size at each extreme (see the " *
    "file header) -- compare the two T_out min/max lines above rather than assuming both improve.")

# Figures
years = t_h ./ steps_per_year
last_year_h = t_h[last_year] .- t_h[last_year[1]]

block_mean(x, n) = [sum(@view x[((i - 1) * n + 1):(i * n)]) / n for i in 1:(length(x) ÷ n)]
t_monthly = range(0, n_years, length=nt ÷ 730) # 8760 / 12
T_out_s_monthly = block_mean(res_stationary.T_out, 730)
T_out_ns_monthly = block_mean(res_ns.T_out, 730)

fig1 = Figure(size=(900, 450))
ax = Axis(fig1[1, 1], xlabel="Time [years]", ylabel="T_out [degC]",
    title="Ground outlet temperature: stationary vs. non-stationary flow")
lines!(ax, years, res_stationary.T_out, color=(:firebrick, 0.2), label="stationary, hourly")
lines!(ax, years, res_ns.T_out, color=(:steelblue, 0.2), label="non-stationary, hourly")
lines!(ax, t_monthly, T_out_s_monthly, color=:firebrick, linewidth=2, label="stationary, monthly")
lines!(ax, t_monthly, T_out_ns_monthly, color=:steelblue, linewidth=2, label="non-stationary, monthly")
axislegend(ax, position=:rb, nbanks=2)
display(fig1)

fig2 = Figure(size=(900, 500))
ax1 = Axis(fig2[1, 1], ylabel="Flow rate [L/min]", title="Active flow rate, last simulated year")
lines!(ax1, last_year_h ./ 24, fill(V_nominal * 60_000, length(last_year)), color=:firebrick,
    label="stationary")
lines!(ax1, last_year_h ./ 24, res_ns.V_t[last_year] .* 60_000, color=:steelblue,
    label="non-stationary")
axislegend(ax1, position=:rb)
ax2 = Axis(fig2[2, 1], xlabel="Time [days]", ylabel="|T_in - T_out| [degC]",
    title="Inlet/outlet split, last simulated year")
lines!(ax2, last_year_h ./ 24, abs.(res_stationary.T_in[last_year] .- res_stationary.T_out[last_year]),
    color=:firebrick, label="stationary")
lines!(ax2, last_year_h ./ 24, abs.(res_ns.T_in[last_year] .- res_ns.T_out[last_year]),
    color=:steelblue, label="non-stationary")
axislegend(ax2, position=:rt)
display(fig2)

fig3 = Figure(size=(600, 400))
ax = Axis(fig3[1, 1], xlabel="Flow-rate state (sorted candidates)", ylabel="Rb,eff [m.K/W]",
    title="Effective borehole thermal resistance per flow-rate state")
barplot!(ax, 1:length(V_candidates), res_ns.Rbe)
ax.xticks = (1:length(V_candidates), string.(round.(sort(V_candidates) .* 60_000, digits=1)) .* " L/min")
display(fig3)
