# Heat pump performance and building-to-ground load conversion, from the simplest description
# (fixed COP) to a more complete one (temperature-dependent COP/CAP). No ground simulation here:
# the source EFT is just swept over a representative range to show how the conversion responds to
# it. See script_gshp_1_fixedCOP.jl / script_gshp_2_variableCOP.jl for the source EFT coming from
# an actual ground simulation, and the README's "Conventions" section for the sign convention.
#
# This script only needs GroundSourceHeatPumpDesign: the heat-pump/load layer has no ground-physics
# dependency (GroundHeatExchanger.jl is only pulled in by the two gshp_* simulation scripts).

using Pkg; Pkg.activate(@__DIR__)
Pkg.instantiate()

using GroundSourceHeatPumpDesign
using CairoMakie

# Manufacturer performance data (1D: source EFT only), single heat pump unit
Tsrc_grid_heat = [-5.0, 0.0, 5.0, 10.0, 15.0]
CAP_heat_vec = [8.0, 9.0, 10.0, 11.0, 12.0] .* 1000 # W
COP_heat_vec = [3.0, 3.4, 3.8, 4.2, 4.6]
Tsrc_grid_cool = [10.0, 20.0, 30.0, 40.0]
CAP_cool_vec = [13.0, 12.0, 11.0, 10.0] .* 1000 # W
COP_cool_vec = [5.8, 5.0, 4.3, 3.7]
perf_heating = (Tsrc_grid_heat, CAP_heat_vec, COP_heat_vec)
perf_cooling = (Tsrc_grid_cool, CAP_cool_vec, COP_cool_vec)
n_HP = 5

# 1. Heat pump curves: CAP and COP interpolated over the manufacturer grid
Tsrc_heat_fine = collect(range(Tsrc_grid_heat[1], Tsrc_grid_heat[end], length=200))
Tsrc_cool_fine = collect(range(Tsrc_grid_cool[1], Tsrc_grid_cool[end], length=200))
CAP_heat_fine, COP_heat_fine = heat_pump_performance(Tsrc_heat_fine, perf_heating...)
CAP_cool_fine, COP_cool_fine = heat_pump_performance(Tsrc_cool_fine, perf_cooling...)

fig1 = Figure(size=(800, 600))
ax1 = Axis(fig1[1, 1], xlabel="Source EFT [degC]", ylabel="CAP [kW]", title="Heating capacity")
lines!(ax1, Tsrc_heat_fine, CAP_heat_fine ./ 1000)
scatter!(ax1, Tsrc_grid_heat, CAP_heat_vec ./ 1000)
ax2 = Axis(fig1[1, 2], xlabel="Source EFT [degC]", ylabel="COP [-]", title="Heating COP")
lines!(ax2, Tsrc_heat_fine, COP_heat_fine)
scatter!(ax2, Tsrc_grid_heat, COP_heat_vec)
ax3 = Axis(fig1[2, 1], xlabel="Source EFT [degC]", ylabel="CAP [kW]", title="Cooling capacity")
lines!(ax3, Tsrc_cool_fine, CAP_cool_fine ./ 1000)
scatter!(ax3, Tsrc_grid_cool, CAP_cool_vec ./ 1000)
ax4 = Axis(fig1[2, 2], xlabel="Source EFT [degC]", ylabel="COP [-]", title="Cooling COP")
lines!(ax4, Tsrc_cool_fine, COP_cool_fine)
scatter!(ax4, Tsrc_grid_cool, COP_cool_vec)
display(fig1)

# A synthetic year of building loads and illustrative source EFTs (in phase: cold source when
# heating, warm source when cooling), used to compare conversion levels below.
t_h = collect(0.0:1.0:(8760.0 - 1))
Qb = 80_000.0 .* cos.(2π .* t_h ./ 8760) # W, + heating (winter) at t=0, - cooling (summer) at t=4380
Tsrc = 12.0 .- 10.0 .* cos.(2π .* t_h ./ 8760) # degC, cold in winter, warm in summer

# 2. Fixed nominal COP: a single design COP, independent of Tsrc (Q_COP, scalar overload)
COP_heating_nominal = sum(COP_heat_vec) / length(COP_heat_vec)
COP_cooling_nominal = sum(COP_cool_vec) / length(COP_cool_vec)
Qg_fixed = Q_COP(Qb, COP_heating_nominal, COP_cooling_nominal)

# 3. Variable COP: COP read off the curve at the (here, illustrative) source EFT at each time step
# (heat_pump_performance), then Q_COP's vector-COP overload applies it point by point.
_, COP_heating_series = heat_pump_performance(Tsrc, Tsrc_grid_heat, CAP_heat_vec, COP_heat_vec)
_, COP_cooling_series = heat_pump_performance(Tsrc, Tsrc_grid_cool, CAP_cool_vec, COP_cool_vec)
Qg_variable = Q_COP(Qb, COP_heating_series, COP_cooling_series)

fig2 = Figure(size=(900, 400))
ax = Axis(fig2[1, 1], xlabel="Time [h]", ylabel="Ground load [kW]",
    title="Fixed nominal COP vs. temperature-dependent COP")
lines!(ax, t_h, Qg_fixed ./ 1000, label="Fixed nominal COP")
lines!(ax, t_h, Qg_variable ./ 1000, label="Variable COP (at illustrative Tsrc)")
axislegend(ax)
display(fig2)

# 4. Capacity capping and heat pump staging: ground_load_capped clamps the building load to what
# the installed heat pumps can deliver, then n_active_heat_pumps reports how many units that takes
# at the (illustrative) source EFT.
CAP_heating_total = n_HP * maximum(CAP_heat_vec)
CAP_cooling_total = n_HP * maximum(CAP_cool_vec)
Qg_capped, Qb_unmet = ground_load_capped(Qb, COP_heating_nominal, COP_cooling_nominal,
    CAP_heating_total, CAP_cooling_total)
Qb_hp = Qb .- Qb_unmet
n_active, COP_actual = n_active_heat_pumps(Qb_hp, Tsrc, n_HP, perf_heating, perf_cooling)

println("Unmet hours: $(count(!=(0.0), Qb_unmet)) / $(length(Qb))")
println("Mean active heat pumps: $(round(sum(n_active) / length(n_active); digits=2)) / $n_HP")

fig3 = Figure(size=(900, 450))
ax1 = Axis(fig3[1, 1], ylabel="Active heat pumps [-]",
    title="Heat pump staging (illustrative source EFT)")
lines!(ax1, t_h, n_active)
ax2 = Axis(fig3[2, 1], xlabel="Time [h]", ylabel="COP [-]")
lines!(ax2, t_h, COP_actual)
display(fig3)
