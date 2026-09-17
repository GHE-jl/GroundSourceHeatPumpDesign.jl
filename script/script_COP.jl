# heat_pump_performance: interpolation of a single heat pump's capacity (CAP) and coefficient of
# performance (COP) from manufacturer data, for both overloads of entering fluid temperature (EFT):
#   - 1D: CAP/COP as a function of source EFT only (Tsrc_grid, CAP_vec, COP_vec)
#   - 2D: CAP/COP as a function of source AND load EFT (Tsrc_grid, Tload_grid, CAP_table, COP_table)
# No building or ground load conversion here (see script_thermal_loads.jl for that), this script is
# only about converting the operating temperature to CAP and COP. See the README's "Conventions"
# section for the sign convention and EFT nomenclature.

using Pkg; Pkg.activate(@__DIR__)
# Pkg.instantiate() # Once per machine, to install dependencies in the local environment

using GroundSourceHeatPumpDesign
using CairoMakie

# 1D overload: heat pump CAP and COP as a function of source EFT only.

# Heat pump performance data from manufacturer, for a single heat pump model.
Tsrc_grid_heat = [-5.0, 0.0, 5.0, 10.0, 15.0]       # [°C] source EFT grid for heating
CAP_heat_vec = [8.0, 9.0, 10.0, 11.0, 12.0] .* 1000 # [W] capacity corresponding to Tsrc_grid_heat
COP_heat_vec = [3.0, 3.4, 3.8, 4.2, 4.6]            # [-] COP corresponding to Tsrc_grid_heat
Tsrc_grid_cool = [10.0, 20.0, 30.0, 40.0]           # [°C] source EFT grid for cooling
CAP_cool_vec = [13.0, 12.0, 11.0, 10.0] .* 1000     # [W] capacity corresponding to Tsrc_grid_cool
COP_cool_vec = [5.8, 5.0, 4.3, 3.7]                 # [-] COP corresponding to Tsrc_grid_cool

# Tuples for heat_pump_performance function
perf_heating_1d = (Tsrc_grid_heat, CAP_heat_vec, COP_heat_vec) # heating
perf_cooling_1d = (Tsrc_grid_cool, CAP_cool_vec, COP_cool_vec) # cooling

# Interpolated (vector overload) over a fine source EFT grid
Tsrc_heat_fine = collect(range(Tsrc_grid_heat[1]-10, Tsrc_grid_heat[end]+10, length=100))
Tsrc_cool_fine = collect(range(Tsrc_grid_cool[1]-10, Tsrc_grid_cool[end]+10, length=100))
CAP_heat_fine, COP_heat_fine = heat_pump_performance(Tsrc_heat_fine, perf_heating_1d...)
CAP_cool_fine, COP_cool_fine = heat_pump_performance(Tsrc_cool_fine, perf_cooling_1d...)

# Plotting the 1D overload results. It shows that CAP and COP fits on the specified curves.
fig1 = Figure()

ax1 = Axis(fig1[1, 1], xlabel="Source EFT [°C]", ylabel="CAP [kW]", title="Heating capacity (1D)")
lines!(ax1, Tsrc_heat_fine, CAP_heat_fine ./ 1000)
scatter!(ax1, Tsrc_grid_heat, CAP_heat_vec ./ 1000)

ax2 = Axis(fig1[1, 2], xlabel="Source EFT [°C]", ylabel="COP [-]", title="Heating COP (1D)")
lines!(ax2, Tsrc_heat_fine, COP_heat_fine)
scatter!(ax2, Tsrc_grid_heat, COP_heat_vec)

ax3 = Axis(fig1[2, 1], xlabel="Source EFT [°C]", ylabel="CAP [kW]", title="Cooling capacity (1D)")
lines!(ax3, Tsrc_cool_fine, CAP_cool_fine ./ 1000)
scatter!(ax3, Tsrc_grid_cool, CAP_cool_vec ./ 1000)

ax4 = Axis(fig1[2, 2], xlabel="Source EFT [°C]", ylabel="COP [-]", title="Cooling COP (1D)")
lines!(ax4, Tsrc_cool_fine, COP_cool_fine)
scatter!(ax4, Tsrc_grid_cool, COP_cool_vec)

display(fig1)

# 2D overload: heat pump CAP and COP as a function of source AND load EFT

Tload_grid_heat = [30.0, 35.0, 40.0, 45.0] # [°C] load EFT grid (e.g. heating loop, hot water)
# CAP and COP tables: rows follow Tsrc_grid_heat, columns follow Tload_grid_heat. CAP and COP
# increase with Tsrc (warmer source) and decrease with Tload (harder to reach a hotter supply).
CAP_table_heat = [7.0 6.6 6.2 5.8
                  8.0 7.5 7.0 6.6
                  9.0 8.4 7.8 7.3
                  10.0 9.3 8.7 8.1
                  11.0 10.2 9.5 8.9] .* 1000 # [W]
COP_table_heat = [2.8 2.5 2.2 2.0
                  3.2 2.9 2.5 2.3
                  3.6 3.2 2.9 2.6
                  4.0 3.6 3.2 2.9
                  4.4 3.9 3.5 3.1] # [-]
perf_heating_2d = (Tsrc_grid_heat, Tload_grid_heat, CAP_table_heat, COP_table_heat)

Tload_grid_cool = [7.0, 12.0, 18.0] # [°C] load EFT grid (e.g. chilled water loop, cooling)
# CAP and COP decrease with Tsrc (warmer source) and decrease with Tload (colder supply required).
CAP_table_cool = [12.0 13.0 14.0
                  11.0 12.0 13.0
                  10.0 11.0 12.0
                  9.0 10.0 11.0] .* 1000 # [W]
COP_table_cool = [5.0 5.6 6.2
                  4.4 4.9 5.4
                  3.8 4.2 4.7
                  3.2 3.6 4.0] # [-]
perf_cooling_2d = (Tsrc_grid_cool, Tload_grid_cool, CAP_table_cool, COP_table_cool)

# Scalar overload: Full surface, evaluated point by point over a fine 2D grid, for plotting
Tsrc_heat_fine2 = collect(range(Tsrc_grid_heat[1], Tsrc_grid_heat[end], length=100))
Tload_heat_fine2 = collect(range(Tload_grid_heat[1], Tload_grid_heat[end], length=100))
Tsrc_cool_fine2 = collect(range(Tsrc_grid_cool[1], Tsrc_grid_cool[end], length=100))
Tload_cool_fine2 = collect(range(Tload_grid_cool[1], Tload_grid_cool[end], length=100))

surf_heat = [heat_pump_performance(Ts, Tl, perf_heating_2d...)
    for Ts in Tsrc_heat_fine2, Tl in Tload_heat_fine2]
CAP_surf_heat = first.(surf_heat)
COP_surf_heat = last.(surf_heat)

surf_cool = [heat_pump_performance(Ts, Tl, perf_cooling_2d...)
    for Ts in Tsrc_cool_fine2, Tl in Tload_cool_fine2]
CAP_surf_cool = first.(surf_cool)
COP_surf_cool = last.(surf_cool)

fig2 = Figure()
ax1 = Axis(fig2[1, 1], xlabel="Source EFT [°C]", ylabel="Load EFT [°C]",
    title="Heating CAP [kW] (2D)")
hm1 = heatmap!(ax1, Tsrc_heat_fine2, Tload_heat_fine2, CAP_surf_heat ./ 1000)
Colorbar(fig2[1, 2], hm1)
ax2 = Axis(fig2[1, 3], xlabel="Source EFT [°C]", ylabel="Load EFT [°C]", title="Heating COP (2D)")
hm2 = heatmap!(ax2, Tsrc_heat_fine2, Tload_heat_fine2, COP_surf_heat)
Colorbar(fig2[1, 4], hm2)
ax3 = Axis(fig2[2, 1], xlabel="Source EFT [°C]", ylabel="Load EFT [°C]",
    title="Cooling CAP [kW] (2D)")
hm3 = heatmap!(ax3, Tsrc_cool_fine2, Tload_cool_fine2, CAP_surf_cool ./ 1000)
Colorbar(fig2[2, 2], hm3)
ax4 = Axis(fig2[2, 3], xlabel="Source EFT [°C]", ylabel="Load EFT [°C]", title="Cooling COP (2D)")
hm4 = heatmap!(ax4, Tsrc_cool_fine2, Tload_cool_fine2, COP_surf_cool)
Colorbar(fig2[2, 4], hm4)
display(fig2)

# Vector overload: a paired (Tsrc, Tload) trajectory, as would occur along a simulated time series,
# swept diagonally across the 2D domain for illustration.
Tsrc_heat_vec = collect(range(Tsrc_grid_heat[1], Tsrc_grid_heat[end], length=100))
Tload_heat_vec = collect(range(Tload_grid_heat[1], Tload_grid_heat[end], length=100))
CAP_vec, COP_vec = heat_pump_performance(Tsrc_heat_vec, Tload_heat_vec, perf_heating_2d...)

Tsrc_cool_vec = collect(range(Tsrc_grid_cool[1], Tsrc_grid_cool[end], length=100))
Tload_cool_vec = collect(range(Tload_grid_cool[1], Tload_grid_cool[end], length=100))
CAP_vec, COP_vec = heat_pump_performance(Tsrc_cool_vec, Tload_cool_vec, perf_cooling_2d...)

fig3 = Figure()
ax1 = Axis(fig3[1, 1], xlabel="Source EFT [°C]", ylabel="CAP [kW]",
    title="Heating capacity along a (Tsrc, Tload) path")
lines!(ax1, Tsrc_heat_vec, CAP_vec ./ 1000)
ax2 = Axis(fig3[1, 2], xlabel="Source EFT [°C]", ylabel="COP [-]",
    title="Heating COP along a (Tsrc, Tload) path")
lines!(ax2, Tsrc_heat_vec, COP_vec)
ax3 = Axis(fig3[2, 1], xlabel="Source EFT [°C]", ylabel="CAP [kW]",
    title="Cooling capacity along a (Tsrc, Tload) path")
lines!(ax3, Tsrc_cool_vec, CAP_vec ./ 1000)
ax4 = Axis(fig3[2, 2], xlabel="Source EFT [°C]", ylabel="COP [-]",
    title="Cooling COP along a (Tsrc, Tload) path")
lines!(ax4, Tsrc_cool_vec, COP_vec)
display(fig3)
