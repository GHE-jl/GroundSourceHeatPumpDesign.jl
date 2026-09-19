"""
    unmet_hours(T, Tlim)

Count the number of time steps where `T` falls outside the operating temperature limits.
# Arguments
    - `T`: Temperature time series [°C]
    - `Tlim`: Operating temperature limits `[T_low, T_high]` [°C]
# Output
    - Number of time steps with `T < T_low` or `T > T_high`
"""
function unmet_hours(T::AbstractVector{<:Real}, Tlim::AbstractVector{<:Real})
    return count(x -> x < Tlim[1] || x > Tlim[2], T)
end

"""
    annual_drift(T, steps_per_year)

Difference between the mean temperature of the last simulated year and the first, an indicator of
the long-term thermal drift of the ground.
# Arguments
    - `T`: Temperature time series [°C], covering an integer number of years
    - `steps_per_year`: Number of time steps per year
        - E.g. `8760` for an hourly simulation
# Output
    - Annual drift, last year's mean minus first year's mean [°C]
"""
function annual_drift(T::AbstractVector{<:Real}, steps_per_year::Integer)
    first_year = @view T[1:steps_per_year]
    last_year = @view T[(end - steps_per_year + 1):end]
    return sum(last_year) / steps_per_year - sum(first_year) / steps_per_year
end

"""
    thermal_energy(Q, dt)

Total thermal energy of a load time series, by rectangular (constant time step) integration.
# Arguments
    - `Q`: Thermal load time series [W, or kW]
    - `dt`: Time step duration [h]
# Output
    - Total energy, in the same power unit times hours 
        - E.g. kWh if `Q` is in kW
"""
function thermal_energy(Q::AbstractVector{<:Real}, dt::Real)
    return sum(Q) * dt
end

"""
    report_gshp_simulation(T_in, T_out, Qb, Qb_hp, Qb_unmet, Qg, Tlim, steps_per_year)

Print a formatted report of a GSHP simulation's temperature and energy results (e.g.
[`gshp_simulation`](@ref)'s output): ground temperature extremes, unmet hours and long-term drift
over the last simulated year, then a table of energy (MWh) and peak power (kW) for the building,
heat pump and ground loads over the whole simulation.
# Arguments
    - `T_in`, `T_out`: Ground temperature entering/leaving the field, time series [°C]
    - `Qb`: Building thermal load, time series [W]
    - `Qb_hp`: Building load actually delivered by the heat pump(s), time series [W]
    - `Qb_unmet`: Building load left unmet (auxiliary), time series [W]
    - `Qg`: Ground thermal load, time series [W]
    - `Tlim`: Operating temperature limits `[T_low, T_high]`, checked against `T_out` [°C]
    - `steps_per_year`: Number of time steps per year (e.g. `8760` for an hourly simulation);
        every time series above must cover an integer number of years
"""
function report_gshp_simulation(T_in::AbstractVector{<:Real}, T_out::AbstractVector{<:Real},
    Qb::AbstractVector{<:Real}, Qb_hp::AbstractVector{<:Real}, Qb_unmet::AbstractVector{<:Real},
    Qg::AbstractVector{<:Real}, Tlim::AbstractVector{<:Real}, steps_per_year::Integer)

    # Check that all time series have the same length
    nt = length(T_in)
    n_years = nt ÷ steps_per_year
    last_year = (nt - steps_per_year + 1):nt

    # Print ground temperature extremes, unmet hours and drift over the last simulated year
    println("Ground temperature, year $n_years (last simulated year):")
    println("  T_in  min/max = $(round(minimum(T_in[last_year]); digits=1)) / " *
        "$(round(maximum(T_in[last_year]); digits=1)) degC")
    println("  T_out min/max = $(round(minimum(T_out[last_year]); digits=1)) / " *
        "$(round(maximum(T_out[last_year]); digits=1)) degC")
    println("  Hours with T_out outside [$(Tlim[1]), $(Tlim[2])] degC: " *
        "$(unmet_hours(T_out[last_year], Tlim))")
    println("  Annual drift (T_out, last year - first year): " *
        "$(round(annual_drift(T_out, steps_per_year); digits=2)) degC")

    # Energy (MWh) and peak power (kW) table
    MWh(Q) = thermal_energy(Q, 1.0) / 1e6
    rows = (
        ("Building heating", MWh(max.(Qb, 0.0)), maximum(Qb) / 1000),
        ("Building cooling", -MWh(min.(Qb, 0.0)), -minimum(Qb) / 1000),
        ("Heat pump heating", MWh(max.(Qb_hp, 0.0)), maximum(Qb_hp) / 1000),
        ("Heat pump cooling", -MWh(min.(Qb_hp, 0.0)), -minimum(Qb_hp) / 1000),
        ("Auxiliary heating", MWh(max.(Qb_unmet, 0.0)), maximum(Qb_unmet) / 1000),
        ("Auxiliary cooling", -MWh(min.(Qb_unmet, 0.0)), -minimum(Qb_unmet) / 1000),
        ("Ground extraction", -MWh(min.(Qg, 0.0)), -minimum(Qg) / 1000),
        ("Ground rejection", MWh(max.(Qg, 0.0)), maximum(Qg) / 1000),
    )

    # Print the energy and peak power table
    label_w = maximum(length(row[1]) for row in rows)
    header = rpad("", label_w) * " | " * lpad("Energy [MWh]", 12) * " | " * lpad("Peak [kW]", 9)
    println("\nEnergy and peak power over $n_years years:")
    println("  " * header)
    println("  " * "-"^length(header))
    for (label, energy, peak) in rows
        println("  " * rpad(label, label_w) * " | " * lpad(round(energy; digits=1), 12) * " | " *
            lpad(round(peak; digits=1), 9))
    end

    return nothing
end
