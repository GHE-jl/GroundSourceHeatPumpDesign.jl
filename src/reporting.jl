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
    - `steps_per_year`: Number of time steps per year (e.g. `8760` for an hourly simulation)
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
    - Total energy, in the same power unit times hours (e.g. kWh if `Q` is in kW)
"""
function thermal_energy(Q::AbstractVector{<:Real}, dt::Real)
    return sum(Q) * dt
end
