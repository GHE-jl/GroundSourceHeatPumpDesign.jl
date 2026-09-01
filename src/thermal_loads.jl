
"""
    Q_COP(Q, COP_heating, COP_cooling)

Convert building thermal load(s) to ground thermal load(s) based on heat pump COP. Inputs can be
scalar or vector values with the following sign convention:
- `Qb > 0`: building heating load -> `Qg < 0`: heat extracted from the ground loop
- `Qb < 0`: building cooling load -> `Qg > 0`: heat rejected to the ground loop
For heating, the conversion is `Qg = -Qb * (1 - 1 / COP_heating)`.
For cooling, the conversion is `Qg = -Qb * (1 + 1 / COP_cooling)`.
# Arguments
    - `Qb`: Building thermal load, scalar or vector [W]
    - `COP_heating`: Heating COP, scalar or vector
    - `COP_cooling`: Cooling COP, scalar or vector
# Output
    - `Qg`: Ground thermal load, scalar or vector with the same shape as `Qb` [W]
"""
function Q_COP(Qb::Real, COP_heating::Real, COP_cooling::Real)
    # Check
    COP_heating > 0 || throw(ArgumentError("COP_heating must be strictly positive."))
    COP_cooling > 0 || throw(ArgumentError("COP_cooling must be strictly positive."))

    # Conversion
    if Qb > 0
        return -Qb * (1 - 1 / COP_heating)
    elseif Qb < 0
        return -Qb * (1 + 1 / COP_cooling)
    end

    # If Qb is zero, return zero with the correct type
    return zero(promote_type(Float64, typeof(Qb), typeof(COP_heating), typeof(COP_cooling)))
end
function Q_COP(Qb::AbstractVector{<:Real}, COP_heating::Real, COP_cooling::Real)
    # Check
    T = promote_type(Float64, eltype(Qb), typeof(COP_heating), typeof(COP_cooling))
    
    # Preallocation
    Qg = similar(Qb, T)

    # Conversion
    for k in eachindex(Qb)
        Qg[k] = Q_COP(Qb[k], COP_heating, COP_cooling)
    end

    return Qg
end
function Q_COP(Qb::AbstractVector{<:Real}, COP_heating::AbstractVector{<:Real},
    COP_cooling::AbstractVector{<:Real})
    # Check
    length(Qb) == length(COP_heating) ||
        throw(DimensionMismatch("Qb and COP_heating must have the same length."))
    length(Qb) == length(COP_cooling) ||
        throw(DimensionMismatch("Qb and COP_cooling must have the same length."))
    T = promote_type(Float64, eltype(Qb), eltype(COP_heating), eltype(COP_cooling))
    
    # Preallocation
    Qg = similar(Qb, T)

    # Conversion
    for k in eachindex(Qb, COP_heating, COP_cooling)
        Qg[k] = Q_COP(Qb[k], COP_heating[k], COP_cooling[k])
    end

    return Qg
end

"""
    ground_load_capped(Qb, COP_heating, COP_cooling, CAP_heating_total, CAP_cooling_total)

Convert a building thermal load to a ground thermal load using a fixed nominal COP (rather than
[`heat_pump_performance`](@ref)'s source-EFT-dependent COP), after clamping the building load to
what the installed heat pumps can actually deliver. This is the direct, non-iterative counterpart
to [`ground_load_from_heat_pump`](@ref): since the source EFT is not yet known before the ground
temperature is simulated, the COP here is a fixed design value instead of being looked up at the
(unknown) operating temperature, so no outer iteration on the ground response is needed. Any load
beyond `CAP_heating_total`/`CAP_cooling_total` is reported as unmet rather than covered.
# Arguments
    - `Qb`: Building thermal load, scalar or vector [W]
    - `COP_heating`: Nominal heating COP used for the conversion (fixed, not EFT-dependent) [-]
    - `COP_cooling`: Nominal cooling COP used for the conversion (fixed, not EFT-dependent) [-]
    - `CAP_heating_total`: Combined heating capacity of the installed heat pumps, e.g.
        `n_HP * maximum(CAP_heat_vec)` [W]
    - `CAP_cooling_total`: Combined cooling capacity of the installed heat pumps [W]
# Output
    - `Qg`: Ground thermal load, same sign convention as [`Q_COP`](@ref) [W]
    - `Qb_unmet`: Building load left unmet by the installed heat pump capacity, same shape as `Qb`
        [W, or kW]
"""
function ground_load_capped(Qb, COP_heating::Real, COP_cooling::Real,
    CAP_heating_total::Real, CAP_cooling_total::Real)
    Qb_capped = clamp.(Qb, -CAP_cooling_total, CAP_heating_total)
    Qb_unmet = Qb .- Qb_capped
    Qg = Q_COP(Qb_capped, COP_heating, COP_cooling)
    return Qg, Qb_unmet
end

"""
    ground_load_from_heat_pump(Q_building, Tsrc, Tload, perf_heating, perf_cooling)

Convert a building thermal load to the corresponding ground load, accounting for the heat pump's
capacity limit at the current source/load entering fluid temperatures (EFTs).
The heat pump's CAP and COP are looked up from `perf_heating` (if `Q_building > 0`, heating) or
`perf_cooling` (if `Q_building < 0`, cooling) via [`heat_pump_performance`](@ref). Each is a tuple
of arguments forwarded to that function (either `(Tsrc_grid, CAP_vec, COP_vec)` for the 1D overload
or `(Tsrc_grid, Tload_grid, CAP_table, COP_table)` for the 2D overload). The heat pump's
contribution is capped at its capacity (any load beyond that is reported as unmet rather than
extrapolated).
# Arguments
    - `Q_building`: Building thermal load, `> 0` for heating, `< 0` for cooling [W]
    - `Tsrc`: Source entering fluid temperature (EFT) [°C]
    - `Tload`: Load entering fluid temperature (EFT) [°C]
    - `perf_heating`: Heat pump heating performance data, given to [`heat_pump_performance`](@ref)
    - `perf_cooling`: Heat pump cooling performance data, given to [`heat_pump_performance`](@ref)
# Output
    - `Q_ground`: Ground load using the same sign convention as [`Q_COP`](@ref) [W, or kW]
    - `Q_unmet`: Building load left unmet because it exceeds the heat pump's capacity [W, or kW]
"""
function ground_load_from_heat_pump(Q_building, Tsrc, Tload, perf_heating, perf_cooling)
    if Q_building > 0
        CAP, COP = heat_pump_performance(Tsrc, Tload, perf_heating...)
        Q_hp = min(Q_building, CAP)
        Qg = -Q_hp * (1 - 1 / COP)
        Qb_unmet = Q_building - Q_hp
    elseif Q_building < 0
        CAP, COP = heat_pump_performance(Tsrc, Tload, perf_cooling...)
        Q_hp = min(-Q_building, CAP)
        Qg = Q_hp * (1 + 1 / COP)
        Qb_unmet = -Q_building - Q_hp
    else
        Qg = 0.0
        Qb_unmet = 0.0
    end

    return Qg, Qb_unmet
end