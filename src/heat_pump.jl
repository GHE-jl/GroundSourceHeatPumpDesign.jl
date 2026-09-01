#TODO Add the degradation coefficient to the heat pump performance function and use it in the sizing methods.

"""
    heat_pump_performance(Tsrc, Tsrc_grid, CAP_vec, COP_vec)
    heat_pump_performance(Tsrc, Tload, Tsrc_grid, Tload_grid, CAP_table, COP_table)

Function that calculates heat pump capacity (CAP) and coefficient of performance (COP) based on the
source and load entering temperatures. The function has two overloads:
- 1D overload: linear interpolation on source entering fluid temperature (EFT). In this case, the
    load entering temperature is not considered and the provided CAP and COP values are only a
    function of the source EFT.
- 2D overload: bilinear interpolation on source and load entering temperatures. In this case, the
    provided CAP and COP values are a function of both source and load EFTs.
For both overloads, values are clamped outside the provided grid bounds.
# Arguments
    - `Tsrc`: Source entering fluid temperature (EFT, or EWT) (Tout of the ground loop) [°C]
    - `Tload`: Load entering fluid temperature (EFT) (Tout of the building loop) [°C]
    - `Tsrc_grid`: Grid of source EFT values corresponding to the provided CAP and COP values [°C]
    - `Tload_grid`: Grid of load EFT values corresponding to the provided CAP and COP values [°C]
    - `CAP_vec`: Vector of heat pump CAPs corresponding to the source EFT grid (1D) [W]
    - `COP_vec`: Vector of heat pump COPs corresponding to the source EFT grid (1D) [W]
    - `CAP_table`: Matrix of heat pump CAPs corresponding to the source and load EFT grids (2D) [W]
    - `COP_table`: Matrix of heat pump COPs corresponding to the source and load EFT grids (2D) [W]
# Output
    - Heat pump CAP and coefficient of performance (COP) corresponding to the provided
        source and load EFTs.
"""
function heat_pump_performance(Tsrc::Real, Tsrc_grid::AbstractVector{<:Real},
    CAP_vec::AbstractVector{<:Real}, COP_vec::AbstractVector{<:Real})

    # Clamp source temperature to table range
    Tsrc_eval = clamp(Tsrc, Tsrc_grid[1], Tsrc_grid[end])

    # Exact lower and upper bound
    if Tsrc_eval == Tsrc_grid[1]
        return CAP_vec[1], COP_vec[1]
    end
    if Tsrc_eval == Tsrc_grid[end]
        return CAP_vec[end], COP_vec[end]
    end

    # Find interval
    i = searchsortedlast(Tsrc_grid, Tsrc_eval)  # Find the largest index Tsrc_grid[i] <= Tsrc_eval
    i = clamp(i, 1, length(Tsrc_grid) - 1) # Ensure i is within valid range for interpolation

    T1, T2 = Tsrc_grid[i], Tsrc_grid[i + 1]
    w = (Tsrc_eval - T1) / (T2 - T1)

    # Linear interpolation
    CAP = (1 - w) * CAP_vec[i] + w * CAP_vec[i + 1]
    COP = (1 - w) * COP_vec[i] + w * COP_vec[i + 1]

    return CAP, COP
end
function heat_pump_performance(Tsrc::AbstractVector{<:Real}, Tsrc_grid::AbstractVector{<:Real},
    CAP_vec::AbstractVector{<:Real}, COP_vec::AbstractVector{<:Real})

    CAP = similar(Tsrc)
    COP = similar(Tsrc)

    for k in eachindex(Tsrc)
        CAP[k], COP[k] = heat_pump_performance(Tsrc[k], Tsrc_grid, CAP_vec, COP_vec)
    end

    return CAP, COP
end
function heat_pump_performance(Tsrc::Real, Tload::Real, Tsrc_grid::AbstractVector{<:Real},
    Tload_grid::AbstractVector{<:Real}, CAP_table::AbstractMatrix{<:Real},
    COP_table::AbstractMatrix{<:Real})

    # Clamp temperatures to table range
    Tsrc_eval = clamp(Tsrc, Tsrc_grid[1], Tsrc_grid[end])
    Tload_eval = clamp(Tload, Tload_grid[1], Tload_grid[end])

    # Exact lower and upper bounds
    if Tsrc_eval == Tsrc_grid[1] && Tload_eval == Tload_grid[1]
        return CAP_table[1, 1], COP_table[1, 1]
    end
    if Tsrc_eval == Tsrc_grid[end] && Tload_eval == Tload_grid[end]
        return CAP_table[end, end], COP_table[end, end]
    end

    # Find bracketing indices
    i = searchsortedlast(Tsrc_grid, Tsrc_eval)
    j = searchsortedlast(Tload_grid, Tload_eval)

    i = clamp(i, 1, length(Tsrc_grid) - 1)
    j = clamp(j, 1, length(Tload_grid) - 1)

    # Grid points
    Tsrc1, Tsrc2 = Tsrc_grid[i], Tsrc_grid[i + 1]
    Tload1, Tload2 = Tload_grid[j], Tload_grid[j + 1]

    # Interpolation weights
    wx = (Tsrc_eval - Tsrc1) / (Tsrc2 - Tsrc1)
    wy = (Tload_eval - Tload1) / (Tload2 - Tload1)

    # Bilinear interpolation helper
    function bilinear_interp(table)
        f11 = table[i, j]
        f12 = table[i, j + 1]
        f21 = table[i + 1, j]
        f22 = table[i + 1, j + 1]

        return (1 - wx) * (1 - wy) * f11 + (1 - wx) * wy * f12 + wx * (1 - wy) * f21 +
            wx * wy * f22
    end

    CAP = bilinear_interp(CAP_table)
    COP = bilinear_interp(COP_table)

    return CAP, COP
end
function heat_pump_performance(Tsrc::AbstractVector{<:Real}, Tload::AbstractVector{<:Real},
    Tsrc_grid::AbstractVector{<:Real}, Tload_grid::AbstractVector{<:Real},
    CAP_table::AbstractMatrix{<:Real}, COP_table::AbstractMatrix{<:Real})

    CAP = similar(Tsrc)
    COP = similar(Tsrc)

    for k in eachindex(Tsrc)
        CAP[k], COP[k] = heat_pump_performance(Tsrc[k], Tload[k], Tsrc_grid, Tload_grid, CAP_table,
            COP_table)
    end

    return CAP, COP
end

"""
    n_active_heat_pumps(Qb_capped, Tsrc, n_HP, perf_heating, perf_cooling)

Number of heat pumps operating at each time step to deliver `Qb_capped`, given the source EFT
`Tsrc` (e.g. from a ground temperature simulation) and a single unit's performance curve. Following
Pasquier et al. (2013), this is a fractional, time-averaged count (between `0` and `n_HP`) rather
than an integer, since a heat pump typically cycles within a time step instead of switching exactly
on its boundary.
# Arguments
    - `Qb_capped`: Building load actually delivered by the heat pumps (post-capacity-capping, e.g.
        from [`ground_load_capped`](@ref)) [W]
    - `Tsrc`: Source entering fluid temperature (EFT) [°C]
    - `n_HP`: Number of installed heat pumps
    - `perf_heating`: Single heat pump's heating performance data, forwarded to
        [`heat_pump_performance`](@ref) (`(Tsrc_grid, CAP_vec, COP_vec)`)
    - `perf_cooling`: Single heat pump's cooling performance data, forwarded to
        [`heat_pump_performance`](@ref)
# Output
    - `n`: Fractional number of heat pumps operating, in `[0, n_HP]`
    - `COP`: Coefficient of performance at each time step, from the same interpolation as `n`
# Reference
    - Pasquier, P., Marcotte, D., Bernier, M., & Kummert, M. (2013). Simulation of ground-coupled
        heat pump systems using a spectral approach. Proceedings of BS2013, 2691–2698.
"""
function n_active_heat_pumps(Qb_capped::AbstractVector{<:Real}, Tsrc::AbstractVector{<:Real},
    n_HP::Real, perf_heating, perf_cooling)

    # Check
    length(Qb_capped) == length(Tsrc) ||
        throw(DimensionMismatch("Qb_capped and Tsrc must have the same length."))

    n = similar(Qb_capped, Float64)
    COP = similar(Qb_capped, Float64)
    for k in eachindex(Qb_capped, Tsrc)
        CAPk, COPk = Qb_capped[k] >= 0 ? heat_pump_performance(Tsrc[k], perf_heating...) :
                     heat_pump_performance(Tsrc[k], perf_cooling...)
        COP[k] = COPk
        n[k] = CAPk > 0 ? clamp(abs(Qb_capped[k]) / CAPk, 0.0, n_HP) : 0.0
    end
    return n, COP
end