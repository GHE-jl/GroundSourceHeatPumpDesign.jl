# Interpolation grid check: at least two points, strictly increasing (required by searchsortedlast)
function _check_grid(grid::AbstractVector{<:Real}, name::String)
    # Check that the provided grid is valid: at least two points and strictly increasing.
    length(grid) >= 2 || throw(ArgumentError("$name must have at least two points."))
    # Check that the grid is strictly increasing
    all(grid[k] < grid[k + 1] for k in firstindex(grid):(lastindex(grid) - 1)) ||
        throw(ArgumentError("$name must be strictly increasing."))
    return nothing
end

"""
    heat_pump_performance(Tsrc, Tsrc_grid, CAP_vec, COP_vec)
    heat_pump_performance(Tsrc, Tload, Tsrc_grid, Tload_grid, CAP_table, COP_table)

Function that calculates heat pump capacity (CAP) and coefficient of performance (COP) based on the
source and load entering temperatures. The function has two overloads:
- 1D overload: linear interpolation on source entering fluid temperature (EFT). In this case, the
    load entering temperature is not considered and the provided CAP and COP values are only a
    function of the source EFT.
- 2D overload: bilinear interpolation on source and load entering temperatures. In this case, the
    provided CAP and COP values are a function of both source and load EFTs. The matrix construction
    is such that rows follow Tsrc_grid, columns follow Tload_grid.
For both overloads, values are clamped outside the provided grid bounds. Grids must be strictly
increasing with at least two points; `CAP_vec`/`COP_vec` must match `length(Tsrc_grid)` (1D), and
`CAP_table`/`COP_table` must have size `(length(Tsrc_grid), length(Tload_grid))` (2D). An
`ArgumentError`/`DimensionMismatch` is thrown otherwise (e.g. for a transposed table).
# Arguments
    - `Tsrc`: Source entering fluid temperature (EFT, or EWT) (Tout of the ground loop) [°C]
    - `Tload`: Load entering fluid temperature (EFT) (Tout of the building loop) [°C]
    - `Tsrc_grid`: Grid of source EFT values corresponding to the provided CAP and COP values [°C]
    - `Tload_grid`: Grid of load EFT values corresponding to the provided CAP and COP values [°C]
    - `CAP_vec`: Vector of heat pump CAPs corresponding to the source EFT grid (1D) [W]
    - `COP_vec`: Vector of heat pump COPs corresponding to the source EFT grid (1D) [-]
    - `CAP_table`: Matrix of heat pump CAPs corresponding to the source and load EFT grids (2D) [W]
    - `COP_table`: Matrix of heat pump COPs corresponding to the source and load EFT grids (2D) [-]
# Output
    - CAP, COP corresponding to the provided source and load EFTs.
"""
function heat_pump_performance(Tsrc::Real, Tsrc_grid::AbstractVector{<:Real},
    CAP_vec::AbstractVector{<:Real}, COP_vec::AbstractVector{<:Real})
    # Method 1: performance-curve-based, 1D overload
    _check_performance(Tsrc_grid, CAP_vec, COP_vec)
    return _interpolate_performance(Tsrc, Tsrc_grid, CAP_vec, COP_vec)
end
function heat_pump_performance(Tsrc::AbstractVector{<:Real}, Tsrc_grid::AbstractVector{<:Real},
    CAP_vec::AbstractVector{<:Real}, COP_vec::AbstractVector{<:Real})
    # Method 1: performance-curve-based, 1D overload (data checked once for the whole series)
    _check_performance(Tsrc_grid, CAP_vec, COP_vec)

    # Preallocation (floating-point output, even for integer inputs)
    T = promote_type(Float64, eltype(Tsrc), eltype(CAP_vec), eltype(COP_vec))
    CAP = similar(Tsrc, T)
    COP = similar(Tsrc, T)

    for k in eachindex(Tsrc)
        CAP[k], COP[k] = _interpolate_performance(Tsrc[k], Tsrc_grid, CAP_vec, COP_vec)
    end

    return CAP, COP
end

# Checks of the 1D performance data
function _check_performance(Tsrc_grid, CAP_vec::AbstractVector, COP_vec::AbstractVector)
    _check_grid(Tsrc_grid, "Tsrc_grid")
    length(CAP_vec) == length(COP_vec) == length(Tsrc_grid) ||
        throw(DimensionMismatch("CAP_vec and COP_vec must have the same length as Tsrc_grid."))
    return nothing
end

# 1D linear interpolation of the performance data (no checks)
function _interpolate_performance(Tsrc::Real, Tsrc_grid, CAP_vec::AbstractVector,
    COP_vec::AbstractVector)
    # Clamp source temperature to table range providing the evaluated temperature for interpolation
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
function heat_pump_performance(Tsrc::Real, Tload::Real, Tsrc_grid::AbstractVector{<:Real},
    Tload_grid::AbstractVector{<:Real}, CAP_table::AbstractMatrix{<:Real},
    COP_table::AbstractMatrix{<:Real})
    # Method 2: performance-curve-based, 2D overload
    _check_performance(Tsrc_grid, Tload_grid, CAP_table, COP_table)
    return _interpolate_performance(Tsrc, Tload, Tsrc_grid, Tload_grid, CAP_table, COP_table)
end

# Checks of the 2D performance data: grids and table orientation (rows follow Tsrc_grid, columns
# follow Tload_grid)
function _check_performance(Tsrc_grid, Tload_grid, CAP_table::AbstractMatrix,
    COP_table::AbstractMatrix)
    _check_grid(Tsrc_grid, "Tsrc_grid")
    _check_grid(Tload_grid, "Tload_grid")
    expected = (length(Tsrc_grid), length(Tload_grid))
    hint = "expected $expected (rows follow Tsrc_grid, columns follow Tload_grid). " *
        "Is it transposed?"
    size(CAP_table) == expected ||
        throw(DimensionMismatch("CAP_table has size $(size(CAP_table)), $hint"))
    size(COP_table) == expected ||
        throw(DimensionMismatch("COP_table has size $(size(COP_table)), $hint"))
    return nothing
end

# 2D bilinear interpolation of the performance data (no checks)
function _interpolate_performance(Tsrc::Real, Tload::Real, Tsrc_grid, Tload_grid,
    CAP_table::AbstractMatrix, COP_table::AbstractMatrix)
    # Clamp temperatures to table range
    Tsrc_eval = clamp(Tsrc, Tsrc_grid[1], Tsrc_grid[end])
    Tload_eval = clamp(Tload, Tload_grid[1], Tload_grid[end])

    # Find bracketing indices
    i = clamp(searchsortedlast(Tsrc_grid, Tsrc_eval), 1, length(Tsrc_grid) - 1)
    j = clamp(searchsortedlast(Tload_grid, Tload_eval), 1, length(Tload_grid) - 1)

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
    # Method 2: performance-curve-based, 2D overload (data checked once for the whole series)
    length(Tsrc) == length(Tload) ||
        throw(DimensionMismatch("Tsrc and Tload must have the same length."))
    _check_performance(Tsrc_grid, Tload_grid, CAP_table, COP_table)

    # Preallocation (floating-point output, even for integer inputs)
    T = promote_type(Float64, eltype(Tsrc), eltype(Tload), eltype(CAP_table), eltype(COP_table))
    CAP = similar(Tsrc, T)
    COP = similar(Tsrc, T)

    for k in eachindex(Tsrc, Tload)
        CAP[k], COP[k] = _interpolate_performance(Tsrc[k], Tload[k], Tsrc_grid, Tload_grid,
            CAP_table, COP_table)
    end

    return CAP, COP
end
"""
    part_load_COP(COP, n_active, Cd)

Effective coefficient of performance of `n_active` heat pumps (fractional, time-averaged count, see
[`n_active_heat_pumps`](@ref)) with a part-load degradation coefficient `Cd`. The `floor(n_active)`
units running the whole time step keep the full-load `COP`; the last unit only runs a fraction
`PLR = n_active - floor(n_active)` of the time step and cycles, so its part-load factor is
`PLF = 1 - Cd·(1 - PLR)`. The electricity of the time step is then proportional to
`floor(n_active) + PLR / PLF`, which gives
`COP_eff = COP · n_active / (floor(n_active) + PLR / PLF)`.
For a single unit (`n_active < 1`), this is the usual `COP_eff = PLF · COP`. With `Cd = 0`, or when
no unit cycles (`PLR = 0`), `COP_eff = COP`. Accepts scalars or vectors of the same length.
# Arguments
    - `COP`: Full-load coefficient of performance [-]
    - `n_active`: Fractional number of heat pumps operating [-]
    - `Cd`: Part-load degradation coefficient, `0 ≤ Cd < 1` (typically 0.1 to 0.25) [-]
# Output
    - `COP_eff`: Effective coefficient of performance over the time step [-]
"""
function part_load_COP(COP::Real, n_active::Real, Cd::Real)
    # Checks
    0 <= Cd < 1 || throw(ArgumentError("Cd must be in [0, 1)."))
    n_active >= 0 || throw(ArgumentError("n_active must be non-negative."))

    # Units running the whole time step, and part-load ratio of the cycling unit
    n_full = floor(n_active)
    PLR = n_active - n_full
    PLR == 0 && return float(COP)

    # Part-load factor of the cycling unit and effective COP of the time step
    PLF = 1 - Cd * (1 - PLR)
    return COP * n_active / (n_full + PLR / PLF)
end
function part_load_COP(COP::AbstractVector{<:Real}, n_active::AbstractVector{<:Real}, Cd::Real)
    length(COP) == length(n_active) ||
        throw(DimensionMismatch("COP and n_active must have the same length."))
    return part_load_COP.(COP, n_active, Cd)
end
