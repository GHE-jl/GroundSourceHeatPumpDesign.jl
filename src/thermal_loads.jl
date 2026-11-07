"""
    ground_load(Qb, COP_heating, COP_cooling)
    ground_load(Qb, COP_heating, COP_cooling, CAP_heating_unit, CAP_cooling_unit; n_HP=1)
    ground_load(Qb, Tsrc, perf_heating, perf_cooling; n_HP=1)
    ground_load(Qb, Tsrc, Tload, perf_heating, perf_cooling; n_HP=1)

Convert a building thermal load to the corresponding ground thermal load, with the sign convention:
- `Qb > 0`: building heating load -> `Qg < 0`: heat extracted from the ground loop
- `Qb < 0`: building cooling load -> `Qg > 0`: heat rejected to the ground loop

# Method 1: fixed COP, no capacity limit: `ground_load(Qb, COP_heating, COP_cooling)`
The base conversion: `Qg = -Qb * (1 - 1 / COP_heating)` for heating, and
`Qg = -Qb * (1 + 1 / COP_cooling)` for cooling. `COP_heating` and `COP_cooling` are a single nominal
value, or, if `Qb` is a vector, can also be vectors of the same length (a per-time-step COP, e.g.
read [`heat_pump_performance`](@ref) beforehand). Returns `Qg` (there is no unmet load without a
capacity limit, so no Qb_unmet or Qb_hp).

# Method 2: fixed COP, capacity-capped: `ground_load(Qb, COP_heating, COP_cooling,
CAP_heating_unit, CAP_cooling_unit; n_HP=1)`
Clamps `Qb` to `[-n_HP * CAP_cooling_unit, n_HP * CAP_heating_unit]`, then applies method 1 to the
clamped load. This is the fixed-COP counterpart to method 3: used when the source entering fluid
temperature (EFT) (and thus the true, temperature-dependent COP) isn't known yet, so a fixed
nominal COP is used instead. Any load beyond capacity is reported as unmet rather than covered.
Returns `(Qg, Qb_unmet, Qb_hp)`.

# Method 3: performance-curve-based: `ground_load(Qb, Tsrc, [Tload,] perf_heating, perf_cooling;
n_HP=1)`
Looks up CAP and COP from `perf_heating` (`Qb >= 0`, heating) or `perf_cooling` (`Qb < 0`, cooling)
via [`heat_pump_performance`](@ref) at the current source (and, for the 2D overload, load) EFT, caps
the heat pump's contribution at `n_HP` times that CAP, then applies method 1 with the resulting COP.
`perf_heating` and `perf_cooling` are tuples of arguments forwarded to `heat_pump_performance`.
Either `(Tsrc_grid, CAP_vec, COP_vec)` for the 1D overload or `(Tsrc_grid, Tload_grid, CAP_table,
COP_table)` for the 2D overload. Returns `(Qg, Qb_unmet, Qb_hp)`.

# Arguments
    - `Qb`: Building thermal load, scalar or vector [W]
    - `COP_heating`, `COP_cooling`: Heating/cooling COP, scalar, or vector matching `Qb` [-]
    - `CAP_heating_unit`, `CAP_cooling_unit`: Heating/cooling capacity of one heat pump unit [W]
    - `n_HP`: Number of installed heat pump units (methods 2-3 only); scales `CAP_heating_unit`/
        `CAP_cooling_unit` (method 2), or the CAP looked up from `perf_heating`/`perf_cooling`
        (method 3), to the total installed capacity. Defaults to `1`, i.e. the CAP argument(s)
        already represent the total installed capacity
    - `Tsrc`, `Tload`: Source/load entering fluid temperature, scalar or vector matching `Qb` [°C]
    - `perf_heating`, `perf_cooling`: Single heat pump unit's heating/cooling performance data,
        forwarded to [`heat_pump_performance`](@ref)
# Output
    - `Qg`: Ground thermal load, same shape as `Qb` [W]
    - `Qb_unmet` (methods 2-3 only): Building load left unmet by the installed heat pump capacity,
        same shape as `Qb` [W]
    - `Qb_hp` (methods 2-3 only): Building load actually delivered by the heat pumps,
        `Qb_hp = Qb - Qb_unmet`, same shape as `Qb` [W]
"""
function ground_load(Qb::Real, COP_heating::Real, COP_cooling::Real)
    # Method 1: fixed COP, no capacity limit
    # Check
    COP_heating > 0 || throw(ArgumentError("COP_heating must be strictly positive."))
    COP_cooling > 0 || throw(ArgumentError("COP_cooling must be strictly positive."))

    # Conversion (same output type for every branch)
    T = promote_type(Float64, typeof(Qb), typeof(COP_heating), typeof(COP_cooling))
    if Qb > 0
        return T(-Qb * (1 - 1 / COP_heating))
    elseif Qb < 0
        return T(-Qb * (1 + 1 / COP_cooling))
    else
        return zero(T)
    end
end
function ground_load(Qb::AbstractVector{<:Real}, COP_heating::Real, COP_cooling::Real)
    # Method 1: fixed COP, no capacity limit
    # Preallocation
    T = promote_type(Float64, eltype(Qb), typeof(COP_heating), typeof(COP_cooling))
    Qg = similar(Qb, T)

    # Conversion
    for k in eachindex(Qb)
        Qg[k] = ground_load(Qb[k], COP_heating, COP_cooling)
    end

    return Qg
end
function ground_load(Qb::AbstractVector{<:Real}, COP_heating::AbstractVector{<:Real},
    COP_cooling::AbstractVector{<:Real})
    # Method 1: fixed COP, no capacity limit
    # Check
    length(Qb) == length(COP_heating) ||
        throw(DimensionMismatch("Qb and COP_heating must have the same length."))
    length(Qb) == length(COP_cooling) ||
        throw(DimensionMismatch("Qb and COP_cooling must have the same length."))

    # Preallocation
    T = promote_type(Float64, eltype(Qb), eltype(COP_heating), eltype(COP_cooling))
    Qg = similar(Qb, T)

    # Conversion
    for k in eachindex(Qb, COP_heating, COP_cooling)
        Qg[k] = ground_load(Qb[k], COP_heating[k], COP_cooling[k])
    end

    return Qg
end
function ground_load(Qb::Real, COP_heating::Real, COP_cooling::Real,
    CAP_heating_unit::Real, CAP_cooling_unit::Real; n_HP::Real=1)
    # Method 2: fixed COP, capacity-capped
    Qb_hp = clamp(Qb, -n_HP * CAP_cooling_unit, n_HP * CAP_heating_unit)
    Qb_unmet = Qb - Qb_hp
    Qg = ground_load(Qb_hp, COP_heating, COP_cooling)
    return Qg, Qb_unmet, Qb_hp
end
function ground_load(Qb::AbstractVector{<:Real}, COP_heating::Real, COP_cooling::Real,
    CAP_heating_unit::Real, CAP_cooling_unit::Real; n_HP::Real=1)
    # Method 2: fixed COP, capacity-capped
    # Preallocation
    T = promote_type(Float64, eltype(Qb), typeof(COP_heating), typeof(COP_cooling),
        typeof(CAP_heating_unit), typeof(CAP_cooling_unit), typeof(n_HP))
    Qg = similar(Qb, T)
    Qb_unmet = similar(Qb, T)
    Qb_hp = similar(Qb, T)

    # Conversion
    for k in eachindex(Qb)
        Qg[k], Qb_unmet[k], Qb_hp[k] = ground_load(Qb[k], COP_heating, COP_cooling,
            CAP_heating_unit, CAP_cooling_unit; n_HP)
    end

    return Qg, Qb_unmet, Qb_hp
end

"""
    _ground_load_at_performance(Qb, CAP, COP)

Shared by both `ground_load` performance-curve methods (method 3) once CAP/COP have been looked up
at the current EFT(s): caps the heat pump's contribution at `CAP`, then hands off to method 1 for
the arithmetic (passing the same `COP` in both slots is safe since `Qb_hp`'s sign always matches
`Qb`'s, so only one of method 1's branches ever fires).
"""
function _ground_load_at_performance(Qb::Real, CAP::Real, COP::Real)
    # Cap the heat pump's contribution at its current CAP
    Qb_hp = Qb >= 0 ? min(Qb, CAP) : max(Qb, -CAP)
    # Compute the ground load and unmet building load
    Qg = ground_load(Qb_hp, COP, COP)
    Qb_unmet = Qb - Qb_hp
    return Qg, Qb_unmet, Qb_hp
end

function ground_load(Qb::Real, Tsrc::Real, perf_heating::Tuple, perf_cooling::Tuple; n_HP::Real=1)
    # Method 3: performance-curve-based, 1D overload
    CAP, COP = Qb >= 0 ? heat_pump_performance(Tsrc, perf_heating...) :
               heat_pump_performance(Tsrc, perf_cooling...)
    return _ground_load_at_performance(Qb, n_HP * CAP, COP)
end
function ground_load(Qb::AbstractVector{<:Real}, Tsrc::AbstractVector{<:Real},
    perf_heating::Tuple, perf_cooling::Tuple; n_HP::Real=1)
    # Method 3: performance-curve-based, 1D overload
    # Check
    length(Qb) == length(Tsrc) || throw(DimensionMismatch("Qb and Tsrc must have the same length."))
    T = promote_type(Float64, eltype(Qb), eltype(Tsrc))

    # Preallocation
    Qg = similar(Qb, T)
    Qb_unmet = similar(Qb, T)
    Qb_hp = similar(Qb, T)

    # Conversion
    for k in eachindex(Qb, Tsrc)
        Qg[k], Qb_unmet[k], Qb_hp[k] = ground_load(Qb[k], Tsrc[k], perf_heating, perf_cooling; n_HP)
    end

    return Qg, Qb_unmet, Qb_hp
end
function ground_load(Qb::Real, Tsrc::Real, Tload::Real, perf_heating::Tuple, perf_cooling::Tuple;
    n_HP::Real=1)
    # Method 3: performance-curve-based, 2D overload
    CAP, COP = Qb >= 0 ? heat_pump_performance(Tsrc, Tload, perf_heating...) :
               heat_pump_performance(Tsrc, Tload, perf_cooling...)
    return _ground_load_at_performance(Qb, n_HP * CAP, COP)
end
function ground_load(Qb::AbstractVector{<:Real}, Tsrc::AbstractVector{<:Real},
    Tload::AbstractVector{<:Real}, perf_heating::Tuple, perf_cooling::Tuple; n_HP::Real=1)
    # Method 3: performance-curve-based, 2D overload
    # Check
    length(Qb) == length(Tsrc) == length(Tload) ||
        throw(DimensionMismatch("Qb, Tsrc and Tload must have the same length."))
    T = promote_type(Float64, eltype(Qb), eltype(Tsrc), eltype(Tload))

    # Preallocation
    Qg = similar(Qb, T)
    Qb_unmet = similar(Qb, T)
    Qb_hp = similar(Qb, T)

    # Conversion
    for k in eachindex(Qb, Tsrc, Tload)
        Qg[k], Qb_unmet[k], Qb_hp[k] = ground_load(Qb[k], Tsrc[k], Tload[k], perf_heating,
            perf_cooling; n_HP)
    end

    return Qg, Qb_unmet, Qb_hp
end

"""
    n_active_heat_pumps(Qb_hp, n_HP, CAP_heating, CAP_cooling, COP_heating, COP_cooling)
    n_active_heat_pumps(Qb_hp, Tsrc, n_HP, perf_heating, perf_cooling)
    n_active_heat_pumps(Qb_hp, Tsrc, Tload, n_HP, perf_heating, perf_cooling)

Number of heat pumps operating at each time step to deliver `Qb_hp`, given a single unit's fixed
capacity and COP, or its performance curve evaluated at the source EFT `Tsrc` (1D) or at `Tsrc` and
the load EFT `Tload` (2D) (e.g. from a ground temperature simulation). Following
Pasquier et al. (2013), this is a fractional, time-averaged count (between `0` and `n_HP`) rather
than an integer, since a heat pump typically cycles within a time step instead of switching exactly
on its boundary.
# Arguments
    - `Qb_hp`: Building load actually delivered by the heat pumps, post-capacity-capping (the
        third output of [`ground_load`](@ref), methods 2-3) [W]
    - `Tsrc`, `Tload`: Source/load entering fluid temperature (EFT) [°C]
    - `n_HP`: Number of installed heat pumps
    - `CAP_heating`, `CAP_cooling`: Fixed heating/cooling capacity of one heat pump unit [W]
    - `COP_heating`, `COP_cooling`: Fixed heating/cooling COP [-]
    - `perf_heating`: Single heat pump's heating performance data, forwarded to
        [`heat_pump_performance`](@ref) (`(Tsrc_grid, CAP_vec, COP_vec)` or
        `(Tsrc_grid, Tload_grid, CAP_table, COP_table)`)
    - `perf_cooling`: Single heat pump's cooling performance data, forwarded to
        [`heat_pump_performance`](@ref)
# Output
    - `n`: Fractional number of heat pumps operating, in `[0, n_HP]`
    - `COP`: Full-load coefficient of performance at each time step (the fixed value, or the
        same interpolation as `n`); see [`part_load_COP`](@ref) for the part-load degradation
# Reference
    - Pasquier, P., Marcotte, D., Bernier, M., & Kummert, M. (2013). Simulation of ground-coupled
        heat pump systems using a spectral approach. Proceedings of BS2013, 2691–2698.
"""
function n_active_heat_pumps(Qb_hp::AbstractVector{<:Real}, Tsrc::AbstractVector{<:Real},
    n_HP::Real, perf_heating::Tuple, perf_cooling::Tuple)

    # Check
    length(Qb_hp) == length(Tsrc) ||
        throw(DimensionMismatch("Qb_hp and Tsrc must have the same length."))

    # Preallocation
    n = similar(Qb_hp, Float64)
    COP = similar(Qb_hp, Float64)

    # Conversion between building load and number of active heat pumps
    for k in eachindex(Qb_hp, Tsrc)
        CAPk, COPk = Qb_hp[k] >= 0 ? heat_pump_performance(Tsrc[k], perf_heating...) :
                     heat_pump_performance(Tsrc[k], perf_cooling...)
        COP[k] = COPk
        n[k] = CAPk > 0 ? clamp(abs(Qb_hp[k]) / CAPk, 0.0, n_HP) : 0.0
    end
    return n, COP
end
function n_active_heat_pumps(Qb_hp::AbstractVector{<:Real}, Tsrc::AbstractVector{<:Real},
    Tload::AbstractVector{<:Real}, n_HP::Real, perf_heating::Tuple, perf_cooling::Tuple)

    # Check
    length(Qb_hp) == length(Tsrc) == length(Tload) ||
        throw(DimensionMismatch("Qb_hp, Tsrc and Tload must have the same length."))

    # Preallocation
    n = similar(Qb_hp, Float64)
    COP = similar(Qb_hp, Float64)

    # Conversion between building load and number of active heat pumps
    for k in eachindex(Qb_hp, Tsrc, Tload)
        CAPk, COPk = Qb_hp[k] >= 0 ? heat_pump_performance(Tsrc[k], Tload[k], perf_heating...) :
                     heat_pump_performance(Tsrc[k], Tload[k], perf_cooling...)
        COP[k] = COPk
        n[k] = CAPk > 0 ? clamp(abs(Qb_hp[k]) / CAPk, 0.0, n_HP) : 0.0
    end
    return n, COP
end
function n_active_heat_pumps(Qb_hp::AbstractVector{<:Real}, n_HP::Real, CAP_heating::Real,
    CAP_cooling::Real, COP_heating::Real, COP_cooling::Real)

    # Preallocation
    n = similar(Qb_hp, Float64)
    COP = similar(Qb_hp, Float64)

    # Conversion between building load and number of active heat pumps
    for k in eachindex(Qb_hp)
        CAPk, COP[k] = Qb_hp[k] >= 0 ? (CAP_heating, COP_heating) : (CAP_cooling, COP_cooling)
        n[k] = CAPk > 0 ? clamp(abs(Qb_hp[k]) / CAPk, 0.0, n_HP) : 0.0
    end
    return n, COP
end
