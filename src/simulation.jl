# Simulation cases, in increasing order of complexity. Only the first two are implemented; see
# each function's docstring and the package tutorial for what each one covers.
#   1. Fixed COP/CAP, stationary flow, no source-side shaving      -> gshp_simulation (Real, Real)
#   2. Variable (temperature-dependent) COP/CAP, stationary flow,
#      no source-side shaving                                      -> gshp_simulation (perf tuples)
#   3. (future) Either COP style, stationary flow, WITH source-side shaving: caps `Qg` further so
#      `T_out` stays within operating limits, following the correction scheme of Pasquier et al.
#      (2013) — needs a couple of extra passes of `fluid_temperature`, still no full outer
#      iteration.
#   4. (future) Either COP style, non-stationary (time-varying) flow rate, with or without
#      shaving: swap the stationary `ground_response`/`fluid_temperature` call for
#      `GroundHeatExchanger.jl`'s `convolution_ns`/`impulse_func_ns`/state-vector machinery, with a
#      flow-rate control law (e.g. a `V(EWT)` map, as in `fNSCiter.m`).

using GroundHeatExchanger: resistance_ULoop_effective, ground_response, FLSModel,
    fluid_temperature, inlet_temperature, outlet_temperature

"""
    _stationary_ground_temperature(t, Qg, xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf,
        V, T0)

Ground temperature entering (`T_in`) and leaving (`T_out`) a borehole field for a constant flow
rate `V`, by stationary (FFT) convolution of the ground load `Qg` with the field's g-function. See
the README's "Conventions" section for the `T_in`/`T_out` nomenclature. Shared by every
[`gshp_simulation`](@ref) method — internal helper, not exported.
"""
function _stationary_ground_temperature(t::AbstractVector{<:Real}, Qg::AbstractVector{<:Real},
    xy::AbstractMatrix{<:Real}, rb::Real, H::Real, D::Real, ks::Real, Cs::Real, s::Real, ro::Real,
    ri::Real, kg::Real, kp::Real, kf::Real, cf::Real, ρf::Real, μf::Real, V::Real, T0::Real)
    L = size(xy, 1) * H
    q = Qg ./ L
    Rbe = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
    g = ground_response(t, rb, xy, FLSModel(H, D, ks, Cs))
    Tf = fluid_temperature(t, q, g, T0, Rbe)
    C = cf * ρf
    T_in = inlet_temperature(Tf, q, H, V, C)
    T_out = outlet_temperature(Tf, q, H, V, C)
    return T_in, T_out, Rbe, g
end

"""
    gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating_unit, CAP_cooling_unit, n_HP,
        xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0)

Case 1, the simplest GSHP simulation: building load to ground load with **fixed**
(temperature-independent) COP and CAP, a **stationary** (constant flow rate) ground response, and
**no source-side load shaving** — any load beyond the installed capacity is simply unmet. See
[`ground_load_capped`](@ref) for the load conversion step and
[`_stationary_ground_temperature`](@ref) for the ground response, and the README's "Conventions"
section for the sign convention and `T_in`/`T_out` nomenclature.

# Arguments
    - `t`: Time vector [s], uniformly spaced, matching `Qb`
    - `Qb`: Building thermal load [W], `> 0` heating, `< 0` cooling, same length as `t`
    - `COP_heating`, `COP_cooling`: Fixed nominal heating/cooling COP of one heat pump unit
    - `CAP_heating_unit`, `CAP_cooling_unit`: Fixed heating/cooling capacity of one heat pump unit
        [W]
    - `n_HP`: Number of installed heat pumps
    - `xy`: Borehole coordinates (nb x 2) [m]
    - `rb`: Borehole radius [m]
    - `H`, `D`: Borehole depth and buried depth [m]
    - `ks`, `Cs`: Ground thermal conductivity [W/mK] and volumetric heat capacity [J/m3K]
    - `s`: Shank spacing [m]
    - `ro`, `ri`: Pipe outer/inner radius [m]
    - `kg`, `kp`: Grout and pipe thermal conductivity [W/mK]
    - `kf`, `cf`, `ρf`, `μf`: Fluid thermal conductivity [W/mK], specific heat [J/kgK], density
        [kg/m3] and dynamic viscosity [kg/(m.s)]
    - `V`: Flow rate per borehole [m3/s]
    - `T0`: Undisturbed ground temperature [°C]

# Output
    Named tuple with fields:
    - `Qg`: Ground thermal load [W]
    - `Qb_hp`: Building load actually delivered by the heat pumps [W]
    - `Qb_unmet`: Building load left unmet (auxiliary) [W]
    - `T_in`, `T_out`: Ground temperature entering/leaving the field [°C]
    - `n_active`: Fractional number of heat pumps operating, in `[0, n_HP]`
    - `Rbe`: Effective borehole thermal resistance [m.K/W]
    - `g`: g-function used for the ground response [°C.m/W]
"""
function gshp_simulation(t::AbstractVector{<:Real}, Qb::AbstractVector{<:Real},
    COP_heating::Real, COP_cooling::Real, CAP_heating_unit::Real, CAP_cooling_unit::Real,
    n_HP::Real, xy::AbstractMatrix{<:Real}, rb::Real, H::Real, D::Real, ks::Real, Cs::Real,
    s::Real, ro::Real, ri::Real, kg::Real, kp::Real, kf::Real, cf::Real, ρf::Real, μf::Real,
    V::Real, T0::Real)

    CAP_heating_total = n_HP * CAP_heating_unit
    CAP_cooling_total = n_HP * CAP_cooling_unit
    Qg, Qb_unmet = ground_load_capped(Qb, COP_heating, COP_cooling, CAP_heating_total,
        CAP_cooling_total)
    Qb_hp = Qb .- Qb_unmet

    T_in, T_out, Rbe, g = _stationary_ground_temperature(t, Qg, xy, rb, H, D, ks, Cs, s, ro, ri,
        kg, kp, kf, cf, ρf, μf, V, T0)

    n_active = similar(Qb_hp, Float64)
    for k in eachindex(Qb_hp)
        CAP_unit = Qb_hp[k] >= 0 ? CAP_heating_unit : CAP_cooling_unit
        n_active[k] = CAP_unit > 0 ? clamp(abs(Qb_hp[k]) / CAP_unit, 0.0, n_HP) : 0.0
    end

    return (; Qg, Qb_hp, Qb_unmet, T_in, T_out, n_active, Rbe, g)
end

"""
    gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP,
        xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0)

Case 2: building load to ground load with **temperature-dependent** COP/CAP
([`heat_pump_performance`](@ref) curves), a **stationary** (constant flow rate) ground response,
and **no source-side load shaving**. The building-to-ground conversion still can't use the curve
directly (the source EFT isn't known until the ground response has been simulated), so it uses the
mean of the curve as a fixed nominal COP for that one pass ([`ground_load_capped`](@ref)); once the
resulting `T_out` is known, [`n_active_heat_pumps`](@ref) reports the number of units operating and
the COP actually achieved — a diagnostic, not fed back into the simulation. See the README's
"Conventions" section for the sign convention and `T_in`/`T_out` nomenclature.

# Arguments
    - `t`: Time vector [s], uniformly spaced, matching `Qb`
    - `Qb`: Building thermal load [W], `> 0` heating, `< 0` cooling, same length as `t`
    - `perf_heating`, `perf_cooling`: Single heat pump's heating/cooling performance data,
        forwarded to [`heat_pump_performance`](@ref) (`(Tsrc_grid, CAP_vec, COP_vec)`)
    - `n_HP`: Number of installed heat pumps
    - remaining arguments: as in the fixed-COP [`gshp_simulation`](@ref) method

# Output
    Named tuple with fields:
    - `Qg`: Ground thermal load [W]
    - `Qb_hp`: Building load actually delivered by the heat pumps [W]
    - `Qb_unmet`: Building load left unmet (auxiliary) [W]
    - `T_in`, `T_out`: Ground temperature entering/leaving the field [°C]
    - `n_active`: Fractional number of heat pumps operating, in `[0, n_HP]`
    - `COP_actual`: COP actually achieved at `T_out`, from the same interpolation as `n_active`
    - `Rbe`: Effective borehole thermal resistance [m.K/W]
    - `g`: g-function used for the ground response [°C.m/W]
# Reference
    - Pasquier, P., Marcotte, D., Bernier, M., & Kummert, M. (2013). Simulation of ground-coupled
        heat pump systems using a spectral approach. Proceedings of BS2013, 2691–2698.
"""
function gshp_simulation(t::AbstractVector{<:Real}, Qb::AbstractVector{<:Real},
    perf_heating, perf_cooling, n_HP::Real, xy::AbstractMatrix{<:Real}, rb::Real, H::Real,
    D::Real, ks::Real, Cs::Real, s::Real, ro::Real, ri::Real, kg::Real, kp::Real, kf::Real,
    cf::Real, ρf::Real, μf::Real, V::Real, T0::Real)

    CAP_heat_vec, COP_heat_vec = perf_heating[2], perf_heating[3]
    CAP_cool_vec, COP_cool_vec = perf_cooling[2], perf_cooling[3]
    COP_heating_nominal = sum(COP_heat_vec) / length(COP_heat_vec)
    COP_cooling_nominal = sum(COP_cool_vec) / length(COP_cool_vec)
    CAP_heating_total = n_HP * maximum(CAP_heat_vec)
    CAP_cooling_total = n_HP * maximum(CAP_cool_vec)

    Qg, Qb_unmet = ground_load_capped(Qb, COP_heating_nominal, COP_cooling_nominal,
        CAP_heating_total, CAP_cooling_total)
    Qb_hp = Qb .- Qb_unmet

    T_in, T_out, Rbe, g = _stationary_ground_temperature(t, Qg, xy, rb, H, D, ks, Cs, s, ro, ri,
        kg, kp, kf, cf, ρf, μf, V, T0)

    n_active, COP_actual = n_active_heat_pumps(Qb_hp, T_out, n_HP, perf_heating, perf_cooling)

    return (; Qg, Qb_hp, Qb_unmet, T_in, T_out, n_active, COP_actual, Rbe, g)
end
