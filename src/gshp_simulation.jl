# Simulation cases, in increasing order of complexity. See each function's docstring and the package
# tutorial for what each one covers.
#   1. Fixed COP/CAP, stationary flow, no source-side shaving, but building side Qb capping
#   2. Variable (temperature-dependent) COP/CAP, stationary flow, no source-side shaving, but
#       building side Qb capping with a small fixed-point loop over the source EFT-dependent COP,
#       typically converging in well under `max_iter` iterations.
#   3. Fixed COP/CAP, non-stationary (time-varying) flow rate, without shaving: the stationary
#      `ground_response`/`fluid_temperature` call is replaced by `GroundHeatExchanger.jl`'s
#      `convolution_ns`/state-vector machinery, with a load-based flow-rate control law. Picks the
#      active flow rate from `Qg` alone, so this stays a single, non-iterative pass (see
#      `_flow_rate_states`). (future: same, but variable COP/CAP) #TODO Shouldn't the control be based on the source-side loads? I don't understand how the controls are made from the explanation.
#   4. (future) Either COP style, stationary flow, WITH source-side shaving: caps `Qg` further so
#      `T_out` stays within operating limits, following the correction scheme of Pasquier et al.
#      (2013).
#   5. (future) All together: Either COP style, non-stationary flow depending on Qg to limit pumping
#       energy, source-side shaving capping.
#
# Each case has two `gshp_simulation` methods: one takes an already-computed g-function `g` and 
# effective borehole resistance `Rbe`, the other takes the field's raw geometry/soil/pipe
# inputs, computes `g`/`Rbe` once, and delegates to the first.

using GroundHeatExchanger: resistance_ULoop_effective, ground_response, FLSModel,
    fluid_temperature, inlet_temperature, outlet_temperature, convolution_ns

"""
    _ground_temperature(t, Qg, nb, H, g, Rbe, cf, ρf, V, T0)

Ground temperature entering (`T_in`) and leaving (`T_out`) a borehole field for a constant flow
rate `V`, from an already-computed g-function `g` and effective borehole resistance `Rbe`, by
stationary (FFT) convolution of the ground load `Qg` with `g`. See the README's "Conventions"
section for the `T_in`/`T_out` nomenclature. Shared by every stationary [`gshp_simulation`](@ref)
method.
"""
function _ground_temperature(t::AbstractVector{<:Real}, Qg::AbstractVector{<:Real}, nb::Integer,
    H::Real, g::AbstractVector{<:Real}, Rbe::Real, cf::Real, ρf::Real, V::Real, T0::Real)
    # Define inputs
    L = nb * H
    q = Qg ./ L
    # Compute the mean fluid temperature and the inlet/outlet split
    Tf = fluid_temperature(t, q, g, T0, Rbe)
    C = cf * ρf # TODO Should this be passed in instead of recomputed here? It's already computed in the calling function, and it's a constant for a given fluid.
    T_in = inlet_temperature(Tf, q, H, V, C)
    T_out = outlet_temperature(Tf, q, H, V, C)
    return T_in, T_out
end

"""
    gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating, CAP_cooling, n_HP,
        nb, H, g, Rbe, cf, ρf, V, T0)
    gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating, CAP_cooling, n_HP,
        xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0)

Case 1, the simplest GSHP simulation: building load to ground load with **fixed**
(temperature-independent) COP and CAP, a **stationary** (constant flow rate) ground response, and
**no source-side load shaving**. Any load on the building side beyond the installed capacity is
simply unmet. See [`ground_load`](@ref) for the load conversion step and
[`_ground_temperature`](@ref) for the ground response, and the README's "Conventions" section for
the sign convention and `T_in`/`T_out` nomenclature.

The first method takes an already-computed g-function `g` and effective borehole resistance `Rbe`.
The second method takes the field's raw geometry/soil/pipe inputs, computes `g` and `Rbe` once via
[`ground_response`](@ref)/[`resistance_ULoop_effective`](@ref), and calls the first method.
# Arguments
    - `t`: Time vector [s], uniformly spaced, matching `Qb`
    - `Qb`: Building thermal load [W], `> 0` heating, `< 0` cooling, same length as `t`
    - `COP_heating`, `COP_cooling`: Fixed nominal heating/cooling COP of one heat pump unit
    - `CAP_heating`, `CAP_cooling`: Fixed heating/cooling capacity of one heat pump unit [W]
    - `n_HP`: Number of installed heat pumps [-]
    - `nb`: Number of boreholes [-] (first method only; `size(xy, 1)` in the second)
    - `g`: Pre-computed g-function [°C.m/W], e.g. from [`ground_response`](@ref) (first method only)
    - `Rbe`: Pre-computed effective borehole thermal resistance [m.K/W], e.g. from
        [`resistance_ULoop_effective`](@ref) (first method only)
    - `xy`: Borehole coordinates (nb x 2) [m] (second method only)
    - `rb`: Borehole radius [m] (second method only)
    - `H`, `D`: Borehole depth and buried depth [m] (`D` second method only)
    - `ks`, `Cs`: Ground thermal conductivity [W/mK] and volumetric heat capacity [J/m3K] (second
        method only)
    - `s`: Shank spacing [m] (second method only)
    - `ro`, `ri`: Pipe outer/inner radius [m] (second method only)
    - `kg`, `kp`: Grout and pipe thermal conductivity [W/mK] (second method only)
    - `kf`, `cf`, `ρf`, `μf`: Fluid thermal conductivity [W/mK], specific heat [J/kgK], density
        [kg/m3] and dynamic viscosity [kg/(m.s)] (`kf`, `μf` second method only, from
        [`fluid_property`](@ref))
    - `V`: Flow rate per borehole [m3/s]
    - `T0`: Undisturbed ground temperature [°C]
# Output
    Named tuple with fields:
    - `Qg`: Ground thermal load [W]
    - `Qb_hp`: Building load actually delivered by the heat pumps (Qb = Qb_hp+Qb_unmet) [W]
    - `Qb_unmet`: Building load left unmet (auxiliary) [W]
    - `T_in`, `T_out`: Ground temperature entering/leaving the field [°C]
    - `n_active`: Fractional number of heat pumps operating, in `[0, n_HP]`
    - `Rbe`: Effective borehole thermal resistance [m.K/W]
    - `g`: g-function used for the ground response [°C.m/W]
"""
function gshp_simulation(t::AbstractVector{<:Real}, Qb::AbstractVector{<:Real},
    COP_heating::Real, COP_cooling::Real, CAP_heating::Real, CAP_cooling::Real,
    n_HP::Real, nb::Integer, H::Real, g::AbstractVector{<:Real}, Rbe::Real, cf::Real, ρf::Real,
    V::Real, T0::Real)

    #TODO Here, the fluid properties should be calculated as a function of the temperature and inputs.
    # So, cf and ρf should be calculated using fluid_property(). The simulation inputs should instead be the variables for the fluid_property calculation.

    # Thermal load conversion
    Qg, Qb_unmet, Qb_hp = ground_load(Qb, COP_heating, COP_cooling, CAP_heating, CAP_cooling; n_HP)

    # Ground loop temperature evaluation
    T_in, T_out = _ground_temperature(t, Qg, nb, H, g, Rbe, cf, ρf, V, T0)

    # Fractional number of heat pumps operating: the fixed-CAP special case of
    # `n_active_heat_pumps` (no performance curve to interpolate CAP from here, just the two fixed
    # scalars `CAP_heating`/`CAP_cooling`).
    n_active = similar(Qb_hp, Float64)
    for k in eachindex(Qb_hp)
        CAP = Qb_hp[k] >= 0 ? CAP_heating : CAP_cooling
        n_active[k] = CAP > 0 ? clamp(abs(Qb_hp[k]) / CAP, 0.0, n_HP) : 0.0
    end #TODO here, why not simply use the n_active_heat_pumps function? It does the same thing, but with more generality. A user not wanting to use CAP could just set it very large. 

    return (; Qg, Qb_hp, Qb_unmet, T_in, T_out, n_active, Rbe, g)
end
function gshp_simulation(t::AbstractVector{<:Real}, Qb::AbstractVector{<:Real},
    COP_heating::Real, COP_cooling::Real, CAP_heating::Real, CAP_cooling::Real,
    n_HP::Real, xy::AbstractMatrix{<:Real}, rb::Real, H::Real, D::Real, ks::Real, Cs::Real,
    s::Real, ro::Real, ri::Real, kg::Real, kp::Real, kf::Real, cf::Real, ρf::Real, μf::Real,
    V::Real, T0::Real)

    Rbe = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
    g = ground_response(t, rb, xy, FLSModel(H, D, ks, Cs))

    return gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating, CAP_cooling,
        n_HP, size(xy, 1), H, g, Rbe, cf, ρf, V, T0)
end

"""
    gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP, nb, H, g, Rbe, cf, ρf, V, T0;
        tol=0.05, max_iter=5)
    gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP,
        xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0; tol=0.05, max_iter=5)

Case 2, building load to ground load with **temperature-dependent** COP/CAP
([`heat_pump_performance`](@ref) curves), a **stationary** (constant flow rate) ground response,
and **no source-side load shaving**. The building-to-ground conversion can't use the curve on its
first pass (the source EFT isn't known until the ground response has been simulated), so it starts
from the mean of the curve as a nominal COP ([`ground_load`](@ref)) to get a first `T_out` estimate;
it then **iterates**, alternating [`ground_load`](@ref)'s curve-based method (using the previous
iteration's `T_out` as the source EFT) with [`_ground_temperature`](@ref), until `T_out` stops
moving by more than `tol` or `max_iter` is reached — a small fixed-point (Picard) loop, not a full
outer iteration over the whole simulation pipeline. This mirrors `fNSCiter.m`'s `while` loop, which
recomputes COP/CAP from the current source EFT every pass and is itself capped at 5 iterations; it
typically converges in well under that. [`n_active_heat_pumps`](@ref) reports the number of units
operating and the COP actually achieved at the converged `T_out`. See the README's "Conventions"
section for the sign convention and `T_in`/`T_out` nomenclature.

Two methods, in increasing order of complexity, exactly as the fixed-COP [`gshp_simulation`](@ref)
(case 1): the first takes precomputed `g`/`Rbe` (reused, unchanged, across every iteration of the
loop below); the second takes the field's raw geometry/soil/pipe inputs, computes them once, and
calls the first.
# Arguments
    - `t`: Time vector [s], uniformly spaced, matching `Qb`
    - `Qb`: Building thermal load [W], `> 0` heating, `< 0` cooling, same length as `t`
    - `perf_heating`, `perf_cooling`: Single heat pump's heating/cooling performance data,
        forwarded to [`heat_pump_performance`](@ref) (`(Tsrc_grid, CAP_vec, COP_vec)`)
    - `n_HP`: Number of installed heat pumps
    - `tol`: Convergence tolerance on `T_out` between iterations [°C]
    - `max_iter`: Maximum number of fixed-point iterations (safety cap; the loop typically stops
        well before this)
    - remaining arguments: as in the fixed-COP [`gshp_simulation`](@ref) method

# Output
    Named tuple with fields:
    - `Qg`: Ground thermal load [W]
    - `Qb_hp`: Building load actually delivered by the heat pumps [W]
    - `Qb_unmet`: Building load left unmet (auxiliary) [W]
    - `T_in`, `T_out`: Ground temperature entering/leaving the field, at convergence [°C]
    - `n_active`: Fractional number of heat pumps operating, in `[0, n_HP]`
    - `COP_actual`: COP actually achieved at the converged `T_out`, from the same interpolation as
        `n_active`
    - `n_iter`: Number of fixed-point iterations actually taken [-]
    - `Rbe`: Effective borehole thermal resistance [m.K/W]
    - `g`: g-function used for the ground response [°C.m/W]
# Reference
    - Pasquier, P., Marcotte, D., Bernier, M., & Kummert, M. (2013). Simulation of ground-coupled
        heat pump systems using a spectral approach. Proceedings of BS2013, 2691–2698.
"""
function gshp_simulation(t::AbstractVector{<:Real}, Qb::AbstractVector{<:Real},
    perf_heating, perf_cooling, n_HP::Real, nb::Integer, H::Real, g::AbstractVector{<:Real},
    Rbe::Real, cf::Real, ρf::Real, V::Real, T0::Real; tol::Real=0.05, max_iter::Integer=5)

    CAP_heat_vec, COP_heat_vec = perf_heating[2], perf_heating[3]
    CAP_cool_vec, COP_cool_vec = perf_cooling[2], perf_cooling[3]
    COP_heating_nominal = sum(COP_heat_vec) / length(COP_heat_vec)
    COP_cooling_nominal = sum(COP_cool_vec) / length(COP_cool_vec)

    # Pass 0: nominal COP (the source EFT isn't known yet) to get a first T_out estimate
    Qg, Qb_unmet = ground_load(Qb, COP_heating_nominal, COP_cooling_nominal, maximum(CAP_heat_vec),
        maximum(CAP_cool_vec); n_HP)
    Qb_hp = Qb .- Qb_unmet
    T_in, T_out = _ground_temperature(t, Qg, nb, H, g, Rbe, cf, ρf, V, T0)

    # Fixed-point loop: re-derive Qg/T_out from the actual, curve-based COP at the previous
    # iteration's T_out, until T_out stops moving (or max_iter is reached, as a safety net).
    n_iter = 0
    for _ in 1:max_iter
        n_iter += 1
        T_out_prev = T_out
        Qg, Qb_unmet = ground_load(Qb, T_out, perf_heating, perf_cooling; n_HP)
        Qb_hp = Qb .- Qb_unmet
        T_in, T_out = _ground_temperature(t, Qg, nb, H, g, Rbe, cf, ρf, V, T0)
        maximum(abs.(T_out .- T_out_prev)) < tol && break
    end

    n_active, COP_actual = n_active_heat_pumps(Qb_hp, T_out, n_HP, perf_heating, perf_cooling)

    return (; Qg, Qb_hp, Qb_unmet, T_in, T_out, n_active, COP_actual, n_iter, Rbe, g)
end
function gshp_simulation(t::AbstractVector{<:Real}, Qb::AbstractVector{<:Real},
    perf_heating, perf_cooling, n_HP::Real, xy::AbstractMatrix{<:Real}, rb::Real, H::Real,
    D::Real, ks::Real, Cs::Real, s::Real, ro::Real, ri::Real, kg::Real, kp::Real, kf::Real,
    cf::Real, ρf::Real, μf::Real, V::Real, T0::Real; tol::Real=0.05, max_iter::Integer=5)

    Rbe = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
    g = ground_response(t, rb, xy, FLSModel(H, D, ks, Cs))

    return gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP, size(xy, 1), H, g, Rbe, cf, ρf,
        V, T0; tol, max_iter)
end

"""
    _flow_rate_states(Qg, nb, V_candidates, ΔT_target, C)

Discrete flow-rate schedule for the non-stationary case: at each time step, the per-borehole flow
rate is picked from `V_candidates` — the controllable set points a real pumping control would
switch between (e.g. `[25, 30, 35] ./ 60000` for 25/30/35 L/min) — as the one closest to the flow
that would hold the field's inlet/outlet split at `ΔT_target` for the instantaneous ground load
`Qg`. This mirrors the flow-rate map (`V.carte`) of `fNSCiter.m`, but skips its continuous
map-then-round step: with only a handful of candidates to choose from, comparing the ideal
continuous flow rate directly against each candidate is simpler and exact. The selection depends on
`Qg` alone (never on a simulated temperature), so it can be computed ahead of, and independently
from, the ground response — this is what keeps the non-stationary
[`gshp_simulation`](@ref) methods a single, non-iterative pass.
# Arguments
    - `Qg`: Ground thermal load [W]
    - `nb`: Number of boreholes [-]
    - `V_candidates`: Controllable per-borehole flow rate set points [m3/s]
    - `ΔT_target`: Target field inlet-outlet temperature difference used to pick the active flow
        rate [K]
    - `C`: Fluid volumetric specific heat [J/m3K] (`cf * ρf`)
# Output
    - `V_t`: Per-borehole flow rate active at each time step, one of `V_candidates` [m3/s]
    - `state_vec`: Index into sorted `V_candidates` active at each time step [-]
"""
function _flow_rate_states(Qg::AbstractVector{<:Real}, nb::Integer,
    V_candidates::AbstractVector{<:Real}, ΔT_target::Real, C::Real)

    V_sorted = sort(V_candidates)
    V_ideal = abs.(Qg) ./ (nb * ΔT_target * C)

    state_vec = similar(Qg, Int)
    V_t = similar(Qg, Float64)
    for k in eachindex(Qg)
        _, i = findmin(abs.(V_sorted .- V_ideal[k]))
        state_vec[k] = i
        V_t[k] = V_sorted[i]
    end

    return V_t, state_vec
end

"""
    _nonstationary_ground_temperature(t, Qg, xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf,
        V_candidates, ΔT_target, T0)

Ground temperature entering (`T_in`) and leaving (`T_out`) a borehole field for a **non-stationary**
(time-varying) flow rate, switching between `V_candidates` under the load-based control law of
[`_flow_rate_states`](@ref). The g-function (pure conduction) does not depend on flow rate and is
shared by every state; only the effective borehole resistance `Rbₑ` — and, downstream, the
inlet/outlet split — changes with the active flow rate. Following the combined transfer function
`h_s(t) = Rbₑ(V_s) + g(t)` per state, the mean fluid temperature is a single non-stationary
convolution (`GroundHeatExchanger.jl`'s `convolution_ns`); this reduces exactly to
[`_ground_temperature`](@ref) (given the single resulting `g`/`Rbe`) when `V_candidates` has one
element. See the README's "Conventions" section for the `T_in`/`T_out` nomenclature. Shared by every
non-stationary [`gshp_simulation`](@ref) method.
# Arguments
    - as the fixed-COP [`gshp_simulation`](@ref) (case 1) raw-input method, except:
    - `V_candidates`: Controllable per-borehole flow rate set points [m3/s]
    - `ΔT_target`: Target field inlet-outlet temperature difference used to pick the active flow
        rate at each time step [K]
# Output
    - `T_in`, `T_out`: Ground temperature entering/leaving the field [°C]
    - `V_t`: Per-borehole flow rate active at each time step [m3/s]
    - `Rbe`: Effective borehole thermal resistance for each flow-rate state, in `sort(V_candidates)`
        order [m.K/W]
    - `g`: g-function used for the ground response, shared by every flow state [°C.m/W]
# Reference
    - Beaudry, G., Pasquier, P., & Nguyen, A. (2024). New formulations and experimental validation
        of non-stationary convolutions for the fast simulation of time-variant flowrates in ground
        heat exchangers. Science and Technology for the Built Environment, 30(3), 208–219.
"""
function _nonstationary_ground_temperature(t::AbstractVector{<:Real}, Qg::AbstractVector{<:Real},
    xy::AbstractMatrix{<:Real}, rb::Real, H::Real, D::Real, ks::Real, Cs::Real, s::Real, ro::Real,
    ri::Real, kg::Real, kp::Real, kf::Real, cf::Real, ρf::Real, μf::Real,
    V_candidates::AbstractVector{<:Real}, ΔT_target::Real, T0::Real)

    nb = size(xy, 1)
    L = nb * H
    q = Qg ./ L
    C = cf * ρf

    V_t, state_vec = _flow_rate_states(Qg, nb, V_candidates, ΔT_target, C)

    g = ground_response(t, rb, xy, FLSModel(H, D, ks, Cs))
    Rbe = [resistance_ULoop_effective(Vs, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
           for Vs in sort(V_candidates)]
    h = hcat([Rbe[i] .+ g for i in eachindex(Rbe)]...)

    Tf = T0 .+ convolution_ns(q, h, state_vec)
    T_in = Tf .+ (q .* H) ./ (2 .* V_t .* C)
    T_out = Tf .- (q .* H) ./ (2 .* V_t .* C)

    return T_in, T_out, V_t, Rbe, g
end

"""
    gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating_unit, CAP_cooling_unit, n_HP,
        xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V_candidates, ΔT_target, T0)

Case 3, fixed-COP GSHP simulation with a **non-stationary** (time-varying) flow rate: the same
building load conversion and fixed COP/CAP as case 1, and still **no source-side load shaving**, but
the ground response now switches between a small set of controllable per-borehole flow rates
`V_candidates` instead of a single constant `V`. The active flow rate at each time step is picked
from `Qg` alone by [`_flow_rate_states`](@ref) — never from the simulated temperature — so this
remains a single, non-iterative pass: a larger flow rate always cuts the `T_in`/`T_out` split for the
same `Qg` (it's proportional to `Qg/V`), which is the guaranteed effect of switching to a higher flow
step. Whether that also lowers the peak `T_out` at a given extreme depends on the load's sign and
size there too (a smaller split pulls `T_out` toward the mean fluid temperature, while the smaller
`Rbₑ` that comes with a bigger flow rate pulls the mean fluid temperature itself down for a
heat-rejecting load and up for a heat-extracting one) — see `script_gshp_3_nonstationaryFlow.jl` for
a worked comparison. Operating limits are still unenforced here. See [`ground_load`](@ref) for the
load conversion and [`_nonstationary_ground_temperature`](@ref) for the ground response, and the
README's "Conventions" section for the sign convention and `T_in`/`T_out` nomenclature.
# Arguments
    - as the fixed-COP [`gshp_simulation`](@ref) (case 1) method, except:
    - `V_candidates`: Controllable per-borehole flow rate set points, e.g. `[25, 30, 35] ./ 60000`
        for 25/30/35 L/min [m3/s]
    - `ΔT_target`: Target field inlet-outlet temperature difference used to pick the active flow
        rate at each time step [K]
# Output
    Named tuple with fields:
    - `Qg`: Ground thermal load [W]
    - `Qb_hp`: Building load actually delivered by the heat pumps [W]
    - `Qb_unmet`: Building load left unmet (auxiliary) [W]
    - `T_in`, `T_out`: Ground temperature entering/leaving the field [°C]
    - `n_active`: Fractional number of heat pumps operating, in `[0, n_HP]`
    - `V_t`: Per-borehole flow rate actually active at each time step [m3/s]
    - `Rbe`: Effective borehole thermal resistance for each flow-rate state [m.K/W]
    - `g`: g-function used for the ground response, shared by every flow state [°C.m/W]
# Reference
    - Beaudry, G., Pasquier, P., & Nguyen, A. (2024). New formulations and experimental validation
        of non-stationary convolutions for the fast simulation of time-variant flowrates in ground
        heat exchangers. Science and Technology for the Built Environment, 30(3), 208–219.
"""
function gshp_simulation(t::AbstractVector{<:Real}, Qb::AbstractVector{<:Real},
    COP_heating::Real, COP_cooling::Real, CAP_heating_unit::Real, CAP_cooling_unit::Real,
    n_HP::Real, xy::AbstractMatrix{<:Real}, rb::Real, H::Real, D::Real, ks::Real, Cs::Real,
    s::Real, ro::Real, ri::Real, kg::Real, kp::Real, kf::Real, cf::Real, ρf::Real, μf::Real,
    V_candidates::AbstractVector{<:Real}, ΔT_target::Real, T0::Real)

    # Thermal load conversion
    Qg, Qb_unmet = ground_load(Qb, COP_heating, COP_cooling, CAP_heating_unit, CAP_cooling_unit; n_HP)
    Qb_hp = Qb .- Qb_unmet # [W] building loads covered by the heat pump

    # Ground loop temperature evaluation, non-stationary flow rate
    T_in, T_out, V_t, Rbe, g = _nonstationary_ground_temperature(t, Qg, xy, rb, H, D, ks, Cs, s, ro,
        ri, kg, kp, kf, cf, ρf, μf, V_candidates, ΔT_target, T0)

    # Fractional number of heat pumps operating
    n_active = similar(Qb_hp, Float64)
    for k in eachindex(Qb_hp)
        CAP_unit = Qb_hp[k] >= 0 ? CAP_heating_unit : CAP_cooling_unit
        n_active[k] = CAP_unit > 0 ? clamp(abs(Qb_hp[k]) / CAP_unit, 0.0, n_HP) : 0.0
    end

    return (; Qg, Qb_hp, Qb_unmet, T_in, T_out, n_active, V_t, Rbe, g)
end