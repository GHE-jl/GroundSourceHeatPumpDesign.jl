# Simulation cases, in increasing order of complexity. See each function's docstring and the
# README's "Simulation cases" section for what each one covers.
#   1.  Fixed COP/CAP, stationary (constant) flow rate.
#   2.  Temperature-dependent COP/CAP curves, stationary flow rate, with a small fixed-point loop on
#       the source EFT.
#   3.  Fixed COP/CAP, non-stationary flow rate picked from the ground load (`_flow_rate_states`).
#   3b. Temperature-dependent COP/CAP curves, non-stationary flow rate (cases 2 and 3 together).
#   4.  Case 1 or 2 with operating limits `Tlim`: source-side mitigation of the ground load
#       (`_mitigation`) to keep the fluid within the limits.
#   5.  Case 3 or 3b with operating limits `Tlim`: largest flow rate near the limits
#       (`_flow_override`) and source-side mitigation.
# Every case takes the part-load degradation coefficient `Cd` (see `part_load_COP`). With the
# `fluid` keyword (raw-input methods), the effective borehole resistance is recalculated at the mean
# fluid temperature (`_fluid_temperature_loop`).
#
# Each case has two `gshp_simulation` methods: one takes an already-computed ground response (the
# g-function `g` and effective borehole resistance `Rbe` for a constant flow rate, or the outlet
# responses `h_out` for a non-stationary flow rate), the other takes the field's raw
# geometry/soil/pipe inputs, computes the ground response once, and delegates to the first.

using GroundHeatExchanger: resistance_ULoop_effective, gfunction, FLSModel,
    fluid_temperature, inlet_temperature, outlet_temperature, convolution_ns,
    outlet_gfunction, DeepANN, fluid_property

"""
    _source_fluid_temperature(t, Qg, nb, H, g, Rbe, C, V, T0; ΔR=nothing)

Source fluid temperature entering (`T_in`) and leaving (`T_out`) a borehole field for a constant
flow rate `V`, from an already-computed g-function `g` and effective borehole resistance `Rbe`, by
stationary (FFT) convolution of the ground load `Qg` with `g`. `ΔR` is the change of the field
resistance `Rbe / (nb·H)` at each time step [K/W], when `Rbe` is recalculated at the fluid
temperature (`fluid` keyword, see [`_fluid_temperature_loop`](@ref)). See the README's
"Conventions" section for the `T_in`/`T_out` nomenclature. Shared by every stationary
[`gshp_simulation`](@ref) methods.
"""
function _source_fluid_temperature(t::AbstractVector{<:Real}, Qg::AbstractVector{<:Real},
    nb::Integer, H::Real, g::AbstractVector{<:Real}, Rbe::Real, C::Real, V::Real, T0::Real;
    ΔR=nothing)
    # Define inputs
    L = nb * H
    q = Qg ./ L

    # Compute the mean fluid temperature and the inlet/outlet split
    Tf = fluid_temperature(t, q, g, T0, Rbe)
    ΔR === nothing || (Tf = Tf .+ Qg .* ΔR)
    T_in = inlet_temperature(Tf, q, H, V, C)
    T_out = outlet_temperature(Tf, q, H, V, C)
    return T_in, T_out
end

"""
    _source_fluid_temperature_ns(Qg, h_out, state_vec, V_candidates, nb, C, T0; ΔR=nothing)

Source fluid temperature entering (`T_in`) and leaving (`T_out`) a borehole field for a
**non-stationary** flow rate. Column `i` of `h_out` is the response of `T_out` to a unit step of
the ground load [K/W] when the flow rate `V_candidates[i]` is active, and `state_vec` gives the
active column at each time step. `T_out` is a single non-stationary convolution
(`GroundHeatExchanger.jl`'s `convolution_ns`), and `T_in = T_out + Qg / (nb · V_t · C)`. `ΔR` is
the change of the field resistance at each time step [K/W], as in
[`_source_fluid_temperature`](@ref), with one column per flow rate. Shared by every
non-stationary [`gshp_simulation`](@ref) method.
# Reference
    - Beaudry, G., Pasquier, P., & Nguyen, A. (2024). New formulations and experimental validation
        of non-stationary convolutions for the fast simulation of time-variant flowrates in ground
        heat exchangers. Science and Technology for the Built Environment, 30(3), 208–219.
"""
function _source_fluid_temperature_ns(Qg::AbstractVector{<:Real}, h_out::AbstractMatrix{<:Real},
    state_vec::AbstractVector{<:Integer}, V_candidates::AbstractVector{<:Real}, nb::Integer,
    C::Real, T0::Real; ΔR=nothing)
    V_t = V_candidates[state_vec]
    T_out = T0 .+ convolution_ns(collect(Float64, Qg), h_out, state_vec)
    ΔR === nothing || (T_out = T_out .+ Qg .* _active(ΔR, state_vec))
    T_in = T_out .+ Qg ./ (nb .* V_t .* C)
    return T_in, T_out
end

"""
    _active(ΔR, state_vec)

Change of the field resistance of the flow rate active at each time step: `ΔR[k, state_vec[k]]`.
"""
_active(ΔR::AbstractMatrix, state_vec::AbstractVector{<:Integer}) =
    ΔR[CartesianIndex.(eachindex(state_vec), state_vec)]

"""
    _outlet_responses(t, xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V_candidates;
        outlet_model=:resistance, Cg=nothing, Cp=1.9e6)

Response of the outlet temperature `T_out` to a unit step of the ground load [K/W], one column per
flow rate of `V_candidates` (sorted increasing), for the non-stationary [`gshp_simulation`](@ref)
methods. Two outlet models are available:
- `:resistance` (default): the mean fluid temperature from the g-function `g` and the effective
    borehole resistance `Rbe(V)`, minus half the loop temperature difference:
    `h_out = (Rbe(V) + g) / (nb·H) - 1 / (2·nb·V·C)`;
- `:gfunction`: the borehole outlet g-function of each flow rate
    (`outlet_gfunction` of `GroundHeatExchanger.jl`, short-term `DeepANN`):
    `h_out = ḡ(V) / (nb·V·C)`. Its short-term part captures the flow-dependent fluid transit and
    borehole thermal capacity that a constant `Rbe(V)` cannot (Beaudry et al. 2024). It needs the
    grout and pipe volumetric heat capacities `Cg` and `Cp`, and the DeepANN validity ranges apply
    (e.g. 15 to 35 L/min per borehole, H from 50 to 250 m).
Returns `h_out`, the effective borehole resistances `Rbe` (one per flow rate) and the g-function
`g`.
"""
function _outlet_responses(t::AbstractVector{<:Real}, xy::AbstractMatrix{<:Real}, rb::Real,
    H::Real, D::Real, ks::Real, Cs::Real, s::Real, ro::Real, ri::Real, kg::Real, kp::Real,
    kf::Real, cf::Real, ρf::Real, μf::Real, V_candidates::AbstractVector{<:Real};
    outlet_model::Symbol=:resistance, Cg=nothing, Cp::Real=1.9e6)
    # Ground response and effective borehole resistance of each flow rate
    nb = size(xy, 1)
    C = cf * ρf
    m = FLSModel(H, D, ks, Cs)
    g = gfunction(t, rb, xy, m)
    Rbe = [resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
           for V in V_candidates]

    # Outlet response of each flow rate
    if outlet_model === :resistance
        h_out = hcat([(Rbe[i] .+ g) ./ (nb * H) .- 1 / (2 * nb * V_candidates[i] * C)
                      for i in eachindex(V_candidates)]...)
    elseif outlet_model === :gfunction
        Cg === nothing &&
            throw(ArgumentError("outlet_model = :gfunction needs the grout heat " *
                "capacity Cg."))
        h_out = hcat([outlet_gfunction(t, ks, Cs, kg, Cg, kp, Cp, C, ri, ro, rb, H,
                          V_candidates[i], s, Rbe[i], m; xy=xy, interp=true, model=DeepANN()) ./
                      (nb * V_candidates[i] * C) for i in eachindex(V_candidates)]...)
    else
        throw(ArgumentError("outlet_model must be :resistance or :gfunction, got " *
            ":$outlet_model"))
    end
    return h_out, Rbe, g
end

"""
    _fluid_temperature_loop(simulate, V, nb, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf;
        fluid=nothing, η_pump=0.5, tol=0.05, max_iter=5)

Recalculation of the effective borehole resistance `Rb*` at the fluid temperature, and pumping
power, for the raw-input [`gshp_simulation`](@ref) methods. `simulate(ΔR)` runs the simulation
with the change of the field resistance `ΔR` [K/W] at each time step (`nothing` for none, one
column per flow rate), and `V` is the flow rate per borehole (a value, or the candidates of the
non-stationary cases).

Without `fluid`, the simulation is run once with `Rb*` at the given fluid properties. With
`fluid = (name, percentage)` (e.g. `(:MPG, 20)`, see [`fluid_property`](@ref)), the fluid
properties and `Rb*` are recalculated at the mean fluid temperature `(T_in + T_out) / 2` of the
previous pass (linear interpolation on a 0.5 °C grid), and the simulation is repeated until the
mean fluid temperature changes by less than `tol` (at most `max_iter` passes). The fluid heat
capacity `ρf·cf` of the inlet/outlet split stays at the given value.

The pumping power ([`pumping_power`](@ref)) uses the viscosity and density at the same mean fluid
temperature, and is zero when the heat pumps are off (`Qb_hp = 0`). Returns the simulation output
merged with `P_pump` and `n_iter_fluid`.
"""
function _fluid_temperature_loop(simulate, V, nb::Integer, H::Real, s::Real, rb::Real, ro::Real,
    ri::Real, ks::Real, kg::Real, kp::Real, kf::Real, cf::Real, ρf::Real, μf::Real;
    fluid=nothing, η_pump::Real=0.5, tol::Real=0.05, max_iter::Integer=5)
    # Fluid properties at a temperature, and Rb* for these properties and a flow rate
    props(T) = fluid === nothing ? (kf, cf, ρf, μf) :
        fluid_property(T, fluid[1]; percentage=fluid[2])
    Rbe(V, p) = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, p...)

    # Linear interpolation of `f(T)` at the temperatures `T`, on a 0.5 °C grid
    function on_grid(f, T)
        T_grid = floor(minimum(T)):0.5:(ceil(maximum(T)) + 0.5)
        f_grid = [f(Tg) for Tg in T_grid]
        k = clamp.(floor.(Int, (T .- T_grid[1]) ./ 0.5) .+ 1, 1, length(T_grid) - 1)
        w = (T .- T_grid[k]) ./ 0.5
        return (1 .- w) .* f_grid[k] .+ w .* f_grid[k .+ 1]
    end

    # First pass with Rb* at the given fluid properties
    res = simulate(nothing)
    V_t = haskey(res, :V_t) ? res.V_t : fill(float(V), length(res.Qg))
    T_mean = (res.T_in .+ res.T_out) ./ 2

    # Passes with Rb* at the mean fluid temperature of the previous pass
    n_iter_fluid = 1
    if fluid !== nothing
        converged, ΔT_mean = false, Inf
        for _ in 2:max_iter
            n_iter_fluid += 1
            ΔR = hcat([(on_grid(T -> Rbe(Vi, props(T)), T_mean) .- Rbe(Vi, (kf, cf, ρf, μf))) ./
                       (nb * H) for Vi in V]...)
            V isa Real && (ΔR = vec(ΔR))

            res = simulate(ΔR)
            V_t = haskey(res, :V_t) ? res.V_t : V_t
            T_mean_prev = T_mean
            T_mean = (res.T_in .+ res.T_out) ./ 2
            ΔT_mean = maximum(abs.(T_mean .- T_mean_prev))
            if ΔT_mean < tol
                converged = true
                break
            end
        end
        converged || @warn "The mean fluid temperature did not converge within $max_iter " *
            "passes of the Rb* recalculation (last change $(round(ΔT_mean; sigdigits=3)) °C)."
    end

    # Pumping power with the fluid properties at the mean fluid temperature, when the heat pumps run
    ρ_t = on_grid(T -> props(T)[3], T_mean)
    μ_t = on_grid(T -> props(T)[4], T_mean)
    P_pump = pumping_power.(V_t, nb, H, ri, ρ_t, μ_t; η=η_pump) .* (res.Qb_hp .!= 0)
    return merge(res, (; P_pump, n_iter_fluid))
end

"""
    _limit_weight(Tlim_at)

Where the operating limits `Tlim` apply, seen from the ground (not the heat pump), as the weight
`w` of `T = w·T_in + (1 - w)·T_out`: `:inlet` (`w = 1`, the heat-pump leaving fluid temperature,
LFT), `:outlet` (`w = 0`, the heat-pump entering fluid temperature, EFT) or `:mean` (`w = 1/2`).
"""
function _limit_weight(Tlim_at::Symbol)
    Tlim_at === :inlet && return 1.0
    Tlim_at === :outlet && return 0.0
    Tlim_at === :mean && return 0.5
    throw(ArgumentError("Tlim_at must be :inlet, :outlet or :mean, got :$Tlim_at"))
end

"""
    _mitigation(Qg_hp, temperature, a, Tlim; tol=1e-3, max_iter=50)

Source-side mitigation (Pasquier et al. 2013): reduce the ground load at the time steps where the
fluid temperature `T = temperature(Qg)` would leave the operating limits `Tlim = [low, high]`, so
that `T` stays at the limit there (`tol` inside it). At each pass, the ground load of those time
steps is corrected by `(Tlim - T) / a`, `a` being the response of `T` to the ground load of the
same time step [K/W].
The correction is repeated until it stops changing (it also accounts for the later time steps that
the earlier corrections already moved inside the limits). The ground load is only reduced towards
zero, never reversed: heat extraction (`Qg_hp < 0`) is shaved when `T < Tlim[1]` and heat rejection
(`Qg_hp > 0`) when `T > Tlim[2]`. The shaved part, `Qs_unmet = Qg_hp - Qg`, must be covered by a
source-side auxiliary system (e.g. a boiler in heating or a dry cooler in cooling).
Returns `(Qg, Qs_unmet, converged)`.
"""
function _mitigation(Qg_hp::AbstractVector{<:Real}, temperature, a::AbstractVector{<:Real},
    Tlim::AbstractVector{<:Real}; tol::Real=1e-3, max_iter::Integer=50)
    Qg = collect(Float64, Qg_hp)
    converged = false
    for _ in 1:max_iter
        T = temperature(Qg)

        # Time steps to correct: outside the limits, or already shaved
        extract = (Qg_hp .< 0) .& ((T .< Tlim[1]) .| (Qg .!= Qg_hp))
        reject = (Qg_hp .> 0) .& ((T .> Tlim[2]) .| (Qg .!= Qg_hp))

        # Corrected ground load, between the heat pump ground load and zero
        Qg_new = copy(Qg)
        Qg_new[extract] .= clamp.(Qg[extract] .+ (Tlim[1] + tol .- T[extract]) ./ a[extract],
            Qg_hp[extract], 0.0)
        Qg_new[reject] .= clamp.(Qg[reject] .+ (Tlim[2] - tol .- T[reject]) ./ a[reject], 0.0,
            Qg_hp[reject])

        # Convergence on the temperature change caused by the correction
        ΔT = maximum(abs.(Qg_new .- Qg) .* a; init=0.0)
        Qg = Qg_new
        if ΔT < tol / 10
            converged = true
            break
        end
    end
    converged || @warn "Source-side mitigation did not converge in $max_iter iterations."
    return Qg, Qg_hp .- Qg, converged
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
`Qg` alone (never on a simulated temperature); near the operating limits, it can be overridden
by [`_flow_override`](@ref).
# Arguments
    - `Qg`: Ground thermal load [W]
    - `nb`: Number of boreholes [-]
    - `V_candidates`: Controllable per-borehole flow rate set points, sorted increasing [m3/s]
    - `ΔT_target`: Target field inlet-outlet temperature difference used to pick the active flow
        rate [K]
    - `C`: Fluid volumetric specific heat [J/m3K] (`cf * ρf`)
# Output
    - `V_t`: Per-borehole flow rate active at each time step, one of `V_candidates` [m3/s]
    - `state_vec`: Index into `V_candidates` active at each time step [-]
"""
function _flow_rate_states(Qg::AbstractVector{<:Real}, nb::Integer,
    V_candidates::AbstractVector{<:Real}, ΔT_target::Real, C::Real)
    V_ideal = abs.(Qg) ./ (nb * ΔT_target * C)
    state_vec = similar(Qg, Int)
    V_t = similar(Qg, Float64)
    for k in eachindex(Qg)
        _, i = findmin(abs.(V_candidates .- V_ideal[k]))
        state_vec[k] = i
        V_t[k] = V_candidates[i]
    end
    return V_t, state_vec
end

"""
    _flow_override(state_vec, temperature, Tlim, T_margin, n_states; max_iter=20)

Flow-rate override of the non-stationary case (as in `fNSCiter.m`): the largest flow rate
(`n_states`) is forced at the time steps where the fluid temperature `T = temperature(state_vec)` is
within `T_margin` of an operating limit `Tlim = [low, high]`. Since the temperature depends on the
flow rates, this is repeated until no more time steps are switched. Returns the new `state_vec`.
"""
function _flow_override(state_vec::AbstractVector{<:Integer}, temperature,
    Tlim::AbstractVector{<:Real}, T_margin::Real, n_states::Integer; max_iter::Integer=20)
    state_vec = copy(state_vec)
    for _ in 1:max_iter
        T = temperature(state_vec)
        near = ((T .< Tlim[1] + T_margin) .| (T .> Tlim[2] - T_margin)) .&
            (state_vec .!= n_states)
        any(near) || break
        state_vec[near] .= n_states
    end
    return state_vec
end

"""
    gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating, CAP_cooling, n_HP,
        nb, H, g, Rbe, cf, ρf, V, T0; Cd=0.0, Tlim=nothing, Tlim_at=:inlet, ΔR=nothing)
    gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating, CAP_cooling, n_HP,
        xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0; fluid=nothing,
        η_pump=0.5, kwargs...)

Case 1: The simplest GSHP simulation: building load to ground load with **fixed**
(temperature-independent) COP and CAP and a **stationary** (constant flow rate) ground response.
Any load on the building side beyond the installed capacity is simply unmet (`Qb_unmet`). With
`CAP_heating = CAP_cooling = Inf`, the heat pumps are uncapped (base case). The steps are:
1. capacity limit of the `n_HP` heat pumps ([`ground_load`](@ref));
2. number of active heat pumps and part-load COP ([`n_active_heat_pumps`](@ref),
    [`part_load_COP`](@ref));
3. ground load from the building load delivered by the heat pumps ([`ground_load`](@ref));
4. (case 4, when `Tlim` is given) source-side mitigation of the ground load to keep the fluid
    temperature at `Tlim_at` within `Tlim` ([`_mitigation`](@ref));
5. ground loop temperatures ([`_source_fluid_temperature`](@ref)).
See the README's "Conventions" section for the sign convention and `T_in`/`T_out` nomenclature.

The first method takes an already-computed g-function `g` and effective borehole resistance `Rbe`.
The second method takes the field's raw geometry/soil/pipe inputs, computes `g` and `Rbe` once via
[`gfunction`](@ref)/[`resistance_ULoop_effective`](@ref), and calls the first method. It also
gives the pumping power and, with `fluid`, recalculates `Rbe` at the mean fluid temperature
([`_fluid_temperature_loop`](@ref)), which makes it a small fixed-point loop.
# Arguments
    - `t`: Time vector [s], uniformly spaced, matching `Qb`
    - `Qb`: Building thermal load [W], `> 0` heating, `< 0` cooling, same length as `t`
    - `COP_heating`, `COP_cooling`: Fixed nominal heating/cooling COP of one heat pump unit
    - `CAP_heating`, `CAP_cooling`: Fixed heating/cooling capacity of one heat pump unit [W]
    - `n_HP`: Number of installed heat pumps [-]
    - `nb`: Number of boreholes [-] (first method only; `size(xy, 1)` in the second)
    - `g`: Pre-computed g-function [°C.m/W], e.g. from [`gfunction`](@ref) (first method only)
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
# Keywords
    - `Cd`: Part-load degradation coefficient of the heat pumps (default `0.0`, no degradation)
    - `Tlim`: Operating temperature limits `[low, high]` [°C]; when given, the ground load is
        mitigated on the source side to respect them (case 4). Default `nothing` (no mitigation)
    - `Tlim_at`: Where `Tlim` applies, seen from the ground: `:inlet` (`T_in`, heat-pump LFT,
        default), `:outlet` (`T_out`, heat-pump EFT) or `:mean`
    - `ΔR`: Change of the field resistance `Rbe / (nb·H)` at each time step [K/W] (first method
        only, used by the `fluid` keyword of the second method)
    - `fluid`: Fluid `(name, percentage)` of [`fluid_property`](@ref), e.g. `(:MPG, 20)`; when
        given, `Rbe` is recalculated at the mean fluid temperature (second method only)
    - `η_pump`: Pump efficiency of the pumping power (default `0.5`, second method only)
# Output
    Named tuple with fields:
    - `Qg`: Ground thermal load, after mitigation [W]
    - `Qb_hp`: Building load actually delivered by the heat pumps (Qb = Qb_hp+Qb_unmet) [W]
    - `Qb_unmet`: Building load left unmet (building-side auxiliary) [W]
    - `Qs_unmet`: Ground load shaved by the mitigation (source-side auxiliary, same sign
        convention as `Qg`), zero without `Tlim` [W]
    - `T_in`, `T_out`: Ground temperature entering/leaving the field [°C]
    - `n_active`: Fractional number of heat pumps operating, in `[0, n_HP]`
    - `COP_actual`: Effective COP at each time step, including the part-load degradation [-]
    - `Rbe`: Effective borehole thermal resistance at the given fluid properties [m.K/W]
    - `g`: g-function used for the ground response [°C.m/W]
    - `P_pump`: Pumping power of the field at each time step ([`pumping_power`](@ref)) [W]
        (second method only)
    - `n_iter_fluid`: Number of passes of the `Rbe` recalculation (`1` without `fluid`, second
        method only)
# Reference
    - Pasquier, P., Marcotte, D., Bernier, M., & Kummert, M. (2013). Simulation of ground-coupled
        heat pump systems using a spectral approach. Proceedings of BS2013, 2691–2698.
"""
function gshp_simulation(t::AbstractVector{<:Real}, Qb::AbstractVector{<:Real},
    COP_heating::Real, COP_cooling::Real, CAP_heating::Real, CAP_cooling::Real,
    n_HP::Real, nb::Integer, H::Real, g::AbstractVector{<:Real}, Rbe::Real, cf::Real, ρf::Real,
    V::Real, T0::Real; Cd::Real=0.0, Tlim=nothing, Tlim_at::Symbol=:inlet, ΔR=nothing)
    C = cf * ρf

    # Fluid temperature where Tlim applies, and its response to the ground load of the same step
    w = _limit_weight(Tlim_at)
    function T_at(Q)
        T_in, T_out = _source_fluid_temperature(t, Q, nb, H, g, Rbe, C, V, T0; ΔR)
        return w .* T_in .+ (1 - w) .* T_out
    end
    a = fill((Rbe + g[1]) / (nb * H) + (w - 0.5) / (nb * V * C), length(Qb))
    ΔR === nothing || (a = a .+ ΔR)

    # Heat pumps: capacity limit, active units, part-load COP and ground load
    _, Qb_unmet, Qb_hp = ground_load(Qb, COP_heating, COP_cooling, CAP_heating, CAP_cooling; n_HP)
    n_active, COP = n_active_heat_pumps(Qb_hp, n_HP, CAP_heating, CAP_cooling, COP_heating,
        COP_cooling)
    COP_actual = part_load_COP(COP, n_active, Cd)
    Qg_hp = ground_load(Qb_hp, COP_actual, COP_actual)

    # Source-side mitigation (case 4)
    Qg, Qs_unmet = Qg_hp, zero(Qg_hp)
    if Tlim !== nothing
        Qg, Qs_unmet, _ = _mitigation(Qg_hp, T_at, a, Tlim)
    end

    # Ground loop temperatures
    T_in, T_out = _source_fluid_temperature(t, Qg, nb, H, g, Rbe, C, V, T0; ΔR)
    return (; Qg, Qb_hp, Qb_unmet, Qs_unmet, T_in, T_out, n_active, COP_actual, Rbe, g)
end
function gshp_simulation(t::AbstractVector{<:Real}, Qb::AbstractVector{<:Real},
    COP_heating::Real, COP_cooling::Real, CAP_heating::Real, CAP_cooling::Real,
    n_HP::Real, xy::AbstractMatrix{<:Real}, rb::Real, H::Real, D::Real, ks::Real, Cs::Real,
    s::Real, ro::Real, ri::Real, kg::Real, kp::Real, kf::Real, cf::Real, ρf::Real, μf::Real,
    V::Real, T0::Real; fluid=nothing, η_pump::Real=0.5, kwargs...)
    nb = size(xy, 1)
    Rbe = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
    g = gfunction(t, rb, xy, FLSModel(H, D, ks, Cs))
    simulate(ΔR) = gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating, CAP_cooling,
        n_HP, nb, H, g, Rbe, cf, ρf, V, T0; ΔR, kwargs...)
    return _fluid_temperature_loop(simulate, V, nb, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf;
        fluid, η_pump)
end

"""
    gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP, nb, H, g, Rbe, cf, ρf, V, T0;
        Tload=nothing, Cd=0.0, Tlim=nothing, Tlim_at=:inlet, tol=0.05, max_iter=5, ΔR=nothing)
    gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP,
        xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0; fluid=nothing,
        η_pump=0.5, kwargs...)

Case 2: Building load to ground load with **temperature-dependent** COP/CAP
([`heat_pump_performance`](@ref) curves) and a **stationary** (constant flow rate) ground response.
The building-to-ground conversion can't use the curve on its first pass (the source EFT isn't
known until the ground response has been simulated), so it starts from the mean of the curve as a
nominal COP to get a first `T_out` estimate. It then **iterates** the steps of case 1 (capacity
limit, active units and part-load COP, ground load, mitigation when `Tlim` is given (case 4),
ground loop temperatures) with the COP/CAP of the curves at the previous iteration's `T_out`, until
`T_out` stops moving by more than `tol` or `max_iter` is reached — a small fixed-point (Picard)
loop, not a full outer iteration over the whole simulation pipeline. This mirrors `fNSCiter.m`'s
`while` loop, which recomputes COP/CAP from the current source EFT every pass and is itself capped
at 5 iterations; it typically converges in well under that. A warning is given otherwise. See the
README's "Conventions" section for the sign convention and `T_in`/`T_out` nomenclature.

Two methods, in increasing order of complexity, exactly as case 1: the first takes precomputed
`g`/`Rbe` (reused, unchanged, across every iteration of the loop); the second takes the field's raw
geometry/soil/pipe inputs, computes them once, and calls the first.
# Arguments
    - `perf_heating`, `perf_cooling`: Single heat pump's heating/cooling performance data,
        forwarded to [`heat_pump_performance`](@ref): `(Tsrc_grid, CAP_vec, COP_vec)` (1D) or
        `(Tsrc_grid, Tload_grid, CAP_table, COP_table)` (2D, needs `Tload`)
    - remaining arguments: as in case 1
# Keywords
    - `Tload`: Load-side entering fluid temperature at each time step [°C], for 2D curves only
    - `tol`: Convergence tolerance on `T_out` between iterations [°C]
    - `max_iter`: Maximum number of fixed-point iterations
    - `Cd`, `Tlim`, `Tlim_at`, `ΔR`, `fluid`, `η_pump`: as in case 1
# Output
    Named tuple with the fields of case 1, plus:
    - `n_iter`: Number of fixed-point iterations actually taken [-]
    - `converged`: Whether `T_out` converged within `max_iter` iterations
# Reference
    - Pasquier, P., Marcotte, D., Bernier, M., & Kummert, M. (2013). Simulation of ground-coupled
        heat pump systems using a spectral approach. Proceedings of BS2013, 2691–2698.
"""
function gshp_simulation(t::AbstractVector{<:Real}, Qb::AbstractVector{<:Real},
    perf_heating::Tuple, perf_cooling::Tuple, n_HP::Real, nb::Integer, H::Real,
    g::AbstractVector{<:Real}, Rbe::Real, cf::Real, ρf::Real, V::Real, T0::Real;
    Tload=nothing, Cd::Real=0.0, Tlim=nothing, Tlim_at::Symbol=:inlet, tol::Real=0.05,
    max_iter::Integer=5, ΔR=nothing)
    C = cf * ρf
    length(perf_heating) == 4 && Tload === nothing &&
        throw(ArgumentError("2D performance curves need the load-side temperature `Tload`."))
    Tl = Tload === nothing ? () : (Tload,)

    # Fluid temperature where Tlim applies, and its response to the ground load of the same step
    w = _limit_weight(Tlim_at)
    function T_at(Q)
        T_in, T_out = _source_fluid_temperature(t, Q, nb, H, g, Rbe, C, V, T0; ΔR)
        return w .* T_in .+ (1 - w) .* T_out
    end
    a = fill((Rbe + g[1]) / (nb * H) + (w - 0.5) / (nb * V * C), length(Qb))
    ΔR === nothing || (a = a .+ ΔR)

    # First source EFT estimate: nominal COP (mean of the curve) and maximum capacity
    COP_heating = sum(perf_heating[end]) / length(perf_heating[end])
    COP_cooling = sum(perf_cooling[end]) / length(perf_cooling[end])
    Qg, _, _ = ground_load(Qb, COP_heating, COP_cooling, maximum(perf_heating[end - 1]),
        maximum(perf_cooling[end - 1]); n_HP)
    T_in, T_out = _source_fluid_temperature(t, Qg, nb, H, g, Rbe, C, V, T0; ΔR)

    # Fixed-point loop on the source EFT T_out
    Qb_hp, Qb_unmet, Qs_unmet, n_active, COP_actual = Qb, zero(Qb), zero(Qb), zero(Qb), zero(Qb)
    n_iter, converged, ΔT_out = 0, false, Inf
    for _ in 1:max_iter
        n_iter += 1
        T_out_prev = T_out

        # Heat pumps at the previous T_out: capacity limit, active units, part-load COP and
        # ground load
        _, Qb_unmet, Qb_hp = ground_load(Qb, T_out, Tl..., perf_heating, perf_cooling; n_HP)
        n_active, COP = n_active_heat_pumps(Qb_hp, T_out, Tl..., n_HP, perf_heating, perf_cooling)
        COP_actual = part_load_COP(COP, n_active, Cd)
        Qg_hp = ground_load(Qb_hp, COP_actual, COP_actual)

        # Source-side mitigation (case 4)
        Qg, Qs_unmet = Qg_hp, zero(Qg_hp)
        if Tlim !== nothing
            Qg, Qs_unmet, _ = _mitigation(Qg_hp, T_at, a, Tlim)
        end

        # Ground loop temperatures and convergence on T_out
        T_in, T_out = _source_fluid_temperature(t, Qg, nb, H, g, Rbe, C, V, T0; ΔR)
        ΔT_out = maximum(abs.(T_out .- T_out_prev))
        if ΔT_out < tol
            converged = true
            break
        end
    end
    converged || @warn "T_out did not converge within $max_iter iterations (tol = $tol °C, last " *
        "change $(round(ΔT_out; sigdigits=3)) °C)."
    return (; Qg, Qb_hp, Qb_unmet, Qs_unmet, T_in, T_out, n_active, COP_actual, n_iter,
        converged, Rbe, g)
end
function gshp_simulation(t::AbstractVector{<:Real}, Qb::AbstractVector{<:Real},
    perf_heating::Tuple, perf_cooling::Tuple, n_HP::Real, xy::AbstractMatrix{<:Real}, rb::Real,
    H::Real, D::Real, ks::Real, Cs::Real, s::Real, ro::Real, ri::Real, kg::Real, kp::Real,
    kf::Real, cf::Real, ρf::Real, μf::Real, V::Real, T0::Real; fluid=nothing,
    η_pump::Real=0.5, kwargs...)
    nb = size(xy, 1)
    Rbe = resistance_ULoop_effective(V, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf)
    g = gfunction(t, rb, xy, FLSModel(H, D, ks, Cs))
    simulate(ΔR) = gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP, nb, H, g, Rbe, cf, ρf,
        V, T0; ΔR, kwargs...)
    return _fluid_temperature_loop(simulate, V, nb, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf, μf;
        fluid, η_pump)
end

"""
    gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating, CAP_cooling, n_HP,
        nb, h_out, cf, ρf, V_candidates, ΔT_target, T0; Cd=0.0, Tlim=nothing, Tlim_at=:inlet,
        T_margin=1.0, ΔR=nothing)
    gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating, CAP_cooling, n_HP,
        xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V_candidates, ΔT_target, T0;
        outlet_model=:resistance, Cg=nothing, Cp=1.9e6, fluid=nothing, η_pump=0.5, kwargs...)

Case 3: the same heat pump steps and fixed COP/CAP as case 1, but the ground response switches
between a small set of controllable per-borehole flow rates `V_candidates` instead of a single
constant `V`. The steps are:
1. capacity limit, active units, part-load COP and ground load, as in case 1;
2. flow rate at each time step picked from the ground load ([`_flow_rate_states`](@ref));
3. (case 5, when `Tlim` is given) largest flow rate where the fluid temperature at `Tlim_at` is
    within `T_margin` of a limit ([`_flow_override`](@ref)), then source-side mitigation of the
    ground load ([`_mitigation`](@ref));
4. ground loop temperatures by non-stationary convolution
    ([`_source_fluid_temperature_ns`](@ref)).
Without `Tlim`, this is a single pass: a larger flow rate always cuts the `T_in`/`T_out` split for
the same `Qg` (it's proportional to `Qg/V`). Whether that also lowers the peak `T_out` at a given
extreme depends on the load's sign and size there too (a smaller split pulls `T_out` toward the
mean fluid temperature, while the smaller `Rbₑ` that comes with a bigger flow rate pulls the mean
fluid temperature itself down for a heat-rejecting load and up for a heat-extracting one) — see
`script_gshp_3_nonstationaryFlow.jl` for a worked comparison.

The first method takes the outlet responses `h_out` (one column per flow rate, see
[`_outlet_responses`](@ref)). The second method takes the field's raw geometry/soil/pipe inputs,
computes them once with the chosen `outlet_model`, and calls the first.
# Arguments
    - `h_out`: Response of `T_out` to a unit step of the ground load, one column per flow rate of
        `V_candidates` [K/W] (first method only)
    - `V_candidates`: Controllable per-borehole flow rate set points, sorted increasing, e.g.
        `[25, 30, 35] ./ 60000` for 25/30/35 L/min [m3/s]
    - `ΔT_target`: Target field inlet-outlet temperature difference used to pick the active flow
        rate at each time step [K]
    - remaining arguments: as in case 1
# Keywords
    - `T_margin`: Temperature margin to the limits below which the largest flow rate is forced
        (case 5, default 1 °C)
    - `outlet_model`: `:resistance` (default) or `:gfunction` (second method only, see
        [`_outlet_responses`](@ref))
    - `Cg`, `Cp`: Grout and pipe volumetric heat capacities [J/m3K], for
        `outlet_model = :gfunction` only (second method only)
    - `ΔR`: as in case 1, one column per flow rate of `V_candidates` (first method only)
    - `Cd`, `Tlim`, `Tlim_at`, `η_pump`: as in case 1
    - `fluid`: as in case 1, with `outlet_model = :resistance` only
# Output
    Named tuple with the fields of case 1 (`Rbe` one per flow rate, `g` second method only), plus:
    - `V_t`: Per-borehole flow rate actually active at each time step [m3/s]
    - `h_out`: Outlet responses used for the ground response [K/W]
# Reference
    - Beaudry, G., Pasquier, P., & Nguyen, A. (2024). New formulations and experimental validation
        of non-stationary convolutions for the fast simulation of time-variant flowrates in ground
        heat exchangers. Science and Technology for the Built Environment, 30(3), 208–219.
"""
function gshp_simulation(t::AbstractVector{<:Real}, Qb::AbstractVector{<:Real},
    COP_heating::Real, COP_cooling::Real, CAP_heating::Real, CAP_cooling::Real,
    n_HP::Real, nb::Integer, h_out::AbstractMatrix{<:Real}, cf::Real, ρf::Real,
    V_candidates::AbstractVector{<:Real}, ΔT_target::Real, T0::Real; Cd::Real=0.0, Tlim=nothing,
    Tlim_at::Symbol=:inlet, T_margin::Real=1.0, ΔR=nothing)
    C = cf * ρf
    issorted(V_candidates) || throw(ArgumentError("V_candidates must be sorted increasing."))

    # Fluid temperature where Tlim applies, for the flow rate states `state`
    w = _limit_weight(Tlim_at)
    function T_at(Q, state)
        T_in, T_out = _source_fluid_temperature_ns(Q, h_out, state, V_candidates, nb, C, T0; ΔR)
        return w .* T_in .+ (1 - w) .* T_out
    end

    # Heat pumps: capacity limit, active units, part-load COP and ground load
    _, Qb_unmet, Qb_hp = ground_load(Qb, COP_heating, COP_cooling, CAP_heating, CAP_cooling; n_HP)
    n_active, COP = n_active_heat_pumps(Qb_hp, n_HP, CAP_heating, CAP_cooling, COP_heating,
        COP_cooling)
    COP_actual = part_load_COP(COP, n_active, Cd)
    Qg_hp = ground_load(Qb_hp, COP_actual, COP_actual)

    # Flow rate at each time step, picked from the ground load
    _, state_vec = _flow_rate_states(Qg_hp, nb, V_candidates, ΔT_target, C)

    # Case 5: largest flow rate near the limits, then source-side mitigation
    Qg, Qs_unmet = Qg_hp, zero(Qg_hp)
    if Tlim !== nothing
        state_vec = _flow_override(state_vec, state -> T_at(Qg_hp, state), Tlim, T_margin,
            length(V_candidates))
        a = h_out[1, state_vec] .+ w ./ (nb .* V_candidates[state_vec] .* C)
        ΔR === nothing || (a = a .+ _active(ΔR, state_vec))
        Qg, Qs_unmet, _ = _mitigation(Qg_hp, Q -> T_at(Q, state_vec), a, Tlim)
    end

    # Ground loop temperatures
    V_t = V_candidates[state_vec]
    T_in, T_out = _source_fluid_temperature_ns(Qg, h_out, state_vec, V_candidates, nb, C, T0;
        ΔR)
    return (; Qg, Qb_hp, Qb_unmet, Qs_unmet, T_in, T_out, n_active, COP_actual, V_t, h_out)
end
function gshp_simulation(t::AbstractVector{<:Real}, Qb::AbstractVector{<:Real},
    COP_heating::Real, COP_cooling::Real, CAP_heating::Real, CAP_cooling::Real,
    n_HP::Real, xy::AbstractMatrix{<:Real}, rb::Real, H::Real, D::Real, ks::Real, Cs::Real,
    s::Real, ro::Real, ri::Real, kg::Real, kp::Real, kf::Real, cf::Real, ρf::Real, μf::Real,
    V_candidates::AbstractVector{<:Real}, ΔT_target::Real, T0::Real;
    outlet_model::Symbol=:resistance, Cg=nothing, Cp::Real=1.9e6, fluid=nothing,
    η_pump::Real=0.5, kwargs...)
    fluid !== nothing && outlet_model !== :resistance &&
        throw(ArgumentError("fluid needs outlet_model = :resistance."))
    nb = size(xy, 1)
    V_sorted = sort(V_candidates)
    h_out, Rbe, g = _outlet_responses(t, xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf,
        V_sorted; outlet_model, Cg, Cp)
    simulate(ΔR) = gshp_simulation(t, Qb, COP_heating, COP_cooling, CAP_heating, CAP_cooling, n_HP,
        nb, h_out, cf, ρf, V_sorted, ΔT_target, T0; ΔR, kwargs...)
    res = _fluid_temperature_loop(simulate, V_sorted, nb, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf,
        μf; fluid, η_pump)
    return merge(res, (; Rbe, g))
end

"""
    gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP, nb, h_out, cf, ρf, V_candidates,
        ΔT_target, T0; Tload=nothing, Cd=0.0, Tlim=nothing, Tlim_at=:inlet, T_margin=1.0,
        tol=0.05, max_iter=5, ΔR=nothing)
    gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP,
        xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V_candidates, ΔT_target, T0;
        outlet_model=:resistance, Cg=nothing, Cp=1.9e6, fluid=nothing, η_pump=0.5, kwargs...)

Case 3b: temperature-dependent COP/CAP curves (case 2) with a non-stationary flow rate (case 3).
The fixed-point loop of case 2 is kept: at each iteration, the heat pump steps use the COP/CAP of
the curves at the previous iteration's `T_out`, the flow rate is picked from the new ground load
and, when `Tlim` is given (case 5), the largest flow rate is forced near the limits and the ground
load is mitigated on the source side, before the ground loop temperatures are updated by
non-stationary convolution.
# Arguments
    - `perf_heating`, `perf_cooling`: as in case 2
    - remaining arguments: as in case 3
# Keywords
    - `Tload`, `tol`, `max_iter`: as in case 2
    - `T_margin`, `outlet_model`, `Cg`, `Cp`, `Cd`, `Tlim`, `Tlim_at`, `ΔR`, `fluid`, `η_pump`: as
        in case 3
# Output
    Named tuple with the fields of case 3, plus `n_iter` and `converged` (as in case 2).
"""
function gshp_simulation(t::AbstractVector{<:Real}, Qb::AbstractVector{<:Real},
    perf_heating::Tuple, perf_cooling::Tuple, n_HP::Real, nb::Integer,
    h_out::AbstractMatrix{<:Real}, cf::Real, ρf::Real, V_candidates::AbstractVector{<:Real},
    ΔT_target::Real, T0::Real; Tload=nothing, Cd::Real=0.0, Tlim=nothing, Tlim_at::Symbol=:inlet,
    T_margin::Real=1.0, tol::Real=0.05, max_iter::Integer=5, ΔR=nothing)
    C = cf * ρf
    issorted(V_candidates) || throw(ArgumentError("V_candidates must be sorted increasing."))
    length(perf_heating) == 4 && Tload === nothing &&
        throw(ArgumentError("2D performance curves need the load-side temperature `Tload`."))
    Tl = Tload === nothing ? () : (Tload,)

    # Fluid temperature where Tlim applies, for the flow rate states `state`
    w = _limit_weight(Tlim_at)
    function T_at(Q, state)
        T_in, T_out = _source_fluid_temperature_ns(Q, h_out, state, V_candidates, nb, C, T0; ΔR)
        return w .* T_in .+ (1 - w) .* T_out
    end

    # First source EFT estimate: nominal COP (mean of the curve), maximum capacity, flow rate
    # picked from the ground load
    COP_heating = sum(perf_heating[end]) / length(perf_heating[end])
    COP_cooling = sum(perf_cooling[end]) / length(perf_cooling[end])
    Qg, _, _ = ground_load(Qb, COP_heating, COP_cooling, maximum(perf_heating[end - 1]),
        maximum(perf_cooling[end - 1]); n_HP)
    _, state_vec = _flow_rate_states(Qg, nb, V_candidates, ΔT_target, C)
    T_in, T_out = _source_fluid_temperature_ns(Qg, h_out, state_vec, V_candidates, nb, C, T0;
        ΔR)

    # Fixed-point loop on the source EFT T_out
    Qb_hp, Qb_unmet, Qs_unmet, n_active, COP_actual = Qb, zero(Qb), zero(Qb), zero(Qb), zero(Qb)
    n_iter, converged, ΔT_out = 0, false, Inf
    for _ in 1:max_iter
        n_iter += 1
        T_out_prev = T_out

        # Heat pumps at the previous T_out: capacity limit, active units, part-load COP and
        # ground load
        _, Qb_unmet, Qb_hp = ground_load(Qb, T_out, Tl..., perf_heating, perf_cooling; n_HP)
        n_active, COP = n_active_heat_pumps(Qb_hp, T_out, Tl..., n_HP, perf_heating, perf_cooling)
        COP_actual = part_load_COP(COP, n_active, Cd)
        Qg_hp = ground_load(Qb_hp, COP_actual, COP_actual)

        # Flow rate at each time step, picked from the ground load. After the first iteration, a
        # flow rate can only increase, so that the switching between the discrete flow rates
        # does not prevent the convergence.
        _, state_new = _flow_rate_states(Qg_hp, nb, V_candidates, ΔT_target, C)
        state_vec = n_iter == 1 ? state_new : max.(state_vec, state_new)

        # Case 5: largest flow rate near the limits, then source-side mitigation
        Qg, Qs_unmet = Qg_hp, zero(Qg_hp)
        if Tlim !== nothing
            state_vec = _flow_override(state_vec, state -> T_at(Qg_hp, state), Tlim, T_margin,
                length(V_candidates))
            a = h_out[1, state_vec] .+ w ./ (nb .* V_candidates[state_vec] .* C)
            ΔR === nothing || (a = a .+ _active(ΔR, state_vec))
            Qg, Qs_unmet, _ = _mitigation(Qg_hp, Q -> T_at(Q, state_vec), a, Tlim)
        end

        # Ground loop temperatures and convergence on T_out
        T_in, T_out = _source_fluid_temperature_ns(Qg, h_out, state_vec, V_candidates, nb, C, T0;
            ΔR)
        ΔT_out = maximum(abs.(T_out .- T_out_prev))
        if ΔT_out < tol
            converged = true
            break
        end
    end
    converged || @warn "T_out did not converge within $max_iter iterations (tol = $tol °C, last " *
        "change $(round(ΔT_out; sigdigits=3)) °C)."
    V_t = V_candidates[state_vec]
    return (; Qg, Qb_hp, Qb_unmet, Qs_unmet, T_in, T_out, n_active, COP_actual, n_iter, converged,
        V_t, h_out)
end
function gshp_simulation(t::AbstractVector{<:Real}, Qb::AbstractVector{<:Real},
    perf_heating::Tuple, perf_cooling::Tuple, n_HP::Real, xy::AbstractMatrix{<:Real}, rb::Real,
    H::Real, D::Real, ks::Real, Cs::Real, s::Real, ro::Real, ri::Real, kg::Real, kp::Real,
    kf::Real, cf::Real, ρf::Real, μf::Real, V_candidates::AbstractVector{<:Real},
    ΔT_target::Real, T0::Real; outlet_model::Symbol=:resistance, Cg=nothing, Cp::Real=1.9e6,
    fluid=nothing, η_pump::Real=0.5, kwargs...)
    fluid !== nothing && outlet_model !== :resistance &&
        throw(ArgumentError("fluid needs outlet_model = :resistance."))
    nb = size(xy, 1)
    V_sorted = sort(V_candidates)
    h_out, Rbe, g = _outlet_responses(t, xy, rb, H, D, ks, Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf,
        V_sorted; outlet_model, Cg, Cp)
    simulate(ΔR) = gshp_simulation(t, Qb, perf_heating, perf_cooling, n_HP, nb, h_out, cf, ρf,
        V_sorted, ΔT_target, T0; ΔR, kwargs...)
    res = _fluid_temperature_loop(simulate, V_sorted, nb, H, s, rb, ro, ri, ks, kg, kp, kf, cf, ρf,
        μf; fluid, η_pump)
    return merge(res, (; Rbe, g))
end
