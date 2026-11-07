using GroundHeatExchanger: Reynolds, friction_factor_Colebrook_White

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
    pumping_power(V, nb, H, ri, ρf, μf; ϵ=5e-6, K=1.5 + 2 * 0.9, η=0.5)

Electrical power of the circulation pump of a borehole field with single U-loops in parallel. The
pressure drop of one U-loop is given by the Darcy-Weisbach equation over the pipe length `2H`, plus
the minor losses of one 180° bend at the bottom and two 90° elbows at the top:
`ΔP = (f·2H/(2ri) + K)·ρf·v²/2`, with `v = V/(π·ri²)` and the Darcy friction factor `f`
([`friction_factor_Colebrook_White`](@ref), `64/Re` in laminar flow). The field power is
`P = ΔP·nb·V/η`, which is the same as the pressure drop over the whole length of the field
(`2H·nb`) for the flow rate of one borehole. The headers and the horizontal pipes are not included.
# Arguments
    - `V`: Flow rate per borehole [m3/s]
    - `nb`: Number of boreholes [-]
    - `H`: Borehole depth [m]
    - `ri`: Pipe inner radius [m]
    - `ρf`, `μf`: Fluid density [kg/m3] and dynamic viscosity [kg/(m.s)]
# Keywords
    - `ϵ`: Pipe roughness [m] (default `5e-6`, HDPE)
    - `K`: Minor loss coefficient of one U-loop [-] (default `1.5` for the 180° bend plus `0.9` for
        each 90° elbow)
    - `η`: Pump efficiency [-] (default `0.5`)
# Output
    - `P`: Electrical pumping power [W]
"""
function pumping_power(V::Real, nb::Integer, H::Real, ri::Real, ρf::Real, μf::Real;
    ϵ::Real=5e-6, K::Real=1.5 + 2 * 0.9, η::Real=0.5)
    v = V / (π * ri^2)                                  # fluid velocity [m/s]
    f = friction_factor_Colebrook_White(Reynolds(v, ri, ρf, μf), ri, ϵ)
    ΔP = (f * 2 * H / (2 * ri) + K) * ρf * v^2 / 2      # pressure drop of one U-loop [Pa]
    return ΔP * nb * V / η
end

"""
    report_gshp_simulation(T_in, T_out, Qb, Qb_hp, Qb_unmet, Qg, Tlim, steps_per_year;
        Tlim_at=:inlet, Qs_unmet=nothing, COP_actual=nothing, P_pump=nothing)
    report_gshp_simulation(res, Qb, Tlim, steps_per_year; Tlim_at=:inlet)

Print a formatted report of a GSHP simulation's temperature and energy results (e.g.
[`gshp_simulation`](@ref)'s output): ground temperature extremes, unmet hours and long-term drift
over the last simulated year, then a table of energy (MWh) and peak power (kW) for the building,
heat pump and ground loads over the whole simulation (and the source-side auxiliary when
`Qs_unmet` is given), and the annual ground load imbalance. With `COP_actual` and `P_pump`, the
electricity of the heat pumps (`|Qb_hp| / COP_actual`) and of the pump is added, with the seasonal
performance factors of the heat pumps (`SPF_hp = Σ|Qb_hp| / ΣW_hp`) and of the system, including
the pumping (`SPF_sys = Σ|Qb_hp| / Σ(W_hp + W_pump)`). The auxiliary systems are not included.

The second method takes the output `res` of [`gshp_simulation`](@ref) and uses every field
available in it.
# Arguments
    - `T_in`, `T_out`: Ground temperature entering/leaving the field, time series [°C]
    - `Qb`: Building thermal load, time series [W]
    - `Qb_hp`: Building load actually delivered by the heat pump(s), time series [W]
    - `Qb_unmet`: Building load left unmet (auxiliary), time series [W]
    - `Qg`: Ground thermal load, time series [W]
    - `Tlim`: Operating temperature limits `[T_low, T_high]`, checked at `Tlim_at` [°C]
    - `steps_per_year`: Number of time steps per year (e.g. `8760` for an hourly simulation);
        every time series above must cover an integer number of years
# Keywords
    - `Tlim_at`: Where `Tlim` is checked, seen from the ground: `:inlet` (`T_in`, heat-pump LFT,
        default), `:outlet` (`T_out`, heat-pump EFT) or `:mean`
    - `Qs_unmet`: Ground load shaved by the source-side mitigation (output of
        [`gshp_simulation`](@ref) with `Tlim`), time series [W]
    - `COP_actual`: Effective COP of the heat pumps, time series [-]
    - `P_pump`: Pumping power, time series [W] (see [`pumping_power`](@ref))
"""
function report_gshp_simulation(T_in::AbstractVector{<:Real}, T_out::AbstractVector{<:Real},
    Qb::AbstractVector{<:Real}, Qb_hp::AbstractVector{<:Real}, Qb_unmet::AbstractVector{<:Real},
    Qg::AbstractVector{<:Real}, Tlim::AbstractVector{<:Real}, steps_per_year::Integer;
    Tlim_at::Symbol=:inlet, Qs_unmet=nothing, COP_actual=nothing, P_pump=nothing)

    # Check that all time series have the same length
    nt = length(T_in)
    n_years = nt ÷ steps_per_year
    last_year = (nt - steps_per_year + 1):nt

    # Fluid temperature where Tlim is checked
    w = _limit_weight(Tlim_at)
    T = w .* T_in .+ (1 - w) .* T_out

    # Print ground temperature extremes, unmet hours and drift over the last simulated year
    println("Ground temperature, year $n_years (last simulated year):")
    println("  T_in  min/max = $(round(minimum(T_in[last_year]); digits=1)) / " *
        "$(round(maximum(T_in[last_year]); digits=1)) degC")
    println("  T_out min/max = $(round(minimum(T_out[last_year]); digits=1)) / " *
        "$(round(maximum(T_out[last_year]); digits=1)) degC")
    println("  Hours with T at the $Tlim_at outside [$(Tlim[1]), $(Tlim[2])] degC: " *
        "$(unmet_hours(T[last_year], Tlim))")
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
    if Qs_unmet !== nothing
        rows = (rows...,
            ("Source auxiliary heating", -MWh(min.(Qs_unmet, 0.0)), -minimum(Qs_unmet) / 1000),
            ("Source auxiliary cooling", MWh(max.(Qs_unmet, 0.0)), maximum(Qs_unmet) / 1000))
    end

    # Electricity of the heat pumps (compressors) and of the pump
    W_hp = COP_actual === nothing ? nothing : abs.(Qb_hp) ./ COP_actual
    if W_hp !== nothing
        rows = (rows..., ("Heat pump electricity", MWh(W_hp), maximum(W_hp) / 1000))
    end
    if P_pump !== nothing
        rows = (rows..., ("Pumping electricity", MWh(P_pump), maximum(P_pump) / 1000))
    end

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

    # Annual ground load imbalance (> 0: more heat rejected than extracted)
    println("\nAnnual ground load imbalance (rejection - extraction): " *
        "$(round(MWh(Qg) / n_years; digits=1)) MWh/year " *
        "($(round(100 * sum(Qg) / sum(abs, Qg); digits=1)) % of the ground load)")

    # Seasonal performance factors
    if W_hp !== nothing
        println("Seasonal performance factor, heat pumps: " *
            "$(round(sum(abs, Qb_hp) / sum(W_hp); digits=2))")
        if P_pump !== nothing
            println("Seasonal performance factor, heat pumps and pumping: " *
                "$(round(sum(abs, Qb_hp) / (sum(W_hp) + sum(P_pump)); digits=2))")
        end
    end

    return nothing
end
function report_gshp_simulation(res::NamedTuple, Qb::AbstractVector{<:Real},
    Tlim::AbstractVector{<:Real}, steps_per_year::Integer; Tlim_at::Symbol=:inlet)
    return report_gshp_simulation(res.T_in, res.T_out, Qb, res.Qb_hp, res.Qb_unmet, res.Qg, Tlim,
        steps_per_year; Tlim_at, Qs_unmet=any(!=(0), res.Qs_unmet) ? res.Qs_unmet : nothing,
        COP_actual=res.COP_actual,
        P_pump=get(res, :P_pump, nothing))
end
