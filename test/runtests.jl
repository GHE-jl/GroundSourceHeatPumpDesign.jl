using Test
using GroundSourceHeatPumpDesign

@testset "GroundSourceHeatPumpDesign.jl" begin

    @testset "heat_pump_performance — 1D" begin
        Tsrc_grid = [-5.0, 0.0, 5.0, 10.0]
        CAP_vec = [8.0, 9.0, 10.0, 11.0]
        COP_vec = [3.0, 3.4, 3.8, 4.2]

        # Exact grid points
        CAP, COP = heat_pump_performance(0.0, Tsrc_grid, CAP_vec, COP_vec)
        @test CAP ≈ 9.0
        @test COP ≈ 3.4

        # Midpoint linear interpolation
        CAP, COP = heat_pump_performance(2.5, Tsrc_grid, CAP_vec, COP_vec)
        @test CAP ≈ 9.5
        @test COP ≈ 3.6

        # Clamping outside the grid
        CAP, COP = heat_pump_performance(-100.0, Tsrc_grid, CAP_vec, COP_vec)
        @test CAP ≈ 8.0 && COP ≈ 3.0
        CAP, COP = heat_pump_performance(100.0, Tsrc_grid, CAP_vec, COP_vec)
        @test CAP ≈ 11.0 && COP ≈ 4.2

        # Vector overload matches scalar calls elementwise
        CAP_v, COP_v = heat_pump_performance([0.0, 2.5], Tsrc_grid, CAP_vec, COP_vec)
        @test CAP_v ≈ [9.0, 9.5]
        @test COP_v ≈ [3.4, 3.6]
    end

    @testset "heat_pump_performance — 2D" begin
        Tsrc_grid = [0.0, 10.0]
        Tload_grid = [30.0, 40.0]
        CAP_table = [8.0 6.0; 10.0 8.0]
        COP_table = [4.0 3.0; 5.0 4.0]

        # Exact grid corners
        CAP, COP = heat_pump_performance(0.0, 30.0, Tsrc_grid, Tload_grid, CAP_table, COP_table)
        @test CAP ≈ 8.0 && COP ≈ 4.0
        CAP, COP = heat_pump_performance(10.0, 40.0, Tsrc_grid, Tload_grid, CAP_table, COP_table)
        @test CAP ≈ 8.0 && COP ≈ 4.0

        # Center point: average of the four corners
        CAP, COP = heat_pump_performance(5.0, 35.0, Tsrc_grid, Tload_grid, CAP_table, COP_table)
        @test CAP ≈ (8.0 + 6.0 + 10.0 + 8.0) / 4
        @test COP ≈ (4.0 + 3.0 + 5.0 + 4.0) / 4

        # Vector overload matches scalar calls elementwise
        CAP_v, COP_v = heat_pump_performance([0.0, 5.0], [30.0, 35.0], Tsrc_grid, Tload_grid,
            CAP_table, COP_table)
        @test CAP_v[1] ≈ 8.0
        @test CAP_v[2] ≈ (8.0 + 6.0 + 10.0 + 8.0) / 4
    end

    @testset "ground_load — method 1 (fixed COP)" begin
        # Heating: Qg = -Qb * (1 - 1/COP_heating)
        @test ground_load(8.0, 4.0, 5.0) ≈ -8.0 * (1 - 1 / 4.0)
        # Cooling: Qg = -Qb * (1 + 1/COP_cooling)
        @test ground_load(-8.0, 4.0, 5.0) ≈ 8.0 * (1 + 1 / 5.0)
        # Zero load
        @test ground_load(0.0, 4.0, 5.0) == 0.0

        # Vector overload, scalar COP
        @test ground_load([8.0, -8.0], 4.0, 5.0) ≈
              [ground_load(8.0, 4.0, 5.0), ground_load(-8.0, 4.0, 5.0)]

        # Vector overload, vector COP
        Qb = [8.0, -8.0]
        @test ground_load(Qb, [4.0, 4.5], [5.0, 5.5]) ≈
              [ground_load(8.0, 4.0, 5.0), ground_load(-8.0, 4.5, 5.5)]

        # Non-positive COP is rejected
        @test_throws ArgumentError ground_load(8.0, 0.0, 5.0)
        @test_throws ArgumentError ground_load(8.0, 4.0, -1.0)
        @test_throws DimensionMismatch ground_load([8.0, -8.0], [4.0], [5.0, 5.5])
    end

    @testset "ground_load — method 2 (fixed COP, capacity-capped)" begin
        # Below capacity: behaves like method 1 with the nominal COP
        Qg, Qb_unmet = ground_load(5.0, 4.0, 5.0, 10.0, 10.0)
        @test Qg ≈ ground_load(5.0, 4.0, 5.0)
        @test Qb_unmet == 0.0

        # Above heating capacity: capped, with unmet remainder
        Qg, Qb_unmet = ground_load(15.0, 4.0, 5.0, 10.0, 10.0)
        @test Qg ≈ ground_load(10.0, 4.0, 5.0)
        @test Qb_unmet ≈ 5.0

        # Above cooling capacity (Qb very negative): capped, with unmet remainder
        Qg, Qb_unmet = ground_load(-15.0, 4.0, 5.0, 10.0, 10.0)
        @test Qg ≈ ground_load(-10.0, 4.0, 5.0)
        @test Qb_unmet ≈ -5.0

        # Vector overload
        Qb = [5.0, 15.0, -15.0]
        Qg, Qb_unmet = ground_load(Qb, 4.0, 5.0, 10.0, 10.0)
        @test Qg ≈ [ground_load(5.0, 4.0, 5.0), ground_load(10.0, 4.0, 5.0),
            ground_load(-10.0, 4.0, 5.0)]
        @test Qb_unmet ≈ [0.0, 5.0, -5.0]

        # n_HP scales the per-unit CAP to the total installed capacity: 3 units of CAP=10.0 behaves
        # like a single CAP_heating_unit/CAP_cooling_unit of 30.0 (n_HP=1)
        Qg3, Qb_unmet3 = ground_load(25.0, 4.0, 5.0, 10.0, 10.0; n_HP=3)
        @test (Qg3, Qb_unmet3) == ground_load(25.0, 4.0, 5.0, 30.0, 30.0)
        @test Qb_unmet3 == 0.0
        Qg3, Qb_unmet3 = ground_load(35.0, 4.0, 5.0, 10.0, 10.0; n_HP=3)
        @test Qg3 ≈ ground_load(30.0, 4.0, 5.0)
        @test Qb_unmet3 ≈ 5.0

        # Vector overload with n_HP
        Qg3_v, Qb_unmet3_v = ground_load(Qb, 4.0, 5.0, 10.0, 10.0; n_HP=3)
        @test Qg3_v ≈ [ground_load(Qb[k], 4.0, 5.0, 10.0, 10.0; n_HP=3)[1] for k in eachindex(Qb)]
        @test all(Qb_unmet3_v .== 0.0) # below the n_HP=3 installed capacity everywhere
    end

    @testset "n_active_heat_pumps" begin
        Tsrc_grid = [0.0, 10.0]
        CAP_vec = [10.0, 10.0]
        COP_vec = [4.0, 4.0]
        perf = (Tsrc_grid, CAP_vec, COP_vec)

        n, COP = n_active_heat_pumps([5.0, 20.0, -5.0], [5.0, 5.0, 5.0], 3, perf, perf)
        @test n ≈ [0.5, 2.0, 0.5]
        @test COP ≈ [4.0, 4.0, 4.0]

        # Capped at n_HP
        n, _ = n_active_heat_pumps([100.0], [5.0], 3, perf, perf)
        @test n ≈ [3.0]

        @test_throws DimensionMismatch n_active_heat_pumps([1.0, 2.0], [1.0], 3, perf, perf)
    end

    @testset "reporting utilities" begin
        T = [0.0, 2.0, -1.0, 5.0, 31.0]
        @test unmet_hours(T, [0.0, 30.0]) == 2   # -1.0 and 31.0 are outside [0, 30]

        T2 = [0.0, 2.0, 4.0, 6.0]   # two years of 2 steps: mean 1.0, then mean 5.0
        @test annual_drift(T2, 2) ≈ 4.0

        @test thermal_energy([1.0, 2.0, 3.0], 1.0) ≈ 6.0
        @test thermal_energy([1.0, 2.0, 3.0], 0.5) ≈ 3.0
    end

    @testset "ground_load — method 3 (performance-curve-based)" begin
        @testset "2D (Tsrc, Tload)" begin
            Tsrc_grid = [0.0, 10.0]
            Tload_grid = [30.0, 40.0]
            CAP_table = [8.0 6.0; 10.0 8.0]
            COP_table = [4.0 3.0; 5.0 4.0]
            perf = (Tsrc_grid, Tload_grid, CAP_table, COP_table)

            # Heating load fully met by the heat pump (below capacity)
            Qg, Qunmet = ground_load(5.0, 5.0, 35.0, perf, perf)
            CAP, COP = heat_pump_performance(5.0, 35.0, Tsrc_grid, Tload_grid, CAP_table, COP_table)
            @test Qg ≈ -5.0 * (1 - 1 / COP)
            @test Qunmet ≈ 0.0

            # Heating load exceeding capacity: capped, with unmet remainder
            Qg, Qunmet = ground_load(100.0, 5.0, 35.0, perf, perf)
            @test Qg ≈ -CAP * (1 - 1 / COP)
            @test Qunmet ≈ 100.0 - CAP

            # Cooling load (negative building load)
            Qg, Qunmet = ground_load(-5.0, 5.0, 35.0, perf, perf)
            @test Qg ≈ 5.0 * (1 + 1 / COP)
            @test Qunmet ≈ 0.0

            # Zero load
            Qg, Qunmet = ground_load(0.0, 5.0, 35.0, perf, perf)
            @test Qg == 0.0 && Qunmet == 0.0

            # Vector overload matches scalar calls elementwise
            Qb = [5.0, 100.0, -5.0]
            Tsrc = [5.0, 5.0, 5.0]
            Tload = [35.0, 35.0, 35.0]
            Qg_v, Qunmet_v = ground_load(Qb, Tsrc, Tload, perf, perf)
            @test Qg_v ≈ [ground_load(Qb[k], Tsrc[k], Tload[k], perf, perf)[1] for k in 1:3]
            @test Qunmet_v ≈ [ground_load(Qb[k], Tsrc[k], Tload[k], perf, perf)[2] for k in 1:3]
        end

        @testset "1D (Tsrc only)" begin
            # Regression test: this form used to raise a MethodError since the 2D
            # heat_pump_performance overload was always called under the hood, regardless of
            # whether perf_heating/perf_cooling were 1D or 2D tuples.
            Tsrc_grid = [0.0, 10.0]
            CAP_vec = [8.0, 10.0]
            COP_vec = [4.0, 5.0]
            perf = (Tsrc_grid, CAP_vec, COP_vec)

            Qg, Qunmet = ground_load(5.0, 5.0, perf, perf)
            CAP, COP = heat_pump_performance(5.0, Tsrc_grid, CAP_vec, COP_vec)
            @test Qg ≈ -5.0 * (1 - 1 / COP)
            @test Qunmet ≈ 0.0

            # Vector overload
            Qb = [5.0, 100.0, -5.0]
            Tsrc = [5.0, 5.0, 5.0]
            Qg_v, Qunmet_v = ground_load(Qb, Tsrc, perf, perf)
            @test Qg_v ≈ [ground_load(Qb[k], Tsrc[k], perf, perf)[1] for k in 1:3]
            @test Qunmet_v ≈ [ground_load(Qb[k], Tsrc[k], perf, perf)[2] for k in 1:3]

            # n_HP scales the CAP looked up from the curve to the total installed capacity: a load
            # beyond a single unit's CAP is fully met once n_HP=2 installed units are considered
            CAP, COP = heat_pump_performance(5.0, Tsrc_grid, CAP_vec, COP_vec)
            Qg1, Qunmet1 = ground_load(15.0, 5.0, perf, perf)
            @test Qunmet1 ≈ 15.0 - CAP # capped at the single unit's CAP
            Qg2, Qunmet2 = ground_load(15.0, 5.0, perf, perf; n_HP=2)
            @test Qunmet2 == 0.0 # below the n_HP=2 installed capacity (2 * CAP)
            @test Qg2 ≈ ground_load(15.0, COP, COP)
        end
    end

    @testset "gshp_simulation" begin
        nt = 48
        t = collect(3600.0:3600.0:(3600.0 * nt))
        Qb = [isodd(k) ? 5000.0 : -3000.0 for k in 1:nt] # W, alternating heating/cooling
        xy = reshape([0.0, 0.0], 1, 2) # single borehole
        rb, H, D = 0.08, 100.0, 2.0
        ks, Cs = 3.0, 2.2e6
        s, ro, ri, kg, kp = 0.05, 0.02, 0.016, 1.5, 0.4
        kf, cf, ρf, μf = 0.6, 4180.0, 1000.0, 0.0013
        V, T0 = 0.0005, 10.0

        @testset "fixed COP" begin
            res = gshp_simulation(t, Qb, 4.0, 5.0, 10_000.0, 12_000.0, 1, xy, rb, H, D, ks, Cs, s,
                ro, ri, kg, kp, kf, cf, ρf, μf, V, T0)
            @test length(res.T_in) == nt
            @test length(res.T_out) == nt
            @test all(0.0 .<= res.n_active .<= 1.0)
            # Below the single unit's capacity everywhere: nothing unmet
            @test all(res.Qb_unmet .== 0.0)
            @test res.Qb_hp == Qb
            @test res.Qg ≈ ground_load(Qb, 4.0, 5.0)

            # The precomputed-g/Rbe method (given the full-input method's own g/Rbe) must reproduce
            # it exactly
            res_low = gshp_simulation(t, Qb, 4.0, 5.0, 10_000.0, 12_000.0, 1, size(xy, 1), H,
                res.g, res.Rbe, cf, ρf, V, T0)
            @test res_low.T_in ≈ res.T_in
            @test res_low.T_out ≈ res.T_out
            @test res_low.Qg == res.Qg
        end

        @testset "variable COP" begin
            perf_heating = ([-5.0, 15.0], [10_000.0, 12_000.0], [3.5, 4.5])
            perf_cooling = ([10.0, 35.0], [11_000.0, 13_000.0], [5.5, 4.0])
            res = gshp_simulation(t, Qb, perf_heating, perf_cooling, 1, xy, rb, H, D, ks, Cs, s,
                ro, ri, kg, kp, kf, cf, ρf, μf, V, T0)
            @test length(res.T_in) == nt
            @test length(res.T_out) == nt
            @test length(res.COP_actual) == nt
            @test all(0.0 .<= res.n_active .<= 1.0)
            @test all(res.Qb_unmet .== 0.0)

            # The fixed-point loop over the source EFT-dependent COP converges well within the
            # default max_iter=5 safety cap (matches fNSCiter.m's identically-capped while loop)
            @test 1 <= res.n_iter < 5

            # At convergence, Qg must be consistent with the curve-based COP at the returned T_out
            # (not the nominal COP used only to bootstrap the first pass) -- this is the actual fix
            # for the "COP never fed back" issue: re-deriving Qg from T_out must reproduce it, up to
            # the loop's own T_out convergence tolerance (Qg was derived from the *previous*
            # iteration's T_out, one step before the returned, converged one).
            Qg_check, _ = ground_load(Qb, res.T_out, perf_heating, perf_cooling)
            @test Qg_check ≈ res.Qg rtol = 1e-3

            # The precomputed-g/Rbe method (given the full-input method's own g/Rbe) must reproduce
            # it exactly, including the iteration count
            res_low = gshp_simulation(t, Qb, perf_heating, perf_cooling, 1, size(xy, 1), H, res.g,
                res.Rbe, cf, ρf, V, T0)
            @test res_low.T_in ≈ res.T_in
            @test res_low.T_out ≈ res.T_out
            @test res_low.n_iter == res.n_iter

            # A tighter tolerance takes at least as many iterations; max_iter=1 forces a single pass
            res_loose = gshp_simulation(t, Qb, perf_heating, perf_cooling, 1, xy, rb, H, D, ks, Cs,
                s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0; max_iter=1)
            @test res_loose.n_iter == 1

            # Regression: the fixed-point loop's ground_load call (method 3) must respect n_HP, not
            # just the pass-0 nominal-COP call. A load beyond a single unit's max CAP (12_000.0) but
            # within n_HP=2 units' combined capacity (>= 2 * 10_000.0, the grid's minimum per-unit
            # CAP) must be unmet with n_HP=1 but fully met with n_HP=2.
            Qb_big = fill(18_000.0, nt)
            res_1unit = gshp_simulation(t, Qb_big, perf_heating, perf_cooling, 1, xy, rb, H, D, ks,
                Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0)
            @test any(!=(0.0), res_1unit.Qb_unmet)
            res_2units = gshp_simulation(t, Qb_big, perf_heating, perf_cooling, 2, xy, rb, H, D, ks,
                Cs, s, ro, ri, kg, kp, kf, cf, ρf, μf, V, T0)
            @test all(res_2units.Qb_unmet .== 0.0)
        end

        @testset "fixed COP, non-stationary flow" begin
            ΔT_target = 4.0 # K

            # A single candidate must reduce exactly to the stationary case
            res_ns = gshp_simulation(t, Qb, 4.0, 5.0, 10_000.0, 12_000.0, 1, xy, rb, H, D, ks, Cs,
                s, ro, ri, kg, kp, kf, cf, ρf, μf, [V], ΔT_target, T0)
            res_s = gshp_simulation(t, Qb, 4.0, 5.0, 10_000.0, 12_000.0, 1, xy, rb, H, D, ks, Cs, s,
                ro, ri, kg, kp, kf, cf, ρf, μf, V, T0)
            @test res_ns.T_in ≈ res_s.T_in
            @test res_ns.T_out ≈ res_s.T_out
            @test all(res_ns.V_t .== V)

            # Multiple candidates: the state actually picked at each step is one of them, and the
            # active flow rate tracks the load (bigger |Qg| -> bigger flow, up to the largest
            # candidate), which is exactly the mechanism that should shrink the T_in/T_out split.
            V_candidates = [0.25, 0.30, 0.35] ./ 1000 # 15/18/21 L/min per borehole
            res = gshp_simulation(t, Qb, 4.0, 5.0, 10_000.0, 12_000.0, 1, xy, rb, H, D, ks, Cs, s,
                ro, ri, kg, kp, kf, cf, ρf, μf, V_candidates, ΔT_target, T0)
            @test length(res.T_in) == nt
            @test length(res.T_out) == nt
            @test all(v -> v in V_candidates, res.V_t)
            @test length(res.Rbe) == length(V_candidates)

            # Same ground loads but forced onto the lowest vs. the highest candidate: the smaller
            # flow rate must show a larger inlet/outlet split for the same load.
            res_low = gshp_simulation(t, Qb, 4.0, 5.0, 10_000.0, 12_000.0, 1, xy, rb, H, D, ks, Cs,
                s, ro, ri, kg, kp, kf, cf, ρf, μf, [V_candidates[1]], ΔT_target, T0)
            res_high = gshp_simulation(t, Qb, 4.0, 5.0, 10_000.0, 12_000.0, 1, xy, rb, H, D, ks, Cs,
                s, ro, ri, kg, kp, kf, cf, ρf, μf, [V_candidates[end]], ΔT_target, T0)
            ΔT_low = abs.(res_low.T_in .- res_low.T_out)
            ΔT_high = abs.(res_high.T_in .- res_high.T_out)
            @test all(ΔT_high .<= ΔT_low .+ 1e-9)
            @test maximum(ΔT_high) < maximum(ΔT_low)
        end
    end
end
