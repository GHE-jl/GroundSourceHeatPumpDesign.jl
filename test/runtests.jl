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

    @testset "Q_COP" begin
        # Heating: Qg = -Qb * (1 - 1/COP_heating)
        @test Q_COP(8.0, 4.0, 5.0) ≈ -8.0 * (1 - 1 / 4.0)
        # Cooling: Qg = -Qb * (1 + 1/COP_cooling)
        @test Q_COP(-8.0, 4.0, 5.0) ≈ 8.0 * (1 + 1 / 5.0)
        # Zero load
        @test Q_COP(0.0, 4.0, 5.0) == 0.0

        # Vector overload, scalar COP
        @test Q_COP([8.0, -8.0], 4.0, 5.0) ≈ [Q_COP(8.0, 4.0, 5.0), Q_COP(-8.0, 4.0, 5.0)]

        # Vector overload, vector COP
        Qb = [8.0, -8.0]
        @test Q_COP(Qb, [4.0, 4.5], [5.0, 5.5]) ≈
              [Q_COP(8.0, 4.0, 5.0), Q_COP(-8.0, 4.5, 5.5)]

        # Non-positive COP is rejected
        @test_throws ArgumentError Q_COP(8.0, 0.0, 5.0)
        @test_throws ArgumentError Q_COP(8.0, 4.0, -1.0)
        @test_throws DimensionMismatch Q_COP([8.0, -8.0], [4.0], [5.0, 5.5])
    end

    @testset "ground_load_capped" begin
        # Below capacity: behaves like Q_COP with the nominal COP
        Qg, Qb_unmet = ground_load_capped(5.0, 4.0, 5.0, 10.0, 10.0)
        @test Qg ≈ Q_COP(5.0, 4.0, 5.0)
        @test Qb_unmet == 0.0

        # Above heating capacity: capped, with unmet remainder
        Qg, Qb_unmet = ground_load_capped(15.0, 4.0, 5.0, 10.0, 10.0)
        @test Qg ≈ Q_COP(10.0, 4.0, 5.0)
        @test Qb_unmet ≈ 5.0

        # Above cooling capacity (Qb very negative): capped, with unmet remainder
        Qg, Qb_unmet = ground_load_capped(-15.0, 4.0, 5.0, 10.0, 10.0)
        @test Qg ≈ Q_COP(-10.0, 4.0, 5.0)
        @test Qb_unmet ≈ -5.0

        # Vector overload
        Qb = [5.0, 15.0, -15.0]
        Qg, Qb_unmet = ground_load_capped(Qb, 4.0, 5.0, 10.0, 10.0)
        @test Qg ≈ [Q_COP(5.0, 4.0, 5.0), Q_COP(10.0, 4.0, 5.0), Q_COP(-10.0, 4.0, 5.0)]
        @test Qb_unmet ≈ [0.0, 5.0, -5.0]
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

    @testset "ground_load_from_heat_pump" begin
        Tsrc_grid = [0.0, 10.0]
        Tload_grid = [30.0, 40.0]
        CAP_table = [8.0 6.0; 10.0 8.0]
        COP_table = [4.0 3.0; 5.0 4.0]
        perf = (Tsrc_grid, Tload_grid, CAP_table, COP_table)

        # Heating load fully met by the heat pump (below capacity)
        Qg, Qunmet = ground_load_from_heat_pump(5.0, 5.0, 35.0, perf, perf)
        CAP, COP = heat_pump_performance(5.0, 35.0, Tsrc_grid, Tload_grid, CAP_table, COP_table)
        @test Qg ≈ -5.0 * (1 - 1 / COP)
        @test Qunmet ≈ 0.0

        # Heating load exceeding capacity: capped, with unmet remainder
        Qg, Qunmet = ground_load_from_heat_pump(100.0, 5.0, 35.0, perf, perf)
        @test Qg ≈ -CAP * (1 - 1 / COP)
        @test Qunmet ≈ 100.0 - CAP

        # Cooling load (negative building load)
        Qg, Qunmet = ground_load_from_heat_pump(-5.0, 5.0, 35.0, perf, perf)
        @test Qg ≈ 5.0 * (1 + 1 / COP)
        @test Qunmet ≈ 0.0

        # Zero load
        Qg, Qunmet = ground_load_from_heat_pump(0.0, 5.0, 35.0, perf, perf)
        @test Qg == 0.0 && Qunmet == 0.0
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
            @test res.Qg ≈ Q_COP(Qb, 4.0, 5.0)
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
        end
    end
end
