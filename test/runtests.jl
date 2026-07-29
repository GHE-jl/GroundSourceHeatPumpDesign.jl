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
end
