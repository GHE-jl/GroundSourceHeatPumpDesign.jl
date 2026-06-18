"""
Script to test the heat pump capacity (CAP) and coefficient of performance (COP) calculations.
"""

includet("../src/thermal_load_analysis.jl")

Tsrc_grid = [0.0, 5.0, 10.0]
CAP_vec   = [140.0, 160.0, 180.0]
COP_vec   = [3.2, 3.8, 4.2]

# Single value
CAP, COP = hp_performance_1D(7.0, Tsrc_grid, CAP_vec, COP_vec)

# Time series
Tsrc_vec = [-2.0, 3.0, 6.0, 12.0]  # includes out-of-range values
CAP_vec_out, COP_vec_out = hp_performance_series_1D(
    Tsrc_vec, Tsrc_grid, CAP_vec, COP_vec
)