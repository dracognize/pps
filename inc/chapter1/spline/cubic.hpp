#pragma once
/// @file cubic.hpp
/// @brief Piecewise cubic spline construction and display.

#include <math_object.hpp>

namespace detail {
	auto build_cubic_spline_history(const RealList &nodes,
									const RealList &values,
									const Eigen::Matrix<double, 2, 5> &boundary_conditions)
		-> SplineHistory;
	auto find_cubic_spline(const RealList &nodes,
						const RealList &values,
						const Eigen::Matrix<double, 2, 5> &boundary_conditions) -> SplineFunction;
	auto render_cubic_spline(const SplineHistory &history) -> void;
} // namespace detail

auto display_cubic_spline(const RealList &nodes,
						  const RealList &values,
						  const Eigen::Matrix<double, 2, 5> &boundary_conditions) -> void;
