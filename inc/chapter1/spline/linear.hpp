#pragma once
/// @file linear.hpp
/// @brief Piecewise linear spline construction.

#include <math_object.hpp>

namespace detail {
	/// @brief Construct the degree-1 spline through the supplied data points.
	/// @details Each interval polynomial is stored in powers of `(x - nodes[i])`.
	/// @param nodes Strictly increasing interpolation nodes.
	/// @param values Function values corresponding to `nodes`.
	/// @return A piecewise linear spline passing through every supplied point.
	/// @throws std::invalid_argument if there are fewer than two points, the input
	///         sizes differ, or the nodes are not strictly increasing.
	auto find_linear_spline(const RealList &nodes, const RealList &values) -> SplineFunction;

	/// @brief Render each interval and its degree-1 polynomial coefficients.
	/// @details Each row represents `p_i(x) = a_i + b_i * (x - nodes[i])`.
	/// @param spline Piecewise linear spline returned by `find_linear_spline`.
	auto render_linear_spline(const SplineFunction &spline) -> void;

	/// @brief Construct and display the degree-1 spline coefficient table.
	/// @param nodes Strictly increasing interpolation nodes.
	/// @param values Function values corresponding to `nodes`.
	/// @throws std::invalid_argument if the input does not satisfy the requirements
	///         of `find_linear_spline`.
}; // namespace detail

auto display_linear_spline(const RealList &nodes, const RealList &values) -> void;