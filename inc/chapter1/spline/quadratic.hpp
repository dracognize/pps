#pragma once
/// @file quadratic.hpp
/// @brief Piecewise quadratic spline construction and display.

#include <math_object.hpp>

namespace detail {
	/// @brief Construct a C1 quadratic spline through the supplied data points.
	/// @details The endpoint derivatives satisfy `alpha*S'_0(x_0) + beta*S'_{n-1}(x_n) = gamma`,
	///          with `beta = (-1)^n * alpha` and nonzero `alpha`.
	/// @param nodes Strictly increasing interpolation nodes.
	/// @param values Function values corresponding to `nodes`.
	/// @return A piecewise quadratic spline passing through every supplied point.
	/// @throws std::invalid_argument if there are fewer than two points, the input
	///         sizes differ, or the nodes are not strictly increasing.
	auto build_quadratic_spline_history(const RealList &nodes,
										const RealList &values,
										Real			alpha,
										Real			beta,
										Real			gamma) -> SplineHistory;
	auto find_quadratic_spline(const RealList &nodes,
							   const RealList &values,
							   Real			   alpha,
							   Real			   beta,
							   Real			   gamma) -> SplineFunction;

	/// @brief Render the interval and coefficients of each quadratic segment.
	/// @param spline Piecewise quadratic spline returned by `find_quadratic_spline`.
	auto render_quadratic_spline(const SplineHistory &history) -> void;
} // namespace detail

/// @brief Construct and display a C1 quadratic spline under an endpoint derivative condition.
/// @param nodes Strictly increasing interpolation nodes.
/// @param values Function values corresponding to `nodes`.
/// @param alpha Coefficient of the derivative at the first endpoint.
/// @param beta Coefficient of the derivative at the last endpoint; must equal `(-1)^n * alpha`.
/// @param gamma Prescribed linear combination of endpoint derivatives.
/// @throws std::invalid_argument if the input does not satisfy the requirements
///         of `detail::find_quadratic_spline`.
auto display_quadratic_spline(
	const RealList &nodes, const RealList &values, Real alpha, Real beta, Real gamma) -> void;
