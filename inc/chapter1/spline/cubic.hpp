#pragma once
/// @file cubic.hpp
/// @brief Piecewise cubic spline construction and display.

#include <math_object.hpp>

namespace detail {
	/// @brief Build the system and coefficient history for a cubic spline.
	/// @details `boundary_conditions` is a 2-by-5 augmented matrix. Each row
	///          contains four boundary coefficients followed by its right-hand side.
	/// @param nodes Spline nodes.
	/// @param values Function values corresponding to `nodes`.
	/// @param boundary_conditions Two augmented boundary-condition rows.
	/// @return History containing the system, solution, widths, and interval polynomials.
	auto build_cubic_spline_history(const RealList					  &nodes,
									const RealList					  &values,
									const Eigen::Matrix<double, 2, 5> &boundary_conditions)
		-> SplineHistory;
	/// @brief Construct the piecewise cubic spline.
	/// @param nodes Spline nodes.
	/// @param values Function values corresponding to `nodes`.
	/// @param boundary_conditions Two augmented boundary-condition rows.
	/// @return Piecewise cubic spline in local powers of `(x - nodes[i])`.
	auto find_cubic_spline(const RealList					 &nodes,
						   const RealList					 &values,
						   const Eigen::Matrix<double, 2, 5> &boundary_conditions)
		-> SplineFunction;
	/// @brief Display the augmented system, nodal second derivatives, and coefficients.
	/// @param history Cubic spline history produced by `build_cubic_spline_history`.
	auto render_cubic_spline(const SplineHistory &history) -> void;
} // namespace detail

/// @brief Construct and display a piecewise cubic spline.
/// @param nodes Spline nodes.
/// @param values Function values corresponding to `nodes`.
/// @param boundary_conditions Two augmented boundary-condition rows.
auto display_cubic_spline(const RealList					&nodes,
						  const RealList					&values,
						  const Eigen::Matrix<double, 2, 5> &boundary_conditions) -> void;
