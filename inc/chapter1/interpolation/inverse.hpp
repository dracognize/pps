#pragma once
/// @file inverse.hpp
/// @brief Inverse interpolation methods.

#include <math_object.hpp>

/// @brief Available methods for constructing the inverse interpolation polynomial.
enum class InterpolationMethod { Lagrange, NewtonForward, NewtonBackward };

namespace detail {

	/// @brief Build inverse-interpolation history using the selected polynomial method.
	auto build_inverse_interpolation_history(const RealList		&nodes,
											 const RealList		&values,
											 Real				 target_value,
											 SizeType			 max_points,
											 InterpolationMethod method)
		-> InverseInterpolationHistory;

	/// @brief Render isolating intervals, each interpolation step, and the roots.
	auto render_inverse_interpolation_history(const RealList					&nodes,
											  const RealList					&values,
											  Real								 target_value,
											  const InverseInterpolationHistory &history) -> void;

	/// @brief Find roots of the tabulated function using inverse interpolation.
	auto inverse_interpolation(const RealList	  &nodes,
							   const RealList	  &values,
							   Real				   target_value,
							   SizeType			   max_points,
							   InterpolationMethod method) -> RealList;

} // namespace detail

/// @brief Display inverse interpolation using the selected polynomial method.
/// @param nodes Independent-variable samples.
/// @param values Function values at the samples.
/// @param target_value Target function value.
/// @param max_points Maximum number of samples used per inverse interpolant.
/// @param method Interpolation method used to build each inverse polynomial.
auto display_inverse_interpolation(const RealList	  &nodes,
								   const RealList	  &values,
								   Real				   target_value,
								   SizeType			   max_points,
								   InterpolationMethod method = InterpolationMethod::NewtonForward)
	-> void;
