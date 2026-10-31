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

/// @brief History of iterative inverse interpolation.
struct InverseIterationHistory {
	IntervalsList		   isolating_intervals;	 ///< Isolating intervals [p_i, p_i + 1].
	std::vector<bool>	   is_forward;			 ///< True if Newton forward was used, false for backward.
	RealTable			   interpolation_nodes;	 ///< Interpolation nodes used for each interval.
	RealTable			   interpolation_values; ///< Interpolation values used for each interval.
	PolynomialList		   polynomials;			 ///< Interpolation polynomials P_i(t).
	PolynomialList		   phi_functions;		 ///< Iteration functions phi(t).
	std::vector<RealTable> iteration_tables;	 ///< Iteration rows [iter, t_0, t, delta] for each interval.
	RealList			   interval_roots;		 ///< Root found in each isolating interval.
	RealList			   roots;				 ///< All distinct roots found.
};

namespace detail {

	/// @brief Find roots using iterative inverse interpolation.
	/// @param nodes Equidistant independent-variable samples.
	/// @param values Function values at the samples.
	/// @param target_value Target function value \bar{y}.
	/// @param tolerance Error tolerance \varepsilon.
	/// @param max_points Maximum number of interpolation nodes m.
	/// @return Roots found.
	auto inverse_iteration(const RealList &nodes,
						   const RealList &values,
						   Real			   target_value,
						   Real			   tolerance,
						   SizeType		   max_points) -> RealList;

	/// @brief Build history for iterative inverse interpolation.
	/// @param nodes Equidistant independent-variable samples.
	/// @param values Function values at the samples.
	/// @param target_value Target function value \bar{y}.
	/// @param tolerance Error tolerance \varepsilon.
	/// @param max_points Maximum number of interpolation nodes m.
	/// @return History of the iterative inverse interpolation.
	auto build_inverse_iteration_history(const RealList &nodes,
										 const RealList &values,
										 Real			 target_value,
										 Real			 tolerance,
										 SizeType		 max_points)
		-> InverseIterationHistory;

	/// @brief Render iterative inverse interpolation history.
	/// @param nodes Independent-variable samples.
	/// @param values Function values at the samples.
	/// @param target_value Target function value \bar{y}.
	/// @param tolerance Error tolerance \varepsilon.
	/// @param history Precomputed iterative history.
	auto render_inverse_iteration_history(
		const RealList				  &nodes,
		const RealList				  &values,
		Real						   target_value,
		const InverseIterationHistory &history) -> void;

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

/// @brief Display iterative inverse interpolation.
/// @param nodes Equidistant independent-variable samples.
/// @param values Function values at the samples.
/// @param target_value Target function value \bar{y}.
/// @param tolerance Error tolerance \varepsilon.
/// @param max_points Maximum number of interpolation nodes m.
auto display_inverse_iteration(const RealList &nodes,
							   const RealList &values,
							   Real			   target_value,
							   Real			   tolerance,
							   SizeType		   max_points) -> void;
