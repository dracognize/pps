#pragma once
/// '@file Lagrange.hpp
/// @brief Lagrange interpolation polynomials.

#include <math_object.hpp>

namespace detail {
	/// @brief Build the full Lagrange interpolation history.
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return Node polynomials, basis polynomials, normalizers, and values.
	auto build_Lagrange_interpolation_history(const RealList &nodes, const RealList &values)
		-> LagrangeInterpolationHistory;
	/// @brief 	Compute the Lagrange interpolation polynomial from a precomputed history.
	/// @param history Lagrange interpolation history.
	/// @return The Lagrange interpolation polynomial.
	auto Lagrange_polynomial(const LagrangeInterpolationHistory &history) -> Polynomial;
	/// @brief 	Compute the Lagrange interpolation polynomial from nodes and values.
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The Lagrange interpolation polynomial.
	auto Lagrange_polynomial(const RealList &nodes, const RealList &values) -> Polynomial;
	/// @brief 	Render the Lagrange interpolation history.
	/// @param nodes Interpolation nodes.
	/// @param history Lagrange interpolation history.
	auto render_Lagrange_interpolation(const RealList					  &nodes,
									   const LagrangeInterpolationHistory &history) -> void;
} // namespace detail

/// @brief 	Display the Lagrange interpolation as plain text and tables.
/// @param nodes Interpolation nodes.
/// @param values Function values at nodes.
auto display_Lagrange_interpolation(const RealList &nodes, const RealList &values) -> void;