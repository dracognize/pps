#pragma once
/// @file Bessel.hpp
/// @brief Bessel interpolation polynomials.

#include <math_object.hpp>

namespace detail {
	/// @brief Build the full Bessel interpolation history.
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation history; empty when invalid.
	auto build_Bessel_history(const RealList &nodes, const RealList &values) -> EquidistantHistory;

	/// @brief Compute the Bessel polynomial for equidistant nodes (interpolation parameter: u = (x
	/// - x(n-1))/h - 0.5).
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation polynomial; empty when invalid.
	auto Bessel_polynomial(const RealList &nodes, const RealList &values) -> Polynomial;

	/// @brief Render the Bessel interpolation history.
	/// @param history Precomputed forward history.
	auto render_Bessel_interpolation(const EquidistantHistory &history) -> void;
} // namespace detail

/// @brief Display the Bessel interpolation as plain text and tables.
/// @param nodes Interpolation nodes.
/// @param values Function values at nodes.
auto display_Bessel_interpolation(const RealList &nodes, const RealList &values) -> void;