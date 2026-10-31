#pragma once
/// @file Stirling.hpp
/// @brief Stirling interpolation polynomials.

#include <math_object.hpp>

namespace detail {
	/// @brief Build the full Stirling interpolation history.
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation history; empty when invalid.
	auto build_Stirling_history(const RealList &nodes, const RealList &values)
		-> EquidistantHistory;

	/// @brief Compute the Stirling polynomial for equidistant nodes.
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation polynomial; empty when invalid.
	auto Stirling_polynomial(const RealList &nodes, const RealList &values) -> Polynomial;

	/// @brief Render the Stirling interpolation history.
	/// @param history Precomputed forward history.
	auto render_Stirling_interpolation(const EquidistantHistory &history) -> void;
} // namespace detail

/// @brief Display the Stirling interpolation as plain text and tables.
/// @param nodes Interpolation nodes.
/// @param values Function values at nodes.
auto display_Stirling_interpolation(const RealList &nodes, const RealList &values) -> void;