#pragma once
/// @file Gauss.hpp
/// @brief Gauss interpolation polynomials.

#include <math_object.hpp>

namespace detail {
	/// @brief Build the full Gauss I interpolation history.
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation history; empty when invalid.
	auto build_GaussI_history(const RealList &nodes, const RealList &values) -> EquidistantHistory;

	/// @brief Compute the Gauss I polynomial for equidistant nodes (interpolation parameter: t = (x
	/// - xn)/h).
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation polynomial; empty when invalid.
	auto GaussI_polynomial(const RealList &nodes, const RealList &values) -> Polynomial;

	/// @brief Render the Gauss I interpolation history.
	/// @param history Precomputed forward history.
	auto render_GaussI_interpolation(const EquidistantHistory &history) -> void;

	/// @brief Build the full Gauss II interpolation history.
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation history; empty when invalid.
	auto build_GaussII_history(const RealList &nodes, const RealList &values) -> EquidistantHistory;

	/// @brief Compute the Gauss II polynomial for equidistant nodes (interpolation parameter: t =
	/// (x - xn)/h).
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation polynomial; empty when invalid.
	auto GaussII_polynomial(const RealList &nodes, const RealList &values) -> Polynomial;

	/// @brief Render the Gauss II interpolation history.
	/// @param history Precomputed forward history.
	auto render_GaussII_interpolation(const EquidistantHistory &history) -> void;
} // namespace detail

/// @brief Display the Gauss I interpolation as plain text and tables.
/// @param nodes Interpolation nodes.
/// @param values Function values at nodes.
auto display_GaussI_interpolation(const RealList &nodes, const RealList &values) -> void;

/// @brief Display the Gauss II interpolation as plain text and tables.
/// @param nodes Interpolation nodes.
/// @param values Function values at nodes.
auto display_GaussII_interpolation(const RealList &nodes, const RealList &values) -> void;