#pragma once
/// @file Newton.hpp
/// @brief Newton interpolation polynomials.

#include <math_object.hpp>

namespace detail {

	/// @brief Build the forward history for general nodes.
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation history; empty when invalid.
	auto build_Newton_general_forward_history(const RealList &nodes, const RealList &values)
		-> NewtonInterpolationHistory;

	/// @brief Compute the Newton forward polynomial for general nodes.
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation polynomial; empty when invalid.
	auto Newton_general_forward_polynomial(const RealList &nodes, const RealList &values)
		-> Polynomial;

	/// @brief Render the forward history for general nodes.
	/// @param nodes Interpolation nodes.
	/// @param history Precomputed forward history.
	auto render_Newton_general_forward(const RealList					&nodes,
									   const NewtonInterpolationHistory &history) -> void;

	/// @brief Build the backward history for general nodes.
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation history; empty when invalid.
	auto build_Newton_general_backward_history(const RealList &nodes, const RealList &values)
		-> NewtonInterpolationHistory;

	/// @brief Compute the Newton backward polynomial for general nodes.
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation polynomial; empty when invalid.
	auto Newton_general_backward_polynomial(const RealList &nodes, const RealList &values)
		-> Polynomial;

	/// @brief Render the backward history for general nodes.
	/// @param nodes Interpolation nodes.
	/// @param history Precomputed backward history.
	auto render_Newton_general_backward(const RealList					 &nodes,
										const NewtonInterpolationHistory &history) -> void;

	/// @brief Build the forward history for equidistant nodes.
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation history; empty when invalid.
	auto build_Newton_equidistant_forward_history(const RealList &nodes, const RealList &values)
		-> EquidistantHistory;

	/// @brief Compute the Newton forward polynomial for equidistant nodes.
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation polynomial; empty when invalid.
	auto Newton_equidistant_forward_polynomial(const RealList &nodes, const RealList &values)
		-> Polynomial;

	/// @brief Render the forward history for equidistant nodes.
	/// @param nodes Interpolation nodes.
	/// @param history Precomputed forward history.
	auto render_Newton_equidistant_forward(const RealList &nodes, const EquidistantHistory &history)
		-> void;

	/// @brief Build the backward history for equidistant nodes.
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation history; empty when invalid.
	auto build_Newton_equidistant_backward_history(const RealList &nodes, const RealList &values)
		-> EquidistantHistory;

	/// @brief Compute the Newton backward polynomial for equidistant nodes.
	/// @param nodes Interpolation nodes.
	/// @param values Function values at nodes.
	/// @return The interpolation polynomial; empty when invalid.
	auto Newton_equidistant_backward_polynomial(const RealList &nodes, const RealList &values)
		-> Polynomial;

	/// @brief Render the backward history for equidistant nodes.
	/// @param nodes Interpolation nodes.
	/// @param history Precomputed backward history.
	auto render_Newton_equidistant_backward(const RealList			 &nodes,
											const EquidistantHistory &history) -> void;
} // namespace detail

/// @brief Display the Newton general forward interpolation as plain text and tables.
/// @param nodes Interpolation nodes.
/// @param values Function values at nodes.
auto display_Newton_general_forward(const RealList &nodes, const RealList &values) -> void;

/// @brief Display the Newton general backward interpolation as plain text and tables.
/// @param nodes Interpolation nodes.
/// @param values Function values at nodes.
auto display_Newton_general_backward(const RealList &nodes, const RealList &values) -> void;

/// @brief Display the Newton equidistant forward interpolation as plain text and tables.
/// @param nodes Interpolation nodes.
/// @param values Function values at nodes.
auto display_Newton_equidistant_forward(const RealList &nodes, const RealList &values) -> void;

/// @brief Display the Newton equidistant backward interpolation as plain text and tables.
/// @param nodes Interpolation nodes.
/// @param values Function values at nodes.
auto display_Newton_equidistant_backward(const RealList &nodes, const RealList &values) -> void;