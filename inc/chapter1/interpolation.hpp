#pragma once
/// @file interpolation.hpp
/// @brief Chebyshev interpolation nodes.

#include <math_object.hpp>

namespace detail {

	/// @brief Generate Chebyshev nodes on an interval.
	/// @param lower Left endpoint.
	/// @param upper Right endpoint.
	/// @param num_nodes Positive node count.
	/// @return Node list of size `num_nodes`.
	auto generate_chebyshev_nodes(Real lower, Real upper, SizeType num_nodes) -> RealList;
	/// @brief Render the `i | x_i` table.
	/// @param nodes Node list to display.
	auto render_nodes_table(const RealList &nodes) -> void;

} // namespace detail

/// @brief Render the `i | x_i` table.
/// @param lower Left endpoint.
/// @param upper Right endpoint.
/// @param num_nodes Positive node count.
auto display_chebyshev_nodes(Real lower, Real upper, SizeType num_nodes) -> void;
