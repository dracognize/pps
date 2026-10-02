#pragma once

#include <math_object.hpp>

namespace detail {

	auto generate_chebyshev_nodes(Real lower, Real upper, SizeType num_nodes) -> RealList;
	auto render_nodes_table(const RealList &nodes) -> void;

} // namespace detail

auto display_chebyshev_nodes(Real lower, Real upper, SizeType num_nodes) -> void;
