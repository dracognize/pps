#pragma once
/// @file interval.hpp
/// @brief Select interpolation samples around a query point.

#include <math_object.hpp>

namespace detail {

	/// @brief Select a requested number of samples from equidistant nodes.
	/// @param nodes Strictly increasing or decreasing nodes; at least two entries.
	/// @param values Sample values, one for each node.
	/// @param x Query point within the node range.
	/// @param point_count Number of consecutive samples to select.
	/// @return Inclusive start and end indices; `{nodes.size(), nodes.size()}` for invalid input.
	auto select_equidistant_interpolation_points(const RealList &nodes,
												 const RealList &values,
												 Real			 x,
												 SizeType		 point_count) -> ExtractInterval;

	/// @brief Select a requested number of samples from monotonic nodes.
	/// @param nodes Strictly increasing or decreasing nodes.
	/// @param values Sample values, one for each node.
	/// @param x Query point within the node range.
	/// @param point_count Number of consecutive samples to select.
	/// @return Inclusive start and end indices; `{nodes.size(), nodes.size()}` for invalid input.
	auto select_interpolation_points(const RealList &nodes,
									 const RealList &values,
									 Real			 x,
									 SizeType		 point_count) -> ExtractInterval;

	/// @brief Render the selected interpolation interval and its points.
	/// @param nodes Full node list.
	/// @param values Corresponding sample values.
	/// @param interval Inclusive range of selected indices.
	auto render_selection_points(const RealList		   &nodes,
								 const RealList		   &values,
								 const ExtractInterval &interval) -> void;

	/// @brief Find index intervals where the sample values are monotonic.
	/// @param values Sample values in node order.
	/// @return Inclusive index intervals of monotonic segments.
	auto find_monotonic_intervals(const RealList &values) -> IntervalsList;

	/// @brief Find monotonic intervals that bracket a target value.
	/// @param values Sample values in node order.
	/// @param target_value Value to bracket.
	/// @return Inclusive index intervals whose endpoint values bracket the target.
	auto find_isolating_intervals(const RealList &values, Real target_value) -> IntervalsList;

	/// @brief Render index intervals with their endpoint values.
	/// @param nodes Sample nodes in node order.
	/// @param intervals Intervals to display.
	/// @param title Table title.
	auto render_value_intervals(const RealList		&nodes,
								const IntervalsList &intervals,
								const char			*title) -> void;

} // namespace detail

/// @brief Display the selected interpolation interval for equidistant nodes.
/// @param nodes Equidistant nodes in increasing or decreasing order.
/// @param values Sample values, one for each node.
/// @param x Query point within the node range.
/// @param point_count Number of consecutive samples to select.
auto display_equidistant_interpolation_points(const RealList &nodes,
											  const RealList &values,
											  Real			  x,
											  SizeType		  point_count) -> void;

/// @brief Display all monotonic intervals in the sample values.
/// @param nodes Samples nodes in node order.
/// @param values Sample values in node order.
auto display_monotonic_intervals(const RealList &nodes, const RealList &values) -> void;

/// @brief Display monotonic intervals that bracket a target value.
/// @param nodes Samples nodes in node order.
/// @param values Sample values in node order.
/// @param target_value Value to bracket.
auto display_isolating_intervals(const RealList &nodes, const RealList &values, Real target_value)
	-> void;
