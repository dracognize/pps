#include <chapter1/interpolation/interval.hpp>
#include <config.hpp>
#include <io/table.hpp>

#include <cmath>
#include <string>
#include <vector>

#include <fmt/format.h>

namespace {

	auto is_valid_input(const RealList &nodes,
						const RealList &values,
						Real			center_point,
						SizeType		point_count) -> bool {
		if (nodes.size() != values.size() || nodes.size() < 2 || point_count == 0 ||
			point_count > nodes.size() || !std::isfinite(center_point))
			return false;

		const bool increasing = nodes.front() < nodes.back();
		for (SizeType i = 1; i < nodes.size(); ++i)
			if (increasing ? nodes[i] <= nodes[i - 1] : nodes[i] >= nodes[i - 1])
				return false;

		return increasing ? center_point >= nodes.front() && center_point <= nodes.back()
						  : center_point <= nodes.front() && center_point >= nodes.back();
	}

	auto select_range(SizeType center_index, SizeType point_count, SizeType node_count)
		-> ExtractInterval {
		auto first_index = static_cast<std::ptrdiff_t>(center_index) -
						   (static_cast<std::ptrdiff_t>(point_count) - 1) / 2;
		auto last_index	 = static_cast<std::ptrdiff_t>(center_index) +
						   static_cast<std::ptrdiff_t>(point_count) / 2;

		if (first_index < 0) {
			last_index += -first_index;
			first_index = 0;
		}
		const auto last_node_index = static_cast<std::ptrdiff_t>(node_count) - 1;
		if (last_index > last_node_index) {
			first_index -= last_index - last_node_index;
			last_index = last_node_index;
		}

		return {static_cast<SizeType>(first_index), static_cast<SizeType>(last_index)};
	}

} // namespace

namespace detail {

	auto select_equidistant_interpolation_points(const RealList &nodes,
												 const RealList &values,
												 Real			 x,
												 SizeType		 point_count) -> ExtractInterval {
		if (!is_valid_input(nodes, values, x, point_count))
			return {nodes.size(), nodes.size()};

		const Real h = nodes[1] - nodes[0];
		const auto p = static_cast<SizeType>((x - nodes[0]) / h);
		return select_range(p, point_count, nodes.size());
	}

	auto select_interpolation_points(const RealList &nodes,
									 const RealList &values,
									 Real			 x,
									 SizeType		 point_count) -> ExtractInterval {
		if (!is_valid_input(nodes, values, x, point_count))
			return {nodes.size(), nodes.size()};

		const SizeType n = nodes.size() - 1;
		SizeType	   p = 0;
		if (nodes[0] < nodes[n]) {
			while (p < n && x > nodes[p + 1])
				++p;
		} else if (nodes[0] > nodes[n]) {
			while (p < n && x < nodes[p + 1])
				++p;
		}

		return select_range(p, point_count, nodes.size());
	}

	auto render_selection_points(const RealList		   &nodes,
								 const RealList		   &values,
								 const ExtractInterval &interval) -> void {
		if (nodes.size() != values.size() || interval.start > interval.end ||
			interval.end >= nodes.size())
			return;

		fmt::print("Selected interval: [{:.{}g}, {:.{}g}]\n",
				   nodes[interval.start],
				   config::kDisplayPrecision,
				   nodes[interval.end],
				   config::kDisplayPrecision);

		std::vector<std::vector<std::string>> rows{{"i", "x_i", "y_i"}};
		for (SizeType i = interval.start; i <= interval.end; ++i) {
			rows.push_back({fmt::format("{}", i),
							fmt::format("{:.{}g}", nodes[i], config::kDisplayPrecision),
							fmt::format("{:.{}g}", values[i], config::kDisplayPrecision)});
		}

		ui::render_table(rows, config::kTableMinWidth);
	}

	auto find_monotonic_intervals(const RealList &values) -> IntervalsList {
		IntervalsList intervals;
		if (values.size() < 2)
			return intervals;

		const SizeType last_index  = values.size() - 1;
		SizeType	   start_index = 0;
		while (start_index < last_index) {
			SizeType end_index = start_index + 1;
			if (values[start_index] == values[end_index]) {
				start_index = end_index;
				continue;
			}

			while (end_index < last_index) {
				const Real previous_difference = values[end_index] - values[end_index - 1];
				const Real next_difference	   = values[end_index + 1] - values[end_index];
				if (previous_difference * next_difference <= 0)
					break;
				++end_index;
			}

			intervals.push_back({start_index, end_index});
			start_index = end_index;
		}

		return intervals;
	}

	auto find_isolating_intervals(const RealList &values, Real target_value) -> IntervalsList {
		IntervalsList isolating_intervals;
		for (const ExtractInterval &interval : find_monotonic_intervals(values)) {
			const Real start_difference = values[interval.start] - target_value;
			const Real end_difference	= values[interval.end] - target_value;
			if (start_difference * end_difference <= 0)
				isolating_intervals.push_back(interval);
		}
		return isolating_intervals;
	}

	auto render_value_intervals(const RealList		&nodes,
								const IntervalsList &intervals,
								const char			*title) -> void {
		SizeType interval_count = 0;

		fmt::println("{}", title);
		std::vector<std::vector<std::string>> rows{
			{"Interval", "p_i", "q_i", "x_(p_i)", "x_(q_i)"}};
		for (const ExtractInterval &interval : intervals) {
			if (interval.start > interval.end || interval.end >= nodes.size())
				continue;
			rows.push_back(
				{fmt::format("i = {}", interval_count++),
				 fmt::format("{}", interval.start),
				 fmt::format("{}", interval.end),
				 fmt::format("{:.{}g}", nodes[interval.start], config::kDisplayPrecision),
				 fmt::format("{:.{}g}", nodes[interval.end], config::kDisplayPrecision)});
		}

		auto transform_table = [](ftxui::Table &table) {
			table.SelectColumn(0).DecorateCells(color(ftxui::Color::Yellow));
			table.SelectColumn(0).DecorateCells(ftxui::bold);
		};

		ui::render_table(rows, config::kTableMinWidth, transform_table);
	}

} // namespace detail

auto display_equidistant_interpolation_points(const RealList &nodes,
											  const RealList &values,
											  Real			  x,
											  SizeType		  point_count) -> void {
	const ExtractInterval interval =
		detail::select_equidistant_interpolation_points(nodes, values, x, point_count);
	detail::render_selection_points(nodes, values, interval);
}

auto display_monotonic_intervals(const RealList &nodes, const RealList &values) -> void {
	const IntervalsList intervals = detail::find_monotonic_intervals(values);
	detail::render_value_intervals(nodes, intervals, "Monotonic intervals");
}

auto display_isolating_intervals(const RealList &nodes, const RealList &values, Real target_value)
	-> void {
	const IntervalsList intervals = detail::find_isolating_intervals(values, target_value);
	detail::render_value_intervals(nodes, intervals, "Isolating intervals");
}
