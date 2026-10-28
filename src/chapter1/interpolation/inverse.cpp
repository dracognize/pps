#include <chapter1/interpolation/horner.hpp>
#include <chapter1/interpolation/interval.hpp>
#include <chapter1/interpolation/inverse.hpp>
#include <chapter1/interpolation/polynomial/Lagrange.hpp>
#include <chapter1/interpolation/polynomial/Newton.hpp>
#include <config.hpp>
#include <io/table.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace detail {

	auto build_inverse_interpolation_history(const RealList		&nodes,
											 const RealList		&values,
											 Real				 target_value,
											 SizeType			 max_points,
											 InterpolationMethod method)
		-> InverseInterpolationHistory {
		InverseInterpolationHistory history;
		const bool has_invalid_input = nodes.size() != values.size() || nodes.size() < 2 ||
									   max_points == 0 || !std::isfinite(target_value);
		if (has_invalid_input)
			return history;

		history.isolating_intervals = find_isolating_intervals(values, target_value);
		for (SizeType i = 0; i < values.size(); ++i) {
			const bool is_inside_isolating_interval =
				std::any_of(history.isolating_intervals.begin(),
							history.isolating_intervals.end(),
							[i](const ExtractInterval &interval) {
								return interval.start <= i && i <= interval.end;
							});
			if (!is_inside_isolating_interval && values[i] == target_value)
				history.roots.push_back(nodes[i]);
		}

		for (const ExtractInterval &isolating_interval : history.isolating_intervals) {
			if (values[isolating_interval.start] == target_value ||
				values[isolating_interval.end] == target_value) {
				if (values[isolating_interval.start] == target_value)
					history.roots.push_back(nodes[isolating_interval.start]);
				if (values[isolating_interval.end] == target_value)
					history.roots.push_back(nodes[isolating_interval.end]);
				continue;
			}

			const SizeType interval_size = isolating_interval.end - isolating_interval.start + 1;
			const SizeType point_count	 = std::min(max_points, interval_size);
			const RealList interval_values(
				values.begin() + static_cast<std::ptrdiff_t>(isolating_interval.start),
				values.begin() + static_cast<std::ptrdiff_t>(isolating_interval.end + 1));
			const RealList interval_nodes(
				nodes.begin() + static_cast<std::ptrdiff_t>(isolating_interval.start),
				nodes.begin() + static_cast<std::ptrdiff_t>(isolating_interval.end + 1));

			const ExtractInterval selected_interval = select_interpolation_points(
				interval_values, interval_nodes, target_value, point_count);
			const bool has_invalid_selection = selected_interval.start >= interval_size ||
											   selected_interval.end >= interval_size ||
											   selected_interval.start > selected_interval.end;
			if (has_invalid_selection)
				continue;

			InverseInterpolationStep step;
			step.isolating_interval = isolating_interval;
			step.selected_interval	= {isolating_interval.start + selected_interval.start,
									   isolating_interval.start + selected_interval.end};
			step.inverse_nodes		= RealList(
				interval_values.begin() + static_cast<std::ptrdiff_t>(selected_interval.start),
				interval_values.begin() + static_cast<std::ptrdiff_t>(selected_interval.end + 1));
			step.inverse_values = RealList(
				interval_nodes.begin() + static_cast<std::ptrdiff_t>(selected_interval.start),
				interval_nodes.begin() + static_cast<std::ptrdiff_t>(selected_interval.end + 1));
			switch (method) {
			case InterpolationMethod::Lagrange:
				step.inverse_polynomial =
					Lagrange_polynomial(step.inverse_nodes, step.inverse_values);
				break;
			case InterpolationMethod::NewtonForward:
				step.inverse_polynomial =
					Newton_general_forward_polynomial(step.inverse_nodes, step.inverse_values);
				break;
			case InterpolationMethod::NewtonBackward:
				step.inverse_polynomial =
					Newton_general_backward_polynomial(step.inverse_nodes, step.inverse_values);
				break;
			}
			step.root = evaluate_polynomial(target_value, step.inverse_polynomial);
			history.roots.push_back(step.root);
			history.steps.push_back(std::move(step));
		}

		std::sort(history.roots.begin(), history.roots.end());
		history.roots.erase(std::unique(history.roots.begin(), history.roots.end()),
							history.roots.end());
		return history;
	}

	auto render_inverse_interpolation_history(const RealList					&nodes,
											  const RealList					&values,
											  Real								 target_value,
											  const InverseInterpolationHistory &history) -> void {
		if (nodes.size() != values.size())
			return;

		fmt::println(
			"Isolating intervals for y = {:.{}g}:", target_value, config::kDisplayPrecision);
		render_value_intervals(values, history.isolating_intervals, "");

		auto transform_table = [](ftxui::Table &table) {
			table.SelectColumn(0).DecorateCells(color(ftxui::Color::Yellow));
			table.SelectColumn(0).DecorateCells(ftxui::bold);
		};

		for (SizeType i = 0; i < history.steps.size(); ++i) {
			const InverseInterpolationStep &step = history.steps[i];
			fmt::println("\n------------------------------------------------------------");
			fmt::println("\nInterval {}: [start..end] = [{}..{}]",
						 i + 1,
						 step.isolating_interval.start,
						 step.isolating_interval.end);
			fmt::println("Selection:     [start..end] = [{}..{}], count = {}",
						 step.selected_interval.start,
						 step.selected_interval.end,
						 step.inverse_nodes.size());

			std::vector<std::vector<std::string>> rows{{"i", "y_i", "x_i"}};
			for (SizeType j = 0; j < step.inverse_nodes.size(); ++j) {
				const SizeType source_index = step.selected_interval.start + j;
				rows.push_back(
					{fmt::format("{}", source_index),
					 fmt::format("{:.{}g}", step.inverse_nodes[j], config::kDisplayPrecision),
					 fmt::format("{:.{}g}", step.inverse_values[j], config::kDisplayPrecision)});
			}
			ui::render_table(rows, config::kTableMinWidth, transform_table);
			fmt::println("P(y) = {}", fmt::format("{:y}", step.inverse_polynomial));
			fmt::println("x_{} = P({:.{}g}) = {:.{}g}",
						 i + 1,
						 target_value,
						 config::kDisplayPrecision,
						 step.root,
						 config::kDisplayPrecision);
		}

		fmt::println("\n------------------------------------------------------------");
		fmt::println("\nRoots for y = {:.{}g}:", target_value, config::kDisplayPrecision);
		fmt::print("S = {{");
		for (SizeType i = 0; i < history.roots.size(); ++i) {
			if (i > 0)
				fmt::print(", ");
			fmt::print("{:.{}g}", history.roots[i], config::kDisplayPrecision);
		}
		fmt::println("}}");
	}

	auto inverse_interpolation(const RealList	  &nodes,
							   const RealList	  &values,
							   Real				   target_value,
							   SizeType			   max_points,
							   InterpolationMethod method) -> RealList {
		return build_inverse_interpolation_history(nodes, values, target_value, max_points, method)
			.roots;
	}

} // namespace detail

auto display_inverse_interpolation(const RealList	  &nodes,
								   const RealList	  &values,
								   Real				   target_value,
								   SizeType			   max_points,
								   InterpolationMethod method) -> void {
	const InverseInterpolationHistory history = detail::build_inverse_interpolation_history(
		nodes, values, target_value, max_points, method);
	detail::render_inverse_interpolation_history(nodes, values, target_value, history);
}
