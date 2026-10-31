#include <chapter1/interpolation/horner.hpp>
#include <chapter1/interpolation/interval.hpp>
#include <chapter1/interpolation/inverse.hpp>
#include <chapter1/interpolation/polynomial/Lagrange.hpp>
#include <chapter1/interpolation/polynomial/Newton.hpp>
#include <config.hpp>
#include <io/output_handling.hpp>

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

	auto inverse_iteration(const RealList &nodes,
						   const RealList &values,
						   Real			   target_value,
						   Real			   tolerance,
						   SizeType		   max_points) -> RealList {
		return build_inverse_iteration_history(nodes, values, target_value, tolerance, max_points)
			.roots;
	}

	auto build_inverse_iteration_history(const RealList &nodes,
										 const RealList &values,
										 Real			 target_value,
										 Real			 tolerance,
										 SizeType		 max_points) -> InverseIterationHistory {
		InverseIterationHistory history;
		const bool invalid = nodes.size() != values.size() || nodes.size() < 2 || max_points < 2 ||
							 tolerance <= 0 || !std::isfinite(target_value) ||
							 !std::isfinite(tolerance);
		if (invalid)
			return history;

		const Real h = std::abs(nodes[1] - nodes[0]);
		if (h < config::kNearZeroThreshold)
			return history;

		for (SizeType i = 1; i < nodes.size() - 1; ++i) {
			if (std::abs(std::abs(nodes[i + 1] - nodes[i]) - h) > 1e-9)
				return history;
		}

		const SizeType n = nodes.size() - 1;
		for (SizeType i = 0; i <= n; ++i) {
			if (values[i] == target_value) {
				history.roots.push_back(nodes[i]);
			}
			if (i < n) {
				const Real d_curr = values[i] - target_value;
				const Real d_next = values[i + 1] - target_value;
				if (d_curr * d_next < 0) {
					history.isolating_intervals.push_back({i, i + 1});
				}
			}
		}

		for (SizeType i = 0; i < history.isolating_intervals.size(); ++i) {
			const SizeType p_i = history.isolating_intervals[i].start;

			const bool is_fwd = (n - p_i > max_points || n - p_i >= p_i);
			history.is_forward.push_back(is_fwd);

			RealList   sel_nodes;
			RealList   sel_values;
			Polynomial poly;

			if (is_fwd) {
				const SizeType k = std::min(n - p_i, max_points - 1);
				sel_nodes  = RealList(nodes.begin() + static_cast<std::ptrdiff_t>(p_i),
									  nodes.begin() + static_cast<std::ptrdiff_t>(p_i + k + 1));
				sel_values = RealList(values.begin() + static_cast<std::ptrdiff_t>(p_i),
									  values.begin() + static_cast<std::ptrdiff_t>(p_i + k + 1));
				poly	   = Newton_equidistant_forward_polynomial(sel_nodes, sel_values);
			} else {
				const SizeType k		 = std::min(p_i + 1, max_points - 1);
				const SizeType start_idx = p_i + 1 - k;
				sel_nodes  = RealList(nodes.begin() + static_cast<std::ptrdiff_t>(start_idx),
									  nodes.begin() + static_cast<std::ptrdiff_t>(p_i + 2));
				sel_values = RealList(values.begin() + static_cast<std::ptrdiff_t>(start_idx),
									  values.begin() + static_cast<std::ptrdiff_t>(p_i + 2));
				poly	   = Newton_equidistant_backward_polynomial(sel_nodes, sel_values);
			}

			if (poly.empty())
				continue;

			const Real dy_pi = values[p_i + 1] - values[p_i];
			if (std::abs(dy_pi) < config::kNearZeroThreshold)
				continue;

			const Real p0 = evaluate_polynomial(0, poly);
			const Real t0 = (target_value - p0) / dy_pi;

			Polynomial phi(std::max(poly.size(), SizeType{2}), 0.0);
			for (SizeType j = 0; j < poly.size(); ++j) {
				phi[j] = -poly[j] / dy_pi;
			}
			phi[0] += target_value / dy_pi;
			phi[1] += 1.0;

			for (Real &c : phi) {
				if (std::abs(c) < config::kNearZeroThreshold)
					c = 0.0;
			}
			while (phi.size() > 1 && phi.back() == 0.0) {
				phi.pop_back();
			}

			Real		   t_prev	= t0;
			Real		   t_curr	= t0;
			Real		   delta	= 0.0;
			const SizeType kMaxIter = 1000;
			RealTable	   iter_rows;

			for (SizeType k = 1; k <= kMaxIter; ++k) {
				t_curr = evaluate_polynomial(t_prev, phi);
				delta  = std::abs(t_curr - t_prev) * h;
				iter_rows.push_back({static_cast<Real>(k), t_prev, t_curr, delta});
				t_prev = t_curr;
				if (delta <= tolerance)
					break;
			}

			const Real root = p0 + t_curr * h;

			history.interpolation_nodes.push_back(std::move(sel_nodes));
			history.interpolation_values.push_back(std::move(sel_values));
			history.polynomials.push_back(std::move(poly));
			history.phi_functions.push_back(std::move(phi));
			history.iteration_tables.push_back(std::move(iter_rows));
			history.interval_roots.push_back(root);
			history.roots.push_back(root);
		}

		std::sort(history.roots.begin(), history.roots.end());
		auto it =
			std::unique(history.roots.begin(), history.roots.end(), [tolerance](Real a, Real b) {
				return std::abs(a - b) <= std::max(tolerance, config::kNearZeroThreshold);
			});
		history.roots.erase(it, history.roots.end());

		return history;
	}

	auto render_inverse_iteration_history(const RealList				&nodes,
										  const RealList				&values,
										  Real							 target_value,
										  const InverseIterationHistory &history) -> void {
		if (nodes.empty() || nodes.size() != values.size())
			return;

		const Real h = std::abs(nodes[1] - nodes[0]);

		auto transform_table = [](ftxui::Table &table) {
			table.SelectColumn(0).DecorateCells(color(ftxui::Color::Yellow));
			table.SelectColumn(0).DecorateCells(ftxui::bold);
		};

		fmt::println(
			"Isolating intervals for ȳ = {:.{}g}:", target_value, config::kDisplayPrecision);
		std::vector<std::vector<std::string>> itv_rows{{"p_i", "[x_i, x_(i+1)]", "[y_i, y_(i+1)]"}};
		for (SizeType i = 0; i < history.isolating_intervals.size(); ++i) {
			const SizeType p_i = history.isolating_intervals[i].start;
			itv_rows.push_back({
				fmt::format("p{} = {}", i, p_i),
				fmt::format("[{:.{}g}, {:.{}g}]",
							nodes[p_i],
							config::kDisplayPrecision,
							nodes[p_i + 1],
							config::kDisplayPrecision),
				fmt::format("[{:.{}g}, {:.{}g}]",
							values[p_i],
							config::kDisplayPrecision,
							values[p_i + 1],
							config::kDisplayPrecision),
			});
		}
		ui::render_table(itv_rows, config::kTableMinWidth, transform_table);

		for (SizeType i = 0; i < history.polynomials.size(); ++i) {
			const SizeType p_i = history.isolating_intervals[i].start;
			fmt::println("\n------------------------------------------------------------");
			fmt::println("\nIsolating interval {}: p{} = {}, [x{}, x{}] = [{:.{}g}, {:.{}g}]",
						 i + 1,
						 i,
						 p_i,
						 p_i,
						 p_i + 1,
						 nodes[p_i],
						 config::kDisplayPrecision,
						 nodes[p_i + 1],
						 config::kDisplayPrecision);

			fmt::println("\nSelected polynomial: {}",
						 history.is_forward[i] ? "Newton equidistant forward interpolation"
											   : "Newton equidistant backward interpolation");

			const auto	  &sel_nodes  = history.interpolation_nodes[i];
			const auto	  &sel_values = history.interpolation_values[i];
			const SizeType start_idx  = history.is_forward[i] ? p_i : (p_i + 2 - sel_nodes.size());
			std::vector<std::vector<std::string>> pts_tbl(3);
			pts_tbl[0].push_back("k");
			pts_tbl[1].push_back("x_k");
			pts_tbl[2].push_back("y_k");
			for (SizeType j = 0; j < sel_nodes.size(); ++j) {
				pts_tbl[0].push_back(fmt::format("{}", start_idx + j));
				pts_tbl[1].push_back(
					fmt::format("{:.{}g}", sel_nodes[j], config::kDisplayPrecision));
				pts_tbl[2].push_back(
					fmt::format("{:.{}g}", sel_values[j], config::kDisplayPrecision));
			}
			ui::render_table(pts_tbl, config::kTableMinWidth, transform_table);

			fmt::println("P_{}(t) = {}", i, fmt::format("{:t}", history.polynomials[i]));
			const Real p0 = evaluate_polynomial(0, history.polynomials[i]);
			fmt::println("P_{}(0) = {:.{}g}", i, p0, config::kDisplayPrecision);

			fmt::println("\nFixed point iterations:");
			fmt::println("ϕ(t) = 1/Δy_p{} * (ȳ - P_{}(t) + Δy_p{}*t) = {}",
						 i,
						 i,
						 i,
						 fmt::format("{:t}", history.phi_functions[i]));
			const Real dy_pi = values[p_i + 1] - values[p_i];
			const Real t0	 = (target_value - p0) / dy_pi;
			fmt::println("t0 = (ȳ - P_{}(0))/Δy_p{} = ({:.{}g} - {:.{}g}) / {:.{}g} = {:.{}g}",
						 i,
						 i,
						 target_value,
						 config::kDisplayPrecision,
						 p0,
						 config::kDisplayPrecision,
						 dy_pi,
						 config::kDisplayPrecision,
						 t0,
						 config::kDisplayPrecision);

			std::vector<std::vector<std::string>> iter_tbl{
				{"iteration", "t0", "t", "Δ = |t - t0|.h"}};
			for (SizeType j = 0; j < history.iteration_tables[i].size(); ++j) {
				const auto &row = history.iteration_tables[i][j];
				iter_tbl.push_back({
					fmt::format("{}", static_cast<SizeType>(row[0])),
					fmt::format("{:.{}g}", row[1], config::kDisplayPrecision),
					fmt::format("{:.{}g}", row[2], config::kDisplayPrecision),
					fmt::format("{:.{}g}", row[3], config::kDisplayPrecision),
				});
			}
			ui::render_table(iter_tbl, config::kTableMinWidth, transform_table);

			const Real t_final =
				history.iteration_tables[i].empty() ? 0.0 : history.iteration_tables[i].back()[2];

			fmt::println("\nRoot: x = P_{}(0) + th = {:.{}g} + ({:.{}g}) * ({:.{}g}) = {:.{}g}",
						 i,
						 p0,
						 config::kDisplayPrecision,
						 t_final,
						 config::kDisplayPrecision,
						 h,
						 config::kDisplayPrecision,
						 history.interval_roots[i],
						 config::kDisplayPrecision);
		}

		fmt::println("\n------------------------------------------------------------");
		fmt::println("\nRoots for ȳ = {:.{}g}:", target_value, config::kDisplayPrecision);
		fmt::print("S = {{");
		for (SizeType i = 0; i < history.roots.size(); ++i) {
			if (i > 0)
				fmt::print(", ");
			fmt::print("{:.{}g}", history.roots[i], config::kDisplayPrecision);
		}
		fmt::println("}}");
	}

} // namespace detail

auto display_inverse_interpolation(const RealList	  &nodes,
								   const RealList	  &values,
								   Real				   target_value,
								   SizeType			   max_points,
								   InterpolationMethod method) -> void {
	ui::DisplayFrame				  display_frame;
	const InverseInterpolationHistory history = detail::build_inverse_interpolation_history(
		nodes, values, target_value, max_points, method);
	detail::render_inverse_interpolation_history(nodes, values, target_value, history);
}

auto display_inverse_iteration(const RealList &nodes,
							   const RealList &values,
							   Real			   target_value,
							   Real			   tolerance,
							   SizeType		   max_points) -> void {
	ui::DisplayFrame display_frame;
	const auto		 history =
		detail::build_inverse_iteration_history(nodes, values, target_value, tolerance, max_points);
	if (history.isolating_intervals.empty() && history.roots.empty()) {
		fmt::println("Invalid Input");
		return;
	}
	detail::render_inverse_iteration_history(nodes, values, target_value, history);
}
