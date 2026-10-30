#include <chapter1/spline/linear.hpp>
#include <io/output_handling.hpp>

#include <config.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace detail {
	auto find_linear_spline(const RealList &nodes, const RealList &values) -> SplineFunction {
		if (nodes.size() < 2)
			throw std::invalid_argument("A linear spline requires at least two nodes");
		if (nodes.size() != values.size())
			throw std::invalid_argument("Linear spline nodes and values must have the same size");

		SplineFunction spline;
		spline.nodes = nodes;
		spline.piecewise_poly.reserve(nodes.size() - 1);

		for (SizeType i = 0; i + 1 < nodes.size(); ++i) {
			const Real width = nodes[i + 1] - nodes[i];
			if (!(width > 0))
				throw std::invalid_argument("Linear spline nodes must be strictly increasing");

			const Real slope = (values[i + 1] - values[i]) / width;
			spline.piecewise_poly.push_back({values[i], slope});
		}

		return spline;
	}

	auto render_linear_spline(const SplineFunction &spline) -> void {
		const SizeType interval_count =
			spline.nodes.size() < 2
				? 0
				: std::min(spline.nodes.size() - 1, spline.piecewise_poly.size());
		if (interval_count == 0)
			return;

		fmt::print("Piecewise polynomials' coefficients:\r\n");
		fmt::print("a_i = (y_(i+1) - y_i)/(x_(i+1) - x_i)\r\n");
		fmt::print("b_i = y_i\r\n");
		std::vector<std::vector<std::string>> rows = {{"i", "Interval", "a_i", "b_i"}};
		rows.reserve(interval_count + 1);

		for (SizeType i = 0; i < interval_count; ++i) {
			const Polynomial &poly		   = spline.piecewise_poly[i];
			const auto		  format_value = [](Real value) {
				return fmt::format("{:.{}g}", value, config::kDisplayPrecision);
			};
			rows.push_back({fmt::format("{}", i),
							fmt::format("[{}, {}]",
										format_value(spline.nodes[i]),
										format_value(spline.nodes[i + 1])),
							poly.empty() ? "0" : format_value(poly[1]),
							poly.size() < 2 ? "0" : format_value(poly[0])});
		}

		auto transform_table = [](ftxui::Table &table) {
			table.SelectColumn(0).DecorateCells(ftxui::color(ftxui::Color::Yellow));
			table.SelectColumn(0).DecorateCells(ftxui::bold);
		};

		ui::render_table(rows, config::kTableMinWidth, transform_table);

		fmt::print("Linear spline function: S(x) = S_i(x - x_i) [x_i <= x <= x_(i+1)]\r\n");
		fmt::print("{}", spline);
	}

} // namespace detail

auto display_linear_spline(const RealList &nodes, const RealList &values) -> void {
	ui::DisplayFrame display_frame;
	const SplineFunction spline = detail::find_linear_spline(nodes, values);
	detail::render_linear_spline(spline);
}
