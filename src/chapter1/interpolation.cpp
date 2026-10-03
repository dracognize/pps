#include <chapter1/interpolation.hpp>
#include <config.hpp>
#include <math_object.hpp>

#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <ftxui/dom/elements.hpp>
#include <ftxui/dom/node.hpp>
#include <ftxui/dom/table.hpp>
#include <ftxui/screen/color.hpp>
#include <ftxui/screen/screen.hpp>

namespace detail {

	auto generate_chebyshev_nodes(Real lower, Real upper, SizeType num_nodes) -> RealList {
		RealList nodes{};
		nodes.reserve(num_nodes);

		for (SizeType i = num_nodes; i > 0; --i) {
			Real node = 0.5 * ((upper - lower) *
								   std::cos((2.0 * i - 1) * std::numbers::pi / (2.0 * num_nodes)) +
							   lower + upper);

			nodes.push_back(node);
		}

		return nodes;
	}

	auto render_nodes_table(const RealList &nodes) -> void {
		using namespace ftxui;

		std::vector<std::string> header = {"i", "x_i"};

		std::vector<std::vector<std::string>> rows;
		rows.push_back(header);

		for (std::size_t i = 0; i < nodes.size(); i++) {
			std::vector<std::string> row;

			row.push_back(fmt::format("{}", i));

			if (std::abs(nodes[i]) < config::kNearZeroThreshold)
				row.push_back("0");
			else
				row.push_back(fmt::format("{:+.{}}", nodes[i], config::kDisplayPrecision));

			rows.push_back(row);
		}

		auto table = Table(rows);

		table.SelectRow(0).Decorate(color(Color::Cyan));

		table.SelectAll().Border(LIGHT);
		table.SelectAll().Separator(LIGHT);

		table.SelectRow(0).DecorateCells(bold);
		table.SelectRow(0).DecorateCells(hcenter);

		table.SelectColumn(0).DecorateCells(hcenter | size(WIDTH, GREATER_THAN, 20));
		table.SelectColumn(1).DecorateCells(hcenter | size(WIDTH, GREATER_THAN, 20));

		auto document = table.Render();
		auto screen	  = Screen::Create(Dimension::Fit(document, true));

		Render(screen, document);
		screen.Print();
		screen.ResetPosition();
		fmt::println("");
	}

} // namespace detail

auto display_chebyshev_nodes(Real lower, Real upper, SizeType num_nodes) -> void {
	auto nodes = detail::generate_chebyshev_nodes(lower, upper, num_nodes);
	detail::render_nodes_table(nodes);
}
