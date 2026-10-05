#include <chapter1/interpolation.hpp>
#include <config.hpp>
#include <math_object.hpp>
#include <ui/table.hpp>

#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <vector>

#include <fmt/format.h>

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
		std::vector<std::string> header = {"i", "x_i"};

		std::vector<std::vector<std::string>> rows;
		rows.reserve(nodes.size() + 1);
		rows.push_back(header);

		for (std::size_t i = 0; i < nodes.size(); i++) {
			std::vector<std::string> row;
			row.reserve(2);

			row.push_back(fmt::format("{}", i));

			if (std::abs(nodes[i]) < config::kNearZeroThreshold)
				row.push_back("0");
			else
				row.push_back(fmt::format("{:+.{}}", nodes[i], config::kDisplayPrecision));

			rows.push_back(row);
		}

		ui::render_table(rows, 20);
	}

} // namespace detail

auto display_chebyshev_nodes(Real lower, Real upper, SizeType num_nodes) -> void {
	const auto nodes = detail::generate_chebyshev_nodes(lower, upper, num_nodes);
	detail::render_nodes_table(nodes);
}
