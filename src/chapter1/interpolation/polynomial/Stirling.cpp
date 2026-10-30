#include <chapter1/interpolation/horner.hpp>
#include <io/output_handling.hpp>
#include <chapter1/interpolation/polynomial/Stirling.hpp>
#include <config.hpp>
#include <math_object.hpp>
#include <math_utils.hpp>

#include <fmt/format.h>

#include <cmath>
#include <string>
#include <vector>

namespace detail {
	namespace {

		auto format_table_value(Real value) -> std::string {
			if (std::abs(value) < config::kNearZeroThreshold)
				return "0";
			return fmt::format("{:.{}g}", value, config::kDisplayPrecision);
		}

	} // namespace

	auto build_Stirling_history(const RealList &nodes, const RealList &values)
		-> EquidistantHistory {
		if (nodes.empty() || nodes.size() != values.size())
			return {};
		if (nodes.size() % 2 == 0)
			return {};
		if (nodes.size() == 1)
			return EquidistantHistory{{{values[0]}}, {{1.0}}};

		SizeType n = nodes.size() / 2;
		Real	 h = nodes[1] - nodes[0];
		for (SizeType i = 1; i < 2 * n; ++i) {
			if (nodes[i + 1] - nodes[i] != h)
				return {};
		}

		EquidistantHistory history{};

		auto &fwd_diffs = history.finite_diffs;
		fwd_diffs.resize(2 * n + 1);
		for (auto &diff_order : fwd_diffs) {
			diff_order.resize(2 * n + 1);
		}
		for (SizeType k = 0; k <= 2 * n; ++k) {
			fwd_diffs[k][0] = values[k];
		}
		for (SizeType i = 1; i <= 2 * n; ++i) {
			for (SizeType j = 0; j <= 2 * n - i; ++j) {
				fwd_diffs[j][i] = fwd_diffs[j + 1][i - 1] - fwd_diffs[j][i - 1];
			}
		}

		RealList t_nodes{};
		for (SizeType i = 1; i <= n - 1; ++i) {
			t_nodes.push_back(std::pow(static_cast<double>(i), 2));
		}
		history.node_polys = detail::build_node_polynomial_history(t_nodes);

		return history;
	}

	auto Stirling_polynomial(const RealList &nodes, const RealList &values) -> Polynomial {
		if (nodes.empty() || nodes.size() != values.size())
			return {};
		if (nodes.size() % 2 == 0)
			return {};
		if (nodes.size() == 1)
			return {values[0]};

		SizeType n = nodes.size() / 2;
		Real	 h = nodes[1] - nodes[0];
		for (SizeType i = 1; i < 2 * n; ++i) {
			if (nodes[i + 1] - nodes[i] != h)
				return {};
		}

		Polynomial result(2 * n + 1, 0.0);

		RealTable fwd_diffs;
		fwd_diffs.resize(2 * n + 1);
		for (auto &diff_order : fwd_diffs) {
			diff_order.resize(2 * n + 1);
		}
		for (SizeType k = 0; k <= 2 * n; ++k) {
			fwd_diffs[k][0] = values[k];
		}
		for (SizeType i = 1; i <= 2 * n; ++i) {
			for (SizeType j = 0; j <= 2 * n - i; ++j) {
				fwd_diffs[j][i] = fwd_diffs[j + 1][i - 1] - fwd_diffs[j][i - 1];
			}
		}

		RealList t_nodes{};
		for (SizeType i = 1; i <= n - 1; ++i) {
			t_nodes.push_back(std::pow(static_cast<double>(i), 2));
		}
		PolynomialList node_polys = detail::build_node_polynomial_history(t_nodes);

		result = result + Polynomial{values[n]};
		for (SizeType i = 0; i <= n - 1; ++i) {
			const auto &node_poly		   = node_polys[i];
			Polynomial	sqrt_argument_poly = utils::compose_power(node_poly, 2);
			Real odd_coeff	= (fwd_diffs[n - i][2 * i + 1] + fwd_diffs[n - i - 1][2 * i + 1]) / 2;
			Real even_coeff = fwd_diffs[n - i - 1][2 * i + 2];

			result = result + odd_coeff / utils::factorial(2 * i + 1) * Polynomial{0, 1} *
								  sqrt_argument_poly;
			result = result + even_coeff / utils::factorial(2 * i + 2) * Polynomial{0, 0, 1} *
								  sqrt_argument_poly;
		}

		return result;
	}

	auto render_Stirling_interpolation(const EquidistantHistory &history) -> void {
		SizeType n = history.finite_diffs.size() / 2;

		fmt::print("Divided differences table Δ[j][i] = Δ[j + 1][i-1] - Δ[j][i-1] (Δ[i][k] = "
				   "Δ^(i)y_k)\r\n");
		std::vector<std::vector<std::string>> diff_table{};
		std::vector<std::string>			  diff_header{};
		diff_header.insert(diff_header.end(), {"i", "y_i"});
		for (SizeType j = 1; j <= 2 * n; ++j) {
			diff_header.push_back(fmt::format("Δ[i][{}]", j));
		}
		diff_table.push_back(diff_header);
		for (SizeType i = 0; i <= 2 * n; ++i) {
			std::vector<std::string> row;
			row.insert(row.end(), {fmt::format("{}", i)});
			for (SizeType j = 0; j <= 2 * n - i; ++j) {
				row.push_back(format_table_value(history.finite_diffs[i][j]));
			}
			diff_table.push_back(row);
		}

		auto number_of_cols = diff_table.front().size();

		auto transform_table = [n, number_of_cols](ftxui::Table &table) {
			table.SelectColumn(0).DecorateCells(color(ftxui::Color::Yellow));
			table.SelectColumn(0).DecorateCells(ftxui::bold);

			ui::Point archor = {0, static_cast<int>(n + 1)};

			for (size_t i = 0; i <= n; ++i) {
				for (size_t j = 0; j < 3; ++j) {
					if (archor.x + 2 * i + j < number_of_cols) {
						table.SelectCell(archor.x + 2 * i + j, archor.y - i)
							.DecorateCells(color(ftxui::Color::Red));
						table.SelectCell(archor.x + 2 * i + j, archor.y - i)
							.DecorateCells(ftxui::bold);
					}
				}
			}
		};

		ui::render_table(diff_table, config::kTableMinWidth, transform_table);

		fmt::print("Forward coefficients (Δ[{}][0] = y{})\r\n", n, n);
		fmt::print("c_(2k-1) = (Δ[n-k+1][2k-1] + Δ[n-k][2k-1]) / 2 (k=1,n), n = {}\r\n", n);
		fmt::print("c_(2k)   = Δ[n-k][2k]            			   (k=0,n), n = {}\r\n", n);
		std::vector<std::vector<std::string>> coeff_table{};
		const auto							  fwd_diffs = history.finite_diffs;
		coeff_table.push_back({"i", "0"});
		coeff_table.push_back({"c_i", format_table_value(fwd_diffs[n][0])});
		coeff_table.push_back({"c_i / i!", format_table_value(fwd_diffs[n][0])});
		for (SizeType i = 0; i <= n - 1; ++i) {
			coeff_table[0].push_back(fmt::format("{}", 2 * i + 1));
			coeff_table[0].push_back(fmt::format("{}", 2 * i + 2));
			Real odd_coeff	= (fwd_diffs[n - i][2 * i + 1] + fwd_diffs[n - i - 1][2 * i + 1]) / 2;
			Real even_coeff = fwd_diffs[n - i - 1][2 * i + 2];
			coeff_table[1].push_back(format_table_value(odd_coeff));
			coeff_table[1].push_back(format_table_value(even_coeff));
			coeff_table[2].push_back(format_table_value(odd_coeff / utils::factorial(2 * i + 1)));
			coeff_table[2].push_back(format_table_value(even_coeff / utils::factorial(2 * i + 2)));
		}
		ui::render_table(coeff_table, config::kTableMinWidth);

		fmt::print("\r\nNode polynomials' coefficients w_(i+1)(t^2) (w0(t^2) = 1)\r\n");
		RealList t_nodes{};
		for (SizeType i = 1; i <= n - 1; ++i) {
			t_nodes.push_back(std::pow(static_cast<double>(i), 2));
		}
		detail::render_nodes_polynomial_table(t_nodes, history.node_polys, "t", 2);
	}
} // namespace detail

namespace {
	auto Stirling_polynomial(const EquidistantHistory &history) -> Polynomial {
		const RealTable		 &fwd_diffs	 = history.finite_diffs;
		const PolynomialList &node_polys = history.node_polys;

		SizeType   n = fwd_diffs.size() / 2;
		Polynomial result(2 * n + 1, 0.0);

		result = result + Polynomial{fwd_diffs[n][0]};
		for (SizeType i = 0; i <= n - 1; ++i) {
			const auto &node_poly		   = node_polys[i];
			Polynomial	sqrt_argument_poly = utils::compose_power(node_poly, 2);
			Real odd_coeff	= (fwd_diffs[n - i][2 * i + 1] + fwd_diffs[n - i - 1][2 * i + 1]) / 2;
			Real even_coeff = fwd_diffs[n - i - 1][2 * i + 2];

			result = result + odd_coeff / utils::factorial(2 * i + 1) * Polynomial{0, 1} *
								  sqrt_argument_poly;
			result = result + even_coeff / utils::factorial(2 * i + 2) * Polynomial{0, 0, 1} *
								  sqrt_argument_poly;
		}

		return result;
	}
} // namespace

auto display_Stirling_interpolation(const RealList &nodes, const RealList &values) -> void {
	ui::DisplayFrame display_frame;
	const auto history = detail::build_Stirling_history(nodes, values);
	if (history.finite_diffs.empty()) {
		fmt::print("Invalid Input\r\n");
		return;
	}
	const SizeType n = nodes.size() / 2;
	const Real	   h = nodes[1] - nodes[0];

	detail::render_Stirling_interpolation(history);

	const auto poly = Stirling_polynomial(history);
	fmt::print(
		"\r\nStirling interpolated polynomials: P(x) = P(x{} + th) = Q(t) = y{} + sum[i=0,{}] "
		"(c_(2i+1)/(2i+1)! * t*w_i(t^2) + (c_(2i+2)/(2i+2)! * t^2*w_i(t^2)",
		n,
		n,
		n - 1);
	fmt::print("\r\n                                           = {}\r\n",
			   fmt::format("{:t}", poly));
	fmt::print("Interpolation parameter: t = (x - x{}) / h\r\n", n);
	fmt::print("x{} = {}\r\n", n, nodes[n]);
	fmt::print("h  = {:.{}g}\r\n", h, config::kDisplayPrecision);
}