#include <chapter1/interpolation/horner.hpp>
#include <io/output_handling.hpp>
#include <chapter1/interpolation/polynomial/Gauss.hpp>
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

	auto build_GaussI_history(const RealList &nodes, const RealList &values) -> EquidistantHistory {
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
		for (SizeType i = 1; i <= n; ++i) {
			t_nodes.push_back(-static_cast<double>(i - 1));
			t_nodes.push_back(static_cast<double>(i));
		}
		history.node_polys = detail::build_node_polynomial_history(t_nodes);

		return history;
	}

	auto GaussI_polynomial(const RealList &nodes, const RealList &values) -> Polynomial {
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
		for (SizeType i = 1; i <= n; ++i) {
			t_nodes.push_back(-static_cast<double>(i - 1));
			t_nodes.push_back(static_cast<double>(i));
		}
		PolynomialList node_polys = detail::build_node_polynomial_history(t_nodes);

		for (SizeType i = 0; i <= 2 * n; ++i) {
			result = result + fwd_diffs[n - (i / 2)][i] / utils::factorial(i) * node_polys[i];
		}

		return result;
	}

	auto render_GaussI_interpolation(const EquidistantHistory &history) -> void {
		SizeType n = history.node_polys.size() / 2;

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
			ui::Point point = {1, static_cast<int>(n + 1)};

			while (point.x < number_of_cols) {
				table.SelectCell(point.x, point.y).DecorateCells(color(ftxui::Color::Red));
				table.SelectCell(point.x, point.y).DecorateCells(ftxui::bold);

				if (point.x % 2 == 0) {
					point.y--;
				}
				point.x++;
			}
		};

		ui::render_table(diff_table, config::kTableMinWidth, transform_table);

		fmt::print("Forward coefficients (Δ[{}][0] = y{})\r\n", n, n);
		fmt::print("c_(2k-1) = Δ[n-(k-1)][2k-1] (k=1,n)\r\n");
		fmt::print("c_(2k)   = Δ[n-k][2k]       (k=0,n)\r\n");
		std::vector<std::vector<std::string>> coeff_table{};
		const auto							  fwd_diffs = history.finite_diffs;
		coeff_table.push_back({"i"});
		coeff_table.push_back({"c_i"});
		coeff_table.push_back({"c_i / i!"});
		for (SizeType i = 0; i <= 2 * n; ++i) {
			coeff_table[0].push_back(fmt::format("{}", i));
			coeff_table[1].push_back(format_table_value(fwd_diffs[n - i / 2][i]));
			coeff_table[2].push_back(
				format_table_value(fwd_diffs[n - i / 2][i] / utils::factorial(i)));
		}
		ui::render_table(coeff_table, config::kTableMinWidth);

		fmt::print("\r\nNode polynomials' coefficients w_(i+1)(t) (w0(t) = 1)\r\n");
		RealList t_nodes{};
		for (SizeType i = 1; i <= n; ++i) {
			t_nodes.push_back(-static_cast<double>(i - 1));
			t_nodes.push_back(static_cast<double>(i));
		}
		detail::render_nodes_polynomial_table(t_nodes, history.node_polys, "t");
	}

	auto build_GaussII_history(const RealList &nodes, const RealList &values)
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
		for (SizeType i = 1; i <= n; ++i) {
			t_nodes.push_back(static_cast<double>(i - 1));
			t_nodes.push_back(-static_cast<double>(i));
		}
		history.node_polys = detail::build_node_polynomial_history(t_nodes);

		return history;
	}

	auto GaussII_polynomial(const RealList &nodes, const RealList &values) -> Polynomial {
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
		for (SizeType i = 1; i <= n; ++i) {
			t_nodes.push_back(static_cast<double>(i - 1));
			t_nodes.push_back(-static_cast<double>(i));
		}
		PolynomialList node_polys = detail::build_node_polynomial_history(t_nodes);

		for (SizeType i = 0; i <= 2 * n; ++i) {
			result = result + fwd_diffs[n - (i + 1) / 2][i] / utils::factorial(i) * node_polys[i];
		}

		return result;
	}

	auto render_GaussII_interpolation(const EquidistantHistory &history) -> void {
		SizeType n = history.node_polys.size() / 2;

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
			ui::Point point = {1, static_cast<int>(n + 1)};

			while (point.x < number_of_cols) {
				table.SelectCell(point.x, point.y).DecorateCells(color(ftxui::Color::Red));
				table.SelectCell(point.x, point.y).DecorateCells(ftxui::bold);

				if (point.x % 2 != 0) {
					point.y--;
				}
				point.x++;
			}
		};

		ui::render_table(diff_table, config::kTableMinWidth, transform_table);

		fmt::print("Forward coefficients (Δ[{}][0] = y{})\r\n", n, n);
		fmt::print("c_(2k-1) = Δ[n-k][2k-1] (k=1,n)\r\n");
		fmt::print("c_(2k)   = Δ[n-k][2k]   (k=0,n)\r\n");
		std::vector<std::vector<std::string>> coeff_table{};
		const auto							  fwd_diffs = history.finite_diffs;
		coeff_table.push_back({"i"});
		coeff_table.push_back({"c_i"});
		coeff_table.push_back({"c_i / i!"});
		for (SizeType i = 0; i <= 2 * n; ++i) {
			coeff_table[0].push_back(fmt::format("{}", i));
			coeff_table[1].push_back(format_table_value(fwd_diffs[n - (i + 1) / 2][i]));
			coeff_table[2].push_back(
				format_table_value(fwd_diffs[n - (i + 1) / 2][i] / utils::factorial(i)));
		}
		ui::render_table(coeff_table, config::kTableMinWidth);

		fmt::print("\r\nNode polynomials' coefficients w_(i+1)(t) (w0(t) = 1)\r\n");
		RealList t_nodes{};
		for (SizeType i = 1; i <= n; ++i) {
			t_nodes.push_back(static_cast<double>(i - 1));
			t_nodes.push_back(-static_cast<double>(i));
		}
		detail::render_nodes_polynomial_table(t_nodes, history.node_polys, "t");
	}
} // namespace detail

namespace {
	auto GaussI_polynomial(const EquidistantHistory &history) -> Polynomial {
		const RealTable		 &fwd_diffs	 = history.finite_diffs;
		const PolynomialList &node_polys = history.node_polys;

		SizeType   n = node_polys.back().size() / 2;
		Polynomial result(2 * n + 1, 0.0);

		for (SizeType i = 0; i <= 2 * n; ++i) {
			result = result + fwd_diffs[n - i / 2][i] / utils::factorial(i) * node_polys[i];
		}

		return result;
	}

	auto GaussII_polynomial(const EquidistantHistory &history) -> Polynomial {
		const RealTable		 &fwd_diffs	 = history.finite_diffs;
		const PolynomialList &node_polys = history.node_polys;

		SizeType   n = node_polys.back().size() / 2;
		Polynomial result(2 * n + 1, 0.0);

		for (SizeType i = 0; i <= 2 * n; ++i) {
			result = result + fwd_diffs[n - (i + 1) / 2][i] / utils::factorial(i) * node_polys[i];
		}

		return result;
	}
} // namespace

auto display_GaussI_interpolation(const RealList &nodes, const RealList &values) -> void {
	ui::DisplayFrame display_frame;
	const auto history = detail::build_GaussI_history(nodes, values);
	if (history.finite_diffs.empty()) {
		fmt::print("Invalid Input\r\n");
		return;
	}
	const SizeType n = nodes.size() / 2;
	const Real	   h = nodes[1] - nodes[0];

	detail::render_GaussI_interpolation(history);

	const auto poly = GaussI_polynomial(history);
	fmt::print("\r\nGauss I interpolated polynomials: P(x) = P(x{} + th) = Q(t) = sum[i=0,{}] "
			   "c_i / i! * w_i(t)",
			   n,
			   2 * n);
	fmt::print("\r\n                                       = {}\r\n", fmt::format("{:t}", poly));
	fmt::print("Interpolation parameter: t = (x - x{}) / h\r\n", n);
	fmt::print("x{} = {}\r\n", n, nodes[n]);
	fmt::print("h  = {:.{}g}\r\n", h, config::kDisplayPrecision);
}

auto display_GaussII_interpolation(const RealList &nodes, const RealList &values) -> void {
	ui::DisplayFrame display_frame;
	const auto history = detail::build_GaussII_history(nodes, values);
	if (history.finite_diffs.empty()) {
		fmt::print("Invalid Input\r\n");
		return;
	}
	const SizeType n = nodes.size() / 2;
	const Real	   h = nodes[1] - nodes[0];

	detail::render_GaussII_interpolation(history);

	const auto poly = GaussII_polynomial(history);
	fmt::print("\r\nGauss II interpolated polynomials: P(x) = P(x{} + th) = Q(t) = sum[i=0,{}] "
			   "c_i / i! * w_i(t)",
			   n,
			   2 * n);
	fmt::print("\r\n                                        = {}\r\n", fmt::format("{:t}", poly));
	fmt::print("Interpolation parameter: t = (x - x{}) / h\r\n", n);
	fmt::print("x{} = {}\r\n", n, nodes[n]);
	fmt::print("h  = {:.{}g}\r\n", h, config::kDisplayPrecision);
}