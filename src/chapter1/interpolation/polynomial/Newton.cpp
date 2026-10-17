#include <chapter1/interpolation/horner.hpp>
#include <chapter1/interpolation/polynomial/Newton.hpp>
#include <config.hpp>
#include <math_object.hpp>
#include <math_utils.hpp>
#include <ui/table.hpp>

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

	auto build_Newton_general_forward_history(const RealList &nodes, const RealList &values)
		-> NewtonInterpolationHistory {
		if (nodes.empty() || nodes.size() != values.size())
			return {};

		SizeType				   n = nodes.size() - 1;
		NewtonInterpolationHistory history{};

		auto &divided_diffs = history.divided_diffs;
		divided_diffs.resize(n + 1);
		for (auto &diff_order : divided_diffs) {
			diff_order.resize(n + 1);
		}
		for (SizeType k = 0; k <= n; ++k) {
			divided_diffs[k][0] = values[k];
		}
		for (SizeType i = 1; i <= n; ++i) {
			for (SizeType j = 0; j <= n - i; ++j) {
				divided_diffs[j][i] = divided_diffs[j + 1][i - 1] - divided_diffs[j][i - 1];
				divided_diffs[j][i] = divided_diffs[j][i] / (nodes[j + i] - nodes[j]);
			}
		}

		RealList sub_nodes(nodes.begin(), nodes.end() - 1);
		history.node_polys = detail::build_node_polynomial_history(sub_nodes);

		return history;
	}

	auto Newton_general_forward_polynomial(const RealList &nodes, const RealList &values)
		-> Polynomial {
		if (nodes.empty() || nodes.size() != values.size())
			return {};

		SizeType   n = nodes.size() - 1;
		Polynomial result(n + 1, 0.0);

		RealTable divided_diffs;
		divided_diffs.resize(n + 1);
		for (auto &diff_order : divided_diffs) {
			diff_order.resize(n + 1);
		}
		for (SizeType k = 0; k <= n; ++k) {
			divided_diffs[k][0] = values[k];
		}
		for (SizeType i = 1; i <= n; ++i) {
			for (SizeType j = 0; j <= n - i; ++j) {
				divided_diffs[j][i] = divided_diffs[j + 1][i - 1] - divided_diffs[j][i - 1];
				divided_diffs[j][i] = divided_diffs[j][i] / (nodes[j + i] - nodes[j]);
			}
		}

		RealList	   sub_nodes(nodes.begin(), nodes.end() - 1);
		PolynomialList node_polys = detail::build_node_polynomial_history(sub_nodes);

		for (SizeType i = 0; i < node_polys.size(); ++i) {
			result = result + divided_diffs[0][i] * node_polys[i];
		}

		return result;
	}

	auto render_Newton_general_forward(const RealList					&nodes,
									   const NewtonInterpolationHistory &history) -> void {
		fmt::print("Divided differences table f[j][i] = (f[j + 1][i-1] - f[j][i-1]) / (x_(j+i) - "
				   "x_j)\r\n");
		std::vector<std::vector<std::string>> diff_table{};
		std::vector<std::string>			  diff_header{};
		diff_header.insert(diff_header.end(), {"i", "x_i", "y_i"});
		for (SizeType j = 1; j < nodes.size(); ++j) {
			diff_header.push_back(fmt::format("Order j = {}", j));
		}
		diff_table.push_back(diff_header);
		for (SizeType i = 0; i < nodes.size(); ++i) {
			std::vector<std::string> row;
			row.insert(row.end(), {fmt::format("{}", i), fmt::format("{}", nodes[i])});
			for (SizeType j = 0; j < nodes.size() - i; ++j) {
				row.push_back(format_table_value(history.divided_diffs[i][j]));
			}
			diff_table.push_back(row);
		}
		ui::render_table(diff_table, config::kTableMinWidth);

		fmt::print("Newton forward coefficients f[0][i] = f[x0,x1,...,x_i] (f[x0] =y0)\r\n");
		std::vector<std::vector<std::string>> coeff_table{};
		coeff_table.push_back({"i"});
		coeff_table.push_back({"f[0][i]"});
		for (SizeType i = 0; i < nodes.size(); ++i) {
			coeff_table[0].push_back(fmt::format("{}", i));
			coeff_table[1].push_back(format_table_value(history.divided_diffs[0][i]));
		}
		ui::render_table(coeff_table, config::kTableMinWidth);

		fmt::print(
			"\r\nNode polynomials' coefficients w_(i+1)(x) = prod[k=0,i](x - x_k) (w0(x) = 1)\r\n");
		RealList sub_nodes(nodes.begin(), nodes.end() - 1);
		detail::render_nodes_polynomial_table(sub_nodes, history.node_polys);
	}

	auto build_Newton_general_backward_history(const RealList &nodes, const RealList &values)
		-> NewtonInterpolationHistory {
		if (nodes.empty() || nodes.size() != values.size())
			return {};

		SizeType				   n = nodes.size() - 1;
		NewtonInterpolationHistory history{};

		auto &divided_diffs = history.divided_diffs;
		divided_diffs.resize(n + 1);
		for (auto &diff_order : divided_diffs) {
			diff_order.resize(n + 1);
		}
		for (SizeType k = 0; k <= n; ++k) {
			divided_diffs[k][0] = values[k];
		}
		for (SizeType i = 1; i <= n; ++i) {
			for (SizeType j = i; j <= n; ++j) {
				divided_diffs[j][i] = divided_diffs[j][i - 1] - divided_diffs[j - 1][i - 1];
				divided_diffs[j][i] = divided_diffs[j][i] / (nodes[j] - nodes[j - i]);
			}
		}

		RealList sub_nodes(nodes.rbegin(), nodes.rend() - 1);
		history.node_polys = detail::build_node_polynomial_history(sub_nodes);

		return history;
	}

	auto Newton_general_backward_polynomial(const RealList &nodes, const RealList &values)
		-> Polynomial {
		if (nodes.empty() || nodes.size() != values.size())
			return {};

		SizeType   n = nodes.size() - 1;
		Polynomial result(n + 1, 0.0);

		RealTable divided_diffs;
		divided_diffs.resize(n + 1);
		for (auto &diff_order : divided_diffs) {
			diff_order.resize(n + 1);
		}
		for (SizeType k = 0; k <= n; ++k) {
			divided_diffs[k][0] = values[k];
		}
		for (SizeType i = 1; i <= n; ++i) {
			for (SizeType j = i; j <= n; ++j) {
				divided_diffs[j][i] = divided_diffs[j][i - 1] - divided_diffs[j - 1][i - 1];
				divided_diffs[j][i] = divided_diffs[j][i] / (nodes[j] - nodes[j - i]);
			}
		}

		RealList	   sub_nodes(nodes.rbegin(), nodes.rend() - 1);
		PolynomialList node_polys = detail::build_node_polynomial_history(sub_nodes);

		for (SizeType i = 0; i < node_polys.size(); ++i) {
			result = result + divided_diffs[n][i] * node_polys[i];
		}

		return result;
	}

	auto render_Newton_general_backward(const RealList					 &nodes,
										const NewtonInterpolationHistory &history) -> void {
		SizeType n = nodes.size() - 1;

		fmt::print("Divided differences table f[j][i] = (f[j][i-1] - f[j-1][i-1]) / (x_j - "
				   "x_(j-i))\r\n");
		std::vector<std::vector<std::string>> diff_table{};
		std::vector<std::string>			  diff_header{};
		diff_header.insert(diff_header.end(), {"i", "x_i", "y_i"});
		for (SizeType j = 1; j <= n; ++j) {
			diff_header.push_back(fmt::format("Order j = {}", j));
		}
		diff_table.push_back(diff_header);
		for (SizeType i = 0; i <= n; ++i) {
			std::vector<std::string> row;
			row.insert(row.end(), {fmt::format("{}", i), fmt::format("{}", nodes[i])});
			for (SizeType j = 0; j <= i; ++j) {
				row.push_back(format_table_value(history.divided_diffs[i][j]));
			}
			diff_table.push_back(row);
		}
		ui::render_table(diff_table, config::kTableMinWidth);

		fmt::print(
			"Newton forward coefficients f[{}][i] = f[x({}-i),x({}-i+1),...,x{}] (f[x{}] =y{})\r\n",
			n,
			n,
			n,
			n,
			n,
			n);
		std::vector<std::vector<std::string>> coeff_table{};
		coeff_table.push_back({"i"});
		coeff_table.push_back({fmt::format("f[{}][i]", n)});
		for (SizeType i = 0; i <= n; ++i) {
			coeff_table[0].push_back(fmt::format("{}", i));
			coeff_table[1].push_back(format_table_value(history.divided_diffs[n][i]));
		}
		ui::render_table(coeff_table, config::kTableMinWidth);

		fmt::print(
			"\r\nNode polynomials' coefficients w_(i+1)(x) = prod[k={}-i,{}](x - x_({}-k)) (w0(x) "
			"= 1)\r\n",
			n,
			n,
			n);
		RealList sub_nodes(nodes.rbegin(), nodes.rend() - 1);
		detail::render_nodes_polynomial_table(sub_nodes, history.node_polys);
	}

	auto build_Newton_equidistant_forward_history(const RealList &nodes, const RealList &values)
		-> EquidistantHistory {
		if (nodes.empty() || nodes.size() != values.size())
			return {};
		if (nodes.size() == 1)
			return EquidistantHistory{{{values[0]}}, {{1.0}}};

		SizeType n = nodes.size() - 1;
		Real	 h = nodes[1] - nodes[0];
		for (SizeType i = 1; i < n; ++i) {
			if (nodes[i + 1] - nodes[i] != h)
				return {};
		}

		EquidistantHistory history{};

		auto &fwd_diffs = history.finite_diffs;
		fwd_diffs.resize(n + 1);
		for (auto &diff_order : fwd_diffs) {
			diff_order.resize(n + 1);
		}
		for (SizeType k = 0; k <= n; ++k) {
			fwd_diffs[k][0] = values[k];
		}
		for (SizeType i = 1; i <= n; ++i) {
			for (SizeType j = 0; j <= n - i; ++j) {
				fwd_diffs[j][i] = fwd_diffs[j + 1][i - 1] - fwd_diffs[j][i - 1];
			}
		}

		RealList t_nodes{};
		for (SizeType i = 0; i < n; ++i) {
			t_nodes.push_back(static_cast<double>(i));
		}
		history.node_polys = detail::build_node_polynomial_history(t_nodes);

		return history;
	}

	auto Newton_equidistant_forward_polynomial(const RealList &nodes, const RealList &values)
		-> Polynomial {
		if (nodes.empty() || nodes.size() != values.size())
			return {};
		if (nodes.size() == 1)
			return {values[0]};

		SizeType n = nodes.size() - 1;
		Real	 h = nodes[1] - nodes[0];
		for (SizeType i = 1; i < n; ++i) {
			if (nodes[i + 1] - nodes[i] != h)
				return {};
		}

		Polynomial result(n + 1, 0.0);

		RealTable fwd_diffs;
		fwd_diffs.resize(n + 1);
		for (auto &diff_order : fwd_diffs) {
			diff_order.resize(n + 1);
		}
		for (SizeType k = 0; k <= n; ++k) {
			fwd_diffs[k][0] = values[k];
		}
		for (SizeType i = 1; i <= n; ++i) {
			for (SizeType j = 0; j <= n - i; ++j) {
				fwd_diffs[j][i] = fwd_diffs[j + 1][i - 1] - fwd_diffs[j][i - 1];
			}
		}

		RealList t_nodes{};
		for (SizeType i = 0; i < n; ++i) {
			t_nodes.push_back(static_cast<double>(i));
		}
		PolynomialList node_polys = detail::build_node_polynomial_history(t_nodes);

		for (SizeType i = 0; i < node_polys.size(); ++i) {
			result = result + fwd_diffs[0][i] / utils::factorial(i) * node_polys[i];
		}

		return result;
	}

	auto render_Newton_equidistant_forward(const RealList &nodes, const EquidistantHistory &history)
		-> void {
		fmt::print("Divided differences table Δ[j][i] = Δ[j + 1][i-1] - Δ[j][i-1] (Δ[i][k] = "
				   "Δ^(i)y_k)\r\n");
		std::vector<std::vector<std::string>> diff_table{};
		std::vector<std::string>			  diff_header{};
		diff_header.insert(diff_header.end(), {"i", "y_i"});
		for (SizeType j = 1; j < nodes.size(); ++j) {
			diff_header.push_back(fmt::format("Δ[i][{}]", j));
		}
		diff_table.push_back(diff_header);
		for (SizeType i = 0; i < nodes.size(); ++i) {
			std::vector<std::string> row;
			row.insert(row.end(), {fmt::format("{}", i)});
			for (SizeType j = 0; j < nodes.size() - i; ++j) {
				row.push_back(format_table_value(history.finite_diffs[i][j]));
			}
			diff_table.push_back(row);
		}
		ui::render_table(diff_table, config::kTableMinWidth);

		fmt::print("Forward coefficients (Δ[0][0] =y0)\r\n");
		std::vector<std::vector<std::string>> coeff_table{};
		coeff_table.push_back({"i"});
		coeff_table.push_back({"Δ[0][i]"});
		coeff_table.push_back({"Δ[0][i] / i!"});
		for (SizeType i = 0; i < nodes.size(); ++i) {
			coeff_table[0].push_back(fmt::format("{}", i));
			coeff_table[1].push_back(format_table_value(history.finite_diffs[0][i]));
			coeff_table[2].push_back(
				format_table_value(history.finite_diffs[0][i] / utils::factorial(i)));
		}
		ui::render_table(coeff_table, config::kTableMinWidth);

		fmt::print(
			"\r\nNode polynomials' coefficients w_(i+1)(t) = prod[k=0,i](t - k) (w0(t) = 1)\r\n");
		RealList sub_nodes(nodes.begin(), nodes.end() - 1);
		detail::render_nodes_polynomial_table(sub_nodes, history.node_polys, "t");
	}

	auto build_Newton_equidistant_backward_history(const RealList &nodes, const RealList &values)
		-> EquidistantHistory {
		if (nodes.empty() || nodes.size() != values.size())
			return {};
		if (nodes.size() == 1)
			return EquidistantHistory{{{values[0]}}, {{1.0}}};

		SizeType n = nodes.size() - 1;
		Real	 h = nodes[1] - nodes[0];
		for (SizeType i = 1; i < n; ++i) {
			if (nodes[i + 1] - nodes[i] != h)
				return {};
		}

		EquidistantHistory history{};

		auto &bwd_diffs = history.finite_diffs;
		bwd_diffs.resize(n + 1);
		for (auto &diff_order : bwd_diffs) {
			diff_order.resize(n + 1);
		}
		for (SizeType k = 0; k <= n; ++k) {
			bwd_diffs[k][0] = values[k];
		}
		for (SizeType i = 1; i <= n; ++i) {
			for (SizeType j = i; j <= n; ++j) {
				bwd_diffs[j][i] = bwd_diffs[j][i - 1] - bwd_diffs[j - 1][i - 1];
			}
		}

		RealList t_nodes{};
		for (SizeType i = 0; i < n; ++i) {
			t_nodes.push_back(-static_cast<double>(i));
		}
		history.node_polys = detail::build_node_polynomial_history(t_nodes);

		return history;
	}

	auto Newton_equidistant_backward_polynomial(const RealList &nodes, const RealList &values)
		-> Polynomial {
		if (nodes.empty() || nodes.size() != values.size())
			return {};
		if (nodes.size() == 1)
			return {values[0]};

		SizeType n = nodes.size() - 1;
		Real	 h = nodes[1] - nodes[0];
		for (SizeType i = 1; i < n; ++i) {
			if (nodes[i + 1] - nodes[i] != h)
				return {};
		}

		Polynomial result(n + 1, 0.0);

		RealTable bwd_diffs;
		bwd_diffs.resize(n + 1);
		for (auto &diff_order : bwd_diffs) {
			diff_order.resize(n + 1);
		}
		for (SizeType k = 0; k <= n; ++k) {
			bwd_diffs[k][0] = values[k];
		}
		for (SizeType i = 1; i <= n; ++i) {
			for (SizeType j = i; j <= n; ++j) {
				bwd_diffs[j][i] = bwd_diffs[j][i - 1] - bwd_diffs[j - 1][i - 1];
			}
		}

		RealList t_nodes{};
		for (SizeType i = 0; i < n; ++i) {
			t_nodes.push_back(-static_cast<double>(i));
		}
		PolynomialList node_polys = detail::build_node_polynomial_history(t_nodes);

		for (SizeType i = 0; i < node_polys.size(); ++i) {
			result = result + bwd_diffs[n][i] / utils::factorial(i) * node_polys[i];
		}

		return result;
	}

	auto render_Newton_equidistant_backward(const RealList			 &nodes,
											const EquidistantHistory &history) -> void {
		SizeType n = nodes.size() - 1;

		fmt::print(
			"Divided differences table ∇[j][i] = ∇[j][i-1] - ∇[j-1][i-1] (∇[i][k] = ∇^(i)y_k)\r\n");
		std::vector<std::vector<std::string>> diff_table{};
		std::vector<std::string>			  diff_header{};
		diff_header.insert(diff_header.end(), {"i", "y_i"});
		for (SizeType j = 1; j <= n; ++j) {
			diff_header.push_back(fmt::format("∇[{}][i]", j));
		}
		diff_table.push_back(diff_header);
		for (SizeType i = 0; i <= n; ++i) {
			std::vector<std::string> row;
			row.insert(row.end(), {fmt::format("{}", i)});
			for (SizeType j = 0; j <= i; ++j) {
				row.push_back(format_table_value(history.finite_diffs[i][j]));
			}
			diff_table.push_back(row);
		}
		ui::render_table(diff_table, config::kTableMinWidth);

		fmt::print("Newton forward coefficients (∇[x{}] =y{})\r\n", n, n);
		std::vector<std::vector<std::string>> coeff_table{};
		coeff_table.push_back({"i"});
		coeff_table.push_back({fmt::format("∇[{}][i]", n)});
		coeff_table.push_back({fmt::format("∇[{}][i] / i!", n)});
		for (SizeType i = 0; i <= n; ++i) {
			coeff_table[0].push_back(fmt::format("{}", i));
			coeff_table[1].push_back(format_table_value(history.finite_diffs[n][i]));
			coeff_table[2].push_back(
				format_table_value(history.finite_diffs[n][i] / utils::factorial(i)));
		}
		ui::render_table(coeff_table, config::kTableMinWidth);

		fmt::print("\r\nNode polynomials' coefficients w_(i+1)(t) = prod[k=0,i](t + k) (w0(t) "
				   "= 1)\r\n");
		RealList t_nodes{};
		for (SizeType i = 0; i < n; ++i) {
			t_nodes.push_back(-static_cast<double>(i));
		}
		detail::render_nodes_polynomial_table(t_nodes, history.node_polys, "t");
	}
} // namespace detail

namespace {
	auto Newton_general_forward_polynomial(const NewtonInterpolationHistory &history)
		-> Polynomial {
		const RealTable		 &divided_diffs = history.divided_diffs;
		const PolynomialList &node_polys	= history.node_polys;

		Polynomial result(node_polys.back().size(), 0.0);

		for (SizeType i = 0; i < node_polys.size(); ++i) {
			result = result + divided_diffs[0][i] * node_polys[i];
		}

		return result;
	}

	auto Newton_general_backward_polynomial(const NewtonInterpolationHistory &history)
		-> Polynomial {
		const RealTable		 &divided_diffs = history.divided_diffs;
		const PolynomialList &node_polys	= history.node_polys;
		const SizeType		  n				= divided_diffs.size() - 1;
		Polynomial			  result(n + 1, 0.0);

		for (SizeType i = 0; i < node_polys.size(); ++i) {
			result = result + divided_diffs[n][i] * node_polys[i];
		}

		return result;
	}

	auto Newton_equidistant_forward_polynomial(const EquidistantHistory &history) -> Polynomial {
		const RealTable		 &fwd_diffs	 = history.finite_diffs;
		const PolynomialList &node_polys = history.node_polys;

		Polynomial result(node_polys.back().size(), 0.0);

		for (SizeType i = 0; i < node_polys.size(); ++i) {
			result = result + fwd_diffs[0][i] / utils::factorial(i) * node_polys[i];
		}

		return result;
	}

	auto Newton_equidistant_backward_polynomial(const EquidistantHistory &history) -> Polynomial {
		const RealTable		 &bwd_diffs	 = history.finite_diffs;
		const PolynomialList &node_polys = history.node_polys;
		const SizeType		  n			 = bwd_diffs.size() - 1;
		Polynomial			  result(n + 1, 0.0);

		for (SizeType i = 0; i < node_polys.size(); ++i) {
			result = result + bwd_diffs[n][i] / utils::factorial(i) * node_polys[i];
		}

		return result;
	}
} // namespace

auto display_Newton_general_forward(const RealList &nodes, const RealList &values) -> void {
	const auto history = detail::build_Newton_general_forward_history(nodes, values);
	if (history.divided_diffs.empty()) {
		fmt::print("Invalid Input\r\n");
		return;
	}

	detail::render_Newton_general_forward(nodes, history);

	const auto poly = Newton_general_forward_polynomial(history);
	fmt::print("\r\nNewton interpolated polynomials: P(x) = sum[i=0,{}] f[0][i] * w_(i)(x)",
			   nodes.size() - 1);
	fmt::print("\r\n                                      = {}\r\n", poly);
}

auto display_Newton_general_backward(const RealList &nodes, const RealList &values) -> void {
	const auto history = detail::build_Newton_general_backward_history(nodes, values);
	if (history.divided_diffs.empty()) {
		fmt::print("Invalid Input\r\n");
		return;
	}

	const SizeType n = nodes.size() - 1;

	detail::render_Newton_general_backward(nodes, history);

	const auto poly = Newton_general_backward_polynomial(history);
	fmt::print("\r\nNewton interpolated polynomials: P(x) = sum[i=0,{}] f[{}][i] * w_(i)(x)", n, n);
	fmt::print("\r\n                                      = {}\r\n", poly);
}

auto display_Newton_equidistant_forward(const RealList &nodes, const RealList &values) -> void {
	const auto history = detail::build_Newton_equidistant_forward_history(nodes, values);
	if (history.finite_diffs.empty()) {
		fmt::print("Invalid Input\r\n");
		return;
	}

	const Real h = nodes[1] - nodes[0];

	detail::render_Newton_equidistant_forward(nodes, history);

	const auto poly = Newton_equidistant_forward_polynomial(history);
	fmt::print("\r\nNewton interpolated polynomials: P(x) = P(x0 + th) = Q(t) = sum[i=0,{}] "
			   "Δ[0][i] / i! * w_i(t)",
			   nodes.size() - 1);
	fmt::print("\r\n                                      = {}\r\n", fmt::format("{:t}", poly));
	fmt::print("Interpolation parameter: t = (x - x0) / h\r\n");
	fmt::print("x0 = {}\r\n", nodes[0]);
	fmt::print("h  = {:.{}g}\r\n", h, config::kDisplayPrecision);
}

auto display_Newton_equidistant_backward(const RealList &nodes, const RealList &values) -> void {
	const auto history = detail::build_Newton_equidistant_backward_history(nodes, values);
	if (history.finite_diffs.empty()) {
		fmt::print("Invalid Input\r\n");
		return;
	}

	const SizeType n = nodes.size() - 1;
	const Real	   h = nodes[1] - nodes[0];

	detail::render_Newton_equidistant_backward(nodes, history);

	const auto poly = Newton_equidistant_backward_polynomial(history);
	fmt::print("\r\nNewton interpolated polynomials: P(x) = P(x{} + th) = Q(t) = sum[i=0,{}] "
			   "∇[{}][i] / i! * w_i(t)",
			   n,
			   n,
			   n);
	fmt::print("\r\n                                      = {}\r\n", fmt::format("{:t}", poly));
	fmt::print("Interpolation parameter: t = (x - x{}) / h\r\n", n);
	fmt::print("x{} = {:.{}g}\r\n", n, nodes[n], config::kDisplayPrecision);
	fmt::print("h  = {:.{}g}\r\n", h, config::kDisplayPrecision);
}