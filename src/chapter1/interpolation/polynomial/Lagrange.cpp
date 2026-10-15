#include <chapter1/interpolation/horner.hpp>
#include <chapter1/interpolation/polynomial/Lagrange.hpp>
#include <config.hpp>
#include <math_object.hpp>
#include <ui/table.hpp>

#include <fmt/format.h>

#include <cmath>
#include <string>
#include <vector>

namespace detail {

	auto build_Lagrange_interpolation_history(const RealList &nodes, const RealList &values)
		-> LagrangeInterpolationHistory {

		if (nodes.empty() || nodes.size() != values.size())
			return {};

		LagrangeInterpolationHistory history{};
		history.normalizers.resize(nodes.size());

		SizeType degree	   = nodes.size() - 1;
		history.node_polys = detail::build_node_polynomial_history(nodes);
		history.values	   = values;

		for (SizeType i = 0; i <= degree; ++i) {
			history.normalizers[i] = 1;
			for (SizeType j = 0; j <= degree; ++j) {
				if (j != i) {
					history.normalizers[i] *= (nodes[i] - nodes[j]);
				}
			}

			const auto division =
				detail::divide_by_linear(history.node_polys.products.back(), {-nodes[i], 1});
			history.basis.push_back(division.quotient);
		}

		return history;
	}

	auto Lagrange_polynomial(const LagrangeInterpolationHistory &history) -> Polynomial {
		const std::vector<Polynomial> &basis	   = history.basis;
		const RealList				  &normalizers = history.normalizers;
		const RealList				  &values	   = history.values;

		Polynomial result(basis.back().size(), 0.0);
		for (SizeType i = 0; i < basis.size(); ++i) {
			result = result + values[i] / normalizers[i] * basis[i];
		}

		return result;
	}

	auto Lagrange_polynomial(const RealList &nodes, const RealList &values) -> Polynomial {
		if (nodes.empty() || nodes.size() != values.size())
			return {};

		SizeType degree = nodes.size() - 1;

		const Polynomial node_poly = detail::node_polynomial(nodes);
		Polynomial		 result(degree + 1, 0.0);

		for (SizeType i = 0; i <= degree; ++i) {
			Real normalizer = 1;
			for (SizeType j = 0; j <= degree; ++j) {
				if (j != i) {
					normalizer *= (nodes[i] - nodes[j]);
				}
			}

			const auto division = detail::divide_by_linear(node_poly, {-nodes[i], 1});
			const auto basis	= division.quotient;

			result = result + values[i] / normalizer * basis;
		}

		return result;
	}

	namespace {

		auto format_table_value(Real value) -> std::string {
			if (std::abs(value) < config::kNearZeroThreshold)
				return "0";
			return fmt::format("{:.{}g}", value, config::kDisplayPrecision);
		}

	} // namespace

	auto render_Lagrange_interpolation(const RealList					  &nodes,
									   const LagrangeInterpolationHistory &history) -> void {
		fmt::print("Node polynomials' coefficients\r\n");
		detail::render_nodes_polynomial_table(nodes, history.node_polys);

		fmt::print("\r\nBasis polynomials' coefficients\r\n");
		std::vector<std::vector<std::string>> basis_table{};
		basis_table.reserve(history.basis.size() + 1);

		std::vector<std::string> basis_header{};
		basis_header.reserve(history.basis.back().size() + 1);
		basis_header.emplace_back("C_i(x)");
		for (SizeType k = history.basis.back().size(); k > 0; --k) {
			basis_header.push_back(fmt::format("c_i{}", k - 1));
		}
		basis_table.push_back(basis_header);
		for (SizeType i = 0; i < history.basis.size(); ++i) {
			std::vector<std::string> row{};
			row.reserve(history.basis[i].size() + 1);
			row.push_back(fmt::format("i = {}", i));

			const auto &basis = history.basis[i];
			for (auto it = basis.rbegin(); it != basis.rend(); ++it) {
				row.push_back(format_table_value(*it));
			}
			basis_table.push_back(row);
		}
		ui::render_table(basis_table, config::kTableMinWidth);

		fmt::print("\r\nBarycentrics' coefficients\r\n");
		std::vector<std::vector<std::string>> bary_table{};
		bary_table.reserve(4);
		bary_table.push_back({"i"});
		bary_table.push_back({"y_i"});
		bary_table.push_back({"D_i"});
		bary_table.push_back({"y_i/D_i"});
		for (SizeType i = 0; i < history.values.size(); ++i) {
			const auto &value_i		 = history.values[i];
			const auto &normalizer_i = history.normalizers[i];

			bary_table[0].push_back(fmt::format("{}", i));
			bary_table[1].push_back(format_table_value(value_i));
			bary_table[2].push_back(format_table_value(normalizer_i));
			bary_table[3].push_back(format_table_value(value_i / normalizer_i));
		}
		ui::render_table(bary_table, config::kTableMinWidth);
	}
} // namespace detail

auto display_Lagrange_interpolation(const RealList &nodes, const RealList &values) -> void {
	const auto history = detail::build_Lagrange_interpolation_history(nodes, values);
	detail::render_Lagrange_interpolation(nodes, history);

	const auto poly = detail::Lagrange_polynomial(history);
	fmt::print("\r\nLagrange interpolated polynomials: L(x) = sum[i=0,n]:y_i / D_i * C_i(x)");
	fmt::print("\r\n                                        = {}\r\n", poly);
}