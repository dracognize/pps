#include <chapter1/spline/quadratic.hpp>
#include <io/output_handling.hpp>

#include <config.hpp>

#include <Eigen/LU>
#include <fmt/format.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace detail {
	auto build_quadratic_spline_history(const RealList &nodes,
										const RealList &values,
										Real			alpha,
										Real			beta,
										Real			gamma) -> SplineHistory {
		if (nodes.size() < 2 || nodes.size() != values.size())
			throw std::invalid_argument(
				"Quadratic spline requires matching nodes and values, with at least two nodes");
		const SizeType n = nodes.size() - 1;
		if (alpha == 0 || beta != (n % 2 == 0 ? alpha : -alpha))
			throw std::invalid_argument(
				"Quadratic spline endpoint coefficients must satisfy beta = (-1)^n * alpha");

		SplineHistory history;
		history.nodes  = nodes;
		history.values = values;
		history.widths.resize(n);

		history.system	  = Eigen::MatrixXd::Zero(nodes.size(), nodes.size());
		history.constants = Eigen::VectorXd::Zero(nodes.size());
		for (SizeType i = 0; i < n; ++i) {
			history.widths[i]	 = nodes[i + 1] - nodes[i];
			history.system(i, i) = history.system(i, i + 1) = 1;
			history.constants[i] = 2 * (values[i + 1] - values[i]) / history.widths[i];
		}
		history.system(n, 0) = alpha;
		history.system(n, n) = beta;
		history.constants[n] = gamma;
		history.derivatives	 = history.system.partialPivLu().solve(history.constants);

		history.piecewise_poly.reserve(n);
		for (SizeType i = 0; i < n; ++i) {
			const Real left	 = history.derivatives[static_cast<Eigen::Index>(i)];
			const Real right = history.derivatives[static_cast<Eigen::Index>(i + 1)];
			history.piecewise_poly.push_back(
				{values[i], left, (right - left) / (2 * history.widths[i])});
		}
		return history;
	}

	namespace {
		auto calculate_quadratic_spline(const SplineHistory &history) -> SplineFunction {
			SplineFunction spline;
			spline.nodes		  = history.nodes;
			spline.piecewise_poly = history.piecewise_poly;
			return spline;
		}
	} // namespace

	auto find_quadratic_spline(const RealList &nodes,
							   const RealList &values,
							   Real			   alpha,
							   Real			   beta,
							   Real			   gamma) -> SplineFunction {
		return calculate_quadratic_spline(
			build_quadratic_spline_history(nodes, values, alpha, beta, gamma));
	}

	auto render_quadratic_spline(const SplineHistory &history) -> void {
		auto transform_table = [](ftxui::Table &table) {
			table.SelectColumn(0).DecorateCells(ftxui::color(ftxui::Color::Yellow));
			table.SelectColumn(0).DecorateCells(ftxui::bold);
		};

		fmt::print("\nNodes differences: h_i = x_(i+1) - x_i\n");
		std::vector<std::vector<std::string>> width_rows = {{"i"}, {"h_i"}};
		for (Eigen::Index i = 0; i < history.widths.size(); ++i) {
			width_rows[0].push_back(fmt::format("{}", i));
			width_rows[1].push_back(
				fmt::format("{:.{}g}", history.widths[i], config::kDisplayPrecision));
		}
		ui::render_table(width_rows, config::kTableMinWidth, transform_table);

		std::vector<SizeType> column_widths(static_cast<SizeType>(history.system.cols()), 0);
		for (Eigen::Index col = 0; col < history.system.cols(); ++col)
			for (Eigen::Index row = 0; row < history.system.rows(); ++row)
				column_widths[static_cast<SizeType>(col)] =
					std::max(column_widths[static_cast<SizeType>(col)],
							 fmt::format("{:.5g}", history.system(row, col)).size());

		SizeType constant_width = 0;
		for (Eigen::Index row = 0; row < history.constants.size(); ++row)
			constant_width =
				std::max(constant_width, fmt::format("{:.5g}", history.constants[row]).size());

		std::vector<std::string> matrix_rows;
		for (Eigen::Index i = 0; i < history.system.rows(); ++i) {
			std::string row = "[ ";
			for (Eigen::Index j = 0; j < history.system.cols(); ++j) {
				if (j != 0)
					row += " ";
				row += fmt::format("{:>{}}",
								   fmt::format("{:.5g}", history.system(i, j)),
								   column_widths[static_cast<SizeType>(j)]);
			}
			row += fmt::format(
				" | {:>{}} ]", fmt::format("{:.5g}", history.constants[i]), constant_width);
			matrix_rows.push_back(std::move(row));
		}

		SizeType row_width = 0;
		for (const auto &row : matrix_rows)
			row_width = std::max(row_width, row.size());

		const std::string border = "\t+" + std::string(row_width + 8, '-') + "+\n";
		fmt::print("Quadratic spline system [A | B]:\n");
		fmt::print("{}", border);
		fmt::print("\t|    {}    |\n", std::string(row_width, ' '));
		for (const auto &row : matrix_rows)
			fmt::print("\t|    {:>{}}    |\n", row, row_width);
		fmt::print("\t|    {}    |\n", std::string(row_width, ' '));
		fmt::print("{}", border);
		fmt::print("\nNodal derivatives (solution of AX = B):\n");

		constexpr Eigen::Index kMaxDerivativeColumns = 10;
		for (Eigen::Index start = 0; start < history.derivatives.size();
			 start += kMaxDerivativeColumns) {
			std::vector<std::vector<std::string>> derivative_rows = {{"i"}, {"m_i = S'_i(x_i)"}};
			const Eigen::Index					  end =
				std::min(start + kMaxDerivativeColumns, history.derivatives.size());
			for (Eigen::Index i = start; i < end; ++i) {
				derivative_rows[0].push_back(fmt::format("{}", i));
				derivative_rows[1].push_back(
					fmt::format("{:.{}g}", history.derivatives[i], config::kDisplayPrecision));
			}
			ui::render_table(derivative_rows, config::kTableMinWidth, transform_table);
		}

		fmt::print("\nPiecewise polynomials' coefficients:\n");
		std::vector<std::vector<std::string>> coefficient_rows = {
			{"i", "Interval", "a_i", "b_i", "c_i"}};
		fmt::print("a_i = (m_(i+1) - m_i)/(2h_i)\r\n");
		fmt::print("b_i = m_i\r\n");
		fmt::print("c_i = y_i\n");
		for (SizeType i = 0; i < history.piecewise_poly.size(); ++i) {
			const auto number = [](Real value) {
				return fmt::format("{:.{}g}", value, config::kDisplayPrecision);
			};
			const auto &poly = history.piecewise_poly[i];
			coefficient_rows.push_back(
				{fmt::format("{}", i),
				 fmt::format("[{}, {}]", number(history.nodes[i]), number(history.nodes[i + 1])),
				 number(poly[2]),
				 number(poly[1]),
				 number(poly[0])});
		}

		ui::render_table(coefficient_rows, config::kTableMinWidth, transform_table);
	}
} // namespace detail

auto display_quadratic_spline(
	const RealList &nodes, const RealList &values, Real alpha, Real beta, Real gamma) -> void {
	ui::DisplayFrame display_frame;

	const auto history = detail::build_quadratic_spline_history(nodes, values, alpha, beta, gamma);
	detail::render_quadratic_spline(history);
	auto spline = detail::calculate_quadratic_spline(history);

	fmt::print("\nQuadratic spline function: S(x) = S_i(x - x_i) [x_i <= x <= x_(i+1)]\r\n");
	fmt::print("{}", spline);
}
