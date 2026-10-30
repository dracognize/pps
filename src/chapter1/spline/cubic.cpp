#include <chapter1/spline/cubic.hpp>
#include <io/output_handling.hpp>

#include <config.hpp>

#include <Eigen/LU>
#include <fmt/format.h>

#include <stdexcept>
#include <string>
#include <vector>

namespace detail {
	auto build_cubic_spline_history(const RealList					  &nodes,
									const RealList					  &values,
									const Eigen::Matrix<double, 2, 5> &boundary_conditions)
		-> SplineHistory {
		if (nodes.size() < 2 || nodes.size() != values.size())
			throw std::invalid_argument(
				"Cubic spline requires matching nodes and values, with at least two nodes");
		const SizeType n = nodes.size() - 1;
		SplineHistory  history;
		history.nodes  = nodes;
		history.values = values;
		history.widths.resize(n);
		history.system	  = Eigen::MatrixXd::Zero(nodes.size(), nodes.size());
		history.constants = Eigen::VectorXd::Zero(nodes.size());
		for (SizeType i = 0; i < n; ++i)
			history.widths[i] = nodes[i + 1] - nodes[i];

		for (SizeType i = 0; i + 1 < n; ++i) {
			const Real h0			 = history.widths[i];
			const Real h1			 = history.widths[i + 1];
			history.system(i, i)	 = h0 / 6;
			history.system(i, i + 1) = (h0 + h1) / 3;
			history.system(i, i + 2) = h1 / 6;
			history.constants[i] =
				(values[i + 1] - values[i]) / h0 - (values[i + 2] - values[i + 1]) / h1;
		}

		for (Eigen::Index r = 0; r < 2; ++r) {
			const Eigen::Index row = static_cast<Eigen::Index>(n - 1) + r;
			const Real		  &a0  = boundary_conditions(r, 0);
			const Real		  &a1  = boundary_conditions(r, 1);
			const Real		  &a2  = boundary_conditions(r, 2);
			const Real		  &a3  = boundary_conditions(r, 3);
			const Real		   h0  = history.widths.front();
			const Real		   hn  = history.widths.back();
			history.system(row, 0) = a0 - h0 * a2 / 3;
			history.system(row, 1) = -h0 * a2 / 6;
			history.system(row, static_cast<Eigen::Index>(n - 1)) += hn * a1 / 6;
			history.system(row, static_cast<Eigen::Index>(n)) += hn * a1 / 3 + a3;
			history.constants[row] = boundary_conditions(r, 4) - a2 * (values[1] - values[0]) / h0 -
									 a1 * (values[n] - values[n - 1]) / hn;
		}

		history.derivatives = history.system.partialPivLu().solve(history.constants);
		for (SizeType i = 0; i < n; ++i) {
			const Real h  = history.widths[i];
			const Real m0 = history.derivatives[static_cast<Eigen::Index>(i)];
			const Real m1 = history.derivatives[static_cast<Eigen::Index>(i + 1)];
			const Real a  = (m1 - m0) / (6 * h);
			const Real b  = m0 / 2;
			const Real c  = (values[i + 1] - values[i]) / h - h * (2 * m0 + m1) / 6;
			history.piecewise_poly.push_back({values[i], c, b, a});
		}
		return history;
	}

	namespace {
		auto calculate_cubic_spline(const SplineHistory &history) -> SplineFunction {
			return {history.nodes, history.piecewise_poly};
		}
	} // namespace

	auto find_cubic_spline(const RealList					 &nodes,
						   const RealList					 &values,
						   const Eigen::Matrix<double, 2, 5> &boundary_conditions)
		-> SplineFunction {
		return calculate_cubic_spline(
			build_cubic_spline_history(nodes, values, boundary_conditions));
	}

	auto render_cubic_spline(const SplineHistory &history) -> void {
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

		fmt::print("\nCubic spline system [A | B]:\n");
		std::vector<SizeType> column_widths(static_cast<SizeType>(history.system.cols()), 0);
		for (Eigen::Index col = 0; col < history.system.cols(); ++col) {
			for (Eigen::Index row = 0; row < history.system.rows(); ++row) {
				column_widths[static_cast<SizeType>(col)] =
					std::max(column_widths[static_cast<SizeType>(col)],
							 fmt::format("{:.5g}", history.system(row, col)).size());
			}
		}
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
		fmt::print("{}", border);
		fmt::print("\t|    {}    |\n", std::string(row_width, ' '));
		for (const auto &row : matrix_rows)
			fmt::print("\t|    {:>{}}    |\n", row, row_width);
		fmt::print("\t|    {}    |\n", std::string(row_width, ' '));
		fmt::print("{}", border);

		fmt::print("\nNodal second derivatives (solution of AX = B):\n");
		constexpr Eigen::Index kMaxDerivativeColumns = 10;
		for (Eigen::Index start = 0; start < history.derivatives.size();
			 start += kMaxDerivativeColumns) {
			std::vector<std::vector<std::string>> derivative_rows = {{"i"}, {"m_i = S''(x_i)"}};
			const Eigen::Index					  end =
				std::min(start + kMaxDerivativeColumns, history.derivatives.size());
			for (Eigen::Index i = start; i < end; ++i) {
				derivative_rows[0].push_back(fmt::format("{}", i));
				derivative_rows[1].push_back(
					fmt::format("{:.{}g}", history.derivatives[i], config::kDisplayPrecision));
			}
			ui::render_table(derivative_rows, config::kTableMinWidth, transform_table);
		}

		fmt::print("\nPiecewise cubic coefficients:\n");
		fmt::print("a_i = (m_(i+1) - m_i)/(6h_i)\r\n");
		fmt::print("b_i = m_i/2h_i\r\n");
		fmt::print("c_i = (y_(i+1) - y_i) / h_i\n");
		fmt::print("d_i = y_i\n");
		std::vector<std::vector<std::string>> rows = {
			{"i", "Interval", "a_i", "b_i", "c_i", "d_i"}};
		for (SizeType i = 0; i < history.piecewise_poly.size(); ++i) {
			const auto &p = history.piecewise_poly[i];
			const auto	f = [](Real v) {
				return fmt::format("{:.{}g}", v, config::kDisplayPrecision);
			};
			rows.push_back({fmt::format("{}", i),
							fmt::format("[{}, {}]", f(history.nodes[i]), f(history.nodes[i + 1])),
							f(p[3]),
							f(p[2]),
							f(p[1]),
							f(p[0])});
		}
		ui::render_table(rows, config::kTableMinWidth, transform_table);
	}
} // namespace detail

auto display_cubic_spline(const RealList					&nodes,
						  const RealList					&values,
						  const Eigen::Matrix<double, 2, 5> &boundary_conditions) -> void {
	ui::DisplayFrame display_frame;
	const auto history = detail::build_cubic_spline_history(nodes, values, boundary_conditions);
	detail::render_cubic_spline(history);
	SplineFunction spline = detail::calculate_cubic_spline(history);

	fmt::print("\nCubic spline function: S(x) = S_i(x - x_i) [x_i <= x <= x_(i+1)]\r\n");
	fmt::print("{}", spline);
}
