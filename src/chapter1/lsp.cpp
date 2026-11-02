#include <chapter1/lsp.hpp>
#include <io/output_handling.hpp>

#include <Eigen/LU>
#include <config.hpp>
#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
	auto format_table_value(Real value) -> std::string {
		if (std::abs(value) < config::kNearZeroThreshold)
			return "0";
		return fmt::format("{:.{}g}", value, config::kDisplayPrecision);
	}

	auto format_lsp_error(LspErrorType type) -> std::string {
		if (type == LspErrorType::SSE)
			return "SSE";
		if (type == LspErrorType::MSE)
			return "MSE";
		if (type == LspErrorType::RMSE)
			return "RMSE";
		return "";
	}

	auto print_matrix_box(const Eigen::MatrixXd &matrix) -> void {
		std::vector<SizeType>				  widths(static_cast<SizeType>(matrix.cols()), 0);
		std::vector<std::vector<std::string>> rows(
			static_cast<SizeType>(matrix.rows()),
			std::vector<std::string>(static_cast<SizeType>(matrix.cols())));
		for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
			for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
				auto &value = rows[static_cast<SizeType>(row)][static_cast<SizeType>(col)];
				value		= fmt::format("{:.{}g}", matrix(row, col), config::kDisplayPrecision);
				widths[static_cast<SizeType>(col)] =
					std::max(widths[static_cast<SizeType>(col)], value.size());
			}
		}

		std::vector<std::string> lines;
		for (const auto &row : rows) {
			std::string line = "[ ";
			for (SizeType col = 0; col < row.size(); ++col) {
				if (col != 0)
					line += "    ";
				line += fmt::format("{:>{}}", row[col], widths[col]);
			}
			line += " ]";
			lines.push_back(std::move(line));
		}
		SizeType width = 0;
		for (const auto &line : lines)
			width = std::max(width, line.size());
		const std::string border = "\t+" + std::string(width + 8, '-') + "+\n";
		fmt::print("{}\t|    {}    |\n", border, std::string(width, ' '));
		for (const auto &line : lines)
			fmt::print("\t|    {:>{}}    |\n", line, width);
		fmt::print("\t|    {}    |\n{}", std::string(width, ' '), border);
	}
} // namespace

namespace detail {
	auto find_design_matrix(const RealList &nodes, const FunctionsList &basis) -> Eigen::MatrixXd {
		Eigen::MatrixXd design(static_cast<Eigen::Index>(nodes.size()),
							   static_cast<Eigen::Index>(basis.size()));
		for (SizeType i = 0; i < nodes.size(); ++i)
			for (SizeType j = 0; j < basis.size(); ++j)
				design(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
					basis[j](nodes[i]);
		return design;
	}

	auto find_approximation_coefficients(const RealList		 &nodes,
										 const RealList		 &values,
										 const FunctionsList &basis) -> RealList {
		if (nodes.size() != values.size())
			throw std::invalid_argument("Least-squares nodes and values must have the same size");
		if (nodes.empty() || basis.empty())
			throw std::invalid_argument(
				"Least-squares approximation requires data and basis functions");

		const Eigen::MatrixXd design = find_design_matrix(nodes, basis);
		const Eigen::VectorXd target = Eigen::Map<const Eigen::VectorXd>(
			values.data(), static_cast<Eigen::Index>(values.size()));
		const Eigen::MatrixXd normal_matrix = design.transpose() * design;
		const Eigen::VectorXd coefficients =
			normal_matrix.partialPivLu().solve(design.transpose() * target);

		return RealList(coefficients.data(), coefficients.data() + coefficients.size());
	}

	auto calculate_sse(const RealList	   &nodes,
					   const RealList	   &values,
					   const FunctionsList &basis,
					   const RealList	   &coeffs) -> Real {
		if (nodes.size() != values.size())
			throw std::invalid_argument("Least-squares nodes and values must have the same size");
		if (nodes.empty())
			throw std::invalid_argument("Least-squares evaluation requires data points");
		if (basis.size() != coeffs.size())
			throw std::invalid_argument("Basis functions and coefficients must have the same size");
		if (basis.empty())
			throw std::invalid_argument(
				"Least-squares evaluation requires basis functions and coefficients");

		Real sse = 0.0;
		for (SizeType i = 0; i < nodes.size(); ++i) {
			Real f_x = 0.0;
			for (SizeType j = 0; j < basis.size(); ++j)
				f_x += coeffs[j] * basis[j](nodes[i]);
			const Real diff = f_x - values[i];
			sse += diff * diff;
		}
		return sse;
	}

	auto calculate_mse(const RealList	   &nodes,
					   const RealList	   &values,
					   const FunctionsList &basis,
					   const RealList	   &coeffs) -> Real {
		return calculate_sse(nodes, values, basis, coeffs) / static_cast<Real>(nodes.size());
	}

	auto calculate_rmse(const RealList		&nodes,
						const RealList		&values,
						const FunctionsList &basis,
						const RealList		&coeffs) -> Real {
		return std::sqrt(calculate_mse(nodes, values, basis, coeffs));
	}

	auto render_lsp_error(const RealList	  &nodes,
						  const RealList	  &values,
						  const FunctionsList &basis,
						  const RealList	  &coeffs,
						  LspErrorType		   error_type) -> void {
		if (nodes.size() != values.size()) {
			throw std::invalid_argument("Least-squares nodes and values must have the same size");
		}
		if (nodes.empty()) {
			throw std::invalid_argument("Least-squares evaluation requires data points");
		}
		if (basis.size() != coeffs.size()) {
			throw std::invalid_argument("Basis functions and coefficients must have the same size");
		}
		if (basis.empty()) {
			throw std::invalid_argument(
				"Least-squares evaluation requires basis functions and coefficients");
		}

		fmt::print("\nCalculate error {}:\n", format_lsp_error(error_type));
		std::vector<std::vector<std::string>> rows = {
			{"i", "x_i", "y_i", "F(x_i)", "(F(x_i) - y_i)^2"}};
		rows.reserve(nodes.size() + 1);

		for (SizeType i = 0; i < nodes.size(); ++i) {
			Real f_xi = 0.0;
			for (SizeType j = 0; j < basis.size(); ++j) {
				f_xi += coeffs[j] * basis[j](nodes[i]);
			}
			const Real diff	   = f_xi - values[i];
			const Real sq_diff = diff * diff;
			rows.push_back({
				fmt::format("{}", i + 1),
				format_table_value(nodes[i]),
				format_table_value(values[i]),
				format_table_value(f_xi),
				format_table_value(sq_diff),
			});
		}

		auto transform_table = [](ftxui::Table &table) {
			table.SelectColumn(0).DecorateCells(ftxui::color(ftxui::Color::Yellow));
			table.SelectColumn(0).DecorateCells(ftxui::bold);
		};

		ui::render_table(rows, config::kTableMinWidth, transform_table);

		std::string formula;
		Real		result = 0.0;
		switch (error_type) {
		case LspErrorType::SSE:
			formula = fmt::format("∑[i=1,{}] (F(x_i) - y_i)^2", nodes.size());
			result	= calculate_sse(nodes, values, basis, coeffs);
			break;
		case LspErrorType::MSE:
			formula = fmt::format("1/n * ∑[i=1,{}] (F(x_i) - y_i)^2", nodes.size());
			result	= calculate_mse(nodes, values, basis, coeffs);
			break;
		case LspErrorType::RMSE:
			formula =
				fmt::format("sqrt(1/{} * ∑[i=1,{}] (F(x_i) - y_i)^2)", nodes.size(), nodes.size());
			result = calculate_rmse(nodes, values, basis, coeffs);
			break;
		}

		fmt::print("\nError: ε = {} = {:.{}g}\n", formula, result, config::kDisplayPrecision);
	}
} // namespace detail

auto display_linear_lsp(const RealList &nodes, const RealList &values, LspErrorType error_type)
	-> void {
	ui::DisplayFrame	display_frame;
	const FunctionsList basis = {
		[](Real x) { return x; },
		[](Real) { return 1.0; },
	};
	const Eigen::MatrixXd design = detail::find_design_matrix(nodes, basis);
	const RealList		  coeffs = detail::find_approximation_coefficients(nodes, values, basis);
	Eigen::MatrixXd		  coeffs_matrix(static_cast<Eigen::Index>(coeffs.size()), 1);
	const Eigen::VectorXd target =
		Eigen::Map<const Eigen::VectorXd>(values.data(), static_cast<Eigen::Index>(values.size()));
	for (SizeType i = 0; i < coeffs.size(); ++i)
		coeffs_matrix(static_cast<Eigen::Index>(i), 0) = coeffs[i];

	fmt::print("Linear least-squares design matrix Φ (order: | x | 1 |):\n");
	print_matrix_box(design);
	fmt::print("Matrix M = Φ^T * Φ:\n");
	print_matrix_box(design.transpose() * design);
	fmt::print("Matrix Φ^T * y:\n");
	print_matrix_box(design.transpose() * target);
	fmt::print("\nLinear least-squares coefficient matrix [a b]^T = (Φ^T * Φ)^(-1) * Φ^T * y\n");
	print_matrix_box(coeffs_matrix);
	fmt::print("\nApproximation function F(x) = {}\n", Polynomial({coeffs[1], coeffs[0]}));
	fmt::print("\n");
	detail::render_lsp_error(nodes, values, basis, coeffs, error_type);
}

auto display_quadratic_lsp(const RealList &nodes, const RealList &values, LspErrorType error_type)
	-> void {
	ui::DisplayFrame	display_frame;
	const FunctionsList basis = {
		[](Real x) { return x * x; },
		[](Real x) { return x; },
		[](Real) { return 1.0; },
	};

	const Eigen::MatrixXd design = detail::find_design_matrix(nodes, basis);
	const RealList		  coeffs = detail::find_approximation_coefficients(nodes, values, basis);
	Eigen::MatrixXd		  coeffs_matrix(static_cast<Eigen::Index>(coeffs.size()), 1);
	const Eigen::VectorXd target =
		Eigen::Map<const Eigen::VectorXd>(values.data(), static_cast<Eigen::Index>(values.size()));
	for (SizeType i = 0; i < coeffs.size(); ++i)
		coeffs_matrix(static_cast<Eigen::Index>(i), 0) = coeffs[i];

	fmt::print("Quadratic least-squares design matrix Φ (order: | x^2 | x | 1 |):\n");
	print_matrix_box(design);
	fmt::print("Matrix M = Φ^T * Φ:\n");
	print_matrix_box(design.transpose() * design);
	fmt::print("Matrix Φ^T * y:\n");
	print_matrix_box(design.transpose() * target);
	fmt::print("\nQuadratic least-squares coefficient matrix [a b c]^T = (Φ^T * Φ)^(-1) * "
			   "Φ^T * y\n:\n");
	print_matrix_box(coeffs_matrix);
	fmt::print("\nApproximation function F(x) = {}\n",
			   Polynomial({coeffs[2], coeffs[1], coeffs[0]}));
	fmt::print("\n");
	detail::render_lsp_error(nodes, values, basis, coeffs, error_type);
}

auto display_general_lsp(const RealList		 &nodes,
						 const RealList		 &values,
						 const FunctionsList &basis,
						 LspErrorType		  error_type) -> void {
	ui::DisplayFrame	  display_frame;
	const Eigen::MatrixXd design = detail::find_design_matrix(nodes, basis);
	const RealList		  coeffs = detail::find_approximation_coefficients(nodes, values, basis);
	Eigen::MatrixXd		  coeffs_matrix(static_cast<Eigen::Index>(coeffs.size()), 1);
	const Eigen::VectorXd target =
		Eigen::Map<const Eigen::VectorXd>(values.data(), static_cast<Eigen::Index>(values.size()));
	for (SizeType i = 0; i < coeffs.size(); ++i)
		coeffs_matrix(static_cast<Eigen::Index>(i), 0) = coeffs[i];

	fmt::print("General least-squares design matrix Φ [φ_(ij) = Φ_j(x_i)]:\n");
	print_matrix_box(design);
	fmt::print("Matrix M = Φ^T * Φ:\n");
	print_matrix_box(design.transpose() * design);
	fmt::print("Matrix Φ^T * y:\n");
	print_matrix_box(design.transpose() * target);
	fmt::print("\nGeneral least-squares coefficient matrix [alpha_i]^T:\n");
	print_matrix_box(coeffs_matrix);

	std::string approximation = "F(x) = ";
	for (SizeType i = 0; i < coeffs.size(); ++i) {
		if (i != 0)
			approximation += " + ";
		approximation += fmt::format(
			"{}*Φ_{}(x)", fmt::format("{:.{}g}", coeffs[i], config::kDisplayPrecision), i);
	}
	fmt::print("\nApproximation function: {}\n", approximation);
	fmt::print("\n");
	detail::render_lsp_error(nodes, values, basis, coeffs, error_type);
}
