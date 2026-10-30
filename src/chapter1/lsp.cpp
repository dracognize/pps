#include <chapter1/lsp.hpp>

#include <Eigen/LU>
#include <config.hpp>
#include <fmt/format.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
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
} // namespace detail

auto display_linear_lsp(const RealList &nodes, const RealList &values) -> void {
	const FunctionsList basis = {
		[](Real x) { return x; },
		[](Real) { return 1.0; },
	};
	const Eigen::MatrixXd design = detail::find_design_matrix(nodes, basis);
	const RealList		  coeffs = detail::find_approximation_coefficients(nodes, values, basis);
	Eigen::MatrixXd		  coeffs_matrix(static_cast<Eigen::Index>(coeffs.size()), 1);
	for (SizeType i = 0; i < coeffs.size(); ++i)
		coeffs_matrix(static_cast<Eigen::Index>(i), 0) = coeffs[i];

	fmt::print("Linear least-squares design matrix Phi:\n");
	print_matrix_box(design);
	fmt::print(
		"\nLinear least-squares coefficient matrix [a b]^T = (Phi^T * Phi)^(-1) * Phi^T * y\n");
	print_matrix_box(coeffs_matrix);
	fmt::print("\nApproximation funciton F(x) = {}\n", Polynomial({coeffs[1], coeffs[0]}));
}

auto display_quadratic_lsp(const RealList &nodes, const RealList &values) -> void {
	const FunctionsList basis = {
		[](Real x) { return x * x; },
		[](Real x) { return x; },
		[](Real) { return 1.0; },
	};

	const Eigen::MatrixXd design = detail::find_design_matrix(nodes, basis);
	const RealList		  coeffs = detail::find_approximation_coefficients(nodes, values, basis);
	Eigen::MatrixXd		  coeffs_matrix(static_cast<Eigen::Index>(coeffs.size()), 1);
	for (SizeType i = 0; i < coeffs.size(); ++i)
		coeffs_matrix(static_cast<Eigen::Index>(i), 0) = coeffs[i];

	fmt::print("Quadratic least-squares design matrix Phi:\n");
	print_matrix_box(design);
	fmt::print("\nQuadratic least-squares coefficient matrix [a b c]^T = (Phi^T * Phi)^(-1) * "
			   "Phi^T * y\n:\n");
	print_matrix_box(coeffs_matrix);
	fmt::print("\nApproximation function F(x) = {}\n",
			   Polynomial({coeffs[2], coeffs[1], coeffs[0]}));
}

auto display_general_lsp(const RealList &nodes, const RealList &values, const FunctionsList &basis)
	-> void {
	const Eigen::MatrixXd design = detail::find_design_matrix(nodes, basis);
	const RealList		  coeffs = detail::find_approximation_coefficients(nodes, values, basis);
	Eigen::MatrixXd		  coeffs_matrix(static_cast<Eigen::Index>(coeffs.size()), 1);
	for (SizeType i = 0; i < coeffs.size(); ++i)
		coeffs_matrix(static_cast<Eigen::Index>(i), 0) = coeffs[i];

	fmt::print("General least-squares design matrix Phi:\n");
	print_matrix_box(design);
	fmt::print("\nGeneral least-squares coefficient matrix [alpha_i]^T:\n");
	print_matrix_box(coeffs_matrix);

	std::string approximation = "F(x) = ";
	for (SizeType i = 0; i < coeffs.size(); ++i) {
		if (i != 0)
			approximation += " + ";
		approximation += fmt::format(
			"{}*phi_{}(x)", fmt::format("{:.{}g}", coeffs[i], config::kDisplayPrecision), i);
	}
	fmt::print("\nApproximation function: {}\n", approximation);
}
