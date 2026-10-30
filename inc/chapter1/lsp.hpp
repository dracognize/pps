#pragma once
/// @file lsp.hpp
/// @brief Least square problem construction and display.

#include <math_object.hpp>

namespace detail {
	/// @brief Build the design matrix for the supplied basis functions.
	/// @details `design(i, j) = basis[j](nodes[i])`.
	/// @param nodes Sample input points.
	/// @param basis Basis functions, one per design-matrix column.
	/// @return Design matrix with one row per sample and one column per basis function.
	auto find_design_matrix(const RealList &nodes, const FunctionsList &basis) -> Eigen::MatrixXd;

	/// @brief Find least-squares approximation coefficients for supplied basis functions.
	/// @param nodes Sample input points.
	/// @param values Sample values corresponding to `nodes`.
	/// @param basis Basis functions mapping a `Real` input to a `Real` value.
	/// @return Coefficients corresponding to `basis`.
	auto find_approximation_coefficients(const RealList		 &nodes,
										 const RealList		 &values,
										 const FunctionsList &basis) -> RealList;
} // namespace detail

/// @brief Display the linear least-squares design matrix and fitted coefficients.
/// @details Uses the basis `{x, 1}` and reports `F(x) = a*x + b`.
/// @param nodes Sample input points.
/// @param values Sample values corresponding to `nodes`.
auto display_linear_lsp(const RealList &nodes, const RealList &values) -> void;

/// @brief Display the quadratic least-squares design matrix and fitted coefficients.
/// @details Uses the basis `{x^2, x, 1}` and reports `F(x) = a*x^2 + b*x + c`.
/// @param nodes Sample input points.
/// @param values Sample values corresponding to `nodes`.
auto display_quadratic_lsp(const RealList &nodes, const RealList &values) -> void;

/// @brief Display a least-squares approximation over arbitrary basis functions.
/// @details Reports the design matrix and `F(x) = sum_i alpha_i * phi_i(x)`.
/// @param nodes Sample input points.
/// @param values Sample values corresponding to `nodes`.
/// @param basis Basis functions in the order used for the returned coefficients.
auto display_general_lsp(const RealList &nodes, const RealList &values, const FunctionsList &basis)
	-> void;
