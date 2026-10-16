#pragma once
/// @file horner.hpp
/// @brief Horner evaluation, linear division, k-th derivatives.

#include <math_object.hpp>

#include <string>

namespace detail {

	/// @brief Evaluate `P(x)` with the Horner scheme.
	/// @param x Point of evaluation.
	/// @param poly Ascending coefficients.
	/// @return `P(x)`, `0` when empty.
	auto evaluate_polynomial(Real x, const Polynomial &poly) -> Real;
	/// @brief Print `P(x)` and `P(x0)` as plain text.
	/// @param x Point of evaluation.
	/// @param poly Ascending coefficients.
	/// @param value Precomputed `P(x)`.
	auto print_polynomial_evaluation(Real x, const Polynomial &poly, Real value) -> void;
	/// @brief Divide `a` by degree-1 `b`.
	/// @param a Dividend.
	/// @param b Divisor with exactly 2 coefficients.
	/// @return Quotient and remainder; zero when `a` is empty or `b` is not degree 1.
	auto divide_by_linear(const Polynomial &a, const Polynomial &b) -> PolynomialDivision;
	/// @brief Print `A = B * Q + R` as plain text.
	/// @param a Dividend.
	/// @param b Degree-1 divisor.
	/// @param division Precomputed division result.
	auto print_horner_division(const Polynomial			&a,
							   const Polynomial			&b,
							   const PolynomialDivision &division) -> void;
	/// @brief Core `P^{(k)}(c)` computation.
	/// @param c Point of evaluation.
	/// @param poly Ascending coefficients.
	/// @param k Derivative order.
	/// @return `P^{(k)}(c)`, `0` when `k` exceeds the degree.
	auto kth_derivative_value(Real c, const Polynomial &poly, SizeType k) -> Real;
	/// @brief Build the full derivative tableau.
	/// @param c Point of evaluation.
	/// @param poly Ascending coefficients.
	/// @param k Derivative order (`k <= deg(P)` required).
	/// @return Quotients `Q_0..Q_k` with evaluations `Q_j(c)`.
	auto build_derivative_history(Real c, const Polynomial &poly, SizeType k) -> DerivativeHistory;
	/// @brief Render the tableau plus the `P^{(k)}(c) = k! * R` line.
	/// @param c Point of evaluation.
	/// @param poly Ascending coefficients.
	/// @param k Derivative order.
	/// @param history Precomputed tableau.
	auto render_derivative_table(Real					  c,
								 const Polynomial		 &poly,
								 SizeType				  k,
								 const DerivativeHistory &history) -> void;
	/// @brief 	Multiply `P(x)` by `(x - c)`.
	/// @param a Ascending coefficients.
	/// @param b Ascending coefficients of degree-1 polynomial (x - c).
	/// @return Ascending coefficients of the product.
	auto multiply_by_linear(const Polynomial &a, const Polynomial &b) -> Polynomial;
	/// @brief Print the result of a Horner multiplication.
	/// @param a Ascending coefficients.
	/// @param b Ascending coefficients of degree-1 polynomial (x - c).
	/// @param product Product of the two polynomials.
	auto print_horner_multiplication(const Polynomial &a,
									 const Polynomial &b,
									 const Polynomial &product) -> void;
	/// @brief Build the polynomial with roots at the given nodes.
	/// @param nodes List of roots.
	/// @return Polynomial with the specified roots.
	auto node_polynomial(const RealList &nodes) -> Polynomial;
	/// @brief Build the history of node polynomials.
	/// @param nodes List of roots.
	/// @return History of node polynomials.
	auto build_node_polynomial_history(const RealList &nodes) -> PolynomialList;
	/// @brief Render the `i | x_i | w_i(x)` table.
	/// @param nodes List of roots.
	/// @param history Precomputed node polynomial history.
	/// @param var Variable label (default "x").
	auto render_nodes_polynomial_table(const RealList		&nodes,
									   const PolynomialList &history,
									   const std::string	&var = "x") -> void;
} // namespace detail

/// @brief Print `P(x)` and `P(x0)` as plain text.
/// @param x Point of evaluation.
/// @param poly Ascending coefficients.
auto display_polynomial_evaluation(Real x, const Polynomial &poly) -> void;

/// @brief Print `A = B * Q + R` as plain text.
/// @param a Dividend.
/// @param b Degree-1 divisor.
auto display_horner_division(const Polynomial &a, const Polynomial &b) -> void;

/// @brief Render the Horner tableau and `P^{(k)}(c) = k! * R`.
/// @param c Point of evaluation.
/// @param poly Ascending coefficients.
/// @param k Derivative order.
auto display_kth_derivative(Real c, const Polynomial &poly, SizeType k) -> void;

/// @brief Multiply `P(x)` by `(x - c)` and print the result.
/// @param a Ascending coefficients.
/// @param b Ascending coefficients of degree-1 polynomial (x - c).
auto display_horner_multiplication(const Polynomial &a, const Polynomial &b) -> void;
/// @brief Render the table of node polynomials.
/// @param nodes List of roots.
/// @param var Variable label (default "x").
auto display_nodes_polynomial_table(const RealList &nodes, const std::string &var = "x") -> void;
