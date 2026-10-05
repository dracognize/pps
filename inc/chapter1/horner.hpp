#pragma once
/// @file horner.hpp
/// @brief Horner evaluation, linear division, k-th derivatives.

#include <math_object.hpp>

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

	/// @brief Factorial as `Real`.
	/// @param k Non-negative order.
	/// @return `k!`, `0! == 1`.
	auto factorial(SizeType k) -> Real;
	/// @brief Highest index with a non-negligible coefficient.
	/// @param poly Ascending coefficients.
	/// @return Effective degree, `0` for the zero polynomial.
	auto effective_degree(const Polynomial &poly) -> SizeType;
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
	/// @param derivative Precomputed `P^{(k)}(c)`.
	auto render_derivative_table(Real					  c,
								 const Polynomial		 &poly,
								 SizeType				  k,
								 const DerivativeHistory &history,
								 Real					  derivative) -> void;

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
