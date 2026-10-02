#pragma once

#include <math_object.hpp>

namespace detail {

	auto evaluate_polynomial(Real x, const Polynomial &poly) -> Real;
	auto print_polynomial_evaluation(Real x, const Polynomial &poly, Real value) -> void;
	auto horner_division(const Polynomial &a, const Polynomial &b) -> PolynomialDivision;
	auto print_horner_division(const Polynomial			&a,
							   const Polynomial			&b,
							   const PolynomialDivision &division) -> void;

} // namespace detail

auto display_polynomial_evaluation(Real x, const Polynomial &poly) -> void;

// Divides a by b, where b has degree 1.
auto display_horner_division(const Polynomial &a, const Polynomial &b) -> void;
