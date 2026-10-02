#include <chapter1/horner.hpp>
#include <config.hpp>
#include <math_object.hpp>

#include <fmt/format.h>

#include <cmath>

namespace detail {

	auto evaluate_polynomial(Real x, const Polynomial &poly) -> Real {
		Real value = 0;

		for (auto it = poly.rbegin(); it != poly.rend(); ++it) {
			value = *it + value * x;
		}

		return value;
	}

	auto print_polynomial_evaluation(Real x, const Polynomial &poly, Real value) -> void {
		fmt::println("P(x) = {}", poly);
		fmt::println(
			"P({:.{}g}) = {:.{}g}", x, config::kDisplayPrecision, value, config::kDisplayPrecision);
	}

	auto horner_division(const Polynomial &a, const Polynomial &b) -> PolynomialDivision {
		PolynomialDivision division{};
		if (a.empty() || b.size() != 2)
			return division;

		division.quotient.resize(a.size() - 1);
		if (division.quotient.empty()) {
			division.remainder = a.front();
			return division;
		}

		const Real lead = b.back();
		const Real c	= -b.front() / lead;

		division.quotient.back() = a.back();

		Real t = a.back();
		for (SizeType i = a.size() - 2; i > 0; --i) {
			t						 = a[i] + t * c;
			division.quotient[i - 1] = t;
		}

		division.remainder = a.front() + t * c;

		if (lead != 1.0)
			for (auto &coeff : division.quotient)
				coeff /= lead;

		return division;
	}

	auto print_horner_division(const Polynomial			&a,
							   const Polynomial			&b,
							   const PolynomialDivision &division) -> void {
		fmt::println("A(x) = {}", a);
		fmt::println("B(x) = {}", b);
		fmt::println("A(x) = B(x) * Q(x) + R");
		fmt::println("Q(x) = {}", division.quotient);

		if (std::abs(division.remainder) < config::kNearZeroThreshold)
			fmt::println("R = 0");
		else
			fmt::println("R = {:.{}g}", division.remainder, config::kDisplayPrecision);
	}

} // namespace detail

auto display_polynomial_evaluation(Real x, const Polynomial &poly) -> void {
	Real value = detail::evaluate_polynomial(x, poly);
	detail::print_polynomial_evaluation(x, poly, value);
}

auto display_horner_division(const Polynomial &a, const Polynomial &b) -> void {
	auto division = detail::horner_division(a, b);
	detail::print_horner_division(a, b, division);
}
