#include <chapter1/horner.hpp>
#include <config.hpp>
#include <math_object.hpp>
#include <ui/table.hpp>

#include <fmt/format.h>

#include <cmath>
#include <string>
#include <vector>

namespace detail {

	auto evaluate_polynomial(Real x, const Polynomial &poly) -> Real {
		Real value = 0;

		for (auto it = poly.rbegin(); it != poly.rend(); ++it) {
			value = *it + value * x;
		}

		return value;
	}

	auto print_polynomial_evaluation(Real x, const Polynomial &poly, Real value) -> void {
		fmt::print("P(x) = {}\r\n", poly);
		fmt::print(
			"P({:.{}g}) = {:.{}g}\r\n", x, config::kDisplayPrecision, value, config::kDisplayPrecision);
	}

	auto divide_by_linear(const Polynomial &a, const Polynomial &b) -> PolynomialDivision {
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
		fmt::print("A(x) = {}\r\n", a);
		fmt::print("B(x) = {}\r\n", b);
		fmt::print("A(x) = B(x) * Q(x) + R\r\n", division.quotient);
		fmt::print("Q(x) = {}\r\n", division.quotient);

		if (std::abs(division.remainder) < config::kNearZeroThreshold)
			fmt::println("R = 0");
		else
			fmt::println("R = {:.{}g}", division.remainder, config::kDisplayPrecision);
	}

	auto factorial(SizeType k) -> Real {
		Real result = 1;
		for (SizeType i = 2; i <= k; ++i)
			result *= static_cast<Real>(i);
		return result;
	}

	auto effective_degree(const Polynomial &poly) -> SizeType {
		for (SizeType i = poly.size(); i > 0; --i) {
			if (std::abs(poly[i - 1]) >= config::kNearZeroThreshold)
				return i - 1;
		}
		return 0;
	}

	auto kth_derivative_value(Real c, const Polynomial &poly, SizeType k) -> Real {
		if (poly.empty())
			return 0;

		if (k > effective_degree(poly))
			return 0;

		Polynomial		 current = poly;
		const Polynomial divisor{-c, 1};
		for (SizeType j = 0; j < k; ++j) {
			auto division = divide_by_linear(current, divisor);
			current		  = division.quotient;
			if (current.empty())
				return 0;
		}

		return factorial(k) * evaluate_polynomial(c, current);
	}

	auto build_derivative_history(Real c, const Polynomial &poly, SizeType k) -> DerivativeHistory {
		DerivativeHistory history{};
		history.quotients.reserve(static_cast<std::size_t>(k) + 1);
		history.remainders.reserve(static_cast<std::size_t>(k) + 1);

		history.quotients.push_back(poly);
		history.remainders.push_back(evaluate_polynomial(c, poly));

		const Polynomial divisor{-c, 1};
		for (SizeType j = 1; j <= k; ++j) {
			auto division = divide_by_linear(history.quotients.back(), divisor);
			history.quotients.push_back(division.quotient);
			history.remainders.push_back(evaluate_polynomial(c, division.quotient));
		}

		return history;
	}

	namespace {

		auto format_table_value(Real value) -> std::string {
			if (std::abs(value) < config::kNearZeroThreshold)
				return "0";
			return fmt::format("{:.{}g}", value, config::kDisplayPrecision);
		}

	} // namespace

	auto render_derivative_table(Real					  c,
								 const Polynomial		 &poly,
								 SizeType				  k,
								 const DerivativeHistory &history,
								 Real					  derivative) -> void {
		fmt::print("P(x) = {}\r\n", poly);
		fmt::print("x = {:.{}g}, k = {}\r\n", c, config::kDisplayPrecision, k);

		std::vector<std::vector<std::string>> rows{};
		rows.reserve(poly.size() + 2);

		std::vector<std::string> header{};
		header.reserve(static_cast<std::size_t>(k) + 2);
		header.emplace_back("Coefficient");
		for (SizeType j = 0; j <= k; ++j)
			header.push_back(fmt::format("k = {}", j));
		rows.push_back(header);

		for (SizeType i = 0; i < poly.size(); ++i) {
			std::vector<std::string> row{};
			row.reserve(static_cast<std::size_t>(k) + 2);
			row.push_back(fmt::format("a{}", i));
			for (SizeType j = 0; j <= k; ++j) {
				const Polynomial &quotient = history.quotients[static_cast<std::size_t>(j)];
				if (i < quotient.size())
					row.push_back(format_table_value(quotient[i]));
				else
					row.emplace_back("");
			}
			rows.push_back(row);
		}

		std::vector<std::string> result_row{};
		result_row.reserve(static_cast<std::size_t>(k) + 2);
		result_row.emplace_back("Result");
		for (SizeType j = 0; j <= k; ++j)
			result_row.push_back(
				format_table_value(history.remainders[static_cast<std::size_t>(j)]));
		rows.push_back(result_row);

		ui::render_table(rows, config::kTableMinWidth);

		const Real remainder_k = history.remainders[static_cast<std::size_t>(k)];
		if (std::abs(remainder_k) < config::kNearZeroThreshold &&
			std::abs(derivative) < config::kNearZeroThreshold) {
			fmt::print("R = 0\r\n");
			fmt::print("P^({})(c) = {}! * R = 0\r\n", k, k);
		} else {
			fmt::print("R = {:.{}g}\r\n", remainder_k, config::kDisplayPrecision);
			fmt::print("P^({})({:.{}g}) = {}! * R = {:.{}g}\r\n",
						 k,
						 c,
						 config::kDisplayPrecision,
						 k,
						 derivative,
						 config::kDisplayPrecision);
		}
	}

} // namespace detail

auto display_polynomial_evaluation(Real x, const Polynomial &poly) -> void {
	const Real value = detail::evaluate_polynomial(x, poly);
	detail::print_polynomial_evaluation(x, poly, value);
}

auto display_horner_division(const Polynomial &a, const Polynomial &b) -> void {
	const auto division = detail::divide_by_linear(a, b);
	detail::print_horner_division(a, b, division);
}

auto display_kth_derivative(Real c, const Polynomial &poly, SizeType k) -> void {
	if (poly.empty()) {
		fmt::print("P(x) = 0\r\n");
		fmt::print("P^({})({:.{}g}) = 0\r\n", k, c, config::kDisplayPrecision);
		return;
	}

	if (k > detail::effective_degree(poly)) {
		fmt::print("P(x) = {}\r\n", poly);
		fmt::print("x = {:.{}g}, k = {}\r\n", c, config::kDisplayPrecision, k);
		fmt::print("k > deg(P), so P^({})(c) = 0\r\n", k);
		return;
	}

	auto history	= detail::build_derivative_history(c, poly, k);
	auto derivative = detail::factorial(k) * history.remainders[static_cast<std::size_t>(k)];
	detail::render_derivative_table(c, poly, k, history, derivative);
}
