#include <chapter1/interpolation/horner.hpp>
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
		fmt::print("P({:.{}g}) = {:.{}g}\r\n",
				   x,
				   config::kDisplayPrecision,
				   value,
				   config::kDisplayPrecision);
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
			history.remainders.push_back(division.remainder);
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
								 const DerivativeHistory &history) -> void {
		fmt::print("P(x) = {}\r\n", poly);
		fmt::print("x = {:.{}g}, k = {}\r\n", c, config::kDisplayPrecision, k);

		std::vector<std::vector<std::string>> rows{};
		rows.reserve(static_cast<size_t>(k) + 2);

		std::vector<std::string> header{};
		header.reserve(poly.size() + 2);
		header.emplace_back("Coefficient");
		for (SizeType j = poly.size(); j > 0; --j)
			header.push_back(fmt::format("a{}", j - 1));
		header.emplace_back("Remainder");
		rows.push_back(header);

		for (SizeType i = 0; i <= k; ++i) {
			std::vector<std::string> row{};
			row.reserve(poly.size() + 2);
			row.push_back(fmt::format("k = {}", i));
			for (SizeType j = 0; j < i; j++) {
				row.emplace_back("");
			}
			const Polynomial &quotient = history.quotients[static_cast<size_t>(i)];
			for (auto it = quotient.rbegin(); it != quotient.rend(); ++it) {
				row.push_back(format_table_value(*it));
			}
			row.push_back(format_table_value(history.remainders[static_cast<size_t>(i)]));
			rows.push_back(row);
		}
		rows[1][poly.size() + 1] = "";

		ui::render_table(rows, config::kTableMinWidth);
	}

	auto multiply_by_linear(const Polynomial &a, const Polynomial &b) -> Polynomial {
		if (a.empty() || b.size() != 2 || b.back() == 0)
			return {};

		Polynomial multiplication(a.size() + 1);
		const Real lead = b.back();
		const Real c	= -b.front() / lead;

		for (SizeType i = 0; i < a.size(); ++i) {
			multiplication[i + 1] = a[i];
			multiplication[i] -= a[i] * c;
		}

		if (lead != 1.0)
			for (auto &coeff : multiplication)
				coeff /= lead;
		return multiplication;
	}

	auto print_horner_multiplication(const Polynomial &a,
									 const Polynomial &b,
									 const Polynomial &product) -> void {
		fmt::print("A(x) = {}\r\n", a);
		fmt::print("B(x) = {}\r\n", b);
		fmt::print("A(x) * B(x) = {}\r\n", product);
	}

	auto node_polynomial(const RealList &nodes) -> Polynomial {
		if (nodes.empty())
			return {1.0};

		Polynomial result{1.0};
		for (const Real c : nodes) {
			const Polynomial divisor{-c, 1.0};
			result = multiply_by_linear(result, divisor);
		}
		return result;
	}

	auto build_node_polynomial_history(const RealList &nodes) -> NodePolynomialHistory {
		NodePolynomialHistory history{};
		history.products.reserve(nodes.size() + 1);

		history.products.push_back({1.0});

		for (const Real c : nodes) {
			const Polynomial multiplier{-c, 1.0};
			auto			 current = multiply_by_linear(history.products.back(), multiplier);
			history.products.push_back(current);
		}
		return history;
	}

	auto render_nodes_polynomial_table(const RealList &nodes, const NodePolynomialHistory &history)
		-> void {
		SizeType k = history.products.back().size() - 1;
		SizeType n = nodes.size();

		std::vector<std::vector<std::string>> rows{};
		rows.reserve(n + 2);

		std::vector<std::string> header{};
		header.reserve(k + 3);

		header.emplace_back("");
		header.emplace_back("");
		for (SizeType j = 0; j <= k; ++j) {
			header.push_back(fmt::format("a{}", k - j));
		}
		rows.push_back(header);

		std::vector<std::string> row{};
		row.reserve(k + 3);

		row.emplace_back("i");
		row.emplace_back("x_i");
		for (SizeType j = 0; j < k; ++j) {
			row.push_back("");
		}
		row.push_back(format_table_value(history.products[0][0]));
		rows.push_back(row);

		for (SizeType i = 0; i < n; ++i) {
			row.clear();
			row.push_back(fmt::format("{}", i));
			row.push_back(format_table_value(nodes[i]));
			for (SizeType j = 0; j < k - i - 1; ++j) {
				row.push_back("");
			}
			const Polynomial &product = history.products[i + 1];
			for (auto it = product.rbegin(); it != product.rend(); ++it) {
				row.push_back(format_table_value(*it));
			}
			rows.push_back(row);
		}

		ui::render_table(rows, config::kTableMinWidth);
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

	auto n		 = detail::effective_degree(poly);
	auto history = detail::build_derivative_history(c, poly, k);
	detail::render_derivative_table(c, poly, std::min(n, k), history);

	for (SizeType i = 0; i <= std::min(n, k); i++) {
		auto value_i	= detail::evaluate_polynomial(c, history.quotients[i]);
		auto derivative = detail::factorial(i) * value_i;

		fmt::print("k = {}\r\n", i);
		if (std::abs(value_i) < config::kNearZeroThreshold &&
			std::abs(derivative) < config::kNearZeroThreshold) {
			fmt::print("P_{}({:.{}g}) = 0\r\n", n - i, c, config::kDisplayPrecision);
			fmt::print("P^({})({:.{}g}) = {}! * P_{}({:.{}g}) = 0\r\n",
					   i,
					   c,
					   config::kDisplayPrecision,
					   i,
					   n - i,
					   c,
					   config::kDisplayPrecision);
		} else {
			fmt::print("P_{}({:.{}g}) = {:.{}g}\r\n",
					   n - i,
					   c,
					   config::kDisplayPrecision,
					   value_i,
					   config::kDisplayPrecision);
			fmt::print("P^({})({:.{}g}) = {}! * P_{}({:.{}g}) = {:.{}g}\r\n\n",
					   i,
					   c,
					   config::kDisplayPrecision,
					   i,
					   n - i,
					   c,
					   config::kDisplayPrecision,
					   derivative,
					   config::kDisplayPrecision);
		}
	}
	if (k > n) {
		fmt::print("k >= {}\r\n", n + 1);
		fmt::print("k > deg(P), so P^(k)({}) = 0\r\n", c);
		return;
	}
}

auto display_horner_multiplication(const Polynomial &a, const Polynomial &b) -> void {
	const auto product = detail::multiply_by_linear(a, b);
	detail::print_horner_multiplication(a, b, product);
}

auto display_nodes_polynomial_table(const RealList &nodes) -> void {
	if (nodes.empty()) {
		fmt::print("No nodes provided.\r\n");
		return;
	}

	fmt::print("Nodes:  ");
	for (const auto &node : nodes) {
		fmt::print("{:.{}g}  ", node, config::kDisplayPrecision);
	}
	fmt::print("\r\n\r\n");
	fmt::print("Node polynomials' coefficients\r\n");
	const auto history = detail::build_node_polynomial_history(nodes);
	detail::render_nodes_polynomial_table(nodes, history);

	fmt::print("Node Polynomial: w_{}(x) = {}\r\n", nodes.size(), history.products.back());
}
