#include <config.hpp>

#include <math_object.hpp>
#include <math_utils.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <fmt/format.h>

/// @brief Polynomial member operators.
auto Polynomial::operator+(const Polynomial &other) const -> Polynomial {
	Polynomial result(std::max(size(), other.size()));
	for (SizeType i = 0; i < size(); ++i)
		result[i] += (*this)[i];
	for (SizeType i = 0; i < other.size(); ++i)
		result[i] += other[i];
	return result;
}
auto Polynomial::operator*(const Real &scalar) const -> Polynomial {
	Polynomial result(*this);
	for (auto &coeff : result)
		coeff *= scalar;
	return result;
}
auto operator*(const Real &scalar, const Polynomial &poly) -> Polynomial {
	return poly * scalar;
}
auto Polynomial::operator*(const Polynomial &other) const -> Polynomial {
	const SizeType deg_a = utils::effective_degree(*this);
	const SizeType deg_b = utils::effective_degree(other);
	Polynomial	   result(deg_a + deg_b + 1, 0.0);
	for (SizeType i = 0; i < other.size(); ++i) {
		if (std::abs(other[i]) < config::kNearZeroThreshold)
			continue;

		for (SizeType k = 0; k < (*this).size(); ++k) {
			if (std::abs((*this)[k]) < config::kNearZeroThreshold)
				continue;
			result[i + k] += (*this)[k] * other[i];
		}
	}
	return result;
}

auto SplineFunction::value(const Real &key) -> Real {
	if (key < nodes.front() || key > nodes.back())
		throw std::out_of_range("SplineFunction evaluation point is outside the node range");

	auto	 interval = std::lower_bound(nodes.begin(), nodes.end(), key);
	SizeType index	  = static_cast<SizeType>(interval - nodes.begin());

	const Polynomial &poly = piecewise_poly[index];

	Real result = 0;
	for (auto it = poly.rbegin(); it != poly.rend(); ++it)
		result = *it + result * (key - nodes[index]);
	return result;
}

namespace {
	auto is_zero_coeff(Real c) -> bool {
		return std::abs(c) < config::kNearZeroThreshold;
	}

	auto is_unit_coeff(Real c) -> bool {
		return std::abs(std::abs(c) - 1.0) < config::kNearZeroThreshold;
	}

	auto format_abs_coeff(Real abs_c) -> std::string {
		return fmt::format("{:.{}g}", abs_c, config::kDisplayPrecision);
	}

	auto monomial_body(Real abs_c, SizeType power, std::string_view var) -> std::string {
		if (power == 0)
			return format_abs_coeff(abs_c);
		if (power == 1)
			return is_unit_coeff(abs_c) ? std::string(var)
										: format_abs_coeff(abs_c) + "*" + std::string(var);
		return is_unit_coeff(abs_c) ? fmt::format("{}^{}", var, power)
									: fmt::format("{}*{}^{}", format_abs_coeff(abs_c), var, power);
	}

	auto
	append_term(std::string &out, bool &is_first, Real coeff, SizeType power, std::string_view var)
		-> void {
		if (is_zero_coeff(coeff))
			return;
		std::string body = monomial_body(std::abs(coeff), power, var);
		if (is_first) {
			out += (coeff < 0 ? "-" : "") + body;
			is_first = false;
		} else {
			out += (coeff < 0 ? " - " : " + ") + body;
		}
	}

} // namespace

/// @cond INTERNAL
auto fmt::formatter<Polynomial>::format(const Polynomial &poly, fmt::format_context &ctx) const
	-> decltype(ctx.out()) {
	std::vector<std::pair<SizeType, Real>> nonzero;
	for (SizeType i = 0; i < poly.size(); ++i) {
		if (!is_zero_coeff(poly[i]))
			nonzero.emplace_back(i, poly[i]);
	}

	std::string_view var = var_;
	while (!var.empty() && std::isspace(static_cast<unsigned char>(var.front())))
		var.remove_prefix(1);
	while (!var.empty() && std::isspace(static_cast<unsigned char>(var.back())))
		var.remove_suffix(1);
	if (var.empty())
		var = "x";

	std::string out;
	if (nonzero.empty()) {
		out = "0";
	} else if (nonzero.size() <= config::kMaxFullTerms) {
		bool is_first = true;
		for (auto it = nonzero.rbegin(); it != nonzero.rend(); ++it)
			append_term(out, is_first, it->second, it->first, var);
	} else {
		bool is_first = true;
		append_term(out, is_first, nonzero.back().second, nonzero.back().first, var);
		append_term(out,
					is_first,
					nonzero[nonzero.size() - 2].second,
					nonzero[nonzero.size() - 2].first,
					var);
		out += " + ...";
		append_term(out, is_first, nonzero.front().second, nonzero.front().first, var);
	}

	return fmt::format_to(ctx.out(), "{}", out);
}

auto fmt::formatter<SplineFunction>::format(const SplineFunction &spline,
											fmt::format_context &ctx) const -> decltype(ctx.out()) {
	const SizeType interval_count =
		spline.nodes.size() < 2 ? 0
								: std::min(spline.nodes.size() - 1, spline.piecewise_poly.size());
	if (interval_count == 0)
		return ctx.out();

	std::vector<std::string> intervals;
	std::vector<std::string> labels;
	std::vector<std::string> polynomials_in_shift;
	std::vector<std::string> polynomials;
	intervals.reserve(interval_count);
	labels.reserve(interval_count);
	polynomials_in_shift.reserve(interval_count);
	polynomials.reserve(interval_count);

	for (SizeType i = 0; i < interval_count; ++i) {
		const Polynomial &poly = spline.piecewise_poly[i];
		const Real		 &node = spline.nodes[i];

		const auto format_node = [](Real node) {
			return fmt::format("{:.{}g}", node, config::kDisplayPrecision);
		};
		intervals.push_back(
			fmt::format("[{}, {}]", format_node(node), format_node(spline.nodes[i + 1])));
		labels.push_back(fmt::format("S_{}(x)", i));
		const Polynomial shifted_poly = utils::shift_polynomial(poly, node);

		std::string shift_variable;
		if (node == 0) {
			shift_variable = "x";
		} else {
			shift_variable =
				fmt::format("(x {} {})", node > 0 ? "-" : "+", format_node(std::abs(node)));
		}

		polynomials_in_shift.push_back(
			fmt::format(fmt::runtime(fmt::format("{{:{}}}", shift_variable)), poly));
		polynomials.push_back(fmt::format("{}", shifted_poly));
	}

	SizeType interval_width			  = 0;
	SizeType label_width			  = 0;
	SizeType shifted_polynomial_width = 0;
	for (SizeType i = 0; i < interval_count; ++i) {
		interval_width = std::max(interval_width, intervals[i].size());
		label_width	   = std::max(label_width, labels[i].size());
		shifted_polynomial_width =
			std::max(shifted_polynomial_width, polynomials_in_shift[i].size());
	}

	auto out = ctx.out();
	for (SizeType i = 0; i < interval_count; ++i) {
		if (i != 0)
			out = fmt::format_to(out, "\n");
		out = fmt::format_to(out,
							 "{:>{}} : {:<{}} = {:<{}} = {}",
							 intervals[i],
							 interval_width,
							 labels[i],
							 label_width,
							 polynomials_in_shift[i],
							 shifted_polynomial_width,
							 polynomials[i]);
	}
	return out;
}

/// @endcond
