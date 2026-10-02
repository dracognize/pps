#include <config.hpp>

#include <math_object.hpp>

#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

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

	auto monomial_body(Real abs_c, SizeType power) -> std::string {
		if (power == 0)
			return format_abs_coeff(abs_c);
		if (power == 1)
			return is_unit_coeff(abs_c) ? "x" : format_abs_coeff(abs_c) + "*x";
		return is_unit_coeff(abs_c) ? fmt::format("x^{}", power)
									: fmt::format("{}*x^{}", format_abs_coeff(abs_c), power);
	}

	auto append_term(std::string &out, bool &is_first, Real coeff, SizeType power) -> void {
		if (is_zero_coeff(coeff))
			return;
		std::string body = monomial_body(std::abs(coeff), power);
		if (is_first) {
			out += (coeff < 0 ? "-" : "") + body;
			is_first = false;
		} else {
			out += (coeff < 0 ? " - " : " + ") + body;
		}
	}

} // namespace

auto fmt::formatter<Polynomial>::format(const Polynomial &poly, fmt::format_context &ctx) const
	-> decltype(ctx.out()) {
	std::vector<std::pair<SizeType, Real>> nonzero;
	for (SizeType i = 0; i < poly.size(); ++i) {
		if (!is_zero_coeff(poly[i]))
			nonzero.emplace_back(i, poly[i]);
	}

	std::string out;
	if (nonzero.empty()) {
		out = "0";
	} else if (nonzero.size() <= config::kMaxFullTerms) {
		bool is_first = true;
		for (auto it = nonzero.rbegin(); it != nonzero.rend(); ++it)
			append_term(out, is_first, it->second, it->first);
	} else {
		bool is_first = true;
		append_term(out, is_first, nonzero.back().second, nonzero.back().first);
		append_term(
			out, is_first, nonzero[nonzero.size() - 2].second, nonzero[nonzero.size() - 2].first);
		out += " + ...";
		append_term(out, is_first, nonzero.front().second, nonzero.front().first);
	}

	return fmt::format_to(ctx.out(), "{}", out);
}
