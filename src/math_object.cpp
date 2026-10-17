#include <config.hpp>

#include <math_object.hpp>

#include <cctype>
#include <cmath>
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
		return is_unit_coeff(abs_c)
				   ? fmt::format("{}^{}", var, power)
				   : fmt::format("{}*{}^{}", format_abs_coeff(abs_c), var, power);
	}

	auto append_term(
		std::string &out, bool &is_first, Real coeff, SizeType power, std::string_view var) -> void {
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
/// @endcond
