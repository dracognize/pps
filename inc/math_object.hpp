#pragma once
/// @file math_object.hpp
/// @brief Core math types: scalars, polynomials, division results.

#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>

/// @brief Floating-point scalar used for every numerical value.
using Real = double;
/// @brief Index and count type.
using SizeType = std::size_t;

/// @brief List of reals (nodes, samples).
using RealList = std::vector<Real>;

/// @brief Polynomial in ascending coefficient order.
/// @details `poly[i]` is the coefficient of `x^i`. Empty means zero.
struct Polynomial : std::vector<Real> {
		using std::vector<Real>::vector;
};

/// @brief Result of dividing by a degree-1 divisor.
struct PolynomialDivision {
		Polynomial quotient;	  ///< Quotient (empty when the dividend was constant).
		Real	   remainder = 0; ///< Scalar remainder.
};

/// @brief Generalized Horner tableau for derivatives at `x = c`.
/// @details Column `j` holds `Q_j`, the quotient after `j` divisions by
/// `(x - c)`. `remainders[j]` is `Q_j(c)`, i.e. `P^{(j)}(c) / j!`.
struct DerivativeHistory {
		std::vector<Polynomial> quotients;	///< Successive quotients `Q_0..Q_k`.
		std::vector<Real>		remainders; ///< Evaluations `Q_j(c)`.
};

/// @cond INTERNAL
template <typename Char>
struct fmt::range_format_kind<Polynomial, Char>
	: std::integral_constant<fmt::range_format, fmt::range_format::disabled> {};

/// @brief Formatter rendering `5*x^3 + 4*x^2 + ...` (`0` when empty).
template <> struct fmt::formatter<Polynomial> {
		constexpr auto parse(fmt::format_parse_context &ctx) {
			return ctx.begin();
		}

		auto format(const Polynomial &poly, fmt::format_context &ctx) const -> decltype(ctx.out());
};
/// @endcond
