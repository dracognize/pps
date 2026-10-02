#pragma once

#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>

using Real	   = double;
using SizeType = std::size_t;

// List of reals.
using RealList = std::vector<Real>;

// poly[i] is coeff of x^i.
struct Polynomial : std::vector<Real> {
		using std::vector<Real>::vector;
};

// poly = quotient * (x - c) + remainder.
struct PolynomialDivision {
		Polynomial quotient;
		Real	   remainder = 0;
};

// Disable fmt range formatting.
template <typename Char>
struct fmt::range_format_kind<Polynomial, Char>
	: std::integral_constant<fmt::range_format, fmt::range_format::disabled> {};

template <> struct fmt::formatter<Polynomial> {
		constexpr auto parse(fmt::format_parse_context &ctx) {
			return ctx.begin();
		}

		auto format(const Polynomial &poly, fmt::format_context &ctx) const -> decltype(ctx.out());
};
