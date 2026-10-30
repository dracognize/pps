#pragma once
/// @file math_object.hpp
/// @brief Core math types: scalars, polynomials, division results.

#include <string_view>
#include <string>
#include <algorithm>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <Eigen/Core>

/// @brief Floating-point scalar used for every numerical value.
using Real = double;
/// @brief Index and count type.
using SizeType = std::size_t;

/// @brief List of reals (nodes, samples).
using RealList = std::vector<Real>;
/// @brief 2D table of reals (list of rows).
using RealTable = std::vector<RealList>;

/// @brief Polynomial in ascending coefficient order.
/// @details `poly[i]` is the coefficient of `x^i`. Empty means zero.
/// @details operators `+` and `*` are defined for polynomial addition, scalar multiplication,
/// polynomial multiplication.
struct Polynomial : std::vector<Real> {
		using std::vector<Real>::vector;

		auto		operator+(const Polynomial &other) const -> Polynomial;
		auto		operator*(const Real &scalar) const -> Polynomial;
		friend auto operator*(const Real &scalar, const Polynomial &poly) -> Polynomial;
		auto		operator*(const Polynomial &other) const -> Polynomial;
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
		RealList				remainders; ///< Evaluations `Q_j(c)`.
};

/// @brief List of polynomials (e.g. node-polynomial history `w_0..w_n`).
using PolynomialList = std::vector<Polynomial>;

/// @brief Inputs and intermediates for Lagrange interpolation.
/// @details For nodes `x_i` with values `y_i`: `basis[i]` is `w(x) / (x - x_i)`
/// and `normalizers[i]` is `prod_{j != i} (x_i - x_j)`, so
/// `P(x) = sum_i values[i] / normalizers[i] * basis[i]`.
struct LagrangeInterpolationHistory {
		PolynomialList			node_polys;	 ///< Node polynomial history `w_0..w_n`.
		std::vector<Polynomial> basis;		 ///< Reduced polynomials `w(x) / (x - x_i)`.
		RealList				normalizers; ///< Barycentric denominators.
		RealList				values;		 ///< Function values `y_i` at the nodes.
};

/// @brief Inputs and intermediates for Newton forward interpolation.
/// @details For nodes `x_i` with values `y_i`: `divided_diffs[j][i]` is the
/// `i`-th order divided difference starting at node `j` (so
/// `divided_diffs[j][0] = y_j`), and `node_polys[k]` is
/// `w_k(x) = prod(x - x_j)` with `w_0 = 1`
struct NewtonInterpolationHistory {
		RealTable	   divided_diffs; ///< Divided-difference table `divided_diffs[j][i]`.
		PolynomialList node_polys;	  ///< Node polynomial history `w_0..w_n`.
};

/// @brief Inputs and intermediates for Newton forward interpolation with equidistant nodes.
/// @details For nodes `x_i` with values `y_i`: `finite_diffs[j][i]` is the
/// `i`-th order finite difference starting at node `j` (so
/// `finite_diffs[j][0] = y_j`), and `node_polys[k]` is
/// `w_k(t) = prod(t + i)` with `w_0 = 1` and `t = (x - x_0) / h`.
struct EquidistantHistory {
		RealTable	   finite_diffs;
		PolynomialList node_polys;
};

/// @brief Inclusive range of selected indices.
struct ExtractInterval {
		SizeType start;
		SizeType end;
};

/// @brief List of inclusive index ranges.
using IntervalsList = std::vector<ExtractInterval>;

/// @brief One inverse-interpolation calculation on an isolating interval.
struct InverseInterpolationStep {
		ExtractInterval isolating_interval;
		ExtractInterval selected_interval;
		RealList		inverse_nodes;
		RealList		inverse_values;
		Polynomial		inverse_polynomial;
		Real			root = 0;
};

/// @brief Full history of inverse interpolation for a target value.
struct InverseInterpolationHistory {
		IntervalsList						  isolating_intervals;
		std::vector<InverseInterpolationStep> steps;
		RealList							  roots;
};

struct SplineFunction {
		RealList	   nodes;
		PolynomialList piecewise_poly;

		auto value(const Real &key) -> Real;
};

/// @cond INTERNAL
template <typename Char>
struct fmt::range_format_kind<Polynomial, Char>
	: std::integral_constant<fmt::range_format, fmt::range_format::disabled> {};

/// @brief Formatter rendering `5*x^3 + 4*x^2 + ...` (`0` when empty).
/// @details The indeterminate is `x` by default. A custom variable name may be
/// passed as the format spec, e.g. `fmt::format("{:t}", poly)` renders in `t`.
/// Blank specs (e.g. `{}`) fall back to `x`.
template <> struct fmt::formatter<Polynomial> {
		constexpr auto parse(fmt::format_parse_context &ctx) -> const char * {
			auto it	 = ctx.begin();
			auto end = ctx.end();
			if (it == end || *it == '}')
				return it; // Keep the default variable name.

			auto begin = it;
			while (it != end && *it != '}')
				++it;
			var_ = std::string_view(begin, static_cast<std::size_t>(it - begin));
			return it;
		}

		auto format(const Polynomial &poly, fmt::format_context &ctx) const -> decltype(ctx.out());

	private:
		/// Indeterminate used when rendering; `x` unless overridden by the spec.
		std::string_view var_{"x"};
};
/// @brief Formatter for piecewise spline polynomials, one interval per line.
template <> struct fmt::formatter<SplineFunction> {
		constexpr auto parse(fmt::format_parse_context &ctx) -> const char * {
			auto it = ctx.begin();
			if (it != ctx.end() && *it != '}')
				throw fmt::format_error("SplineFunction formatter has no format specifiers");
			return it;
		}

		auto format(const SplineFunction &spline, fmt::format_context &ctx) const
			-> decltype(ctx.out());
};

/// @brief Formatter for Eigen matrices using a square-bracket matrix layout.
template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct fmt::formatter<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>> {
		constexpr auto parse(fmt::format_parse_context &ctx) -> const char * {
			auto it = ctx.begin();
			if (it != ctx.end() && *it != '}')
				throw fmt::format_error("Eigen matrix formatter has no format specifiers");
			return it;
		}

		auto format(const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> &matrix,
				fmt::format_context &ctx) const -> decltype(ctx.out()) {
			if (matrix.size() == 0)
				return fmt::format_to(ctx.out(), "[]");

			std::vector<std::vector<std::string>> values(
				static_cast<SizeType>(matrix.rows()),
				std::vector<std::string>(static_cast<SizeType>(matrix.cols())));
			std::vector<SizeType> widths(static_cast<SizeType>(matrix.cols()), 0);
			for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
				for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
					auto &value = values[static_cast<SizeType>(row)][static_cast<SizeType>(col)];
					value = fmt::format("{}", matrix(row, col));
					widths[static_cast<SizeType>(col)] =
						std::max(widths[static_cast<SizeType>(col)], value.size());
				}
			}

			auto out = ctx.out();
			for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
				if (row != 0)
					out = fmt::format_to(out, "\n");
				out = fmt::format_to(out, "[");
				for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
					if (col != 0)
						out = fmt::format_to(out, " ");
					out = fmt::format_to(out,
										 "{:>{}}",
									 values[static_cast<SizeType>(row)][static_cast<SizeType>(col)],
									 widths[static_cast<SizeType>(col)]);
				}
				out = fmt::format_to(out, "]");
			}
			return out;
		}
};
/// @endcond
