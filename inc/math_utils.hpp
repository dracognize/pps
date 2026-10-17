#pragma once
/// @file math_utils.hpp
/// @brief Math utility functions.

#include <math_object.hpp>

namespace utils {
	/// @brief Factorial as `Real`.
	/// @param k Non-negative order.
	/// @return `k!`, `0! == 1`.
	auto factorial(SizeType k) -> Real;
	/// @brief Highest index with a non-negligible coefficient.
	/// @param poly Ascending coefficients.
	/// @return Effective degree, `0` for the zero polynomial.
	auto effective_degree(const Polynomial &poly) -> SizeType;
	/// @brief Compose a polynomial with itself `k` times.
	/// @details Computes `poly(x^k)`.
	/// @param poly The polynomial to compose.
	/// @param pow The number of times to compose.
	/// @return The composed polynomial.
	auto compose_power(const Polynomial &poly, SizeType pow) -> Polynomial;
} // namespace utils