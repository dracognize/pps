#include <config.hpp>
#include <math_object.hpp>

#include <cmath>
#include <utility>

namespace utils {
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

	auto compose_power(const Polynomial &poly, SizeType pow) -> Polynomial {
		if (poly.empty())
			return {};
		if (pow == 1)
			return poly;

		Polynomial result((poly.size() - 1) * pow + 1, 0.0);

		if (pow == 0) {
			for (auto &coeff : poly)
				result[0] += coeff;
			return result;
		}

		for (SizeType i = 0; i < poly.size(); ++i)
			result[i * pow] = poly[i];

		return result;
	}

	auto shift_polynomial(const Polynomial &poly, Real shift) -> Polynomial {
		if (poly.empty())
			return {};

		// Horner composition with the linear polynomial (x - xk).
		Polynomial result;
		for (auto it = poly.rbegin(); it != poly.rend(); ++it) {
			Polynomial next(result.size() + 1, 0.0);
			for (SizeType i = 0; i < result.size(); ++i) {
			next[i] -= shift * result[i];
				next[i + 1] += result[i];
			}
			next[0] += *it;
			result = std::move(next);
		}
		return result;
	}
} // namespace utils
