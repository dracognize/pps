#include <config.hpp>
#include <math_object.hpp>

#include <cmath>

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
} // namespace utils