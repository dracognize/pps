#pragma once

#include <math_object.hpp>

namespace config {

	constexpr double kNearZeroThreshold = 1e-12; // Radius of zero neighborhood that reduced to zero
	constexpr int	 kDisplayPrecision	= 7;	 // Precision for number display

	constexpr SizeType kMaxFullTerms = 99; // How many coefficients that still printing full terms

} // namespace config
