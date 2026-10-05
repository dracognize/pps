#pragma once
/// @file config.hpp
/// @brief Display tolerances shared by printers and tables.

#include <math_object.hpp>

namespace config {

	/// @brief Magnitude below which a value counts as zero.
	constexpr double kNearZeroThreshold = 1e-12;
	/// @brief Significant digits used by every display format.
	constexpr int kDisplayPrecision = 7;
	/// @brief Term count above which the printer collapses to `lead + ...`.
	constexpr SizeType kMaxFullTerms = 99;
	/// @brief Minimum FTXUI column width for tableau tables.
	constexpr int kTableMinWidth = 12;

} // namespace config
