#include <chapter1/horner.hpp>
#include <chapter1/interpolation.hpp>
#include <math_object.hpp>

#include <fmt/base.h>

// C++ entry point: code against pps_lib here.
// Alternative: code in Python in main.py against the pps module.
auto main() -> int {
	Polynomial poly = {0, 0, 0, 0, 0, 0, 1};
	Real	   x	= 2;
	Polynomial b	= {-x, 1};

	fmt::println("Chebyshev Nodes:");
	display_chebyshev_nodes(-1.0, 1.0, 5);

	fmt::println("\nPolynomial Evaluation");
	display_polynomial_evaluation(x, poly);

	fmt::println("\nHorner Division");
	display_horner_division(poly, b);
}
