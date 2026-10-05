#include <chapter1/horner.hpp>
#include <chapter1/interpolation.hpp>
#include <math_object.hpp>

#include <fmt/base.h>

// C++ entry point: demos the API in dependency order.
//   - Pure compute (`detail::evaluate_polynomial`,
//     `detail::divide_by_linear`, `detail::generate_chebyshev_nodes`,
//     `detail::kth_derivative_value`) returns values, never prints.
//   - `display_*` computes via `detail::` then renders (plain text or an
//     FTXUI table through `ui::render_table`).
// Alternative: code in Python in main.py against pps.detail + display_*.
auto main() -> int {
	// P(x) = 2 + 3x + 4x^2 + 5x^3, evaluated at c = 2.
	const Polynomial poly = {2, 3, 4, 5};
	const Real		 c	  = 2;

	fmt::println("Chebyshev Nodes:");
	display_chebyshev_nodes(-1.0, 1.0, 5);

	fmt::println("\nPolynomial Evaluation");
	display_polynomial_evaluation(c, poly);

	fmt::println("\nHorner Division");
	const Polynomial b = {-c, 1}; // (x - c)
	display_horner_division(poly, b);

	fmt::println("\nk-th Derivative at x = c (generalized Horner tableau)");
	for (SizeType k = 0; k <= 4; ++k) {
		fmt::println("");
		display_kth_derivative(c, poly, k);
	}
}
