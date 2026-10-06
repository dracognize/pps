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

	fmt::print("Chebyshev Nodes:\r\n");
	display_chebyshev_nodes(-1.0, 1.0, 5);

	fmt::print("\r\nPolynomial Evaluation\r\n");
	display_polynomial_evaluation(c, poly);

	fmt::print("\r\nHorner Division\r\n");
	const Polynomial b = {-c, 1}; // (x - c)
	display_horner_division(poly, b);

	fmt::print("\r\nk-th Derivative at x = c (generalized Horner tableau)\r\n");
	for (SizeType k = 0; k <= 4; ++k) {
		fmt::print("\r\nk = {}\r\n", k);
		display_kth_derivative(c, poly, k);
	}
}
