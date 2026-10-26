#include <chapter1/interpolation/horner.hpp>
#include <chapter1/interpolation/optimal.hpp>
#include <chapter1/interpolation/polynomial/Bessel.hpp>
#include <chapter1/interpolation/polynomial/Gauss.hpp>
#include <chapter1/interpolation/polynomial/Lagrange.hpp>
#include <chapter1/interpolation/polynomial/Newton.hpp>
#include <chapter1/interpolation/polynomial/Stirling.hpp>
#include <io/input_handling.hpp>
#include <math_object.hpp>

#include <fmt/base.h>

// C++ entry point: exercises every feature implemented in chapter 1.
//   - `detail::` calls print the pure compute result for cross-checking.
//   - `display_*` calls render the full plain-text/FTXUI-table output.
// The interpolation samples come from f(x) = x^4 - 2x + 1, so the
// degree-4 interpolants should reproduce f exactly.
auto main() -> int {
	// P(x) = 1 - 2x + x^4, evaluated at c = 2.
	const Polynomial poly = {1, -2, 0, 0, 1};
	const Real		 c	  = 2;

	// ------------------------------------------------------------------
	// Horner: evaluation, division, derivatives, multiplication, nodes.
	// ------------------------------------------------------------------
	fmt::print("=== Horner: polynomial evaluation ===\r\n");
	display_polynomial_evaluation(c, poly);
	fmt::print("pure: P({}) = {}\r\n", c, detail::evaluate_polynomial(c, poly));

	fmt::print("\r\n=== Horner: division by (x - c) ===\r\n");
	const Polynomial linear = {-c, 1};
	display_horner_division(poly, linear);
	const auto division = detail::divide_by_linear(poly, linear);
	fmt::print("pure: Q(x) = {}, R = {}\r\n", division.quotient, division.remainder);

	fmt::print("\r\n=== Horner: k-th derivative at x = c ===\r\n");
	for (SizeType k = 0; k <= 5; ++k) {
		fmt::print("\r\nk = {}\r\n", k);
		display_kth_derivative(c, poly, k);
		fmt::print("pure: P^({})({}) = {}\r\n", k, c, detail::kth_derivative_value(c, poly, k));
	}

	const auto derivative_history = detail::build_derivative_history(c, poly, 4);
	fmt::print("\r\nderivative history: {} quotients, {} remainders\r\n",
			   derivative_history.quotients.size(),
			   derivative_history.remainders.size());

	fmt::print("\r\n=== Horner: multiplication by (x - c) ===\r\n");
	display_horner_multiplication(poly, linear);
	fmt::print("pure: product = {}\r\n", detail::multiply_by_linear(poly, linear));

	const RealList nodes = {0, 1, 2, 3, 4};

	fmt::print("\r\n=== Horner: node polynomials ===\r\n");
	display_nodes_polynomial_table(nodes);
	fmt::print("pure: node polynomial = {}\r\n", detail::node_polynomial(nodes));

	const auto node_history = detail::build_node_polynomial_history(nodes);
	fmt::print("pure: history size = {}, final node polynomial = {}\r\n",
			   node_history.size(),
			   node_history.back());

	fmt::print("\r\n=== Horner: node polynomials with variable t and pow 2 ===\r\n");
	const RealList squared_nodes = {1, 4, 9};
	display_nodes_polynomial_table(squared_nodes, "t", 2);

	// ------------------------------------------------------------------
	// Polynomial arithmetic.
	// ------------------------------------------------------------------
	fmt::print("\r\n=== Polynomial arithmetic ===\r\n");
	const Polynomial a = {1, 2, 3};
	const Polynomial b = {4, 5};
	fmt::print("({}) + ({}) = {}\r\n", a, b, a + b);
	fmt::print("({}) * ({}) = {}\r\n", a, b, a * b);
	fmt::print("({}) * 2.5 = {}\r\n", a, a * 2.5);
	fmt::print("2.5 * ({}) = {}\r\n", a, 2.5 * a);

	// ------------------------------------------------------------------
	// Chebyshev nodes.
	// ------------------------------------------------------------------
	fmt::print("\r\n=== Chebyshev nodes on [-1, 1] ===\r\n");
	display_chebyshev_nodes(-1.0, 1.0, 5);
	const auto chebyshev_nodes = detail::generate_chebyshev_nodes(-1.0, 1.0, 5);
	fmt::print("pure: {}\r\n", chebyshev_nodes);

	// ------------------------------------------------------------------
	// Interpolation samples: f(x) = x^4 - 2x + 1 on x = 0..4.
	// ------------------------------------------------------------------
	const RealList values = {1, 0, 13, 76, 249};

	fmt::print("\r\n=== Lagrange interpolation ===\r\n");
	display_Lagrange_interpolation(nodes, values);
	const auto lagrange_history = detail::build_Lagrange_interpolation_history(nodes, values);
	const auto lagrange_poly	= detail::Lagrange_polynomial(nodes, values);
	fmt::print("pure: L(x) = {}\r\n", lagrange_poly);
	fmt::print("pure: L(2.5) = {}, P(2.5) = {}\r\n",
			   detail::evaluate_polynomial(2.5, lagrange_poly),
			   detail::evaluate_polynomial(2.5, poly));
	fmt::print("pure: {} basis polynomials, {} normalizers\r\n",
			   lagrange_history.basis.size(),
			   lagrange_history.normalizers.size());

	// Non-equidistant nodes for the general Newton formulas.
	const RealList general_nodes  = {0, 1, 2, 4, 6};
	const RealList general_values = {1, 0, 13, 249, 1285};

	fmt::print("\r\n=== Newton general forward interpolation ===\r\n");
	display_Newton_general_forward(general_nodes, general_values);
	const auto newton_general_forward =
		detail::build_Newton_general_forward_history(general_nodes, general_values);
	const auto newton_general_forward_poly =
		detail::Newton_general_forward_polynomial(general_nodes, general_values);
	fmt::print("pure: P(x) = {}\r\n", newton_general_forward_poly);
	fmt::print("pure: P(2.5) = {}, f(2.5) = {}\r\n",
			   detail::evaluate_polynomial(2.5, newton_general_forward_poly),
			   detail::evaluate_polynomial(2.5, poly));
	fmt::print("pure: {} divided-difference rows, {} node polynomials\r\n",
			   newton_general_forward.divided_diffs.size(),
			   newton_general_forward.node_polys.size());

	fmt::print("\r\n=== Newton general backward interpolation ===\r\n");
	display_Newton_general_backward(general_nodes, general_values);
	const auto newton_general_backward =
		detail::build_Newton_general_backward_history(general_nodes, general_values);
	const auto newton_general_backward_poly =
		detail::Newton_general_backward_polynomial(general_nodes, general_values);
	fmt::print("pure: P(x) = {}\r\n", newton_general_backward_poly);
	fmt::print("pure: P(2.5) = {}, f(2.5) = {}\r\n",
			   detail::evaluate_polynomial(2.5, newton_general_backward_poly),
			   detail::evaluate_polynomial(2.5, poly));
	fmt::print("pure: {} divided-difference rows, {} node polynomials\r\n",
			   newton_general_backward.divided_diffs.size(),
			   newton_general_backward.node_polys.size());

	fmt::print("\r\n=== Newton equidistant forward interpolation ===\r\n");
	display_Newton_equidistant_forward(nodes, values);
	const auto newton_equidistant_forward =
		detail::build_Newton_equidistant_forward_history(nodes, values);
	const auto newton_equidistant_forward_poly =
		detail::Newton_equidistant_forward_polynomial(nodes, values);
	fmt::print("pure: Q(t) = {}\r\n", fmt::format("{:t}", newton_equidistant_forward_poly));
	fmt::print("pure: Q(t = 2.5) = {}, f(2.5) = {}\r\n",
			   detail::evaluate_polynomial(2.5, newton_equidistant_forward_poly),
			   detail::evaluate_polynomial(2.5, poly));
	fmt::print("pure: {} finite-difference rows, {} node polynomials\r\n",
			   newton_equidistant_forward.finite_diffs.size(),
			   newton_equidistant_forward.node_polys.size());

	fmt::print("\r\n=== Newton equidistant backward interpolation ===\r\n");
	display_Newton_equidistant_backward(nodes, values);
	const auto newton_equidistant_backward =
		detail::build_Newton_equidistant_backward_history(nodes, values);
	const auto newton_equidistant_backward_poly =
		detail::Newton_equidistant_backward_polynomial(nodes, values);
	fmt::print("pure: Q(t) = {}\r\n", fmt::format("{:t}", newton_equidistant_backward_poly));
	// x = 2.5, x_n = 4, h = 1 => t = -1.5
	fmt::print("pure: Q(t = -1.5) = {}, f(2.5) = {}\r\n",
			   detail::evaluate_polynomial(-1.5, newton_equidistant_backward_poly),
			   detail::evaluate_polynomial(2.5, poly));
	fmt::print("pure: {} finite-difference rows, {} node polynomials\r\n",
			   newton_equidistant_backward.finite_diffs.size(),
			   newton_equidistant_backward.node_polys.size());

	fmt::print("\r\n=== Gauss I interpolation ===\r\n");
	display_GaussI_interpolation(nodes, values);
	const auto gauss_i_history = detail::build_GaussI_history(nodes, values);
	const auto gauss_i_poly	   = detail::GaussI_polynomial(nodes, values);
	fmt::print("pure: Q(t) = {}\r\n", fmt::format("{:t}", gauss_i_poly));
	// x = 2.5, x_n = 2, h = 1 => t = 0.5
	fmt::print("pure: Q(t = 0.5) = {}, f(2.5) = {}\r\n",

			   detail::evaluate_polynomial(0.5, gauss_i_poly),
			   detail::evaluate_polynomial(2.5, poly));
	fmt::print("pure: {} finite-difference rows, {} node polynomials\r\n",
			   gauss_i_history.finite_diffs.size(),
			   gauss_i_history.node_polys.size());

	fmt::print("\r\n=== Gauss II interpolation ===\r\n");
	display_GaussII_interpolation(nodes, values);
	const auto gauss_ii_history = detail::build_GaussII_history(nodes, values);
	const auto gauss_ii_poly	= detail::GaussII_polynomial(nodes, values);
	fmt::print("pure: Q(t) = {}\r\n", fmt::format("{:t}", gauss_ii_poly));
	// x = 2.5, x_n = 2, h = 1 => t = 0.5
	fmt::print("pure: Q(t = 0.5) = {}, f(2.5) = {}\r\n",
			   detail::evaluate_polynomial(0.5, gauss_ii_poly),
			   detail::evaluate_polynomial(2.5, poly));
	fmt::print("pure: {} finite-difference rows, {} node polynomials\r\n",
			   gauss_ii_history.finite_diffs.size(),
			   gauss_ii_history.node_polys.size());

	fmt::print("\r\n=== Stirling interpolation ===\r\n");
	display_Stirling_interpolation(nodes, values);
	const auto stirling_history = detail::build_Stirling_history(nodes, values);
	const auto stirling_poly	= detail::Stirling_polynomial(nodes, values);
	fmt::print("pure: Q(t) = {}\r\n", fmt::format("{:t}", stirling_poly));
	// x = 2.5, x_n = 2, h = 1 => t = 0.5
	fmt::print("pure: Q(t = 0.5) = {}, f(2.5) = {}\r\n",
			   detail::evaluate_polynomial(0.5, stirling_poly),
			   detail::evaluate_polynomial(2.5, poly));
	fmt::print("pure: {} finite-difference rows, {} node polynomials\r\n",
			   stirling_history.finite_diffs.size(),
			   stirling_history.node_polys.size());

	// Bessel needs an even number of equidistant nodes.
	const RealList bessel_nodes	 = {0, 1, 2, 3, 4, 5};
	const RealList bessel_values = {1, 0, 13, 76, 249, 616};

	fmt::print("\r\n=== Bessel interpolation ===\r\n");
	display_Bessel_interpolation(bessel_nodes, bessel_values);
	const auto bessel_history = detail::build_Bessel_history(bessel_nodes, bessel_values);
	const auto bessel_poly	  = detail::Bessel_polynomial(bessel_nodes, bessel_values);
	fmt::print("pure: Q(u) = {}\r\n", fmt::format("{:u}", bessel_poly));
	// x = 2.5, x_(n-1) = 2, h = 1 => u = 0
	fmt::print("pure: Q(u = 0) = {}, f(2.5) = {}\r\n",
			   detail::evaluate_polynomial(0.0, bessel_poly),
			   detail::evaluate_polynomial(2.5, poly));
	fmt::print("pure: {} finite-difference rows, {} node polynomials\r\n",
			   bessel_history.finite_diffs.size(),
			   bessel_history.node_polys.size());
}
