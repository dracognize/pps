#include <io/input_handling.hpp>
#include <math_object.hpp>

#include <chapter1/interpolation/horner.hpp>
#include <chapter1/interpolation/interval.hpp>
#include <chapter1/interpolation/inverse.hpp>
#include <chapter1/interpolation/optimal.hpp>
#include <chapter1/interpolation/polynomial/Bessel.hpp>
#include <chapter1/interpolation/polynomial/Gauss.hpp>
#include <chapter1/interpolation/polynomial/Lagrange.hpp>
#include <chapter1/interpolation/polynomial/Newton.hpp>
#include <chapter1/interpolation/polynomial/Stirling.hpp>
#include <chapter1/lsp.hpp>
#include <chapter1/spline/cubic.hpp>
#include <chapter1/spline/linear.hpp>
#include <chapter1/spline/quadratic.hpp>

#include <fmt/base.h>

// C++ entry point: exercises every feature implemented.
//   - `detail::` calls print the pure compute result for cross-checking.
//   - `display_*` calls render the full plain-text/FTXUI-table output.
// The interpolation samples come from f(x) = x^4 - 2x + 1, so the
// degree-4 interpolants should reproduce f exactly.
auto main() -> int {
	/////////////////////////////////////////////////////////////////////////
	/*
	 * Chapter1: Approximate Function
	 */
	/////////////////////////////////////////////////////////////////////////

	{
		/// @brief Find Chebyshev nodes form
		{
			Real	 a;
			Real	 b;
			SizeType n;

			display_chebyshev_nodes(a, b, n);
		}

		/// @brief Find polynomial value
		{
			Real	   c;
			Polynomial P = {};

			display_polynomial_evaluation(c, P);
		}

		/// @brief Divide polynomial by linear polynomial
		{
			Polynomial P = {};
			Real	   c;

			display_horner_division(P, {-c, 1});
		}

		/// @brief Find k-th derivative value
		{
			Real	   c;
			Polynomial P = {};
			SizeType   k;

			display_kth_derivative(c, P, k);
		}

		/// @brief Multiply polynomial by linear polynomial
		{
			Polynomial P = {};
			Real	   c;

			display_horner_multiplication(P, {-c, 1});
		}

		/// @brief Expand node polynomial omega_{n+1}(x)
		{
			RealList x = {};
			ih::nodes_input(x);

			display_nodes_polynomial_table(x);
		}

		/// @brief Extract equidistant interpolation points
		{
			RealList x = {};
			RealList y = {};
			Real	 x_target;
			SizeType k;

			ih::data_point_input(x, y);

			display_equidistant_interpolation_points(x, y, x_target, k);
		}

		/// @brief Lagrange interpolation
		{
			RealList x = {};
			RealList y = {};
			ih::data_point_input(x, y);

			display_Lagrange_interpolation(x, y);
		}

		/// @brief Newton general forward interpolation
		{
			RealList x = {};
			RealList y = {};
			ih::data_point_input(x, y);

			display_Newton_general_forward(x, y);
		}

		/// @brief Newton equidistant forward interpolation
		{
			RealList x = {};
			RealList y = {};
			ih::data_point_input(x, y);

			display_Newton_equidistant_forward(x, y);
		}

		/// @brief Newton general backward interpolation
		{
			RealList x = {};
			RealList y = {};
			ih::data_point_input(x, y);

			display_Newton_general_backward(x, y);
		}

		/// @brief Newton equidistant backward interpolation
		{
			RealList x = {};
			RealList y = {};
			ih::data_point_input(x, y);

			display_Newton_equidistant_backward(x, y);
		}

		/// @brief Gauss I interpolation
		{
			RealList x = {};
			RealList y = {};
			ih::data_point_input(x, y);

			display_GaussI_interpolation(x, y);
		}

		/// @brief Gauss II interpolation
		{
			RealList x = {};
			RealList y = {};
			ih::data_point_input(x, y);

			display_GaussII_interpolation(x, y);
		}

		/// @brief Stirling interpolation
		{
			RealList x = {};
			RealList y = {};
			ih::data_point_input(x, y);

			display_Stirling_interpolation(x, y);
		}

		/// @brief Bessel interpolation
		{
			RealList x = {};
			RealList y = {};
			ih::data_point_input(x, y);

			display_Bessel_interpolation(x, y);
		}

		/// @brief Find monotonic intervals
		{
			RealList x = {};
			RealList y = {};
			ih::data_point_input(x, y);

			display_monotonic_intervals(x, y);
		}

		/// @brief Find isolating intervals
		{
			RealList x = {};
			RealList y = {};
			Real	 y_target;

			ih::data_point_input(x, y);

			display_isolating_intervals(x, y, y_target);
		}

		/// @brief Inverse interpolation method
		{
			RealList x = {};
			RealList y = {};
			Real	 y_target;
			SizeType m;

			ih::data_point_input(x, y);

			display_inverse_interpolation(x, y, y_target, m);
		}

		/// @brief Iterative inverse interpolation method
		{
			RealList x = {};
			RealList y = {};
			Real	 y_target;
			Real	 epsilon;
			SizeType m;

			ih::data_point_input(x, y);

			display_inverse_iteration(x, y, y_target, epsilon, m);
		}

		/// @brief Linear spline
		{
			RealList x = {};
			RealList y = {};
			ih::data_point_input(x, y);

			display_linear_spline(x, y);
		}

		/// @brief Quadratic spline
		{
			/// Boundary condition:
			/// weight[0]*S'_0(x0)  + weight[1]*S'_{n-1}(xn) = bias

			RealList x		= {};
			RealList y		= {};
			RealList weight = {0, 0};
			Real	 bias;

			ih::data_point_input(x, y);

			display_quadratic_spline(x, y, weight[0], weight[1], bias);
		}

		/// @brief Cubic spline
		{
			/// @brief Boundary condition 1:
			/// weight[0][0]*S''_0(x0) + weight[0][1]*S''_{n-1}(xn)
			///    + weight[0][2]*S'_0(x0)  + weight[0][3]*S'_{n-1}(xn) = bias[0]
			/// @brief Boundary condition 2:
			/// weight[1][0]*S''_0(x0) + weight[1][1]*S''_{n-1}(xn)
			///	   + weight[1][2]*S'_0(x0) + weight[1][3]*S'_{n-1}(xn) = bias[1]

			RealList  x		 = {};
			RealList  y		 = {};
			RealTable weight = {{0, 0, 0, 0}, {0, 0, 0, 0}};
			RealList  bias	 = {0, 0};

			ih::data_point_input(x, y);
			Eigen::Matrix<double, 2, 5> boundary_conditions;
			boundary_conditions << 
				weight[0][0], weight[0][1], weight[0][2], weight[0][3], bias[0],
				weight[1][0], weight[1][1], weight[1][2], weight[1][3], bias[1];

			display_cubic_spline(x, y, boundary_conditions);
		}

		/// @brief Linear least squares approximation
		{
			RealList x = {};
			RealList y = {};
			ih::data_point_input(x, y);

			display_linear_lsp(x, y);
		}

		/// @brief Quadratic least squares approximation
		{
			RealList x = {};
			RealList y = {};
			ih::data_point_input(x, y);

			display_quadratic_lsp(x, y);
		}

		/// @brief General least squares approximation
		{
			using namespace std;

			RealList	  x = {};
			RealList	  y = {};
			FunctionsList phi;

			ih::data_point_input(x, y);
			phi = {[](Real x) { return 1; },
				   [](Real x) { return x; },
				   [](Real x) { return pow(x, 2); },
				   [](Real x) { return exp(x); },
				   [](Real x) { return log(x) / log(2); }};

			display_general_lsp(x, y, phi);
		}
	}
}
