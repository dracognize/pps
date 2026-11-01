"""Python entry point: exercises every feature implemented, mirroring main.cpp."""

import math
from pps import *


class ih:
    """Input helpers mirroring C++ namespace ih."""

    @staticmethod
    def nodes_input(nodes: RealList) -> None:
        nodes.clear()
        with open("datas/data_points/nodes.txt", "r", encoding="utf-8") as f:
            nodes.extend(float(token) for token in f.read().split())

    @staticmethod
    def data_point_input(nodes: RealList, values: RealList) -> None:
        nodes.clear()
        values.clear()
        with open("datas/data_points/nodes.txt", "r", encoding="utf-8") as f:
            nodes.extend(float(token) for token in f.read().split())
        with open("datas/data_points/values.txt", "r", encoding="utf-8") as f:
            values.extend(float(token) for token in f.read().split())
        if len(nodes) != len(values):
            raise ValueError(
                "The nodes and values files must have the same number of entries"
            )


def main() -> int:
    #########################################################################
    # Chapter1: Approximate Function
    #########################################################################

    # Find Chebyshev nodes form
    if True:
        a: Real = 0.0
        b: Real = 1.0
        n: SizeType = 5

        display_chebyshev_nodes(a, b, n)

    # Find polynomial value
    if True:
        c: Real = 0.0
        P: Polynomial = Polynomial([])

        display_polynomial_evaluation(c, P)

    # Divide polynomial by linear polynomial
    if True:
        P: Polynomial = Polynomial([])
        c: Real = 0.0

        display_horner_division(P, [-c, 1.0])

    # Find k-th derivative value
    if True:
        c: Real = 0.0
        P: Polynomial = Polynomial([])
        k: SizeType = 1

        display_kth_derivative(c, P, k)

    # Multiply polynomial by linear polynomial
    if True:
        P: Polynomial = Polynomial([])
        c: Real = 0.0

        display_horner_multiplication(P, [-c, 1.0])

    # Expand node polynomial omega_{n+1}(x)
    if True:
        x: RealList = []
        ih.nodes_input(x)

        display_nodes_polynomial_table(x)

    # Extract equidistant interpolation points
    if True:
        x: RealList = []
        y: RealList = []
        x_target: Real = 0.0
        k: SizeType = 3
        ih.data_point_input(x, y)

        display_equidistant_interpolation_points(x, y, x_target, k)

    # Lagrange interpolation
    if True:
        x: RealList = []
        y: RealList = []
        ih.data_point_input(x, y)

        display_Lagrange_interpolation(x, y)

    # Newton general forward interpolation
    if True:
        x: RealList = []
        y: RealList = []
        ih.data_point_input(x, y)

        display_Newton_general_forward(x, y)

    # Newton equidistant forward interpolation
    if True:
        x: RealList = []
        y: RealList = []
        ih.data_point_input(x, y)

        display_Newton_equidistant_forward(x, y)

    # Newton general backward interpolation
    if True:
        x: RealList = []
        y: RealList = []
        ih.data_point_input(x, y)

        display_Newton_general_backward(x, y)

    # Newton equidistant backward interpolation
    if True:
        x: RealList = []
        y: RealList = []
        ih.data_point_input(x, y)

        display_Newton_equidistant_backward(x, y)

    # Gauss I interpolation
    if True:
        x: RealList = []
        y: RealList = []
        ih.data_point_input(x, y)

        display_GaussI_interpolation(x, y)

    # Gauss II interpolation
    if True:
        x: RealList = []
        y: RealList = []
        ih.data_point_input(x, y)

        display_GaussII_interpolation(x, y)

    # Stirling interpolation
    if True:
        x: RealList = []
        y: RealList = []
        ih.data_point_input(x, y)

        display_Stirling_interpolation(x, y)

    # Bessel interpolation
    if True:
        x: RealList = []
        y: RealList = []
        ih.data_point_input(x, y)

        display_Bessel_interpolation(x, y)

    # Find monotonic intervals
    if True:
        x: RealList = []
        y: RealList = []
        ih.data_point_input(x, y)

        display_monotonic_intervals(x, y)

    # Find isolating intervals
    if True:
        x: RealList = []
        y: RealList = []
        y_target: Real = 0.0
        ih.data_point_input(x, y)

        display_isolating_intervals(x, y, y_target)

    # Inverse interpolation method
    if True:
        x: RealList = []
        y: RealList = []
        y_target: Real = 0.0
        m: SizeType = 4
        ih.data_point_input(x, y)

        display_inverse_interpolation(x, y, y_target, m)

    # Iterative inverse interpolation method
    if True:
        x: RealList = []
        y: RealList = []
        y_target: Real = 0.0
        epsilon: Real = 1e-4
        m: SizeType = 4
        ih.data_point_input(x, y)

        display_inverse_iteration(x, y, y_target, epsilon, m)

    # Linear spline
    if True:
        x: RealList = []
        y: RealList = []
        ih.data_point_input(x, y)

        display_linear_spline(x, y)

    # Quadratic spline
    if True:
        # Boundary condition:
        # weight[0]*S'_0(x0)  + weight[1]*S'_{n-1}(xn) = bias

        x: RealList = []
        y: RealList = []
        quadratic_weight: RealList = [0.0, 0.0]
        quadratic_bias: Real = 0.0
        ih.data_point_input(x, y)

        display_quadratic_spline(
            x, y, quadratic_weight[0], quadratic_weight[1], quadratic_bias
        )

    # Cubic spline
    if True:
        # Boundary condition 1:
        # weight[0][0]*S''_0(x0) + weight[0][1]*S''_{n-1}(xn)
        #    + weight[0][2]*S'_0(x0)  + weight[0][3]*S'_{n-1}(xn) = bias[0]
        # Boundary condition 2:
        # weight[1][0]*S''_0(x0) + weight[1][1]*S''_{n-1}(xn)
        #	   + weight[1][2]*S'_0(x0) + weight[1][3]*S'_{n-1}(xn) = bias[1]

        x: RealList = []
        y: RealList = []
        weight: list[list[float]] = [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]]
        bias: list[float] = [0.0, 0.0]

        ih.data_point_input(x, y)

        boundary_conditions: list[list[float]] = [
            [weight[0][0], weight[0][1], weight[0][2], weight[0][3], bias[0]],
            [weight[1][0], weight[1][1], weight[1][2], weight[1][3], bias[1]],
        ]

        display_cubic_spline(x, y, boundary_conditions)

    # Linear least squares approximation
    if True:
        x: RealList = []
        y: RealList = []
        error_type: LspErrorType = LspErrorType.RMSE
        ih.data_point_input(x, y)

        display_linear_lsp(x, y, error_type)

    # Quadratic least squares approximation
    if True:
        x: RealList = []
        y: RealList = []
        error_type: LspErrorType = LspErrorType.RMSE
        ih.data_point_input(x, y)

        display_quadratic_lsp(x, y, error_type)

    # General least squares approximation
    if True:
        x: RealList = []
        y: RealList = []
        error_type: LspErrorType = LspErrorType.RMSE

        ih.data_point_input(x, y)
        phi = [
            lambda x: 1.0,
            lambda x: x,
            lambda x: math.pow(x, 2),
            lambda x: math.exp(x),
            lambda x: math.log(x) / math.log(2) if x > 0 else 0.0,
        ]

        display_general_lsp(x, y, phi, error_type)

    return 0


if __name__ == "__main__":
    main()
