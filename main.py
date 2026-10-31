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
        x_bar: Real = 0.0
        k: SizeType = 3
        ih.data_point_input(x, y)

        display_equidistant_interpolation_points(x, y, x_bar, k)

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
        y_bar: Real = 0.0
        ih.data_point_input(x, y)

        display_isolating_intervals(x, y, y_bar)

    # Inverse interpolation method
    if True:
        x: RealList = []
        y: RealList = []
        y_bar: Real = 0.0
        m: SizeType = 4
        ih.data_point_input(x, y)

        display_inverse_interpolation(x, y, y_bar, m)

    # Iterative inverse interpolation method
    if True:
        x: RealList = []
        y: RealList = []
        y_bar: Real = 0.0
        epsilon: Real = 1e-4
        m: SizeType = 4
        ih.data_point_input(x, y)

        display_inverse_iteration(x, y, y_bar, epsilon, m)

    # Linear spline
    if True:
        x: RealList = []
        y: RealList = []
        ih.data_point_input(x, y)

        display_linear_spline(x, y)

    # Quadratic spline
    if True:
        x: RealList = []
        y: RealList = []
        weight: RealList = [0,0]
        bias: Real = 1.0
        ih.data_point_input(x, y)

        display_quadratic_spline(x, y, weight[0], weight[1], bias)

    # Cubic spline
    if True:
        x: RealList = []
        y: RealList = []
        # Boundary condition 1:
        # a0_0*S''_0(x0) + a0_1*S''_{n-1}(xn) + a0_2*S'_0(x0) + a0_3*S'_{n-1}(xn) = gamma0
        # Boundary condition 2:
        # a1_0*S''_0(x0) + a1_1*S''_{n-1}(xn) + a1_2*S'_0(x0) + a1_3*S'_{n-1}(xn) = gamma1

        a0_0: Real = 1.0
        a0_1: Real = 0.0
        a0_2: Real = 0.0
        a0_3: Real = 0.0
        gamma0: Real = 0.0

        a1_0: Real = 0.0
        a1_1: Real = 1.0
        a1_2: Real = 0.0
        a1_3: Real = 0.0
        gamma1: Real = 0.0

        ih.data_point_input(x, y)

        boundary_conditions: list[list[float]] = [
            [a0_0, a0_1, a0_2, a0_3, gamma0],
            [a1_0, a1_1, a1_2, a1_3, gamma1],
        ]

        display_cubic_spline(x, y, boundary_conditions)

    # Linear least squares approximation
    if True:
        x: RealList = []
        y: RealList = []
        ih.data_point_input(x, y)

        display_linear_lsp(x, y)

    # Quadratic least squares approximation
    if True:
        x: RealList = []
        y: RealList = []
        ih.data_point_input(x, y)

        display_quadratic_lsp(x, y)

    # General least squares approximation
    if True:
        x: RealList = []
        y: RealList = []

        ih.data_point_input(x, y)
        phi = [
            lambda x: 1.0,
            lambda x: x,
            lambda x: x**2,
            lambda x: math.exp(x),
            lambda x: math.log(x) / math.log(2) if x > 0 else 0.0,
        ]

        display_general_lsp(x, y, phi)

    return 0


if __name__ == "__main__":
    main()
