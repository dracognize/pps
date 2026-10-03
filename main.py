from pps import Polynomial, Real, display_chebyshev_nodes, display_polynomial_evaluation

poly_a: Polynomial = Polynomial([0, 0, 0, 0, 0, 0, 1])
poly_b: Polynomial = Polynomial([3, 2, 0, 4, 0, 0, 1])

x: Real = 10

poly: Polynomial = poly_a + poly_b

display_polynomial_evaluation(x, poly)
display_chebyshev_nodes(-1, 1, 11)
