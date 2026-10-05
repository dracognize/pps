"""Python demo mirroring main.cpp (detail compute + display_*)."""

from pps import *

# P(x) = 2 + 3x + 4x^2 + 5x^3, evaluated at c = 2.
poly: Polynomial = Polynomial([2, 3, 4, 5])
c: Real = 2

print("Chebyshev Nodes:")
display_chebyshev_nodes(-1, 1, 5)
print("pure:", detail.generate_chebyshev_nodes(-1, 1, 5))

print("\nPolynomial Evaluation")
display_polynomial_evaluation(c, poly)
print("pure:", detail.evaluate_polynomial(c, poly))

print("\nHorner Division")
display_horner_division(poly, Polynomial([-c, 1]))

print("\nk-th Derivative at x = c (generalized Horner tableau)")
for k in range(5):
    print("\n", end="")
    display_kth_derivative(c, poly, k)
    print("pure:", detail.kth_derivative_value(c, poly, k))
