"""Type stubs for the compiled pps extension module.

Mirrors python/bindings.cpp. Keep in sync when the bindings change.
A .pyi file is never executed, so this does not affect runtime.

Conventions (same as C++):
  - Pure compute lives in ``pps.detail`` (mirrors ``detail::``).
  - Point first: ``detail.evaluate_polynomial(x, poly)``,
    ``detail.kth_derivative_value(c, poly, k)``.
  - Dividend first: ``detail.divide_by_linear(a, b)``.
  - ``poly`` accepts either a ``Polynomial`` or a plain ascending-coeff list.
  - Pure ``detail.*`` return values; ``display_*`` render to stdout
    and return ``None``.
"""

from collections.abc import Callable, Iterator
from enum import Enum
from typing import overload

Real = float
RealList = list
SizeType = int

class Polynomial:
    """Ascending-coeff polynomial: ``coeffs[i]`` is the ``x**i`` coefficient."""

    def __init__(self, coeffs: list[float] = ...) -> None: ...
    def __len__(self) -> int: ...
    def __getitem__(self, i: int) -> float: ...
    def __setitem__(self, i: int, v: float) -> None: ...
    def __iter__(self) -> Iterator[float]: ...
    def to_list(self) -> list[float]: ...
    def __add__(self, other: Polynomial) -> Polynomial: ...
    def __sub__(self, other: Polynomial) -> Polynomial: ...
    @overload
    def __mul__(self, other: Polynomial) -> Polynomial: ...
    @overload
    def __mul__(self, other: float) -> Polynomial: ...
    def __rmul__(self, other: float) -> Polynomial: ...
    def __truediv__(self, other: float) -> Polynomial: ...
    def __neg__(self) -> Polynomial: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class PolynomialDivision:
    """``a = b * quotient + remainder`` for degree-1 ``b``."""

    quotient: Polynomial
    remainder: float

class SplineFunction:
    """Piecewise polynomial spline; each segment stored in local powers of ``(x - nodes[i])``."""

    nodes: list[float]
    def value(self, x: float) -> float:
        """Evaluate the spline at ``x``."""
        ...
    def __repr__(self) -> str: ...

class InterpolationMethod(Enum):
    """Polynomial method used to build each inverse interpolant."""

    Lagrange = ...
    NewtonForward = ...
    NewtonBackward = ...

class _DetailModule:
    """Pure compute (mirrors C++ ``detail::``); no I/O."""

    # ---- Horner / Chebyshev ----
    def evaluate_polynomial(self, x: float, poly: Polynomial | list[float]) -> float:
        """Horner evaluation of ``P(x)``; empty poly yields ``0``."""
        ...
    def divide_by_linear(
        self,
        a: Polynomial | list[float],
        b: Polynomial | list[float],
    ) -> PolynomialDivision:
        """Divide by degree-1 ``b``; ``b`` must have exactly 2 coeffs."""
        ...
    def generate_chebyshev_nodes(
        self, lower: float, upper: float, num_nodes: int
    ) -> list[float]:
        """``num_nodes`` Chebyshev nodes on ``[lower, upper]`` (``num_nodes >= 1``)."""
        ...
    def kth_derivative_value(
        self, c: float, poly: Polynomial | list[float], k: int
    ) -> float:
        """``P^{(k)}(c)`` via repeated Horner division; ``0`` when ``k > deg(P)``."""
        ...

    # ---- Lagrange ----
    def Lagrange_polynomial(
        self, nodes: list[float], values: list[float]
    ) -> Polynomial:
        """Lagrange interpolation polynomial through ``(nodes, values)``."""
        ...

    # ---- Newton general ----
    def Newton_general_forward_polynomial(
        self, nodes: list[float], values: list[float]
    ) -> Polynomial:
        """Newton forward interpolation polynomial (general nodes)."""
        ...
    def Newton_general_backward_polynomial(
        self, nodes: list[float], values: list[float]
    ) -> Polynomial:
        """Newton backward interpolation polynomial (general nodes)."""
        ...

    # ---- Newton equidistant ----
    def Newton_equidistant_forward_polynomial(
        self, nodes: list[float], values: list[float]
    ) -> Polynomial:
        """Newton forward polynomial (equidistant nodes, param ``t = (x - x0)/h``)."""
        ...
    def Newton_equidistant_backward_polynomial(
        self, nodes: list[float], values: list[float]
    ) -> Polynomial:
        """Newton backward polynomial (equidistant nodes, param ``t = (x - xn)/h``)."""
        ...

    # ---- Gauss ----
    def GaussI_polynomial(
        self, nodes: list[float], values: list[float]
    ) -> Polynomial:
        """Gauss I polynomial (equidistant nodes, param ``t = (x - xn)/h``)."""
        ...
    def GaussII_polynomial(
        self, nodes: list[float], values: list[float]
    ) -> Polynomial:
        """Gauss II polynomial (equidistant nodes, param ``t = (x - xn)/h``)."""
        ...

    # ---- Stirling ----
    def Stirling_polynomial(
        self, nodes: list[float], values: list[float]
    ) -> Polynomial:
        """Stirling interpolation polynomial (equidistant nodes)."""
        ...

    # ---- Bessel ----
    def Bessel_polynomial(
        self, nodes: list[float], values: list[float]
    ) -> Polynomial:
        """Bessel polynomial (equidistant nodes, param ``u = (x - x(n-1))/h - 0.5``)."""
        ...

    # ---- Interval selection ----
    def select_equidistant_interpolation_points(
        self,
        nodes: list[float],
        values: list[float],
        x: float,
        point_count: int,
    ) -> tuple[int, int]:
        """Select ``point_count`` equidistant samples around ``x``. Returns ``(start, end)`` inclusive."""
        ...
    def select_interpolation_points(
        self,
        nodes: list[float],
        values: list[float],
        x: float,
        point_count: int,
    ) -> tuple[int, int]:
        """Select ``point_count`` samples from monotonic nodes around ``x``. Returns ``(start, end)`` inclusive."""
        ...
    def find_monotonic_intervals(
        self, values: list[float]
    ) -> list[tuple[int, int]]:
        """Monotonic segments in ``values``. Returns list of ``(start, end)`` pairs."""
        ...
    def find_isolating_intervals(
        self, values: list[float], target_value: float
    ) -> list[tuple[int, int]]:
        """Monotonic intervals bracketing ``target_value``. Returns list of ``(start, end)`` pairs."""
        ...

    # ---- Inverse interpolation ----
    def inverse_interpolation(
        self,
        nodes: list[float],
        values: list[float],
        target_value: float,
        max_points: int,
        method: InterpolationMethod = InterpolationMethod.NewtonForward,
    ) -> list[float]:
        """Roots of ``f(x) = target_value`` via inverse interpolation."""
        ...
    def inverse_iteration(
        self,
        nodes: list[float],
        values: list[float],
        target_value: float,
        tolerance: float,
        max_points: int,
    ) -> list[float]:
        """Roots of ``f(x) = target_value`` via iterative inverse interpolation."""
        ...

    # ---- Splines ----
    def find_linear_spline(
        self, nodes: list[float], values: list[float]
    ) -> SplineFunction:
        """Piecewise linear spline through ``(nodes, values)``."""
        ...
    def find_quadratic_spline(
        self,
        nodes: list[float],
        values: list[float],
        alpha: float,
        beta: float,
        gamma: float,
    ) -> SplineFunction:
        """C1 quadratic spline with endpoint condition ``alpha*S'(x0) + beta*S'(xn) = gamma``."""
        ...
    def find_cubic_spline(
        self,
        nodes: list[float],
        values: list[float],
        boundary_conditions: list[list[float]],
    ) -> SplineFunction:
        """Cubic spline; ``boundary_conditions`` is a 2×5 list-of-lists (augmented rows)."""
        ...

    # ---- LSP ----
    def find_approximation_coefficients(
        self,
        nodes: list[float],
        values: list[float],
        basis: list[Callable[[float], float]],
    ) -> list[float]:
        """Least-squares coefficients for arbitrary basis functions."""
        ...

detail: _DetailModule

# ---- Horner display ----
def display_polynomial_evaluation(x: float, poly: Polynomial | list[float]) -> None:
    """Print ``P(x)`` and ``P(x0)`` as plain text."""
    ...

def display_horner_division(
    a: Polynomial | list[float],
    b: Polynomial | list[float],
) -> None:
    """Print ``A = B * Q + R`` for degree-1 ``B``."""
    ...

def display_chebyshev_nodes(lower: float, upper: float, num_nodes: int) -> None:
    """Render the ``i | x_i`` FTXUI table."""
    ...

def display_kth_derivative(c: float, poly: Polynomial | list[float], k: int) -> None:
    """Render the generalized Horner tableau plus ``P^{(k)}(c) = k! * R``."""
    ...

# ---- Interpolation display ----
def display_Lagrange_interpolation(nodes: list[float], values: list[float]) -> None:
    """Render Lagrange interpolation tables and polynomial."""
    ...

def display_Newton_general_forward(nodes: list[float], values: list[float]) -> None:
    """Render Newton general forward interpolation tables and polynomial."""
    ...

def display_Newton_general_backward(nodes: list[float], values: list[float]) -> None:
    """Render Newton general backward interpolation tables and polynomial."""
    ...

def display_Newton_equidistant_forward(nodes: list[float], values: list[float]) -> None:
    """Render Newton equidistant forward interpolation tables and polynomial."""
    ...

def display_Newton_equidistant_backward(nodes: list[float], values: list[float]) -> None:
    """Render Newton equidistant backward interpolation tables and polynomial."""
    ...

def display_GaussI_interpolation(nodes: list[float], values: list[float]) -> None:
    """Render Gauss I interpolation tables and polynomial."""
    ...

def display_GaussII_interpolation(nodes: list[float], values: list[float]) -> None:
    """Render Gauss II interpolation tables and polynomial."""
    ...

def display_Stirling_interpolation(nodes: list[float], values: list[float]) -> None:
    """Render Stirling interpolation tables and polynomial."""
    ...

def display_Bessel_interpolation(nodes: list[float], values: list[float]) -> None:
    """Render Bessel interpolation tables and polynomial."""
    ...

# ---- Interval display ----
def display_equidistant_interpolation_points(
    nodes: list[float], values: list[float], x: float, point_count: int
) -> None:
    """Render the selected equidistant interpolation interval."""
    ...

def display_monotonic_intervals(nodes: list[float], values: list[float]) -> None:
    """Render all monotonic segments in the sample values."""
    ...

def display_isolating_intervals(
    nodes: list[float], values: list[float], target_value: float
) -> None:
    """Render monotonic intervals bracketing ``target_value``."""
    ...

# ---- Inverse interpolation display ----
def display_inverse_interpolation(
    nodes: list[float],
    values: list[float],
    target_value: float,
    max_points: int,
    method: InterpolationMethod = InterpolationMethod.NewtonForward,
) -> None:
    """Render full inverse interpolation: isolating intervals, each step, and roots."""
    ...

def display_inverse_iteration(
    nodes: list[float],
    values: list[float],
    target_value: float,
    tolerance: float,
    max_points: int,
) -> None:
    """Render iterative inverse interpolation: isolating intervals, each iteration step, and roots."""
    ...

# ---- Spline display ----
def display_linear_spline(nodes: list[float], values: list[float]) -> None:
    """Render the piecewise linear spline coefficient table."""
    ...

def display_quadratic_spline(
    nodes: list[float],
    values: list[float],
    alpha: float,
    beta: float,
    gamma: float,
) -> None:
    """Render the C1 quadratic spline coefficient table."""
    ...

def display_cubic_spline(
    nodes: list[float],
    values: list[float],
    boundary_conditions: list[list[float]],
) -> None:
    """Render the cubic spline system and coefficient table.
    ``boundary_conditions`` is a 2×5 list-of-lists (augmented rows).
    """
    ...

# ---- LSP display ----
def display_linear_lsp(nodes: list[float], values: list[float]) -> None:
    """Render the linear least-squares design matrix and ``F(x) = a*x + b``."""
    ...

def display_quadratic_lsp(nodes: list[float], values: list[float]) -> None:
    """Render the quadratic least-squares design matrix and ``F(x) = a*x^2 + b*x + c``."""
    ...

def display_general_lsp(
    nodes: list[float],
    values: list[float],
    basis: list[Callable[[float], float]],
) -> None:
    """Render least-squares over arbitrary basis functions."""
    ...
