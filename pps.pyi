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

from collections.abc import Iterator
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

class _DetailModule:
    """Pure compute (mirrors C++ ``detail::``); no I/O."""

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

detail: _DetailModule

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
