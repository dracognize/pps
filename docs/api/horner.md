# Horner scheme

Defined in `inc/chapter1/horner.hpp`. Pure compute lives in `detail::`
(exposed to Python as `pps.detail`); `display_*` functions render to
stdout.

## detail::evaluate_polynomial

```cpp
Real detail::evaluate_polynomial(Real x, const Polynomial &poly);
```

Horner evaluation of `P(x)`. Returns `0` for an empty polynomial.

## detail::divide_by_linear

```cpp
PolynomialDivision detail::divide_by_linear(const Polynomial &a, const Polynomial &b);
```

Division of `a` by degree-1 `b` (e.g. `{-c, 1}` for
`(x - c)`). Returns a zero division when `a` is empty or `b` does
not have exactly 2 coefficients.

## detail::kth_derivative_value

```cpp
Real detail::kth_derivative_value(Real c, const Polynomial &poly, SizeType k);
```

Value of `P^{(k)}(c)` via `k` divisions by `(x - c)`
followed by `k! * Q_k(c)`. Returns `0` when `k` exceeds the
polynomial degree. `k = 0` is plain evaluation.

## display_polynomial_evaluation

```cpp
void display_polynomial_evaluation(Real x, const Polynomial &poly);
```

Prints `P(x)` and `P(x0)` as plain text (no table).

## display_horner_division

```cpp
void display_horner_division(const Polynomial &a, const Polynomial &b);
```

Prints `A = B * Q + R` as plain text (no table).

## display_kth_derivative

```cpp
void display_kth_derivative(Real c, const Polynomial &poly, SizeType k);
```

Renders the generalized Horner tableau (columns `k = 0..k` of
triangular quotients plus a `Result` row of Taylor coefficients)
and the final `P^{(k)}(c) = k! * R` line. Prints a short
`k > deg(P)` note instead of a table when the order is too high.
