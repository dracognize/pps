# Core types

Defined in `inc/math_object.hpp`. Tolerances in `inc/config.hpp`.

## Scalars

```cpp
using Real = double;
using SizeType = std::size_t;
using RealList = std::vector<Real>;
```

## Polynomial

```cpp
struct Polynomial : std::vector<Real> { ... };
```

Ascending coefficient order: `poly[i]` is the coefficient of `x^i`.
An empty polynomial is the zero polynomial. Degree queries ignore
trailing coefficients below `config::kNearZeroThreshold` (`1e-12`).
Pretty-prints as `5*x^3 + 4*x^2 + ...` (`0` when empty).

## PolynomialDivision

```cpp
struct PolynomialDivision {
    Polynomial quotient;
    Real remainder = 0;
};
```

Result of dividing by a degree-1 divisor:
`a = b * quotient + remainder`. The quotient is empty when the
dividend was constant.

## DerivativeHistory

```cpp
struct DerivativeHistory {
    std::vector<Polynomial> quotients;
    std::vector<Real> remainders;
};
```

Generalized Horner tableau for derivatives at `x = c`. Column `j`
holds `Q_j`, the quotient after `j` divisions by `(x - c)`;
`remainders[j]` is `Q_j(c)`, i.e. the Taylor coefficient
`P^{(j)}(c) / j!`.

## Config

| Constant | Value | Meaning |
|---|---|---|
| `kNearZeroThreshold` | `1e-12` | Magnitude below which a value counts as zero |
| `kDisplayPrecision` | `7` | Significant digits in every display format |
| `kMaxFullTerms` | `99` | Term count above which printing collapses |
| `kTableMinWidth` | `12` | Minimum FTXUI column width |
