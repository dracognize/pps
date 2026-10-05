# Introduction

`pps` is a numerical methods playground: a C++23 core (`pps_lib`) with an
FTXUI command-line demo and a Python binding of the same API.

| Symbol | Meaning |
|---|---|
| `Polynomial` | Ascending coefficients (`poly[i]` is `x^i`) |
| `detail::evaluate_polynomial(x, poly)` | Horner evaluation of `P(x)` |
| `detail::divide_by_linear(a, b)` | Division by degree-1 `b` |
| `detail::kth_derivative_value(c, poly, k)` | `P^{(k)}(c)` via repeated Horner division |
| `detail::generate_chebyshev_nodes(lo, hi, n)` | Chebyshev nodes on `[lo, hi]` |
| `display_*` | Same computation rendered to stdout |

Argument order is fixed everywhere (C++ and Python): point first
(`detail::evaluate_polynomial(x, poly)`,
`detail::kth_derivative_value(c, poly, k)`), dividend first
(`detail::divide_by_linear(a, b)`), interval first
(`detail::generate_chebyshev_nodes(lo, hi, n)`). Python exposes
the same names under `pps.detail`.

```sh
mise run dev      # build everything, run the C++ demo
mise run py:run   # run the Python demo against the compiled module
mise run docs     # build this book
```
