# pps

Numerical methods playground (C++23 core with an FTXUI CLI and a Python
binding). Each chapter lives under `inc/chapter<N>/` + `src/chapter<N>/`
with shared scalars in `inc/math_object.hpp`.

| Symbol | Meaning |
|---|---|
| `Polynomial` | Ascending coefficients (`poly[i]` is `x^i`) |
| `detail::evaluate_polynomial(x, poly)` | Horner evaluation of `P(x)` |
| `detail::divide_by_linear(a, b)` | Division by degree-1 `b` |
| `detail::kth_derivative_value(c, poly, k)` | `P^{(k)}(c)` via repeated Horner division |
| `detail::generate_chebyshev_nodes(lo, hi, n)` | Chebyshev nodes on `[lo, hi]` |
| `display_*` | Same computation rendered to stdout |

Python mirrors this: pure compute lives in `pps.detail`
(`detail.evaluate_polynomial`, `detail.divide_by_linear`, ...).

## API reference (mdBook)

```sh
mise run docs   # build build/docs/index.html, then open it
```

The book sources live in `docs/` (`book.toml` at the repo root) and reuse
the Yue light/dark palette (`docs/theme/pps.css`). `mdbook` itself is
mise-managed — `mise run bootstrap-docs` installs it (normal
build/run via `mise run bootstrap` only needs cmake/ninja, so docs
never slow down regular users; `mise run bootstrap-all` gets
everything for CI/fresh setups) — so no system package is needed.
