# Interpolation

Defined in `inc/chapter1/interpolation.hpp`. Pure compute lives in
`detail::` (exposed to Python as `pps.detail`).

## detail::generate_chebyshev_nodes

```cpp
RealList detail::generate_chebyshev_nodes(Real lower, Real upper, SizeType num_nodes);
```

`num_nodes` Chebyshev nodes on `[lower, upper]`, ordered from the
left endpoint toward the right.

## display_chebyshev_nodes

```cpp
void display_chebyshev_nodes(Real lower, Real upper, SizeType num_nodes);
```

Renders the `i | x_i` table. Near-zero nodes print as `0`.
