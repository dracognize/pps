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

## detail::inverse_iteration

```cpp
RealList detail::inverse_iteration(
    const RealList &nodes,
    const RealList &values,
    Real target_value,
    Real tolerance,
    SizeType max_points);
```

Finds roots of the tabulated function for `target_value` using iterative inverse interpolation.

## display_inverse_iteration

```cpp
void display_inverse_iteration(
    const RealList &nodes,
    const RealList &values,
    Real target_value,
    Real tolerance,
    SizeType max_points);
```

Renders isolating intervals, selection tables, polynomials, iterative values table, and the resulting roots.
