# Piecewise linear spline

Defined in `inc/chapter1/spline/linear.hpp`.

## `detail::render_linear_spline`

Renders a table with one row per interval. The columns show the interval and
the coefficients `a_i` and `b_i` of
`p_i(x) = a_i + b_i * (x - x_i)`, matching the local polynomial stored in
`SplineFunction::piecewise_poly`.

## `detail::display_linear_spline`

Builds the degree-1 spline from the supplied nodes and values, then displays its
coefficient table. Throws `std::invalid_argument` when the inputs do not satisfy
the requirements of `detail::find_linear_spline`.

## Quadratic spline

Defined in `inc/chapter1/spline/quadratic.hpp`. Each interval polynomial is
`p_i(x) = c_i + b_i * (x - x_i) + a_i * (x - x_i)^2`. The spline is C1 and
uses the endpoint condition `alpha*S'_0(x_0) + beta*S'_{n-1}(x_n) = gamma`,
where `beta = (-1)^n * alpha` and `alpha != 0`. The implementation forms and
solves the textbook system for the nodal derivatives `m_i`, then computes
`a_i = (m_{i+1} - m_i) / (2*h_i)`, `b_i = m_i`, and `c_i = y_i`.

`detail::find_quadratic_spline` constructs the spline,
`detail::render_quadratic_spline` displays its coefficient table, and
`display_quadratic_spline(nodes, values, alpha, beta, gamma)` performs both steps.
