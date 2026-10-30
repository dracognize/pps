#include <chapter1/interpolation/optimal.hpp>
#include <chapter1/interpolation/horner.hpp>
#include <chapter1/interpolation/polynomial/Lagrange.hpp>
#include <chapter1/interpolation/polynomial/Newton.hpp>
#include <chapter1/interpolation/polynomial/Gauss.hpp>
#include <chapter1/interpolation/polynomial/Stirling.hpp>
#include <chapter1/interpolation/polynomial/Bessel.hpp>
#include <chapter1/interpolation/interval.hpp>
#include <chapter1/interpolation/inverse.hpp>
#include <chapter1/spline/linear.hpp>
#include <chapter1/spline/quadratic.hpp>
#include <chapter1/spline/cubic.hpp>
#include <chapter1/lsp.hpp>
#include <math_object.hpp>

// Python bindings: thin adapters over the C++ API.
//
// Rules:
//   - Pure compute lives in `detail::` and is exposed as `pps.detail`;
//     call `detail::` directly, never via a forwarding wrapper.
//   - `display_*` stays top-level: compute via `detail::` then render.
//   - Argument order mirrors C++: point first `(x, poly)`, dividend first
//     `(a, b)`, interval first `(lower, upper, n)`.
//   - `Polynomial | list[float]` inputs arrive as `vector<Real>` and convert
//     via `to_poly`; async validation (`check_*`) maps C++ preconditions to
//     Python `ValueError`s.

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

#include <Eigen/Core>
#include <fmt/format.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace {

	/// Convert a Python list of ascending coeffs into the C++ type.
	auto to_poly(const std::vector<Real> &coeffs) -> Polynomial {
		return Polynomial(coeffs.begin(), coeffs.end());
	}

	/// Element-wise add (pads the shorter operand with zeros).
	auto add_poly(const Polynomial &a, const Polynomial &b) -> Polynomial {
		Polynomial out(std::max(a.size(), b.size()), 0.0);
		for (std::size_t i = 0; i < a.size(); ++i)
			out[i] += a[i];
		for (std::size_t i = 0; i < b.size(); ++i)
			out[i] += b[i];
		return out;
	}

	/// Element-wise subtract (pads the shorter operand with zeros).
	auto sub_poly(const Polynomial &a, const Polynomial &b) -> Polynomial {
		Polynomial out(std::max(a.size(), b.size()), 0.0);
		for (std::size_t i = 0; i < a.size(); ++i)
			out[i] += a[i];
		for (std::size_t i = 0; i < b.size(); ++i)
			out[i] -= b[i];
		return out;
	}

	/// Naive O(n*m) convolution; empty input yields the zero polynomial.
	auto mul_poly(const Polynomial &a, const Polynomial &b) -> Polynomial {
		if (a.empty() || b.empty())
			return {};
		Polynomial out(a.size() + b.size() - 1, 0.0);
		for (std::size_t i = 0; i < a.size(); ++i)
			for (std::size_t j = 0; j < b.size(); ++j)
				out[i + j] += a[i] * b[j];
		return out;
	}

	/// Enforce the Horner precondition: divisor must be degree 1 (2 coeffs).
	auto check_degree_one(const Polynomial &b) -> void {
		if (b.size() != 2) {
			throw std::invalid_argument(fmt::format(
				"horner divisor must have degree 1 (2 coefficients), got {}", b.size()));
		}
	}

	/// Enforce `num_nodes >= 1` before narrowing `long long` -> `SizeType`.
	auto check_node_count(long long n) -> SizeType {
		if (n < 1) {
			throw std::invalid_argument(
				fmt::format("num_nodes must be a positive integer, got {}", n));
		}
		return static_cast<SizeType>(n);
	}

	/// Enforce `k >= 0` before narrowing `long long` -> `SizeType`.
	auto check_derivative_order(long long k) -> SizeType {
		if (k < 0) {
			throw std::invalid_argument(fmt::format("k must be a non-negative integer, got {}", k));
		}
		return static_cast<SizeType>(k);
	}

	/// Enforce `point_count >= 1` before narrowing `long long` -> `SizeType`.
	auto check_point_count(long long n) -> SizeType {
		if (n < 1) {
			throw std::invalid_argument(
				fmt::format("point_count must be a positive integer, got {}", n));
		}
		return static_cast<SizeType>(n);
	}

	/// Enforce `max_points >= 1` before narrowing `long long` -> `SizeType`.
	auto check_max_points(long long n) -> SizeType {
		if (n < 1) {
			throw std::invalid_argument(
				fmt::format("max_points must be a positive integer, got {}", n));
		}
		return static_cast<SizeType>(n);
	}

	/// Convert a 2x5 Python nested list into the Eigen fixed-size matrix used by
	/// cubic spline routines. Raises `ValueError` if the shape is wrong.
	auto to_bc_matrix(const std::vector<std::vector<double>> &rows)
		-> Eigen::Matrix<double, 2, 5> {
		if (rows.size() != 2) {
			throw std::invalid_argument(fmt::format(
				"boundary_conditions must have exactly 2 rows, got {}", rows.size()));
		}
		for (std::size_t r = 0; r < 2; ++r) {
			if (rows[r].size() != 5) {
				throw std::invalid_argument(fmt::format(
					"boundary_conditions row {} must have exactly 5 columns, got {}", r,
					rows[r].size()));
			}
		}
		Eigen::Matrix<double, 2, 5> m;
		for (std::size_t r = 0; r < 2; ++r)
			for (std::size_t c = 0; c < 5; ++c)
				m(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c)) = rows[r][c];
		return m;
	}

} // namespace

PYBIND11_MODULE(pps, m) {
	m.doc() = "pps: Numeric Method(C++ core exposed to Python)";

	// Mirror the C++ math types as Python names. Real IS float (no separate
	// extended type), RealList IS list, SizeType IS int — so `x: Real = 10`
	// annotates exactly like C++ `Real x = 10` with zero overhead.
	const py::object builtins = py::module::import("builtins");
	m.attr("Real")       = builtins.attr("float");
	m.attr("RealList")   = builtins.attr("list");
	m.attr("SizeType")   = builtins.attr("int");

	// ---- Math types (mirror `inc/math_object.hpp`) ----
	py::class_<Polynomial>(m, "Polynomial")
		.def(py::init([](const std::vector<Real> &coeffs) { return to_poly(coeffs); }),
			 py::arg("coeffs") = std::vector<Real>{},
			 "Coeffs in ascending order: coeffs[i] is the x^i coefficient.")
		.def("__len__", [](const Polynomial &p) { return p.size(); })
		.def("__getitem__",
			 [](const Polynomial &p, std::size_t i) -> Real {
				 if (i >= p.size())
					 throw py::index_error();
				 return p[i];
			 })
		.def("__setitem__",
			 [](Polynomial &p, std::size_t i, Real v) {
				 if (i >= p.size())
					 throw py::index_error();
				 p[i] = v;
			 })
		.def(
			"__iter__",
			[](const Polynomial &p) { return py::make_iterator(p.begin(), p.end()); },
			py::keep_alive<0, 1>())
		.def("to_list", [](const Polynomial &p) { return std::vector<Real>(p.begin(), p.end()); })
		.def("__add__", [](const Polynomial &a, const Polynomial &b) { return add_poly(a, b); })
		.def("__sub__", [](const Polynomial &a, const Polynomial &b) { return sub_poly(a, b); })
		.def("__mul__", [](const Polynomial &a, const Polynomial &b) { return mul_poly(a, b); })
		.def("__neg__",
			 [](const Polynomial &a) {
				 Polynomial out = a;
				 for (auto &c : out)
					 c = -c;
				 return out;
			 })
		.def("__rmul__",
			 [](const Polynomial &a, Real s) {
				 Polynomial out = a;
				 for (auto &c : out)
					 c *= s;
				 return out;
			 })
		.def("__mul__",
			 [](const Polynomial &a, Real s) {
				 Polynomial out = a;
				 for (auto &c : out)
					 c *= s;
				 return out;
			 })
		.def("__truediv__",
			 [](const Polynomial &a, Real s) {
				 if (s == 0.0)
					 throw std::invalid_argument("division by zero");
				 Polynomial out = a;
				 for (auto &c : out)
					 c /= s;
				 return out;
			 })
		.def("__repr__", [](const Polynomial &p) { return fmt::format("Polynomial({})", p); })
		.def("__str__", [](const Polynomial &p) { return fmt::format("{}", p); });

	// `DerivativeHistory` is intentionally *not* bound: it is a display-only
	// tableau (`quotients` + `remainders`) consumed by `display_kth_derivative`.
	// Python users needing the numbers should call `kth_derivative` per order.
	py::class_<PolynomialDivision>(m, "PolynomialDivision")
		.def_readonly("quotient", &PolynomialDivision::quotient)
		.def_readonly("remainder", &PolynomialDivision::remainder)
		.def("__repr__", [](const PolynomialDivision &d) {
			return fmt::format(
				"PolynomialDivision(quotient={}, remainder={})", d.quotient, d.remainder);
		});

	// SplineFunction is returned by detail compute functions and also exposes
	// `value(x)` so Python can evaluate the spline at any point.
	py::class_<SplineFunction>(m, "SplineFunction")
		.def("value",
			 [](const SplineFunction &s, Real x) { return s.value(x); },
			 py::arg("x"),
			 "Evaluate the spline at x.")
		.def_readonly("nodes", &SplineFunction::nodes)
		.def("__repr__", [](const SplineFunction &s) {
			return fmt::format("SplineFunction({})", s);
		});

	// InterpolationMethod enum for inverse interpolation.
	py::enum_<InterpolationMethod>(m, "InterpolationMethod")
		.value("Lagrange", InterpolationMethod::Lagrange)
		.value("NewtonForward", InterpolationMethod::NewtonForward)
		.value("NewtonBackward", InterpolationMethod::NewtonBackward)
		.export_values();

	// ---- Pure compute (mirror C++ detail::; no I/O) ----
	py::module_ detail = m.def_submodule("detail", "Core compute API (mirrors C++ detail::)");
	detail.def(
		"evaluate_polynomial",
		[](Real x, const std::vector<Real> &coeffs) {
			return detail::evaluate_polynomial(x, to_poly(coeffs));
		},
		py::arg("x"),
		py::arg("poly"),
		"Evaluate P(x) with Horner scheme. poly may be a Polynomial or a plain list.");

	detail.def(
		"divide_by_linear",
		[](const std::vector<Real> &a, const std::vector<Real> &b) {
			Polynomial pa = to_poly(a);
			Polynomial pb = to_poly(b);
			check_degree_one(pb);
			return detail::divide_by_linear(pa, pb);
		},
		py::arg("a"),
		py::arg("b"),
		"Divide a by degree-1 b. Returns PolynomialDivision(quotient, remainder).");

	detail.def(
		"generate_chebyshev_nodes",
		[](Real lower, Real upper, long long num_nodes) {
			return detail::generate_chebyshev_nodes(lower, upper, check_node_count(num_nodes));
		},
		py::arg("lower"),
		py::arg("upper"),
		py::arg("num_nodes"),
		"Chebyshev nodes on [lower, upper] as a list of floats.");

	detail.def(
		"kth_derivative_value",
		[](Real c, const std::vector<Real> &coeffs, long long k) {
			return detail::kth_derivative_value(c, to_poly(coeffs), check_derivative_order(k));
		},
		py::arg("c"),
		py::arg("poly"),
		py::arg("k"),
		"Value of P^{(k)}(c) via repeated Horner division. Returns 0 when k > deg(P).");

	// ---- Lagrange ----
	detail.def(
		"Lagrange_polynomial",
		[](const RealList &nodes, const RealList &values) {
			return detail::Lagrange_polynomial(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"),
		"Lagrange interpolation polynomial through (nodes, values).");

	// ---- Newton general ----
	detail.def(
		"Newton_general_forward_polynomial",
		[](const RealList &nodes, const RealList &values) {
			return detail::Newton_general_forward_polynomial(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"),
		"Newton forward interpolation polynomial (general nodes).");

	detail.def(
		"Newton_general_backward_polynomial",
		[](const RealList &nodes, const RealList &values) {
			return detail::Newton_general_backward_polynomial(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"),
		"Newton backward interpolation polynomial (general nodes).");

	// ---- Newton equidistant ----
	detail.def(
		"Newton_equidistant_forward_polynomial",
		[](const RealList &nodes, const RealList &values) {
			return detail::Newton_equidistant_forward_polynomial(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"),
		"Newton forward interpolation polynomial (equidistant nodes, param t = (x - x0)/h).");

	detail.def(
		"Newton_equidistant_backward_polynomial",
		[](const RealList &nodes, const RealList &values) {
			return detail::Newton_equidistant_backward_polynomial(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"),
		"Newton backward interpolation polynomial (equidistant nodes, param t = (x - xn)/h).");

	// ---- Gauss ----
	detail.def(
		"GaussI_polynomial",
		[](const RealList &nodes, const RealList &values) {
			return detail::GaussI_polynomial(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"),
		"Gauss I interpolation polynomial (equidistant nodes, param t = (x - xn)/h).");

	detail.def(
		"GaussII_polynomial",
		[](const RealList &nodes, const RealList &values) {
			return detail::GaussII_polynomial(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"),
		"Gauss II interpolation polynomial (equidistant nodes, param t = (x - xn)/h).");

	// ---- Stirling ----
	detail.def(
		"Stirling_polynomial",
		[](const RealList &nodes, const RealList &values) {
			return detail::Stirling_polynomial(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"),
		"Stirling interpolation polynomial (equidistant nodes).");

	// ---- Bessel ----
	detail.def(
		"Bessel_polynomial",
		[](const RealList &nodes, const RealList &values) {
			return detail::Bessel_polynomial(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"),
		"Bessel interpolation polynomial (equidistant nodes, param u = (x - x(n-1))/h - 0.5).");

	// ---- Interval selection ----
	detail.def(
		"select_equidistant_interpolation_points",
		[](const RealList &nodes, const RealList &values, Real x, long long point_count) {
			auto iv = detail::select_equidistant_interpolation_points(
				nodes, values, x, check_point_count(point_count));
			return py::make_tuple(iv.start, iv.end);
		},
		py::arg("nodes"),
		py::arg("values"),
		py::arg("x"),
		py::arg("point_count"),
		"Select point_count equidistant samples around x. Returns (start, end) inclusive indices.");

	detail.def(
		"select_interpolation_points",
		[](const RealList &nodes, const RealList &values, Real x, long long point_count) {
			auto iv = detail::select_interpolation_points(
				nodes, values, x, check_point_count(point_count));
			return py::make_tuple(iv.start, iv.end);
		},
		py::arg("nodes"),
		py::arg("values"),
		py::arg("x"),
		py::arg("point_count"),
		"Select point_count samples from monotonic nodes around x. Returns (start, end) inclusive indices.");

	detail.def(
		"find_monotonic_intervals",
		[](const RealList &values) {
			auto ivs = detail::find_monotonic_intervals(values);
			std::vector<py::tuple> out;
			out.reserve(ivs.size());
			for (const auto &iv : ivs)
				out.push_back(py::make_tuple(iv.start, iv.end));
			return out;
		},
		py::arg("values"),
		"Find monotonic segments in values. Returns list of (start, end) inclusive index pairs.");

	detail.def(
		"find_isolating_intervals",
		[](const RealList &values, Real target_value) {
			auto ivs = detail::find_isolating_intervals(values, target_value);
			std::vector<py::tuple> out;
			out.reserve(ivs.size());
			for (const auto &iv : ivs)
				out.push_back(py::make_tuple(iv.start, iv.end));
			return out;
		},
		py::arg("values"),
		py::arg("target_value"),
		"Find monotonic intervals bracketing target_value. Returns list of (start, end) pairs.");

	// ---- Inverse interpolation ----
	detail.def(
		"inverse_interpolation",
		[](const RealList &nodes,
		   const RealList &values,
		   Real           target_value,
		   long long      max_points,
		   InterpolationMethod method) {
			return detail::inverse_interpolation(
				nodes, values, target_value, check_max_points(max_points), method);
		},
		py::arg("nodes"),
		py::arg("values"),
		py::arg("target_value"),
		py::arg("max_points"),
		py::arg("method") = InterpolationMethod::NewtonForward,
		"Find roots of f(x) = target_value by inverse interpolation. Returns list of roots.");

	// ---- Splines ----
	detail.def(
		"find_linear_spline",
		[](const RealList &nodes, const RealList &values) {
			return detail::find_linear_spline(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"),
		"Build a piecewise linear spline. Returns SplineFunction.");

	detail.def(
		"find_quadratic_spline",
		[](const RealList &nodes, const RealList &values, Real alpha, Real beta, Real gamma) {
			return detail::find_quadratic_spline(nodes, values, alpha, beta, gamma);
		},
		py::arg("nodes"),
		py::arg("values"),
		py::arg("alpha"),
		py::arg("beta"),
		py::arg("gamma"),
		"Build a C1 quadratic spline with endpoint derivative condition alpha*S'(x0) + beta*S'(xn) = gamma. Returns SplineFunction.");

	detail.def(
		"find_cubic_spline",
		[](const RealList                        &nodes,
		   const RealList                        &values,
		   const std::vector<std::vector<double>> &boundary_conditions) {
			return detail::find_cubic_spline(nodes, values, to_bc_matrix(boundary_conditions));
		},
		py::arg("nodes"),
		py::arg("values"),
		py::arg("boundary_conditions"),
		"Build a cubic spline. boundary_conditions is a 2x5 list-of-lists (augmented rows). Returns SplineFunction.");

	// ---- LSP compute ----
	detail.def(
		"find_approximation_coefficients",
		[](const RealList &nodes, const RealList &values, const FunctionsList &basis) {
			return detail::find_approximation_coefficients(nodes, values, basis);
		},
		py::arg("nodes"),
		py::arg("values"),
		py::arg("basis"),
		"Least-squares coefficients for arbitrary basis functions. Returns list of coefficients.");

	// ---- Display (compute + render to stdout; return None) ----
	m.def(
		"display_polynomial_evaluation",
		[](Real x, const std::vector<Real> &coeffs) {
			display_polynomial_evaluation(x, to_poly(coeffs));
		},
		py::arg("x"),
		py::arg("poly"));

	m.def(
		"display_horner_division",
		[](const std::vector<Real> &a, const std::vector<Real> &b) {
			Polynomial pa = to_poly(a);
			Polynomial pb = to_poly(b);
			check_degree_one(pb);
			display_horner_division(pa, pb);
		},
		py::arg("a"),
		py::arg("b"));

	m.def(
		"display_chebyshev_nodes",
		[](Real lower, Real upper, long long num_nodes) {
			display_chebyshev_nodes(lower, upper, check_node_count(num_nodes));
		},
		py::arg("lower"),
		py::arg("upper"),
		py::arg("num_nodes"));

	m.def(
		"display_kth_derivative",
		[](Real c, const std::vector<Real> &coeffs, long long k) {
			display_kth_derivative(c, to_poly(coeffs), check_derivative_order(k));
		},
		py::arg("c"),
		py::arg("poly"),
		py::arg("k"));

	// ---- Interpolation display ----
	m.def(
		"display_Lagrange_interpolation",
		[](const RealList &nodes, const RealList &values) {
			display_Lagrange_interpolation(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"));

	m.def(
		"display_Newton_general_forward",
		[](const RealList &nodes, const RealList &values) {
			display_Newton_general_forward(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"));

	m.def(
		"display_Newton_general_backward",
		[](const RealList &nodes, const RealList &values) {
			display_Newton_general_backward(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"));

	m.def(
		"display_Newton_equidistant_forward",
		[](const RealList &nodes, const RealList &values) {
			display_Newton_equidistant_forward(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"));

	m.def(
		"display_Newton_equidistant_backward",
		[](const RealList &nodes, const RealList &values) {
			display_Newton_equidistant_backward(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"));

	m.def(
		"display_GaussI_interpolation",
		[](const RealList &nodes, const RealList &values) {
			display_GaussI_interpolation(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"));

	m.def(
		"display_GaussII_interpolation",
		[](const RealList &nodes, const RealList &values) {
			display_GaussII_interpolation(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"));

	m.def(
		"display_Stirling_interpolation",
		[](const RealList &nodes, const RealList &values) {
			display_Stirling_interpolation(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"));

	m.def(
		"display_Bessel_interpolation",
		[](const RealList &nodes, const RealList &values) {
			display_Bessel_interpolation(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"));

	// ---- Interval display ----
	m.def(
		"display_equidistant_interpolation_points",
		[](const RealList &nodes, const RealList &values, Real x, long long point_count) {
			display_equidistant_interpolation_points(
				nodes, values, x, check_point_count(point_count));
		},
		py::arg("nodes"),
		py::arg("values"),
		py::arg("x"),
		py::arg("point_count"));

	m.def(
		"display_monotonic_intervals",
		[](const RealList &nodes, const RealList &values) {
			display_monotonic_intervals(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"));

	m.def(
		"display_isolating_intervals",
		[](const RealList &nodes, const RealList &values, Real target_value) {
			display_isolating_intervals(nodes, values, target_value);
		},
		py::arg("nodes"),
		py::arg("values"),
		py::arg("target_value"));

	// ---- Inverse interpolation display ----
	m.def(
		"display_inverse_interpolation",
		[](const RealList        &nodes,
		   const RealList        &values,
		   Real                   target_value,
		   long long              max_points,
		   InterpolationMethod    method) {
			display_inverse_interpolation(
				nodes, values, target_value, check_max_points(max_points), method);
		},
		py::arg("nodes"),
		py::arg("values"),
		py::arg("target_value"),
		py::arg("max_points"),
		py::arg("method") = InterpolationMethod::NewtonForward);

	// ---- Spline display ----
	m.def(
		"display_linear_spline",
		[](const RealList &nodes, const RealList &values) {
			display_linear_spline(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"));

	m.def(
		"display_quadratic_spline",
		[](const RealList &nodes, const RealList &values, Real alpha, Real beta, Real gamma) {
			display_quadratic_spline(nodes, values, alpha, beta, gamma);
		},
		py::arg("nodes"),
		py::arg("values"),
		py::arg("alpha"),
		py::arg("beta"),
		py::arg("gamma"));

	m.def(
		"display_cubic_spline",
		[](const RealList                        &nodes,
		   const RealList                        &values,
		   const std::vector<std::vector<double>> &boundary_conditions) {
			display_cubic_spline(nodes, values, to_bc_matrix(boundary_conditions));
		},
		py::arg("nodes"),
		py::arg("values"),
		py::arg("boundary_conditions"),
		"boundary_conditions is a 2x5 list-of-lists (augmented rows).");

	// ---- LSP display ----
	m.def(
		"display_linear_lsp",
		[](const RealList &nodes, const RealList &values) {
			display_linear_lsp(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"));

	m.def(
		"display_quadratic_lsp",
		[](const RealList &nodes, const RealList &values) {
			display_quadratic_lsp(nodes, values);
		},
		py::arg("nodes"),
		py::arg("values"));

	m.def(
		"display_general_lsp",
		[](const RealList &nodes, const RealList &values, const FunctionsList &basis) {
			display_general_lsp(nodes, values, basis);
		},
		py::arg("nodes"),
		py::arg("values"),
		py::arg("basis"));
}
