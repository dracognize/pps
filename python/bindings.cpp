#include <chapter1/horner.hpp>
#include <chapter1/interpolation.hpp>
#include <math_object.hpp>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <fmt/format.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace {

	auto to_poly(const std::vector<Real> &coeffs) -> Polynomial {
		return Polynomial(coeffs.begin(), coeffs.end());
	}

	auto add_poly(const Polynomial &a, const Polynomial &b) -> Polynomial {
		Polynomial out(std::max(a.size(), b.size()), 0.0);
		for (std::size_t i = 0; i < a.size(); ++i)
			out[i] += a[i];
		for (std::size_t i = 0; i < b.size(); ++i)
			out[i] += b[i];
		return out;
	}

	auto sub_poly(const Polynomial &a, const Polynomial &b) -> Polynomial {
		Polynomial out(std::max(a.size(), b.size()), 0.0);
		for (std::size_t i = 0; i < a.size(); ++i)
			out[i] += a[i];
		for (std::size_t i = 0; i < b.size(); ++i)
			out[i] -= b[i];
		return out;
	}

	auto mul_poly(const Polynomial &a, const Polynomial &b) -> Polynomial {
		if (a.empty() || b.empty())
			return {};
		Polynomial out(a.size() + b.size() - 1, 0.0);
		for (std::size_t i = 0; i < a.size(); ++i)
			for (std::size_t j = 0; j < b.size(); ++j)
				out[i + j] += a[i] * b[j];
		return out;
	}

	auto check_degree_one(const Polynomial &b) -> void {
		if (b.size() != 2) {
			throw std::invalid_argument(fmt::format(
				"horner divisor must have degree 1 (2 coefficients), got {}", b.size()));
		}
	}

	auto check_node_count(long long n) -> SizeType {
		if (n < 1) {
			throw std::invalid_argument(
				fmt::format("num_nodes must be a positive integer, got {}", n));
		}
		return static_cast<SizeType>(n);
	}

} // namespace

PYBIND11_MODULE(pps, m) {
	m.doc() = "pps: Numeric Method(C++ core exposed to Python)";

	// Mirror the C++ math types as Python names. Real IS float (no separate
	// extended type), RealList IS list, SizeType IS int — so `x: Real = 10`
	// annotates exactly like C++ `Real x = 10` with zero overhead.
	const py::object builtins = py::module::import("builtins");
	m.attr("Real")			  = builtins.attr("float");
	m.attr("RealList")		  = builtins.attr("list");
	m.attr("SizeType")		  = builtins.attr("int");

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

	py::class_<PolynomialDivision>(m, "PolynomialDivision")
		.def_readonly("quotient", &PolynomialDivision::quotient)
		.def_readonly("remainder", &PolynomialDivision::remainder)
		.def("__repr__", [](const PolynomialDivision &d) {
			return fmt::format(
				"PolynomialDivision(quotient={}, remainder={})", d.quotient, d.remainder);
		});

	m.def(
		"evaluate",
		[](Real x, const std::vector<Real> &coeffs) {
			return detail::evaluate_polynomial(x, to_poly(coeffs));
		},
		py::arg("x"),
		py::arg("poly"),
		"Evaluate P(x) with Horner scheme. poly may be a Polynomial or a plain list.");

	m.def(
		"horner_division",
		[](const std::vector<Real> &a, const std::vector<Real> &b) {
			Polynomial pa = to_poly(a);
			Polynomial pb = to_poly(b);
			check_degree_one(pb);
			return detail::horner_division(pa, pb);
		},
		py::arg("a"),
		py::arg("b"),
		"Divide a by degree-1 b. Returns PolynomialDivision(quotient, remainder).");

	m.def(
		"chebyshev_nodes",
		[](Real lower, Real upper, long long num_nodes) {
			return detail::generate_chebyshev_nodes(lower, upper, check_node_count(num_nodes));
		},
		py::arg("lower"),
		py::arg("upper"),
		py::arg("num_nodes"),
		"Chebyshev nodes on [lower, upper] as a list of floats.");

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
}
