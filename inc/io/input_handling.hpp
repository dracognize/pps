#pragma once
/// @file input_handling.hpp
/// @brief Input helpers for numerical data.

#include <math_object.hpp>

namespace ih {
	/// @brief Read sample nodes and values from the data-point files.
	/// @details Reads whitespace-separated real numbers from
	///          `datas/data_points/nodes.txt` and `datas/data_points/values.txt`.
	/// @param nodes Receives the values from `nodes.txt`.
	/// @param values Receives the values from `values.txt`.
	/// @throws std::runtime_error if either file cannot be opened, contains invalid
	///         data, or the number of nodes and values differs.
	auto data_point_input(RealList &nodes, RealList &values) -> void;
} // namespace ih
