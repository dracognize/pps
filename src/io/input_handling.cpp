#include <io/input_handling.hpp>

#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace ih {
	namespace {
		auto read_real_list(const std::string &path) -> RealList {
			std::ifstream input(path);
			if (!input)
				throw std::runtime_error("Unable to open data file: " + path);

			RealList result;
			std::string token;
			while (input >> token) {
				try {
					size_t parsed_chars = 0;
					const Real value = std::stod(token, &parsed_chars);
					if (parsed_chars != token.size())
						throw std::runtime_error("Invalid numeric token");
					result.push_back(value);
				} catch (const std::exception &) {
					throw std::runtime_error("Invalid numeric data in file: " + path);
				}
			}
			if (input.bad())
				throw std::runtime_error("Unable to read data file: " + path);

			return result;
		}
	} // namespace

	auto nodes_input(RealList &nodes) -> void {
		RealList loaded_nodes  = read_real_list("datas/data_points/nodes.txt");
		nodes  = std::move(loaded_nodes);
	}

	auto data_point_input(RealList &nodes, RealList &values) -> void {
		RealList loaded_nodes  = read_real_list("datas/data_points/nodes.txt");
		RealList loaded_values = read_real_list("datas/data_points/values.txt");
		if (loaded_nodes.size() != loaded_values.size())
			throw std::runtime_error(
				"The nodes and values files must have the same number of entries");

		nodes  = std::move(loaded_nodes);
		values = std::move(loaded_values);
	}
} // namespace ih
