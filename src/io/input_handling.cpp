#include <io/input_handling.hpp>

#include <fstream>

namespace ih {
	auto input_text(std::string rel_file_dir, RealList &keys, RealList &values) -> void {
		std::ifstream fin{rel_file_dir};
		double		  temp_key, temp_value;

		while (fin >> temp_key >> temp_value) {
			keys.push_back(temp_key);
			values.push_back(temp_value);
		}
	}
} // namespace ih
