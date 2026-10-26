#pragma once

#include <math_object.hpp>

#include <string>

namespace ih {
	auto input_text(std::string rel_file_dir, RealList &keys, RealList &values) -> void;
}