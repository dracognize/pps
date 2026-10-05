#pragma once
/// @file table.hpp
/// @brief Shared FTXUI table renderer.

#include <string>
#include <vector>

namespace ui {

	/// @brief Render a string grid with the repo-wide table style.
	/// @param rows Grid where `rows[0]` is the header row. Use `""` for blank cells.
	/// @param min_width Minimum cell width applied to every column.
	auto render_table(const std::vector<std::vector<std::string>> &rows, int min_width) -> void;

} // namespace ui
