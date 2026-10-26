#pragma once
/// @file table.hpp
/// @brief Shared FTXUI table renderer.

#include <functional>
#include <string>
#include <vector>

#include <ftxui/dom/table.hpp>

namespace ui {

	/// @brief Render a string grid with the repo-wide table style.
	/// @param rows Grid where `rows[0]` is the header row. Use `""` for blank cells.
	/// @param min_width Minimum cell width applied to every column.
	auto render_table(
		const std::vector<std::vector<std::string>> &rows,
		int											 min_width,
		std::function<void(ftxui::Table &)> transform = [](ftxui::Table &table) -> void { return; })
		-> void;

	struct Point {
			int x;
			int y;
	};

} // namespace ui
