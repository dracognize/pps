#pragma once
/// @file output_handling.hpp
/// @brief Shared output helpers.

#include <functional>
#include <string>
#include <vector>

#include <ftxui/dom/table.hpp>

namespace ui {
	/// @brief Render a string grid with the shared FTXUI table style.
	/// @param rows Grid where `rows[0]` is the header row. Use empty strings for blank cells.
	/// @param min_width Minimum cell width applied to each column.
	/// @param transform Optional table decorations applied before rendering.
	auto render_table(
		const std::vector<std::vector<std::string>> &rows,
		int min_width,
		std::function<void(ftxui::Table &)> transform = [](ftxui::Table &) -> void {}) -> void;

	/// @brief Coordinate in an FTXUI table, expressed as column and row.
	struct Point {
			int x;
			int y;
	};

	/// @brief Prints a long separator before a display operation.
	struct DisplayFrame {
		DisplayFrame();
		~DisplayFrame();
	};
} // namespace ui
