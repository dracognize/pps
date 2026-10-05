#include <ui/table.hpp>

#include <string>
#include <vector>

#include <fmt/format.h>
#include <ftxui/dom/elements.hpp>
#include <ftxui/dom/table.hpp>
#include <ftxui/screen/screen.hpp>

namespace ui {

	auto render_table(const std::vector<std::vector<std::string>> &rows, int min_width) -> void {
		using namespace ftxui;

		if (rows.empty())
			return;

		auto table = Table(rows);

		table.SelectRow(0).Decorate(color(Color::Cyan));
		table.SelectRow(0).DecorateCells(bold);
		table.SelectRow(0).DecorateCells(hcenter);

		table.SelectAll().Border(LIGHT);
		table.SelectAll().Separator(LIGHT);

		const int num_cols = static_cast<int>(rows.front().size());
		for (int col = 0; col < num_cols; ++col)
			table.SelectColumn(col).DecorateCells(hcenter | size(WIDTH, GREATER_THAN, min_width));

		auto document = table.Render();
		auto screen	  = Screen::Create(Dimension::Fit(document, true));

		Render(screen, document);
		screen.Print();
		screen.ResetPosition();
		fmt::println("");
	}

} // namespace ui
