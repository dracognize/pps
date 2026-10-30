#include <io/output_handling.hpp>

#include <ftxui/dom/elements.hpp>
#include <ftxui/screen/screen.hpp>

#include <fmt/core.h>
#include <fmt/format.h>

#include <string>

#ifdef _WIN32
#include <windows.h>
#endif

namespace ui {
	DisplayFrame::DisplayFrame() {
		fmt::print("\n{}\n", std::string(80, '='));
	}
	DisplayFrame::~DisplayFrame() {
		fmt::print("\n{}\n\n", std::string(80, '='));
	}

	auto render_table(const std::vector<std::vector<std::string>> &rows,
					  int min_width,
					  std::function<void(ftxui::Table &)> transform) -> void {
		using namespace ftxui;
		if (rows.empty())
			return;

#ifdef _WIN32
		HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
		DWORD original_mode;
		GetConsoleMode(hOut, &original_mode);
#endif
		auto table = Table(rows);
		table.SelectRow(0).Decorate(color(Color::Cyan));
		table.SelectRow(0).DecorateCells(bold);
		table.SelectRow(0).DecorateCells(hcenter);
		table.SelectAll().Border(LIGHT);
		table.SelectAll().Separator(LIGHT);
		for (int col = 0; col < static_cast<int>(rows.front().size()); ++col)
			table.SelectColumn(col).DecorateCells(hcenter | size(WIDTH, GREATER_THAN, min_width));
		transform(table);

		auto document = table.Render();
		auto screen = Screen::Create(Dimension::Fit(document, true));
		Render(screen, document);
		screen.Print();
#ifdef _WIN32
		SetConsoleMode(hOut, original_mode);
#endif
		fmt::println("");
	}
} // namespace ui
