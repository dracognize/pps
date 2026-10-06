#include <ui/table.hpp>

#include <string>
#include <vector>

#include <ftxui/dom/elements.hpp>
#include <ftxui/dom/table.hpp>
#include <ftxui/screen/screen.hpp>

#include <fmt/core.h>

#ifdef _WIN32
#include <windows.h>
#endif


namespace ui {

	auto render_table(const std::vector<std::vector<std::string>> &rows, int min_width) -> void {
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

		const int num_cols = static_cast<int>(rows.front().size());
		for (int col = 0; col < num_cols; ++col)
			table.SelectColumn(col).DecorateCells(hcenter | size(WIDTH, GREATER_THAN, min_width));

		auto document = table.Render();
		auto screen	  = Screen::Create(Dimension::Fit(document, true));

		Render(screen, document);
		screen.Print();
		// NOTE: `ResetPosition()` only *returns* the cursor-up sequence for
		// animation loops — it prints nothing, so calling it here is a no-op.
		// End the one-shot table on the same stream instead (upstream pattern
		// from FTXUI's `examples/dom/table.cpp`): on Windows consoles (VT
		// mode, which FTXUI enables) a bare `\n` moves down without returning
		// the carriage, so anything printed after `Print()` would start
		// mid-line and misalign every following row.
		#ifdef _WIN32
		SetConsoleMode(hOut, original_mode);
		#endif
		fmt::println("");
	}

} // namespace ui
