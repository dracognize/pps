# Tables

Defined in `inc/ui/table.hpp`.

## render_table

```cpp
void ui::render_table(
    const std::vector<std::vector<std::string>> &rows,
    int min_width);
```

Renders a pre-built string grid with the repo-wide FTXUI style
(cyan header row, light borders, centered cells). `rows[0]` is the
header; use `""` for intentionally blank (triangular) cells.
Every chapter funnel its tables through here so style changes
happen in one place.
