cmake_minimum_required(VERSION 4.0)
project(
	pps
	LANGUAGES CXX
	VERSION 1.0.0
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

set(CMAKE_COLOR_DIAGNOSTICS ON)

include(FetchContent)

find_package(ftxui QUIET)
if(NOT ftxui_FOUND)
	FetchContent_Declare(
		ftxui
		GIT_REPOSITORY https://github.com/ArthurSonzogni/FTXUI.git
		GIT_TAG v7.0.3
	)
	FetchContent_MakeAvailable(ftxui)
endif()

find_package(fmt QUIET)
if(NOT fmt_FOUND)
	FetchContent_Declare(
		fmt
		GIT_REPOSITORY https://github.com/fmtlib/fmt.git
		GIT_TAG 12.2.0
	)
	FetchContent_MakeAvailable(fmt)
endif()

add_executable(
	pps
	main.cpp
	src/math_object.cpp
	src/chapter1/interpolation.cpp
	src/chapter1/horner.cpp
)

target_include_directories(
	pps
	PRIVATE inc/
)

target_link_libraries(
	pps
	PRIVATE ftxui::ftxui
	PRIVATE fmt::fmt
)
