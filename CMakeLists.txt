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

# Needed: pps_lib links into the Python shared module.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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

find_package(Eigen3 3.4 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
	set(EIGEN_BUILD_TESTING    OFF CACHE BOOL "" FORCE)
	set(EIGEN_BUILD_DOC        OFF CACHE BOOL "" FORCE)
	set(EIGEN_BUILD_PKGCONFIG  OFF CACHE BOOL "" FORCE)
	set(BUILD_TESTING          OFF CACHE BOOL "" FORCE)
	FetchContent_Declare(
		eigen
		GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
		GIT_TAG 3.4.0
		GIT_SHALLOW TRUE
	)
	FetchContent_MakeAvailable(eigen)
	if(NOT TARGET Eigen3::Eigen)
		add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
		target_include_directories(Eigen3::Eigen INTERFACE
			"${eigen_SOURCE_DIR}"
		)
	endif()
endif()

# Core math library: shared by the CLI and the Python module.
add_library(
	pps_lib
	STATIC
	src/math_object.cpp
	src/math_utils.cpp
	src/io/input_handling.cpp
	src/io/output_handling.cpp
	src/chapter1/interpolation/optimal.cpp
	src/chapter1/interpolation/horner.cpp
	src/chapter1/interpolation/polynomial/Lagrange.cpp
	src/chapter1/interpolation/polynomial/Newton.cpp
	src/chapter1/interpolation/polynomial/Gauss.cpp
	src/chapter1/interpolation/polynomial/Stirling.cpp
	src/chapter1/interpolation/polynomial/Bessel.cpp
	src/chapter1/interpolation/interval.cpp
	src/chapter1/interpolation/inverse.cpp
	src/chapter1/spline/linear.cpp
	src/chapter1/spline/quadratic.cpp
	src/chapter1/spline/cubic.cpp
	src/chapter1/lsp.cpp
)

target_include_directories(
	pps_lib
	PUBLIC inc/
)

target_link_libraries(
	pps_lib
	PUBLIC ftxui::ftxui
	PUBLIC fmt::fmt
	PUBLIC Eigen3::Eigen
)

add_executable(
	pps
	main.cpp
)

target_link_libraries(
	pps
	PRIVATE pps_lib
)

# Python module (import pps). Prefers a system-installed pybind11,
# falls back to FetchContent only if none is found.
find_package(Python COMPONENTS Interpreter Development QUIET)
if(Python_FOUND)
	find_package(pybind11 QUIET)
	if(NOT pybind11_FOUND)
		FetchContent_Declare(
			pybind11
			GIT_REPOSITORY https://github.com/pybind/pybind11.git
			GIT_TAG v3.0.1
		)
		FetchContent_MakeAvailable(pybind11)
	endif()

	pybind11_add_module(
		pps_py
		python/bindings.cpp
	)

	target_link_libraries(
		pps_py
		PRIVATE pps_lib
	)

	set_target_properties(
		pps_py
		PROPERTIES OUTPUT_NAME pps
	)
else()
	message(STATUS "Python Development not found, skipping pps Python module")
endif()
