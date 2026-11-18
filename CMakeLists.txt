# CMakeList.txt : Top-level CMake project file,
# Do global configuration and include sub-projects here.
#
cmake_minimum_required (VERSION 3.24)

##########################################################################
#
#   PROJECT
#

# VERSION default "1.0.0":
if(NOT DEFINED APP_VERSION)

  if(NOT DEFINED APP_VERSION_MAJOR)
    set(APP_VERSION_MAJOR 1)
  endif()
  if(NOT DEFINED APP_VERSION_MINOR)
    set(APP_VERSION_MINOR 0)
  endif()
  if(NOT DEFINED APP_VERSION_PATCH)
    set(APP_VERSION_PATCH 0)
  endif()

  set(APP_VERSION "${APP_VERSION_MAJOR}.${APP_VERSION_MINOR}.${APP_VERSION_PATCH}")

endif()

project(
    HYDROGEN_ORBITALS
  VERSION "${APP_VERSION}"
  LANGUAGES C CXX
  DESCRIPTION "Visualization of the hydrogen atom orbitals using the analytical solution of the time independent Schrodinger equation for the hydrogen atom"
)


# C++ language standard to compile for
set(CMAKE_CXX_STANDARD 20)

# To avoid portability issues, language extensions (enabled by default)
# are disabled before any target is created. This allows to stick to the
# language standard such as "-std=c++20" instead of "-std=gnu++20" (which
# is c++20 with GNU extensions):
set(CMAKE_CXX_EXTENSIONS OFF)

# Apps & libraries versioning:
set(APP_VER "${APP_VERSION_MAJOR}${APP_VERSION_MINOR}")
set(APP_VER_CFG "${APP_VER}$<$<CONFIG:Debug>:d>")

##########################################################################
#
#   OPTIONS
#

#########
# GNU standard directories for installation
include(GNUInstallDirs)

#########
# Check whether the compiler supports IPO/LTO (interprocedural optimization):
include(CheckIPOSupported OPTIONAL RESULT_VARIABLE check_ipo_module_found)
if(check_ipo_module_found)
  check_ipo_supported(
    RESULT LTO_OPTIMIZATION_SUPPORTED
    OUTPUT LTO_OPTIMIZATION_ERROR
    LANGUAGES C CXX
  )
else()
  set(LTO_OPTIMIZATION_SUPPORTED NO)
  set(LTO_OPTIMIZATION_ERROR "CheckIPOSupported module not found.")
endif()

#########
# USE_ASAN option:
option(USE_ASAN "Use Address Sanitizer: OFF/ON" OFF)
message(DEBUG "Use Address Sanitizer: ${USE_ASAN}")

#########
# MSVC_RUNTIME option (for compilers targeting the MSVC ABI - other compilers ignore it):
# Address Sanitizer requires a static runtime library, but on Windows it defaults to MSVC's dynamic library.
if(USE_ASAN)
  # Select the static MSVC runtime (for compilers targeting the MSVC ABI):
  set(
    CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>"
    CACHE STRING "MSVC Runtime Library for compilers targeting the MSVC ABI"
    FORCE
  )
endif()

##########################################################################
#
#   TESTS
#
#   Use the following command (in the project's root directory) to
#   run all tests defined in the project:
#
#       ctest --test-dir ./build/ [-V[V]] 
#

# enable_testing()

##########################################################################
#
#   SOURCE PATHS
#

cmake_path(SET ROOT_DIR "${PROJECT_SOURCE_DIR}")

##########################################################################
#
#   TARGETS to BUILD
#

# include CLI11 (command line parser for C++11)
add_subdirectory ("extern/CLI11")

# hydrogen orbitals app
add_subdirectory ("orbitals")
