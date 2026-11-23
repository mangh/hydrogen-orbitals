#
#   Original app by Derek Harrison
#

set(target "orbitals")

set(_source_files
  main.cpp
)

##########################################################################
#
#   HEADER-ONLY CODE (shared by the application and the tests)
#

add_library(orbitals_headers INTERFACE)
target_include_directories(orbitals_headers INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

# Warning options set differently for MSVC-like compilers (/W4) and all others (-Wall -Wextra ...):
target_compile_options(orbitals_headers INTERFACE
  "$<IF:$<STREQUAL:${CMAKE_CXX_COMPILER_FRONTEND_VARIANT},MSVC>,/W4,-Wall;-Wextra;-Wpedantic;-Wshadow;-Wconversion>"
)

##########################################################################
#
#   APPLICATION
#

# Add executable target to the project using the specified source files:
add_executable(
  ${target}
  ${_source_files}
)

##########################################################################
#
#   ADDRESS SANITIZER
#

if(USE_ASAN)
  # The following is required to avoid: lld-link: error: /INFERASANLIBS is not allowed in .drectve
  target_compile_definitions(${target} PRIVATE _DISABLE_STRING_ANNOTATION _DISABLE_VECTOR_ANNOTATION)
  target_compile_options(${target} PRIVATE -fsanitize=address)
  target_link_options(${target} PRIVATE -fsanitize=address)
endif()

##########################################################################
#
#   LINK OPTIONS & LIBRARIES
#

target_link_libraries(${target} PRIVATE
  orbitals_headers
  CLI11::CLI11
)

# IPO/LTO (interprocedural optimization), if the compiler supports it:
if(LTO_OPTIMIZATION_SUPPORTED)
  set_property(TARGET ${target} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
  message(WARNING "IPO/LTO is not supported: ${LTO_OPTIMIZATION_ERROR}")
endif()

set_target_properties(${target} PROPERTIES
  LINK_WHAT_YOU_USE ON
)
