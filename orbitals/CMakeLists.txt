#
#   Original app by Derek Harrison
#

set(target "orbitals")

set(_source_files
  main.cpp
)

# Add executable target to the project using the specified source files:
add_executable(
  ${target}
  ${_source_files}
)

##########################################################################
#
#   COMPILE DEFINITIONS & OPTIONS
#

# Warning options set differently for MSVC-like compilers (/W4) and all others (-Wall -Wextra ...):
target_compile_options(${target} PRIVATE
  "$<IF:$<STREQUAL:${CMAKE_CXX_COMPILER_FRONTEND_VARIANT},MSVC>,/W4,-Wall;-Wextra;-Wpedantic;-Wshadow;-Wconversion>"
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
