#
#   Original app by Derek Harrison
#

set(target "ultimate")

set(_source_files
  main.cpp
)

# # Set explicit dependency of source files on METROLOGY headers
# # (for simplicity, it is assumed here that all source files depend on those headers):
# set_property(
#   SOURCE ${_source_files}
#   PROPERTY OBJECT_DEPENDS "${METROLOGY_CHANGE_TIP}"
# )

# Add executable target to the project using the specified source files:
add_executable(
  ${target}
  ${_source_files}
)

##########################################################################
#
#   COMPILE DEFINITIONS & OPTIONS
#

target_compile_features(${target} PRIVATE cxx_std_17)
# Warning options set differently for MSVC-like compiler (/W3) and all others (-Wall -Wextra):
target_compile_options(${target} PRIVATE $<IF:$<STREQUAL:${CMAKE_CXX_COMPILER_FRONTEND_VARIANT},MSVC>,/W3,-Wall -Wextra>)

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

# IPO/LTO optimization makes the application just as efficient
# whether it uses units of measure or simple numbers.
if(LTO_OPTIMIZATION_SUPPORTED)
  set_property(TARGET ${target} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
  message(WARNING "IPO/LTO is not supported: ${LTO_OPTIMIZATION_ERROR}")
endif()

# # target_link_libraries() command below imports the so called Usage Requirements of
# # the METROLOGY_LIBRARY, which include (among other things) its include directories
# # (so there is no need to use the target_include_directories() command):
# target_link_libraries(${target}
#   PRIVATE METROLOGY_LIBRARY
# )

set_target_properties(${target} PROPERTIES
  LINK_WHAT_YOU_USE ON
)

# set_target_properties(${target} PROPERTIES
#   OUTPUT_NAME "${target}${APP_VER_CFG}"
# )
