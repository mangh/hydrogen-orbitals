#
#   Tests: run them with `ctest --test-dir build` (or `ctest --preset release`)
#

##########################################################################
#
#   C++ unit tests (no framework: test_orbitals <name> runs one test)
#

add_executable(test_orbitals test_orbitals.cpp)
target_link_libraries(test_orbitals PRIVATE orbitals_headers)

foreach(_test IN ITEMS laguerre legendre spherical_harmonics radial_extent mesh_probability)
  add_test(NAME ${_test} COMMAND test_orbitals ${_test})
endforeach()

##########################################################################
#
#   Cross-check against SciPy (only when Python with NumPy and SciPy is available):
#   orbitals.exe writes a few states, check_against_scipy.py compares psi with an independent formula.
#

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_FOUND)
  execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import numpy, scipy"
    RESULT_VARIABLE _scipy_missing
    OUTPUT_QUIET ERROR_QUIET
  )
endif()

if(Python3_FOUND AND NOT _scipy_missing)
  foreach(_state IN ITEMS "1;0;0" "3;2;-2" "4;3;1" "20;19;19" "20;5;-3")
    list(GET _state 0 _n)
    list(GET _state 1 _l)
    list(GET _state 2 _m)
    set(_name "n${_n}_l${_l}_m${_m}")
    set(_dir "${CMAKE_CURRENT_BINARY_DIR}/data/${_name}")

    add_test(NAME generate_${_name}
      COMMAND orbitals -n ${_n} -l ${_l} -m ${_m} -o "${_dir}")
    set_tests_properties(generate_${_name} PROPERTIES FIXTURES_SETUP data_${_name})

    add_test(NAME scipy_${_name}
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_against_scipy.py" "${_dir}")
    set_tests_properties(scipy_${_name} PROPERTIES FIXTURES_REQUIRED data_${_name})
  endforeach()
else()
  message(STATUS "Python 3 with NumPy and SciPy not found: the SciPy cross-check tests are skipped")
endif()
