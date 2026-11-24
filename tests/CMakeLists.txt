#
#   Tests: run them with `ctest --test-dir build` (or `ctest --preset release`)
#

##########################################################################
#
#   C++ unit tests (no framework: test_orbitals <name> runs one test)
#

add_executable(test_orbitals test_orbitals.cpp)
target_link_libraries(test_orbitals PRIVATE orbitals_headers)

foreach(_test IN ITEMS laguerre legendre spherical_harmonics radial_extent mesh_probability
                       expectation_values gregory_weights radial_convergence resolution orbital_names)
  add_test(NAME ${_test} COMMAND test_orbitals ${_test})
endforeach()

##########################################################################
#
#   Command line: invalid options must be rejected (non-zero exit code)
#

set(_bad_options
  "bad_polar_points|-n 3 -l 2 -m 1 --polar-points 20"
  "bad_azimuthal_points|-n 3 -l 2 -m 2 --azimuthal-points 5"
  "bad_radial_step|-n 1 -l 0 -m 0 --radial-step 0"
  "bad_basis|-n 2 -l 1 -m 0 --basis imaginary"
  "bad_l|-n 2 -l 2 -m 0"
)
foreach(_case IN LISTS _bad_options)
  string(REPLACE "|" ";" _parts "${_case}")
  list(GET _parts 0 _name)
  list(GET _parts 1 _args)
  separate_arguments(_args UNIX_COMMAND "${_args}")
  add_test(NAME cli_${_name}
    COMMAND orbitals ${_args} -o "${CMAKE_CURRENT_BINARY_DIR}/data/rejected")
  set_tests_properties(cli_${_name} PROPERTIES WILL_FAIL TRUE)
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
  foreach(_state IN ITEMS "1;0;0;complex" "3;2;-2;complex" "4;3;1;complex" "20;19;19;complex" "20;5;-3;complex"
                          "2;1;-1;real" "3;2;-2;real" "4;3;3;real" "20;19;-19;real" "20;5;-3;real")
    list(GET _state 0 _n)
    list(GET _state 1 _l)
    list(GET _state 2 _m)
    list(GET _state 3 _basis)
    set(_name "n${_n}_l${_l}_m${_m}_${_basis}")
    set(_dir "${CMAKE_CURRENT_BINARY_DIR}/data/${_name}")

    add_test(NAME generate_${_name}
      COMMAND orbitals -n ${_n} -l ${_l} -m ${_m} --basis ${_basis} -o "${_dir}")
    set_tests_properties(generate_${_name} PROPERTIES FIXTURES_SETUP data_${_name})

    add_test(NAME scipy_${_name}
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_against_scipy.py" "${_dir}")
    set_tests_properties(scipy_${_name} PROPERTIES FIXTURES_REQUIRED data_${_name})
  endforeach()
else()
  message(STATUS "Python 3 with NumPy and SciPy not found: the SciPy cross-check tests are skipped")
endif()
