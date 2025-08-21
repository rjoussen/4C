# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# Declare an option with given name, description and default value (ON or OFF).
# This function is almost equivalent to the builtin CMake option() command,
# except that it also prints info on whether the option is set.
# Options need to start with "FOUR_C_".
function(four_c_process_global_option option_name description default)
  if(NOT option_name MATCHES "FOUR_C_.*")
    message(FATAL_ERROR "Disallowed option '${option_name}'. Option needs to start with 'FOUR_C_'.")
  endif()
  option("${option_name}" "${description} (default: ${default})" "${default}")
  if(${option_name})
    message(STATUS "Option ${option_name} = ON")
  else()
    message(STATUS "Option ${option_name} = OFF")
  endif()
endfunction()

# Initialize a cache variable. This function is almost equivalent to the builtin CMake
# set() command, except that it also prints info on whether the variable is set and allows to
# ensure a consistent style. Note that you need to use four_c_process_global_option() for the
# common case of a BOOL variable.
#
# Usage:
#   four_c_initialize_cache_variable(variable_name
#     TYPE <STRING|FILEPATH|PATH>
#     DESCRIPTION "Description of the variable"
#     DEFAULT "Default value of the variable")
#
function(four_c_process_cache_variable variable_name)
  if(NOT variable_name MATCHES "FOUR_C_.*")
    message(
      FATAL_ERROR
        "Disallowed variable name '${variable_name}'. Variable name needs to start with 'FOUR_C_'."
      )
  endif()

  set(options "")
  set(oneValueArgs TYPE DESCRIPTION DEFAULT)
  set(multiValueArgs "")
  cmake_parse_arguments(
    _parsed
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )

  if(DEFINED _parsed_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "There are unparsed arguments: ${_parsed_UNPARSED_ARGUMENTS}")
  endif()

  # Ensure that only one of the options STRING or (FILE)PATH for the cache variable type is set.
  if(NOT _parsed_TYPE MATCHES "STRING|FILEPATH|PATH")
    message(
      FATAL_ERROR
        "Invalid type '${_parsed_TYPE}' for cache variable '${variable_name}'. Allowed types are STRING, FILEPATH and PATH."
      )
  endif()

  set(${variable_name}
      "${_parsed_DEFAULT}"
      CACHE ${_parsed_TYPE} "${_parsed_DESCRIPTION} (default: ${_parsed_DEFAULT})"
      )
  if(${variable_name} STREQUAL "")
    message(STATUS "Cache variable ${variable_name} not set")
  else()
    message(STATUS "Cache variable ${variable_name} = ${${variable_name}}")
  endif()
endfunction()
