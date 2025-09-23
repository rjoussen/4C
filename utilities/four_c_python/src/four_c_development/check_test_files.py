# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Check input test files for errors."""

import os
import argparse
import re
from four_c_common_utils import common_utils as utils


def check_inputtests(filenames: list) -> list:
    """Check if all input test files are listed in
    tests/list_of_tests.cmake.

    Args:
        filenames: List of input test files to check.

    Returns:
        List of error messages if any input test files are missing
        in tests/list_of_tests.cmake.
    """

    with open("tests/list_of_tests.cmake", "r") as cmakefile:
        list_of_test_lines = "\n".join(cmakefile.readlines())

    errors = []

    for input_test in filenames:
        # check, whether this input file is in tests/list_of_tests.cmake
        expected_test_name = os.path.splitext(os.path.basename(input_test))[0]
        if (
            re.search(r"\b" + re.escape(expected_test_name) + r"\b", list_of_test_lines)
            is None
        ):
            errors.append(input_test)

    if errors:
        errors = [
            "The following input files are missing in tests/list_of_tests.cmake:",
            "",
        ] + errors

    return errors


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs="*")
    parser.add_argument(
        "--out",
        type=str,
        default=None,
        help="Add this tag if the error message should be written to a file.",
    )
    args = parser.parse_args()

    # error file (None for sys.stderr)
    errfile = args.out

    # check input file tests
    errors = check_inputtests(args.filenames)

    utils.pretty_print_error_report("", errors, errfile)

    if errors:
        raise RuntimeError("\n".join(errors))


if __name__ == "__main__":

    main()
