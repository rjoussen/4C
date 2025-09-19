# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Check header guards in C/C++ files."""

import argparse
import os
from four_c_common_utils import common_utils as utils


def check_header_guards(filenames, allerrors):
    num_wrong_header_guards = 0
    for file_path in filenames:
        file_name = os.path.splitext(os.path.basename(file_path))[0]
        define_name = file_name.replace("-", "_").replace(".", "_").upper() + "_HPP"
        define_name = define_name.replace("4C", "FOUR_C")

        template_header_if = "#ifndef " + define_name + "\n"
        template_header_define = "#define " + define_name + "\n"

        def has_wrong_header_guards():
            all_lines = utils.file_contents(file_path)
            for line_num, line in enumerate(all_lines):
                if (
                    line == template_header_if
                    and all_lines[line_num + 1] == template_header_define
                ):
                    return False
            return True

        if has_wrong_header_guards():
            allerrors.append("Wrong header guard in " + file_path)
            num_wrong_header_guards += 1

    return num_wrong_header_guards


def main():
    # build command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "filenames", nargs="*", help="List of files to be checked for header guards."
    )

    parser.add_argument(
        "--out",
        type=str,
        default=None,
        help="Add this tag if the error message should be written to a file.",
    )
    args = parser.parse_args()

    # error file (None for sys.stderr)
    errfile = args.out
    errors = 0
    allerrors = []
    try:
        errors += check_header_guards(args.filenames, allerrors)
    except ValueError:
        print("Something went wrong! Check the error functions in this script again!")
        errors += 1

    utils.pretty_print_error_report(
        "Wrong header guards in the following files:", allerrors, errfile
    )
    return errors


if __name__ == "__main__":
    import sys

    sys.exit(main())
