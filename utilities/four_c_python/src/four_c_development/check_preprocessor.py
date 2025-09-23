# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import sys
import re

pattern = re.compile(r"#if\s+[0-9]+")


def check_if_directives(file_path):
    errors = 0
    with open(file_path, "r") as file:
        lines = file.readlines()
        for line_num, line in enumerate(lines, 1):
            if pattern.search(line):
                print(f"{file_path}:{line_num}: Found '{line.strip()}'")
                errors += 1
    return errors


def main():
    files = sys.argv[1:]
    errors = 0
    for file in files:
        errors += check_if_directives(file)
    sys.exit(errors)


if __name__ == "__main__":
    main()
