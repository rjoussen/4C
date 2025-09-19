# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""
Check that our header files do not have include cycles.

This script was inspired by and adapted from the one in deal.II
https://github.com/dealii/dealii/blob/master/contrib/utilities/detect_include_cycles.py
"""

from glob import glob
import networkx as nx
import re

match_includes = re.compile(r"# *include *\"(4C_.*.hpp)\"")


def add_includes_for_file(header_file_path, graph):
    """
    Add the includes of a header file to the file's node as a directed edge.
    """
    with open(header_file_path) as f:
        lines = f.readlines()

    header_file_name = header_file_path.split("/")[-1]

    for line in lines:
        m = match_includes.match(line)
        if m:
            included_file = m.group(1)
            graph.add_edge(header_file_name, included_file)


def main():
    # Create a list of all header files in the include folder.
    filelist = glob("src/**/*.hpp", recursive=True)
    assert filelist, "Please call the script from the top-level directory."

    # For each header file, add the includes as the edges of a directed graph.
    graph = nx.DiGraph()
    for header_file_name in filelist:
        add_includes_for_file(header_file_name, graph)

    # Then figure out whether there are cycles and if so, print them:
    cycles = nx.simple_cycles(graph)
    cycles_as_list = list(cycles)
    if len(cycles_as_list) > 0:
        print(f"Cycles in the include graph detected! Please fix them manually.")
        for cycle in cycles_as_list:
            print(cycle)
        exit(1)


if __name__ == "__main__":
    main()
