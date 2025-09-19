# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Compare two vtk data structures with accounting for different ordering of the
points and or cells.
"""

# Import python modules.
import numpy as np
from vtk.util import numpy_support as VN
import vtk


def compare_arrays(
    array_1_vtk,
    array_2_vtk,
    *,
    reorder_1=None,
    reorder_2=None,
    name=None,
    tol_float=None
):
    """
    Compare two vtk arrays.
    """

    array_1 = VN.vtk_to_numpy(array_1_vtk)
    array_2 = VN.vtk_to_numpy(array_2_vtk)

    # In some cases (e.g. vtp files) it is possible that we get empty arrays
    # here. Therefore, if both arrays have no entries we consider them equal
    # here.
    if len(array_1) == 0 and len(array_2) == 0:
        return

    if reorder_1 is None:
        reorder_1 = [i for i in range(len(array_1))]
    if reorder_2 is None:
        reorder_2 = [i for i in range(len(array_2))]

    if len(array_1.shape) == 1:
        array_1_sorted = array_1[reorder_1]
        array_2_sorted = array_2[reorder_2]
    else:
        array_1_sorted = array_1[reorder_1, :]
        array_2_sorted = array_2[reorder_2, :]

    np.testing.assert_allclose(
        array_1_sorted,
        array_2_sorted,
        atol=tol_float,
        rtol=0,
        err_msg="assertion failed for {0}".format(name),
    )


def compare_data_sets(data1, data2, name, **kwargs):
    """
    Compare data sets obtained from vtk objects.
    """

    # Both data sets need to have the same number of arrays.
    if not data1.GetNumberOfArrays() == data2.GetNumberOfArrays():
        data1_names = [data1.GetArrayName(i) for i in range(data1.GetNumberOfArrays())]
        data2_names = [data2.GetArrayName(i) for i in range(data2.GetNumberOfArrays())]
        raise ValueError(
            "Length of vtk data objects for {} does not "
            "match!\nNames in data1: {}\nNames in data2: {}".format(
                name, data1_names, data2_names
            )
        )

    # Compare each array.
    numerrors = 0
    for i in range(data1.GetNumberOfArrays()):

        # Get the arrays with the same name.
        name = data1.GetArrayName(i)
        print("Checking quantity {0} for equality".format(name))
        if name is not None and name.startswith("uid_"):
            # We do not compare the unique ID arrays, as they won't match in
            # general.
            continue
        array1 = data1.GetArray(name)
        array2 = data2.GetArray(name)
        try:
            compare_arrays(array1, array2, name=name, **kwargs)
        except AssertionError as e:
            print(e)
            numerrors += 1

    if numerrors > 0:
        raise AssertionError()


def get_unique_reordering_map(vtk_data):
    """
    Return a reordering map for vtk_data. The data will be ordered according
    unique ID arrays with the names 'uid_0_<name>' then 'uid_1_<name>', ... .
    The unique IDs will be converted to integers, so be very careful when using
    float data as unique IDs.
    """

    # Get all names and indices of unique ID arrays.
    uid_names_to_indices = {}
    for i in range(vtk_data.GetNumberOfArrays()):
        if vtk_data.GetArrayName(i).startswith("uid_"):
            uid_names_to_indices[vtk_data.GetArrayName(i)] = i

    if len(uid_names_to_indices.keys()) == 0:
        # If no UID data is given, return the trivial (i.e. no reordering)
        # reorder map.
        return None, None

    keys_sorted = list(uid_names_to_indices.keys())
    keys_sorted.sort()

    sort_data = []
    for key in keys_sorted:
        np_array_double = VN.vtk_to_numpy(vtk_data.GetArray(uid_names_to_indices[key]))
        np_array_int = np_array_double.astype(np.int32)
        sort_data.append(np_array_int)

    # Sort the UID vales according to the naming.
    sorted_indices = np.lexsort(sort_data)

    n_data = vtk_data.GetNumberOfTuples()
    sorted_indices_reverse = [i for i in range(n_data)]
    for i in range(n_data):
        sorted_indices_reverse[sorted_indices[i]] = i

    return sorted_indices, sorted_indices_reverse


def compare_vtk_data(vtk_1, vtk_2, **kwargs):
    """
    Compare two vtk data structures with accounting for a different ordering of
    the points and or cells.

    If the data to compare comes in different order, e.g. was computed with a
    different number of ranks, it can still be compared here by giving cell
    and or point data arrays with the name 'uid_0_<name>', 'uid_1_<name>' and
    so on. The data will be ordered according to those unique ID arrays (first
    'uid_0_<name>' then 'uid_1_<name>', ...). The unique IDs will be converted
    to integers, so be very careful when using float data as unique IDs.

    Args
    ----
    vtk_1, vtk_2: vtk object
        The vtk objects to compare.
    """

    # Check global data size.
    if not vtk_1.GetNumberOfPoints() == vtk_2.GetNumberOfPoints():
        raise ValueError(
            "Number of points {}=={} not equal".format(
                vtk_1.GetNumberOfPoints(), vtk_2.GetNumberOfPoints()
            )
        )
    if not vtk_1.GetNumberOfCells() == vtk_2.GetNumberOfCells():
        raise ValueError(
            "Number of cells {}=={} not equal".format(
                vtk_1.GetNumberOfCells(), vtk_2.GetNumberOfCells()
            )
        )
    n_points = vtk_1.GetNumberOfPoints()
    n_cells = vtk_1.GetNumberOfCells()

    # Get the UID ordering map for both vtk data values.
    uid_ordering_point_1, uid_ordering_point_1_reverse = get_unique_reordering_map(
        vtk_1.GetPointData()
    )
    uid_ordering_cell_1, _ = get_unique_reordering_map(vtk_1.GetCellData())
    uid_ordering_point_2, uid_ordering_point_2_reverse = get_unique_reordering_map(
        vtk_2.GetPointData()
    )
    uid_ordering_cell_2, _ = get_unique_reordering_map(vtk_2.GetCellData())

    # Check the coordinates of the points.
    compare_arrays(
        vtk_1.GetPoints().GetData(),
        vtk_2.GetPoints().GetData(),
        reorder_1=uid_ordering_point_1,
        reorder_2=uid_ordering_point_2,
        name="point_positions",
        **kwargs
    )

    for i in range(n_cells):

        # Get the i-th cell with or without reordering.
        if uid_ordering_cell_1 is not None:
            cell_1_id = uid_ordering_cell_1[i]
        else:
            cell_1_id = i
        if uid_ordering_cell_2 is not None:
            cell_2_id = uid_ordering_cell_2[i]
        else:
            cell_2_id = i
        cell_1 = vtk_1.GetCell(cell_1_id)
        cell_2 = vtk_2.GetCell(cell_2_id)

        # Compare cell types.
        if not cell_1.GetCellType() == cell_2.GetCellType():
            raise ValueError("Cell types do not match!")

        # Compare connectivity, in the case of reordering we have to do a
        # reverse search here.
        if uid_ordering_point_1_reverse is not None:
            connectivity_1 = [
                uid_ordering_point_1_reverse[cell_1.GetPointId(j)]
                for j in range(cell_1.GetNumberOfPoints())
            ]
        else:
            connectivity_1 = [
                cell_1.GetPointId(j) for j in range(cell_1.GetNumberOfPoints())
            ]
        if uid_ordering_point_2_reverse is not None:
            connectivity_2 = [
                uid_ordering_point_2_reverse[cell_2.GetPointId(j)]
                for j in range(cell_2.GetNumberOfPoints())
            ]
        else:
            connectivity_2 = [
                cell_2.GetPointId(j) for j in range(cell_2.GetNumberOfPoints())
            ]

        if not connectivity_1 == connectivity_2:
            raise ValueError("Wrong connectivity!")

        if cell_1.GetCellType() == 42:
            # In the case of polyhedrons we have to check that the face connectivity is correct.

            def get_polyhedron_id_list_sorted(
                vtk_object, cell_id, uid_ordering_point_reverse
            ):
                """Return the face connectivity list of the polyhedron with the point ids sorted"""

                id_vtk_list = vtk.vtkIdList()
                vtk_object.GetFaceStream(cell_id, id_vtk_list)
                id_list = [
                    id_vtk_list.GetId(i) for i in range(id_vtk_list.GetNumberOfIds())
                ]

                index = 1
                id_list_ordered = [id_list[0]]
                n_faces = id_list[0]
                for i_face in range(n_faces):
                    n_points = id_list[index]
                    id_list_ordered.append(n_points)
                    index += 1
                    for i_point in range(n_points):
                        point_id = id_list[index]
                        if uid_ordering_point_reverse is not None:
                            id_list_ordered.append(uid_ordering_point_reverse[point_id])
                        else:
                            id_list_ordered.append(point_id)
                        index += 1
                return id_list_ordered

            face_id_list_1 = get_polyhedron_id_list_sorted(
                vtk_1, cell_1_id, uid_ordering_point_1_reverse
            )
            face_id_list_2 = get_polyhedron_id_list_sorted(
                vtk_2, cell_2_id, uid_ordering_point_2_reverse
            )
            if not face_id_list_1 == face_id_list_2:
                raise ValueError("Wrong face connectivity!")

    num_errors = 0
    if n_cells > 0:
        try:
            compare_data_sets(
                vtk_1.GetCellData(),
                vtk_2.GetCellData(),
                "cell data",
                reorder_1=uid_ordering_cell_1,
                reorder_2=uid_ordering_cell_2,
                **kwargs
            )
        except AssertionError as e:
            num_errors += 1

    try:
        compare_data_sets(
            vtk_1.GetPointData(),
            vtk_2.GetPointData(),
            "point data",
            reorder_1=uid_ordering_point_1,
            reorder_2=uid_ordering_point_2,
            **kwargs
        )
    except AssertionError as e:
        num_errors += 1

    try:
        compare_data_sets(
            vtk_1.GetFieldData(), vtk_2.GetFieldData(), "field data", **kwargs
        )
    except AssertionError as e:
        num_errors += 1

    if num_errors > 0:
        raise AssertionError()
