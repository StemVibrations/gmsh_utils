from gmsh_utils.gmsh_IO import GmshIO
from utils import TestUtils

import numpy as np

class TestGmshIO:
    """
    Tests for the GmshIO class.
    """

    def test_generate_mesh_2D(self):
        """
        Checks whether mesh data generated for 2D geometries is not empty.

        """
        # define the points of the surface as a list of tuples
        input_points = [(0, 0, 0), (1, 0, 0), (1, 3, 0), (0, 3, 0), (-1, 1.5, 0)]
        # define the mesh size
        mesh_size = 0.1
        # set a name label for the surface
        name_label = "Soil Layer"
        # if True, saves mesh data to separate mdpa files
        save_file = False
        # if True, opens gmsh interface
        gmsh_interface = False
        # set a name for mesh output file
        mesh_output_name = "test_2D"
        # set output directory of the mesh
        mesh_output_dir = "."

        # test 2D geometry
        # define geometry dimension; input "2" for 2D
        dims = 2
        # input depth of geometry if 3D
        extrusion_length = [0, 0, 0]

        gmsh_io = GmshIO()
        gmsh_io.generate_gmsh_mesh(input_points, extrusion_length, mesh_size, dims, name_label,
                                   mesh_output_name, mesh_output_dir, save_file,
                                   gmsh_interface)

        mesh_data = gmsh_io.mesh_data

        assert mesh_data["nodes"]["coordinates"].size > 0  # check if node_coords is not empty
        assert mesh_data["nodes"]["ids"].size > 0  # check if node_tags is not empty
        assert list(mesh_data["elements"].keys()) == ["LINE_2N", "TRIANGLE_3N",
                                                      "POINT_1N"]  # check if correct elements are present

        # check each element type contains ids and nodes
        for value in mesh_data["elements"].values():
            assert value["element_ids"].size > 0
            assert value["connectivities"].size > 0

    def test_generate_mesh_3D(self):
        """
        Checks whether mesh data generated for 3D geometries is not empty.

        """

        # define the points of the surface as a list of tuples
        input_points = [(0, 0, 0), (1, 0, 0), (1, 3, 0), (0, 3, 0), (-1, 1.5, 0)]
        # define the mesh size
        element_size = 0.1
        # set a name label for the surface
        name_label = "Soil Layer"
        # if True, saves mesh data to separate mdpa files
        save_file = False
        # if True, opens gmsh interface
        gmsh_interface = False
        # test 3D geometry
        # define geometry dimension; input "3" for 3D to extrude the 2D surface
        dims = 3
        # input depth of geometry if 3D
        extrusion_length = [0, 0, 1]
        # set a name for mesh output file
        mesh_output_name = "test_3D"
        # set output directory of the mesh
        mesh_output_dir = "."

        gmsh_io = GmshIO()

        gmsh_io.generate_gmsh_mesh(input_points, extrusion_length, element_size, dims, name_label,
                                   mesh_output_name, mesh_output_dir, save_file,
                                   gmsh_interface)

        mesh_data = gmsh_io.mesh_data

        assert mesh_data["nodes"]["coordinates"].size > 0  # check if node_coords is not empty
        assert mesh_data["nodes"]["ids"].size > 0  # check if node_tags is not empty
        assert list(mesh_data["elements"].keys()) == ["LINE_2N", "TRIANGLE_3N", "TETRAHEDRON_4N",
                                                      "POINT_1N"]  # check if correct elements are present

        # check each element type contains ids and nodes
        for value in mesh_data["elements"].values():
            assert value["element_ids"].size > 0
            assert value["connectivities"].size > 0

    def test_read_gmsh_geo_2D(self):
        """
        Checks whether a gmsh .geo file is read correctly. For a 2d geometry
        """
        geo_file = r"test_data/column_2D.geo"

        gmsh_io = GmshIO()
        gmsh_io.read_gmsh_geo(geo_file)

        geo_data = gmsh_io.geo_data

        # check if the coordinates of the points are correct
        expected_points = {1: [0, 0, 0], 2: [0.5, 0, 0], 3: [0.5, 1, 0], 4: [0, 1, 0], 11: [0, 2, 0], 12: [0.5, 2.0, 0]}
        expected_lines = {5: [1, 2], 6: [2, 3], 7: [3, 4], 8: [1, 4], 13: [4, 11], 14: [11, 12], 15: [3, 12]}
        expected_surfaces = {10: [5, 6, 7, 8], 17: [-7, -13, -14, -15]}  # negative sign means reversed orientation

        expected_physical_groups = {'group_1': {'geometry_id': 10, 'id': 1, 'ndim': 2},
                                    'group_2': {'geometry_id': 17, 'id': 2, 'ndim': 2}}

        expected_geo_data = {"points": expected_points,
                             "lines": expected_lines,
                             "surfaces": expected_surfaces,
                             "volumes": {},
                             "physical_groups": expected_physical_groups}

        # check if expected and actual geo data are equal
        TestUtils.assert_dictionary_almost_equal(expected_geo_data, geo_data)

    def test_read_gmsh_geo_3D(self):
        """
        Checks whether a gmsh .geo file is read correctly. For a 3d geometry
        """
        geo_file = r"test_data/column_3D_tetra4.geo"

        gmsh_io = GmshIO()
        gmsh_io.read_gmsh_geo(geo_file)

        geo_data = gmsh_io.geo_data
        expected_points = {1: [0., 0., 0.], 2: [0.5, 0., 0.], 3: [0.5, 1., 0.], 4: [0., 1., 0.], 11: [0., 2., 0.],
                           12: [0.5, 2., 0.], 13: [0., 0., -0.5], 14: [0.5, 0., -0.5], 18: [0.5, 1., -0.5],
                           22: [0., 1., -0.5], 23: [0., 2., -0.5], 32: [0.5, 2., -0.5]}
        expected_lines = {5: [1, 2], 6: [2, 3], 7: [3, 4], 8: [1, 4], 13: [4, 11], 14: [11, 12], 15: [3, 12],
                          19: [13, 14], 20: [14, 18], 21: [18, 22], 22: [13, 22], 24: [1, 13], 25: [2, 14], 29: [3, 18],
                          33: [4, 22], 41: [22, 23], 43: [18, 32], 44: [23, 32], 46: [11, 23], 55: [12, 32]}
        expected_surfaces = {10: [5, 6, 7, 8], 17: [-7, -13, -14, -15], 26: [5, -19, -24, 25], 30: [6, -20, -25, 29],
                             34: [7, -21, -29, 33], 38: [8, -22, 24, -33], 39: [19, 20, 21, 22],
                             48: [-13, 33, -41, -46], 56: [-15, -29, -43, 55], 60: [-14, -44, 46, -55],
                             61: [-21, 41, 43, 44]}
        expected_volumes = {1: [-10, 26, 30, 34, 38, 39], 2: [-17, -34, -48, -56, -60, 61]}
        expected_physical_groups = {'group_1': {'geometry_id': 1, 'id': 1, 'ndim': 3},
                                    'group_2': {'geometry_id': 2, 'id': 2, 'ndim': 3}}

        expected_geo_data = {"points": expected_points,
                             "lines": expected_lines,
                             "surfaces": expected_surfaces,
                             "volumes": expected_volumes,
                             "physical_groups": expected_physical_groups}

        # check if expected and actual geo data are equal
        TestUtils.assert_dictionary_almost_equal(expected_geo_data, geo_data)

    def test_read_gmsh_msh_2D(self):
        """
        Checks whether a gmsh .msh file is read correctly in 2D.

        """

        msh_file = r"test_data/block_2D.msh"

        gmsh_io = GmshIO()
        gmsh_io.read_gmsh_msh(msh_file)

        mesh_data = gmsh_io.mesh_data

        expected_mesh_data = {'elements': {'TRIANGLE_3N': {'connectivities': [[1, 2, 4], [4, 2, 3]],
                                                           'element_ids': [1, 2]}},
                              'nodes': {'coordinates': [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]],
                                        'ids': [1, 2, 3, 4]}}
        # check if the coordinates of the points are correct
        TestUtils.assert_dictionary_almost_equal(expected_mesh_data, mesh_data)

    def test_read_gmsh_msh_3D(self):
        """
        Checks whether a gmsh .msh file is read correctly in 3D.

        """

        msh_file = r"test_data/block_3D.msh"

        gmsh_io = GmshIO()

        # read mesh file
        gmsh_io.read_gmsh_msh(msh_file)

        # Get mesh data
        mesh_data = gmsh_io.mesh_data

        # set expected mesh data
        expected_mesh_data = {'nodes': {'coordinates': np.array([[0., 0., 0.],
                                                                 [1., 0., 0.],
                                                                 [1., 1., 0.],
                                                                 [0., 1., 0.],
                                                                 [0., 0., -1.],
                                                                 [1., 0., -1.],
                                                                 [1., 1., -1.],
                                                                 [0., 1., -1.]]),
                                        'ids': [1, 2, 3, 4, 5, 6, 7, 8]},
                             'elements': {'TETRAHEDRON_4N': {'element_ids': [1, 2, 3, 4, 5, 6],
                                                             'connectivities': [[2, 1, 4, 8],
                                                                                [5, 6, 8, 2],
                                                                                [5, 2, 8, 1],
                                                                                [2, 4, 3, 7],
                                                                                [8, 6, 7, 2],
                                                                                [8, 2, 7, 4]]}}}


        # check if the coordinates of the points are correct
        TestUtils.assert_dictionary_almost_equal(expected_mesh_data, mesh_data)
