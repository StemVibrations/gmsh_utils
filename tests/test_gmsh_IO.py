from gmsh_utils.gmsh_IO import GmshIO
from utils import TestUtils

import gmsh
import numpy as np
import pytest


class TestGmshIO:
    """
    Tests for the GmshIO class.
    """

    @pytest.fixture(autouse=True)
    def close_gmsh(self):
        """
        Closes the gmsh interface after each test.
        """
        if gmsh.is_initialized():
            gmsh.finalize()

    @pytest.fixture
    def expected_geo_data_3D(self):
        """
        Expected geometry data for a 3D geometry. The geometry is 2 stacked blocks, where the top and bottom blocks
        are in different groups.
        """
        expected_points = {1: [0., 0., 0.], 2: [0.5, 0., 0.], 3: [0.5, 1., 0.], 4: [0., 1., 0.], 11: [0., 2., 0.],
                           12: [0.5, 2., 0.], 13: [0., 0., -0.5], 14: [0.5, 0., -0.5], 18: [0.5, 1., -0.5],
                           22: [0., 1., -0.5], 23: [0., 2., -0.5], 32: [0.5, 2., -0.5]}
        expected_lines = {5: [1, 2], 6: [2, 3], 7: [3, 4], 8: [4, 1], 13: [4, 11], 14: [11, 12], 15: [12, 3],
                          19: [13, 14], 20: [14, 18], 21: [18, 22], 22: [22, 13], 24: [1, 13], 25: [2, 14],
                          29: [3, 18], 33: [4, 22], 41: [23, 22], 43: [18, 32], 44: [32, 23], 46: [11, 23],
                          55: [12, 32]}
        expected_surfaces = {10: [5, 6, 7, 8], 17: [-13, -7, -15, -14], 26: [5, 25, -19, -24], 30: [6, 29, -20, -25],
                             34: [7, 33, -21, -29], 38: [8, 24, -22, -33], 39: [19, 20, 21, 22],
                             48: [-13, 33, -41, -46], 56: [-15, 55, -43, -29], 60: [-14, 46, -44, -55],
                             61: [41, -21, 43, 44]}
        expected_volumes = {1: [-10, 39, 26, 30, 34, 38], 2: [-17, 61, -48, -34, -56, -60]}
        expected_physical_groups = {'group_1': {'geometry_ids': [1], 'id': 1, 'ndim': 3},
                                    'group_2': {'geometry_ids': [2], 'id': 2, 'ndim': 3}}

        return {"points": expected_points,
                "lines": expected_lines,
                "surfaces": expected_surfaces,
                "volumes": expected_volumes,
                "physical_groups": expected_physical_groups}

    @pytest.fixture
    def expected_geo_data_3D_with_shared_group(self):
        """
        Expected geometry data for a 3D geometry. The geometry is 2 stacked blocks, where the top and bottom blocks
        are in different groups.
        """
        expected_points = {1: [0., 0., 0.], 2: [0.5, 0., 0.], 3: [0.5, 1., 0.], 4: [0., 1., 0.], 11: [0., 2., 0.],
                           12: [0.5, 2., 0.], 13: [0., 0., -0.5], 14: [0.5, 0., -0.5], 18: [0.5, 1., -0.5],
                           22: [0., 1., -0.5], 23: [0., 2., -0.5], 32: [0.5, 2., -0.5]}
        expected_lines = {5: [1, 2], 6: [2, 3], 7: [3, 4], 8: [4, 1], 13: [4, 11], 14: [11, 12], 15: [12, 3],
                          19: [13, 14], 20: [14, 18], 21: [18, 22], 22: [22, 13], 24: [1, 13], 25: [2, 14],
                          29: [3, 18], 33: [4, 22], 41: [23, 22], 43: [18, 32], 44: [32, 23], 46: [11, 23],
                          55: [12, 32]}
        expected_surfaces = {10: [5, 6, 7, 8], 17: [-13, -7, -15, -14], 26: [5, 25, -19, -24], 30: [6, 29, -20, -25],
                             34: [7, 33, -21, -29], 38: [8, 24, -22, -33], 39: [19, 20, 21, 22],
                             48: [-13, 33, -41, -46], 56: [-15, 55, -43, -29], 60: [-14, 46, -44, -55],
                             61: [41, -21, 43, 44]}
        expected_volumes = {1: [-10, 39, 26, 30, 34, 38], 2: [-17, 61, -48, -34, -56, -60]}
        expected_physical_groups = {'group_1': {'geometry_ids': [1], 'id': 1, 'ndim': 3},
                                    'group_2': {'geometry_ids': [2], 'id': 2, 'ndim': 3},
                                    'gravity': {'geometry_ids': [1, 2], 'id': 3, 'ndim': 3}}

        return {"points": expected_points,
                "lines": expected_lines,
                "surfaces": expected_surfaces,
                "volumes": expected_volumes,
                "physical_groups": expected_physical_groups}

    def test_generate_mesh_2D(self):
        """
        Checks whether mesh data generated for 2D geometries is not empty.

        """
        # if "True", mesh size is defined by the user in the dictionary below; if "False" default_mesh_size is used
        arbitrary_mesh_size: bool = False
        # if "arbitrary_mesh_size = False", define a mesh size, if set to -1 mesh size is logically chosen by Gmsh
        # itself based on the geometry
        default_mesh_size: float = -1
        # define the points of the surface as a list of tuples
        input_points_list = [[(0, 0, 0), (3, 0, 0), (3, 1, 0), (0, 1, 0)],
                             [(3, 0, 0), (5, 0, 0), (5, 1, 0), (4, 1.5, 0), (3, 1, 0)],
                             [(0, 1, 0), (2, 1, 0), (2, 3, 0), (0, 3, 0)],
                             [(2, 1, 0), (3, 1, 0), (4, 1.5, 0), (5, 1, 0), (5, 3, 0), (2, 3, 0)],
                             [(0, 3, 0), (2.5, 3, 0), (2, 4, 0), (0, 4, 0)],
                             [(0.8, 4, 0), (1.2, 4, 0), (1.2, 4.1, 0), (0.8, 4.1, 0)]]
        # define the name labels for the surfaces
        name_label_list = ["First Soil Layer", "FSL", "Second Soil Layer", "SSL", "Embankment", "Soil Ballast"]
        # define the mesh size for each surface
        mesh_size_list = []
        # define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
        dims = 2
        # if 3D, input depth of geometry to be extruded from 2D surface
        extrusion_length = [0, 0, 0]
        # if "True", saves mesh data to separate mdpa files; otherwise "False"
        save_file = False
        # if "True", opens gmsh interface; otherwise "False"
        open_gmsh_gui = False
        # set a name for mesh output file
        mesh_output_name = "test_2D"
        # set output directory
        mesh_output_dir = "."

        gmsh_io = GmshIO()

        gmsh_io.generate_geometry(input_points_list, extrusion_length, dims, mesh_output_name,
                                  name_label_list, default_mesh_size)
        gmsh_io.generate_extract_mesh(dims, mesh_output_name, mesh_output_dir, mesh_size_list,
                                      save_file, open_gmsh_gui, arbitrary_mesh_size)

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
        # if "True", mesh size is defined by the user in the dictionary below; if "False" default_mesh_size is used
        arbitrary_mesh_size: bool = False
        # if "arbitrary_mesh_size = False", define a mesh size, if set to -1 mesh size is logically chosen by Gmsh
        # itself based on the geometry
        default_mesh_size: float = -1
        # define the points of the surface as a list of tuples
        input_points_list = [[(0, 0, 0), (3, 0, 0), (3, 1, 0), (0, 1, 0)],
                             [(3, 0, 0), (5, 0, 0), (5, 1, 0), (4, 1.5, 0), (3, 1, 0)],
                             [(0, 1, 0), (2, 1, 0), (2, 3, 0), (0, 3, 0)],
                             [(2, 1, 0), (3, 1, 0), (4, 1.5, 0), (5, 1, 0), (5, 3, 0), (2, 3, 0)],
                             [(0, 3, 0), (2.5, 3, 0), (2, 4, 0), (0, 4, 0)],
                             [(0.8, 4, 0), (1.2, 4, 0), (1.2, 4.1, 0), (0.8, 4.1, 0)]]
        # define the name labels for the surfaces
        name_label_list = ["First Soil Layer", "FSL", "Second Soil Layer", "SSL", "Embankment", "Soil Ballast"]
        # define the mesh size for each surface
        mesh_size_list = []
        # define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
        dims = 3
        # if 3D, input depth of geometry to be extruded from 2D surface
        extrusion_length = [0, 0, 3]
        # if "True", saves mesh data to separate mdpa files; otherwise "False"
        save_file = False
        # if "True", opens gmsh interface; otherwise "False"
        open_gmsh_gui = False
        # set a name for mesh output file
        mesh_output_name = "test_3D"
        # set output directory
        mesh_output_dir = "."

        gmsh_io = GmshIO()

        gmsh_io.generate_geometry(input_points_list, extrusion_length, dims, mesh_output_name,
                                  name_label_list, default_mesh_size)
        gmsh_io.generate_extract_mesh(dims, mesh_output_name, mesh_output_dir, mesh_size_list,
                                      save_file, open_gmsh_gui, arbitrary_mesh_size)

        mesh_data = gmsh_io.mesh_data

        assert mesh_data["nodes"]["coordinates"].size > 0  # check if node_coords is not empty
        assert mesh_data["nodes"]["ids"].size > 0  # check if node_tags is not empty
        assert list(mesh_data["elements"].keys()) == ["LINE_2N", "TRIANGLE_3N", "TETRAHEDRON_4N",
                                                      "POINT_1N"]  # check if correct elements are present

        # check each element type contains ids and nodes
        for value in mesh_data["elements"].values():
            assert value["element_ids"].size > 0
            assert value["connectivities"].size > 0

    def test_generate_different_mesh_sizes_2D(self):
        """
        Checks whether mesh data generated for 2D geometries is not empty.

        """
        # if "True", mesh size is defined by the user in the dictionary below; if "False" default_mesh_size is used
        arbitrary_mesh_size: bool = True
        # if "arbitrary_mesh_size = False", define a mesh size, if set to -1 mesh size is logically chosen by Gmsh
        # itself based on the geometry
        default_mesh_size: float = -1
        # define the points of the surface as a list of tuples
        input_points_list = [[(0, 0, 0), (3, 0, 0), (5, 1.5, 0), (2, 1, 0), (0, 1, 0)],
                             [(3, 0, 0), (5, 0, 0), (5, 1.5, 0)],
                             [(0, 1, 0), (0, 3, 0), (2, 3, 0), (2, 1, 0)]]
        # define the name labels for the surfaces
        name_label_list = ["First left Soil Layer", "Second right Soil Layer", "Third top Soil Layer"]
        # define the mesh size for each surface
        mesh_size_list = [10, 0.1, 0.2]
        # define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
        dims = 2
        # if 3D, input depth of geometry to be extruded from 2D surface
        extrusion_length = [0, 0, 0]
        # if "True", saves mesh data to separate mdpa files; otherwise "False"
        save_file = False
        # if "True", opens gmsh interface; otherwise "False"
        open_gmsh_gui = False
        # set a name for mesh output file
        mesh_output_name = "test_2D"
        # set output directory
        mesh_output_dir = "."

        gmsh_io = GmshIO()

        gmsh_io.generate_geometry(input_points_list, extrusion_length, dims, mesh_output_name,
                                  name_label_list, default_mesh_size)
        gmsh_io.generate_extract_mesh(dims, mesh_output_name, mesh_output_dir, mesh_size_list,
                                      save_file, open_gmsh_gui, arbitrary_mesh_size)

        mesh_data = gmsh_io.mesh_data

        assert mesh_data["nodes"]["coordinates"].size > 0  # check if node_coords is not empty
        assert mesh_data["nodes"]["ids"].size > 0  # check if node_tags is not empty
        assert list(mesh_data["elements"].keys()) == ["LINE_2N", "TRIANGLE_3N",
                                                      "POINT_1N"]  # check if correct elements are present

        # check each element type contains ids and nodes
        for value in mesh_data["elements"].values():
            assert value["element_ids"].size > 0
            assert value["connectivities"].size > 0

    def test_generate_different_mesh_sizes_3D(self):
        """
        Checks whether mesh data generated for 3D geometries is not empty.

        """
        # if "True", mesh size is defined by the user in the dictionary below; if "False" default_mesh_size is used
        arbitrary_mesh_size: bool = True
        # if "arbitrary_mesh_size = False", define a mesh size, if set to -1 mesh size is logically chosen by Gmsh
        # itself based on the geometry
        default_mesh_size: float = -1
        # define the points of the surface as a list of tuples
        input_points_list = [[(0, 0, 0), (3, 0, 0), (5, 1.5, 0), (2, 1, 0), (0, 1, 0)],
                             [(3, 0, 0), (5, 0, 0), (5, 1.5, 0)],
                             [(0, 1, 0), (0, 3, 0), (2, 3, 0), (2, 1, 0)]]
        # define the name labels for the surfaces
        name_label_list = ["First left Soil Layer", "Second right Soil Layer", "Third top Soil Layer"]
        # define the mesh size for each surface
        mesh_size_list = [10, 0.1, 0.2]
        # define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
        dims = 3
        # if 3D, input depth of geometry to be extruded from 2D surface
        extrusion_length = [0, 0, 3]
        # if "True", saves mesh data to separate mdpa files; otherwise "False"
        save_file = False
        # if "True", opens gmsh interface; otherwise "False"
        open_gmsh_gui = False
        # set a name for mesh output file
        mesh_output_name = "test_3D"
        # set output directory
        mesh_output_dir = "."

        gmsh_io = GmshIO()

        gmsh_io.generate_geometry(input_points_list, extrusion_length, dims, mesh_output_name,
                                  name_label_list, default_mesh_size)
        gmsh_io.generate_extract_mesh(dims, mesh_output_name, mesh_output_dir, mesh_size_list,
                                      save_file, open_gmsh_gui, arbitrary_mesh_size)

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
        geo_file = r"tests/test_data/column_2D.geo"

        gmsh_io = GmshIO()
        gmsh_io.read_gmsh_geo(geo_file)

        geo_data = gmsh_io.geo_data

        # check if the coordinates of the points are correct
        expected_points = {1: [0, 0, 0], 2: [0.5, 0, 0], 3: [0.5, 1, 0], 4: [0, 1, 0], 11: [0, 2, 0], 12: [0.5, 2.0, 0]}
        expected_lines = {5: [1, 2], 6: [2, 3], 7: [3, 4], 8: [4, 1], 13: [4, 11], 14: [11, 12], 15: [12, 3]}
        expected_surfaces = {10: [5, 6, 7, 8], 17: [-13, -7, -15, -14]}  # negative sign means reversed orientation

        expected_physical_groups = {'group_1': {'geometry_ids': [10], 'id': 1, 'ndim': 2},
                                    'group_2': {'geometry_ids': [17], 'id': 2, 'ndim': 2}}

        expected_geo_data = {"points": expected_points,
                             "lines": expected_lines,
                             "surfaces": expected_surfaces,
                             "volumes": {},
                             "physical_groups": expected_physical_groups}

        # check if expected and actual geo data are equal
        TestUtils.assert_dictionary_almost_equal(expected_geo_data, geo_data)

    def test_read_gmsh_geo_3D(self, expected_geo_data_3D):
        """
        Checks whether a gmsh .geo file is read correctly. For a 3d geometry
        """
        geo_file = r"tests/test_data/column_3D_tetra4.geo"

        gmsh_io = GmshIO()
        gmsh_io.read_gmsh_geo(geo_file)

        geo_data = gmsh_io.geo_data

        # check if expected and actual geo data are equal
        TestUtils.assert_dictionary_almost_equal(expected_geo_data_3D, geo_data)

    def test_read_gmsh_msh_2D(self):
        """
        Checks whether a gmsh .msh file is read correctly in 2D.

        """

        msh_file = r"tests/test_data/block_2D.msh"

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

        msh_file = r"tests/test_data/block_3D.msh"

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

    def test_generate_geo_from_geo_data(self, expected_geo_data_3D):
        """
        Checks if the gmsh geometry is correctly generated from the geo data dictionary.

        This test sets the geo_data dictionary manually, this data is then used to generate the gmsh geometry input.
        The generated geometry input is then read and the geo data is extracted from it. The input and output should be
        equal, besides reoriented lines within surfaces.
        """

        geo_data = expected_geo_data_3D

        gmsh_io = GmshIO()

        # set private attribute geo_data
        gmsh_io._GmshIO__geo_data = geo_data

        # generate gmsh geo input from geo data dictionary
        gmsh_io.generate_geo_from_geo_data()

        # retrieve the generated geo input
        gmsh_io.extract_geo_data()

        new_geo_data = gmsh_io.geo_data

        # only check absolute and sorted values in surfaces, because the values can be reoriented by occ
        for surface_id in geo_data["surfaces"].keys():
            geo_data["surfaces"][surface_id] = np.sort(np.abs(geo_data["surfaces"][surface_id]).astype(int)).tolist()

        for surface_id in new_geo_data["surfaces"].keys():
            new_geo_data["surfaces"][surface_id] = np.sort(np.abs(new_geo_data["surfaces"][surface_id]).astype(int))

        # only check absolute and sorted values in volumes, because the values can be reoriented by occ
        for volume_id in new_geo_data["volumes"].keys():
            new_geo_data["volumes"][volume_id] = np.sort(np.abs(new_geo_data["volumes"][volume_id]).astype(int))

        for volume_id in geo_data["volumes"].keys():
            geo_data["volumes"][volume_id] = np.sort(np.abs(geo_data["volumes"][volume_id]).astype(int))

        # check if expected and actual geo data are equal
        TestUtils.assert_dictionary_almost_equal(geo_data, new_geo_data)

    def test_generate_geo_from_geo_data_with_shared_group(self, expected_geo_data_3D_with_shared_group):
        """
        Checks if the gmsh geometry is correctly generated from the geo data dictionary. In this test, two volumes
        are added to a separate group. And the same volumes are added to a shared group.

        This test sets the geo_data dictionary manually, this data is then used to generate the gmsh geometry input.
        The generated geometry input is then read and the geo data is extracted from it. The input and output should be
        equal, besides reoriented lines within surfaces.
        """

        geo_data = expected_geo_data_3D_with_shared_group

        gmsh_io = GmshIO()

        # set private attribute geo_data
        gmsh_io._GmshIO__geo_data = geo_data

        # generate gmsh geo input from geo data dictionary
        gmsh_io.generate_geo_from_geo_data()

        # retrieve the generated geo input
        gmsh_io.extract_geo_data()

        new_geo_data = gmsh_io.geo_data

        # only check absolute and sorted values in surfaces, because the values can be reoriented by occ
        for surface_id in geo_data["surfaces"].keys():
            geo_data["surfaces"][surface_id] = np.sort(np.abs(geo_data["surfaces"][surface_id]).astype(int)).tolist()

        for surface_id in new_geo_data["surfaces"].keys():
            new_geo_data["surfaces"][surface_id] = np.sort(np.abs(new_geo_data["surfaces"][surface_id]).astype(int))

        # only check absolute and sorted values in volumes, because the values can be reoriented by occ
        for volume_id in new_geo_data["volumes"].keys():
            new_geo_data["volumes"][volume_id] = np.sort(np.abs(new_geo_data["volumes"][volume_id]).astype(int))

        for volume_id in geo_data["volumes"].keys():
            geo_data["volumes"][volume_id] = np.sort(np.abs(geo_data["volumes"][volume_id]).astype(int))

        # check if expected and actual geo data are equal
        TestUtils.assert_dictionary_almost_equal(geo_data, new_geo_data)

    def test_generate_mesh(self):
        """
        Checks whether a mesh is generated correctly from a gmsh .geo file. A 2D block mesh is generated.

        """

        geo_file = r"tests/test_data/block_2D.geo"

        gmsh_io = GmshIO()

        # read geo file
        gmsh_io.read_gmsh_geo(geo_file)

        # generate mesh
        gmsh_io.generate_mesh(2, 1)

        # get mesh data
        mesh_data = gmsh_io.mesh_data

        # set expected mesh data
        expected_mesh_data = {'nodes': {'coordinates': np.array([[0., 0., 0.],
                                                                 [1., 0., 0.],
                                                                 [1., 1., 0.],
                                                                 [0., 1., 0.],
                                                                 [0.5, 0.5, 0.]]),
                                        'ids': np.array([1, 2, 3, 4, 5])},
                              'elements': {'LINE_2N': {'element_ids': np.array([9, 10, 11, 12]),
                                                       'connectivities': np.array([[1, 2],
                                                                                   [2, 3],
                                                                                   [3, 4],
                                                                                   [4, 1]])},
                                           'TRIANGLE_3N': {'element_ids': np.array([1, 2, 3, 4]),
                                                           'connectivities': np.array([[2, 5, 1],
                                                                                       [1, 5, 4],
                                                                                       [3, 5, 2],
                                                                                       [4, 5, 3]])},
                                           'POINT_1N': {'element_ids': np.array([5, 6, 7, 8]),
                                                        'connectivities': np.array([[1],
                                                                                    [2],
                                                                                    [3],
                                                                                    [4]])}}}

        # check if the coordinates of the points are correct
        TestUtils.assert_dictionary_almost_equal(expected_mesh_data, mesh_data)

    def test_physical_groups_in_geometry_data_2D(self):
        """
        Checks whether geometry data in 2D geometry has physical groups
    """
        # define the default mesh size
        default_mesh_size = -1
        # define the points of the surface as a list of tuples
        input_points_list = [[(0, 0, 0), (3, 0, 0), (3, 1, 0), (0, 1, 0)],
                             [(0, 1, 0), (3, 1, 0), (3, 2, 0), (0, 2, 0)],
                             [(1, 2, 0), (2, 2, 0), (2, 2.5, 0), (1, 2.5, 0)]]
        # define the name labels for the surfaces
        name_label_list = ["Soil Layer", "Soil Embankment", "Soil Ballast"]

        # define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
        dims = 2
        # if 3D, input depth of geometry to be extruded from 2D surface
        extrusion_length = [0, 0, 0]
        # set a name for mesh output file
        mesh_output_name = "test_2D"
        # set output directory
        mesh_output_dir = "."

        gmsh_io = GmshIO()

        gmsh_io.generate_geometry(input_points_list, extrusion_length, dims, mesh_output_name,
                                  name_label_list, default_mesh_size)

        geo_data = gmsh_io.geo_data

        expected_physical_groups = {'Soil Layer': {'ndim': 2, 'id': 1, 'geometry_ids': [1]},
                                    'Soil Embankment': {'ndim': 2, 'id': 2, 'geometry_ids': [2]},
                                    'Soil Ballast': {'ndim': 2, 'id': 3, 'geometry_ids': [3]}}

        # check if expected and actual geo data are equal
        TestUtils.assert_dictionary_almost_equal(expected_physical_groups, geo_data["physical_groups"])

    def test_physical_groups_in_geometry_data_3D(self):
        """
        Checks whether geometry data in 3D geometry has physical groups

        """

        # define the default mesh size
        default_mesh_size = 1
        # define the points of the surface as a list of tuples
        input_points_list = [[(0, 0, 0), (3, 0, 0), (3, 1, 0), (0, 1, 0)],
                             [(0, 1, 0), (3, 1, 0), (3, 2, 0), (0, 2, 0)],
                             [(1, 2, 0), (2, 2, 0), (2, 2.5, 0), (1, 2.5, 0)]]
        # define the name labels for the surfaces
        name_label_list = ["Soil Layer", "Soil Embankment", "Soil Ballast"]

        # define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
        dims = 3
        # if 3D, input depth of geometry to be extruded from 2D surface
        extrusion_length = [0, 0, 3]
        # set a name for mesh output file
        mesh_output_name = "test_3D"
        # set output directory
        mesh_output_dir = "."

        gmsh_io = GmshIO()

        gmsh_io.generate_geometry(input_points_list, extrusion_length, dims, mesh_output_name,
                                  name_label_list, default_mesh_size)

        geo_data = gmsh_io.geo_data

        expected_physical_groups = {'Soil Layer': {'ndim': 3, 'id': 1, 'geometry_ids': [1]},
                                    'Soil Embankment': {'ndim': 3, 'id': 2, 'geometry_ids': [2]},
                                    'Soil Ballast': {'ndim': 3, 'id': 3, 'geometry_ids': [3]}}

        # check if expected and actual geo data are equal
        TestUtils.assert_dictionary_almost_equal(expected_physical_groups, geo_data["physical_groups"])
