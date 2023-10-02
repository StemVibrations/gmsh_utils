from pathlib import Path
import re

import gmsh
import numpy as np
import pytest
from _pytest.capture import CaptureFixture

from gmsh_utils.gmsh_IO import GmshIO
from utils import TestUtils


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

        # expected_points = {1: [0.0, 0.0, 0.0], 2: [0.5, 0.0, 0.0], 3: [0.5, 1.0, 0.0], 4: [0.0, 1.0, 0.0],
        #                    5: [0.5, 0.0, -0.5], 6: [0.0, 0.0, -0.5], 7: [0.5, 1.0, -0.5], 8: [0.0, 1.0, -0.5],
        #                    9: [0.0, 2.0, 0.0], 10: [0.5, 2.0, 0.0], 11: [0.0, 2.0, -0.5], 12: [0.5, 2.0, -0.5]}
        # expected_lines = {1: [1, 2], 2: [2, 3], 3: [3, 4], 4: [4, 1], 5: [2, 5], 6: [6, 5], 7: [1, 6], 8: [5, 7],
        #                   9: [7, 8], 10: [8, 6], 11: [3, 7], 12: [4, 8], 13: [4, 9], 14: [9, 10], 15: [10, 3],
        #                   16: [9, 11], 17: [11, 8], 18: [7, 12], 19: [12, 11], 20: [10, 12]}
        # expected_surfaces = {1: [1, 2, 3, 4],
        #              2: [1, 5, -6, -7],
        #              3: [6, 8, 9, 10],
        #              4: [2, 11, -8, -5],
        #              5: [3, 12, -9, -11],
        #              6: [4, 7, -10, -12],
        #              7: [13, 14, 15, 3],
        #              8: [13, 16, 17, -12],
        #              9: [17, -9, 18, 19],
        #              10: [15, 11, 18, -20],
        #              11: [14, 20, 19, -16]}
        # expected_volumes={1: [1, -2, -3, -4, -5, -6], 2: [-7, 8, -9, 5, 10, 11]}
        # expected_physical_groups = {'group_1': {'geometry_ids': [1], 'id': 1, 'ndim': 3},
        #                             'group_2': {'geometry_ids': [2], 'id': 2, 'ndim': 3}}

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
        # define the default mesh size
        default_mesh_size = -1

        # define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
        dims = 2

        input_dict = {'First left Soil Layer': {"element_size": default_mesh_size,
                                                "coordinates": [(0, 0, 0), (3, 0, 0), (3, 1, 0), (0, 1, 0)],
                                                "ndim": dims},
                      'FSL': {"element_size": default_mesh_size,
                              "coordinates": [(3, 0, 0), (5, 0, 0), (5, 1, 0), (4, 1.5, 0), (3, 1, 0)],
                              "ndim": dims},
                      'Third top Soil Layer': {"element_size": default_mesh_size,
                                               "coordinates": [(0, 1, 0), (2, 1, 0), (2, 3, 0), (0, 3, 0)],
                                               "ndim": dims},
                      'SSL': {"element_size": default_mesh_size,
                              "coordinates": [(2, 1, 0), (3, 1, 0), (4, 1.5, 0), (5, 1, 0), (5, 3, 0), (2, 3, 0)],
                              "ndim": dims},
                      'Soil Ballast': {"element_size": default_mesh_size,
                                       "coordinates": [(0, 3, 0), (2.5, 3, 0), (2, 4, 0), (0, 4, 0)],
                                       "ndim": dims},
                      'Line Track': {"element_size": default_mesh_size,
                                     "coordinates": [(0.8, 4, 0), (1.2, 4, 0), (1.2, 4.1, 0), (0.8, 4.1, 0)],
                                     "ndim": dims}
                      }

        # if "True", saves mesh data to separate mdpa files; otherwise "False"
        save_file = False
        # if "True", opens gmsh interface; otherwise "False"
        open_gmsh_gui = False
        # set a name for mesh output file
        mesh_output_name = "test_2D"
        # set output directory
        mesh_output_dir = "."

        gmsh_io = GmshIO()

        gmsh_io.generate_geometry(input_dict, mesh_output_name)
        gmsh_io.generate_extract_mesh(dims, mesh_output_name, mesh_output_dir, save_file, open_gmsh_gui)

        mesh_data = gmsh_io.mesh_data

        # check if nodes are not empty
        assert len(mesh_data["nodes"]) > 0

        # check if correct element types are present
        assert list(mesh_data["elements"].keys()) == ["LINE_2N", "TRIANGLE_3N",
                                                      "POINT_1N"]

        # check elements are not empty
        for value in mesh_data["elements"].values():
            assert len(value) > 0

    def test_generate_mesh_3D(self):
        """
        Checks whether mesh data generated for 3D geometries is not empty.

        """

        # define the default mesh size
        default_mesh_size = 1

        # define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
        dims = 3
        # if 3D, input depth of geometry to be extruded from 2D surface
        extrusion_length = [0, 0, 3]

        input_dict = {'First left Soil Layer': {"element_size": default_mesh_size,
                                                "coordinates": [(0, 0, 0), (3, 0, 0), (3, 1, 0), (0, 1, 0)],
                                                "ndim": dims,
                                                "extrusion_length": extrusion_length},
                      'FSL': {"element_size": default_mesh_size,
                              "coordinates": [(3, 0, 0), (5, 0, 0), (5, 1, 0), (4, 1.5, 0), (3, 1, 0)],
                              "ndim": dims,
                              "extrusion_length": extrusion_length},
                      'Third top Soil Layer': {"element_size": default_mesh_size,
                                               "coordinates": [(0, 1, 0), (2, 1, 0), (2, 3, 0), (0, 3, 0)],
                                               "ndim": dims,
                                               "extrusion_length": extrusion_length},
                      'SSL': {"element_size": default_mesh_size,
                              "coordinates": [(2, 1, 0), (3, 1, 0), (4, 1.5, 0), (5, 1, 0), (5, 3, 0), (2, 3, 0)],
                              "ndim": dims,
                              "extrusion_length": extrusion_length},
                      'Soil Ballast': {"element_size": default_mesh_size,
                                       "coordinates": [(0, 3, 0), (2.5, 3, 0), (2, 4, 0), (0, 4, 0)],
                                       "ndim": dims,
                                       "extrusion_length": extrusion_length},
                      'Line Track': {"element_size": default_mesh_size,
                                     "coordinates": [(0.8, 4, 0), (1.2, 4, 0), (1.2, 4.1, 0), (0.8, 4.1, 0)],
                                     "ndim": dims,
                                     "extrusion_length": extrusion_length}
                      }

        # if "True", saves mesh data to separate mdpa files; otherwise "False"
        save_file = False
        # if "True", opens gmsh interface; otherwise "False"
        open_gmsh_gui = False
        # set a name for mesh output file
        mesh_output_name = "test_3D"
        # set output directory
        mesh_output_dir = "."

        gmsh_io = GmshIO()

        gmsh_io.generate_geometry(input_dict, mesh_output_name)
        gmsh_io.generate_extract_mesh(dims, mesh_output_name, mesh_output_dir, save_file, open_gmsh_gui)

        mesh_data = gmsh_io.mesh_data

        # check if nodes are not empty
        assert len(mesh_data["nodes"]) > 0

        # check if correct element types are present
        assert list(mesh_data["elements"].keys()) == ["LINE_2N", "TRIANGLE_3N", "TETRAHEDRON_4N",
                                                      "POINT_1N"]

        # check elements are not empty
        for value in mesh_data["elements"].values():
            assert len(value) > 0

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

        expected_mesh_data = {"ndim": 2,
                              'nodes': {1: [0., 0., 0.],
                                        2: [1., 0., 0.],
                                        3: [1., 1., 0.],
                                        4: [0., 1., 0.]},
                              'elements': {'TRIANGLE_3N': {1: [1, 2, 4],
                                                           2: [4, 2, 3]}}}

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
        expected_mesh_data = {'nodes': {1: [0.0, 0.0, 0.0], 2: [1.0, 0.0, 0.0],
                                        3: [1.0, 1.0, 0.0], 4: [0.0, 1.0, 0.0],
                                        5: [0.0, 0.0, -1.0], 6: [1.0, 0.0, -1.0],
                                        7: [1.0, 1.0, -1.0], 8: [0.0, 1.0, -1.0]},
                              'elements': {'TETRAHEDRON_4N': {1: [2, 1, 4, 8], 2: [5, 6, 8, 2], 3: [5, 2, 8, 1],
                                                              4: [2, 4, 3, 7], 5: [8, 6, 7, 2], 6: [8, 2, 7, 4]}}}


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
        expected_mesh_data = {"ndim": 2,
                              'nodes': {1: [0., 0., 0.],
                                        2: [1., 0., 0.],
                                        3: [1., 1., 0.],
                                        4: [0., 1., 0.],
                                        5: [0.5, 0.5, 0.]},
                              'elements': {'LINE_2N': {9: [1, 2],
                                                       10: [2, 3],
                                                       11: [3, 4],
                                                       12: [4, 1]},
                                           'TRIANGLE_3N': {1: [2, 5, 1],
                                                           2: [1, 5, 4],
                                                           3: [3, 5, 2],
                                                           4: [4, 5, 3]},
                                           'POINT_1N': {5: [1],
                                                        6: [2],
                                                        7: [3],
                                                        8: [4]}},
                              'physical_groups': {'group_1': {"ndim": 2,
                                                              "node_ids": [1, 2, 3, 4, 5],
                                                              "element_ids": [1, 2, 3, 4],
                                                              "element_type": "TRIANGLE_3N"}}}

        # check if the coordinates of the points are correct
        TestUtils.assert_dictionary_almost_equal(expected_mesh_data, mesh_data)

    def test_extract_mesh_with_multiple_physical_groups(self):
        """
        Checks whether a mesh is generated correctly from a gmsh .geo file. A 2D column mesh is generated. Consisting of
        2 blocks. Each block has its own physical group. The two blocks are also combined in a physical group. Also
        lines and points are added to the physical groups.

        """
        geo_file = r"tests/test_data/column_2D_more_groups.geo"

        gmsh_io = GmshIO()

        # read geo file
        gmsh_io.read_gmsh_geo(geo_file)

        # generate mesh
        gmsh_io.generate_mesh(2, 1)

        # get mesh data
        mesh_data = gmsh_io.mesh_data

        # set expected mesh group data
        expected_groups_in_mesh_data = {'group_1': {"ndim": 2,
                                                    "node_ids": [1, 2, 3, 4, 7],
                                                    "element_ids": [5, 6, 7, 8],
                                                    "element_type": "TRIANGLE_3N"},
                                        'group_2': {"ndim": 2,
                                                    "node_ids": [3, 4, 5, 6, 8],
                                                    "element_ids": [9, 10, 11, 12],
                                                    "element_type": "TRIANGLE_3N"},
                                        "combined_group": {"ndim": 2,
                                                           "node_ids": [1, 2, 3, 4, 5, 6, 7, 8],
                                                           "element_ids": [5, 6, 7, 8, 9, 10, 11, 12],
                                                           "element_type": "TRIANGLE_3N"},
                                        "line_group": {"ndim": 1,
                                                       "node_ids": [1, 2, 3, 6],
                                                       "element_ids": [3, 4],
                                                       "element_type": "LINE_2N"},
                                        "point_group": {"ndim": 0,
                                                        "node_ids": [1, 2],
                                                        "element_ids": [1, 2],
                                                        "element_type": "POINT_1N"}}

        # check if the coordinates of the points are correct
        TestUtils.assert_dictionary_almost_equal(expected_groups_in_mesh_data, mesh_data["physical_groups"])

    def test_physical_groups_in_geometry_data_2D(self):
        """
        Checks whether geometry data in 2D geometry has physical groups
        """

        # define the default mesh size
        default_mesh_size = -1

        # define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
        dims = 2

        # set a name for mesh output file
        mesh_output_name = "test_2D"

        # set input dictionary
        input_dict = {'Soil Layer': {"element_size": default_mesh_size,
                                     "coordinates": [(0, 0, 0), (3, 0, 0), (3, 1, 0), (0, 1, 0)],
                                     "ndim": dims},
                      'Soil Embankment': {"element_size": default_mesh_size,
                                          "coordinates": [(0, 1, 0), (3, 1, 0), (3, 2, 0), (0, 2, 0)],
                                          "ndim": dims},
                      'Soil Ballast': {"element_size": default_mesh_size,
                                       "coordinates": [(1, 2, 0), (2, 2, 0), (2, 2.5, 0), (1, 2.5, 0)],
                                       "ndim": dims}}

        gmsh_io = GmshIO()

        gmsh_io.generate_geometry(input_dict, mesh_output_name)

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

        # define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
        dims = 3
        # if 3D, input depth of geometry to be extruded from 2D surface
        extrusion_length = [0, 0, 3]
        # set a name for mesh output file
        mesh_output_name = "test_3D"

        # set input dictionary
        input_dict = {'Soil Layer': {"element_size": default_mesh_size,
                                     "coordinates": [(0, 0, 0), (3, 0, 0), (3, 1, 0), (0, 1, 0)],
                                     "ndim": dims,
                                     "extrusion_length": extrusion_length},
                      'Soil Embankment': {"element_size": default_mesh_size,
                                          "coordinates": [(0, 1, 0), (3, 1, 0), (3, 2, 0), (0, 2, 0)],
                                          "ndim": dims,
                                          "extrusion_length": extrusion_length},
                      'Soil Ballast': {"element_size": default_mesh_size,
                                       "coordinates": [(1, 2, 0), (2, 2, 0), (2, 2.5, 0), (1, 2.5, 0)],
                                       "ndim": dims,
                                       "extrusion_length": extrusion_length}}

        gmsh_io = GmshIO()
        gmsh_io.generate_geometry(input_dict, mesh_output_name)

        geo_data = gmsh_io.geo_data

        expected_physical_groups = {'Soil Layer': {'ndim': 3, 'id': 1, 'geometry_ids': [1]},
                                    'Soil Embankment': {'ndim': 3, 'id': 2, 'geometry_ids': [2]},
                                    'Soil Ballast': {'ndim': 3, 'id': 3, 'geometry_ids': [3]}}

        # check if expected and actual geo data are equal
        TestUtils.assert_dictionary_almost_equal(expected_physical_groups, geo_data["physical_groups"])

    def test_finalize_gmsh(self):
        """
        Checks whether gmsh is finalized after calling finalize_gmsh

        """

        # initialize gmsh
        gmsh_io = GmshIO()
        gmsh.initialize()

        # check if gmsh is initialized
        assert gmsh.isInitialized()

        # finalize gmsh
        gmsh_io.finalize_gmsh()

        # check if gmsh is finalized
        assert not gmsh.isInitialized()

    def test_synchronize_gmsh(self):
        """
        Checks whether gmsh is synchronized after calling synchronize_gmsh. This test checks whether the geo data
        is updated with the created points and the physical group after calling synchronize_gmsh.

        """

        # initialize gmsh
        gmsh_io = GmshIO()
        gmsh.initialize()

        # create two points and add to physical group
        point_id1 = gmsh.model.occ.addPoint(0, 0, 0)
        point_id2 = gmsh.model.occ.addPoint(1, 0, 0)

        gmsh.model.addPhysicalGroup(0, [point_id1, point_id2], name="test")

        # extract geo data
        gmsh_io.extract_geo_data()
        empty_geo_data = gmsh_io.geo_data

        # check if geo data is empty before synchronizing
        expected_empty_geo_data = {'lines': {}, 'physical_groups': {}, 'points': {}, 'surfaces': {}, 'volumes': {}}
        TestUtils.assert_dictionary_almost_equal(empty_geo_data, expected_empty_geo_data)

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()

        # extract geo data
        gmsh_io.extract_geo_data()
        filled_geo_data = gmsh_io.geo_data

        # check if geo data is filled after synchronizing
        expected_filled_geo_data = {'points': {1: [0., 0., 0.], 2: [1., 0., 0.]},
                                    'lines': {},
                                    'surfaces': {},
                                    'volumes': {},
                                    'physical_groups': {'test': {'geometry_ids': [1, 2], 'id': 1, 'ndim': 0}}}

        TestUtils.assert_dictionary_almost_equal(filled_geo_data, expected_filled_geo_data)

    def test_synchronize_gmsh_with_intersection_point_on_line(self):
        """
        Checks whether gmsh is synchronized after calling synchronize_gmsh. This test checks whether the geo data
        is updated with the created points and the physical group after calling synchronize_gmsh.

        This test is for the case where the intersection points are on two lines.

        """

        # initialize gmsh
        gmsh_io = GmshIO()
        gmsh.initialize()

        # create two points and add to physical group
        point_id1 = gmsh.model.occ.addPoint(0, 0, 0)
        point_id2 = gmsh.model.occ.addPoint(1, 0, 0)
        point_id5 = gmsh.model.occ.addPoint(2, 0, 0)

        # create a line
        line_id = gmsh.model.occ.addLine(point_id1, point_id2)
        line_id2 = gmsh.model.occ.addLine(point_id2, point_id5)
        gmsh.model.addPhysicalGroup(1, [line_id, line_id2], name="line")

        point_id3 = gmsh.model.occ.addPoint(0.5, 0, 0)
        gmsh.model.addPhysicalGroup(0, [point_id3], name="first_new_point")

        point_id4 = gmsh.model.occ.addPoint(1.5, 0, 0)
        gmsh.model.addPhysicalGroup(0, [point_id4], name="second_new_point")

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()

        # extract geo data
        gmsh_io.extract_geo_data()
        filled_geo_data = gmsh_io.geo_data

        expected_geo_data = {'points': {1: [0., 0., 0.], 2: [1., 0., 0.], 3: [2.0, 0.0,0.0],
                                        4: [0.5, 0., 0.], 5: [1.5, 0., 0.]},
                             'lines': {1: [1, 4], 2: [4, 2], 3: [2, 5], 4: [5, 3]},
                             'surfaces': {},
                             'volumes': {},
                             'physical_groups': {'line': {'geometry_ids': [1, 2, 3, 4], 'id': 1, 'ndim': 1},
                                                 'first_new_point': {'geometry_ids': [4], 'id': 2, 'ndim': 0},
                                                 'second_new_point': {'geometry_ids': [5], 'id': 3, 'ndim': 0}}}

        TestUtils.assert_dictionary_almost_equal(filled_geo_data, expected_geo_data)

        # synchronize gmsh again to check if the points are not added again
        gmsh_io.synchronize_gmsh()

        # extract geo data
        gmsh_io.extract_geo_data()
        filled_geo_data = gmsh_io.geo_data

        # check if geo data hasn't changed after re-synchronizing
        TestUtils.assert_dictionary_almost_equal(filled_geo_data, expected_geo_data)

    def test_synchronize_gmsh_with_intersection_line_on_surface(self):
        """
        Checks whether gmsh is synchronized after calling synchronize_gmsh. This test checks whether the geo data
        is updated with the created points and the physical group after calling synchronize_gmsh.

        This test is for the case where the intersection line is on a surface boundary.
        """

        # initialize gmsh
        gmsh_io = GmshIO()
        gmsh.initialize()

        # create surface points
        point_id1 = gmsh.model.occ.addPoint(0, 0, 0)
        point_id2 = gmsh.model.occ.addPoint(1, 0, 0)
        point_id3 = gmsh.model.occ.addPoint(1, 1, 0)
        point_id4 = gmsh.model.occ.addPoint(0, 1, 0)

        # create a surface
        line_id1 = gmsh.model.occ.addLine(point_id1, point_id2)
        line_id2 = gmsh.model.occ.addLine(point_id2, point_id3)
        line_id3 = gmsh.model.occ.addLine(point_id3, point_id4)
        line_id4 = gmsh.model.occ.addLine(point_id4, point_id1)

        curve_loop_id = gmsh.model.occ.addCurveLoop([line_id1, line_id2, line_id3, line_id4])
        surface_id: int = gmsh.model.occ.addPlaneSurface([curve_loop_id])
        gmsh.model.addPhysicalGroup(2, [surface_id], tag=-1, name="surface")

        # create an intersection line
        point_id5 = gmsh.model.occ.addPoint(0.25, 0.0, 0)
        point_id6 = gmsh.model.occ.addPoint(0.75, 0.0, 0)
        line_id5 = gmsh.model.occ.addLine(point_id5, point_id6)
        gmsh.model.addPhysicalGroup(1, [line_id5], name="new_line")

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()

        # extract geo data
        gmsh_io.extract_geo_data()
        filled_geo_data = gmsh_io.geo_data

        expected_geo_data = {'points': {1: [0., 0., 0.], 2: [1., 0., 0.], 3: [1, 1., 0.],
                                        4: [0., 1., 0.], 5: [0.25, 0., 0.], 6: [0.75, 0., 0.]},
                             'lines': {2: [2, 3], 3: [3, 4], 4: [4, 1], 5: [5, 6], 6: [1, 5], 7: [6, 2]},
                             'surfaces': {1: [6, 5, 7, 2, 3, 4]},
                             'volumes': {},
                             'physical_groups': {'surface': {'geometry_ids': [1], 'id': 1, 'ndim': 2},
                                                 'new_line': {'geometry_ids': [5], 'id': 2, 'ndim': 1}}}

        # check if geo data is as expected
        TestUtils.assert_dictionary_almost_equal(filled_geo_data, expected_geo_data)

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()

        # extract geo data
        gmsh_io.extract_geo_data()
        filled_geo_data = gmsh_io.geo_data

        # check if geo data hasn't changed after re-synchronizing
        TestUtils.assert_dictionary_almost_equal(filled_geo_data, expected_geo_data)

    def test_synchronize_gmsh_with_intersection_line_through_surface(self):
        """
        Checks whether gmsh is synchronized after calling synchronize_gmsh. This test checks whether the geo data
        is updated with the created points and the physical group after calling synchronize_gmsh.

        This test is for the case where the intersection line is on the surface.
        """

        # initialize gmsh
        gmsh_io = GmshIO()
        gmsh.initialize()

        # create surface points
        point_id1 = gmsh.model.occ.addPoint(0, 0, 0)
        point_id2 = gmsh.model.occ.addPoint(1, 0, 0)
        point_id3 = gmsh.model.occ.addPoint(1, 1, 0)
        point_id4 = gmsh.model.occ.addPoint(0, 1, 0)

        # create a surface
        line_id1 = gmsh.model.occ.addLine(point_id1, point_id2)
        line_id2 = gmsh.model.occ.addLine(point_id2, point_id3)
        line_id3 = gmsh.model.occ.addLine(point_id3, point_id4)
        line_id4 = gmsh.model.occ.addLine(point_id4, point_id1)

        curve_loop_id = gmsh.model.occ.addCurveLoop([line_id1, line_id2, line_id3, line_id4])
        surface_id: int = gmsh.model.occ.addPlaneSurface([curve_loop_id])
        gmsh.model.addPhysicalGroup(2, [surface_id], tag=-1, name="surface")

        # create an intersection line
        point_id5 = gmsh.model.occ.addPoint(0.0, 0.0, 0)
        point_id6 = gmsh.model.occ.addPoint(1, 1.0, 0)
        line_id5 = gmsh.model.occ.addLine(point_id5, point_id6)
        gmsh.model.addPhysicalGroup(1, [line_id5], name="new_line")

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()

        # extract geo data
        gmsh_io.extract_geo_data()
        filled_geo_data = gmsh_io.geo_data

        expected_geo_data = {'points': {1: [0.0, 0.0, 0.0], 2: [1.0, 0.0, 0.0], 3: [1.0, 1.0, 0.0], 4: [0.0, 1.0, 0.0]},
                             'lines': {1: [1, 2], 2: [2, 3], 3: [3, 4], 4: [4, 1], 5: [1, 3]},
                             'surfaces': {1: [-5, 1, 2], 2: [4, 5, 3]},
                             'volumes': {},
                             'physical_groups': {'surface': {'geometry_ids': [1, 2], 'id': 1, 'ndim': 2},
                                                 'new_line': {'geometry_ids': [5], 'id': 2, 'ndim': 1}}}

        # check if geo data is as expected
        TestUtils.assert_dictionary_almost_equal(filled_geo_data, expected_geo_data)

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()

        # extract geo data
        gmsh_io.extract_geo_data()
        filled_geo_data = gmsh_io.geo_data

        # check if geo data hasn't changed after re-synchronizing
        TestUtils.assert_dictionary_almost_equal(filled_geo_data, expected_geo_data)

    def test_synchronize_gmsh_with_intersection_line_through_surface_with_multiple_groups(self):
        """
        Checks whether gmsh is synchronized after calling synchronize_gmsh. This test checks whether the geo data
        is updated with the created points and the physical group after calling synchronize_gmsh.

        This test is for the case where the intersection line is through a surface which contains multiple groups.
        """

        # initialize gmsh
        gmsh_io = GmshIO()
        gmsh.initialize()

        # create surface points
        point_id1 = gmsh.model.occ.addPoint(0, 0, 0)
        point_id2 = gmsh.model.occ.addPoint(1, 0, 0)
        point_id3 = gmsh.model.occ.addPoint(1, 1, 0)
        point_id4 = gmsh.model.occ.addPoint(0, 1, 0)

        # create a surface
        line_id1 = gmsh.model.occ.addLine(point_id1, point_id2)
        line_id2 = gmsh.model.occ.addLine(point_id2, point_id3)
        line_id3 = gmsh.model.occ.addLine(point_id3, point_id4)
        line_id4 = gmsh.model.occ.addLine(point_id4, point_id1)

        curve_loop_id = gmsh.model.occ.addCurveLoop([line_id1, line_id2, line_id3, line_id4])
        surface_id: int = gmsh.model.occ.addPlaneSurface([curve_loop_id])
        gmsh.model.addPhysicalGroup(2, [surface_id], tag=-1, name="surface")

        gmsh.model.addPhysicalGroup(2, [surface_id], tag=-1, name="group2_surface")

        # create an intersection line
        point_id5 = gmsh.model.occ.addPoint(0.0, 0.0, 0)
        point_id6 = gmsh.model.occ.addPoint(1, 1.0, 0)
        line_id5 = gmsh.model.occ.addLine(point_id5, point_id6)
        gmsh.model.addPhysicalGroup(1, [line_id5], name="new_line")

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()

        # extract geo data
        gmsh_io.extract_geo_data()
        filled_geo_data = gmsh_io.geo_data

        expected_geo_data = {'points': {1: [0.0, 0.0, 0.0], 2: [1.0, 0.0, 0.0], 3: [1.0, 1.0, 0.0], 4: [0.0, 1.0, 0.0]},
                             'lines': {1: [1, 2], 2: [2, 3], 3: [3, 4], 4: [4, 1], 5: [1, 3]},
                             'surfaces': {1: [-5, 1, 2], 2: [4, 5, 3]},
                             'volumes': {},
                             'physical_groups': {'surface': {'geometry_ids': [1, 2], 'id': 1, 'ndim': 2},
                                                 'group2_surface': {'geometry_ids': [1, 2], 'id': 2, 'ndim': 2},
                                                 'new_line': {'geometry_ids': [5], 'id': 3, 'ndim': 1}}}

        # check if geo data is as expected
        TestUtils.assert_dictionary_almost_equal(filled_geo_data, expected_geo_data)

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()

        # extract geo data
        gmsh_io.extract_geo_data()
        filled_geo_data = gmsh_io.geo_data

        # check if geo data hasn't changed after re-synchronizing
        TestUtils.assert_dictionary_almost_equal(filled_geo_data, expected_geo_data)

    def test_synchronize_gmsh_with_intersection_line_through_multiple_surfaces(self):
        """
        Checks whether gmsh is synchronized after calling synchronize_gmsh. This test checks whether the geo data
        is updated with the created points and the physical group after calling synchronize_gmsh.

        This test is for the case where the intersection line is going through multiple surfaces, where each surface
        has its own physical group.
        """

        # initialize gmsh
        gmsh_io = GmshIO()
        gmsh.initialize()

        # create first surface
        # create surface points
        point_id1 = gmsh.model.occ.addPoint(0, 0, 0)
        point_id2 = gmsh.model.occ.addPoint(1, 0, 0)
        point_id3 = gmsh.model.occ.addPoint(1, 1, 0)
        point_id4 = gmsh.model.occ.addPoint(0, 1, 0)

        # create a surface
        line_id1 = gmsh.model.occ.addLine(point_id1, point_id2)
        line_id2 = gmsh.model.occ.addLine(point_id2, point_id3)
        line_id3 = gmsh.model.occ.addLine(point_id3, point_id4)
        line_id4 = gmsh.model.occ.addLine(point_id4, point_id1)

        curve_loop_id = gmsh.model.occ.addCurveLoop([line_id1, line_id2, line_id3, line_id4])
        surface_id: int = gmsh.model.occ.addPlaneSurface([curve_loop_id])
        gmsh.model.addPhysicalGroup(2, [surface_id], tag=-1, name="surface")

        # create second surface
        # create surface points
        point_id5 = gmsh.model.occ.addPoint(0, 2, 0)
        point_id6 = gmsh.model.occ.addPoint(1, 2, 0)

        # create a surface
        line_id5 = gmsh.model.occ.addLine(point_id4, point_id5)
        line_id6 = gmsh.model.occ.addLine(point_id5, point_id6)
        line_id7 = gmsh.model.occ.addLine(point_id6, point_id3)

        curve_loop_id2 = gmsh.model.occ.addCurveLoop([line_id3, line_id5, line_id6, line_id7])
        surface_id2: int = gmsh.model.occ.addPlaneSurface([curve_loop_id2])
        gmsh.model.addPhysicalGroup(2, [surface_id2], tag=-1, name="surface2")

        # create an intersection line
        point_id7 = gmsh.model.occ.addPoint(0.0, 0.0, 0)
        point_id8 = gmsh.model.occ.addPoint(1, 2.0, 0)
        line_id8 = gmsh.model.occ.addLine(point_id7, point_id8)
        gmsh.model.addPhysicalGroup(1, [line_id8], name="new_line")

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()

        # extract geo data
        gmsh_io.extract_geo_data()
        filled_geo_data = gmsh_io.geo_data

        expected_geo_data = {'points': {1: [0.0, 0.0, 0.0], 2: [1.0, 0.0, 0.0], 3: [1.0, 1.0, 0.0], 4: [0.0, 1.0, 0.0],
                                        5: [0.0, 2.0, 0.0], 6: [1.0, 2.0, 0.0], 7: [0.5, 1.0, 0.0]},
                             'lines': {1: [1, 2], 2: [2, 3], 4: [4, 1], 5: [4, 5], 6: [5, 6],
                                       7: [6, 3], 8: [1, 7], 9: [3, 7], 10: [7, 4], 11: [7, 6]},
                             'surfaces': {1: [-8, 1, 2, 9], 2: [4, 8, 10], 3: [7, 9, 11], 4: [-11, 10, 5, 6]},
                             'volumes': {},
                             'physical_groups': {'surface': {'geometry_ids': [1, 2], 'id': 1, 'ndim': 2},
                                                 'surface2': {'geometry_ids': [3, 4], 'id': 2, 'ndim': 2},
                                                 'new_line': {'geometry_ids': [8, 11], 'id': 3, 'ndim': 1}}}

        # check if geo data is as expected
        TestUtils.assert_dictionary_almost_equal(filled_geo_data, expected_geo_data)

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()

        # extract geo data
        gmsh_io.extract_geo_data()
        filled_geo_data = gmsh_io.geo_data

        # check if geo data hasn't changed after re-synchronizing
        TestUtils.assert_dictionary_almost_equal(filled_geo_data, expected_geo_data)

    def test_synchronize_gmsh_with_intersection_surface_through_volume(self):
        """
        Checks whether gmsh is synchronized after calling synchronize_gmsh. This test checks whether the geo data
        is updated with the created points and the physical group after calling synchronize_gmsh.

        This test is for the case where the intersection surface is going through a volume.
        """

        # initialize gmsh
        gmsh_io = GmshIO()
        gmsh.initialize()

        # create first surface
        # create surface points
        point_id1 = gmsh.model.occ.addPoint(0, 0, 0)
        point_id2 = gmsh.model.occ.addPoint(1, 0, 0)
        point_id3 = gmsh.model.occ.addPoint(1, 1, 0)
        point_id4 = gmsh.model.occ.addPoint(0, 1, 0)

        point_id5 = gmsh.model.occ.addPoint(0, 0, 1)
        point_id6 = gmsh.model.occ.addPoint(1, 0, 1)
        point_id7 = gmsh.model.occ.addPoint(1, 1, 1)
        point_id8 = gmsh.model.occ.addPoint(0, 1, 1)


        # create a volume
        line_id1 = gmsh.model.occ.addLine(point_id1, point_id2)
        line_id2 = gmsh.model.occ.addLine(point_id2, point_id3)
        line_id3 = gmsh.model.occ.addLine(point_id3, point_id4)
        line_id4 = gmsh.model.occ.addLine(point_id4, point_id1)

        line_id5 = gmsh.model.occ.addLine(point_id5, point_id6)
        line_id6 = gmsh.model.occ.addLine(point_id6, point_id7)
        line_id7 = gmsh.model.occ.addLine(point_id7, point_id8)
        line_id8 = gmsh.model.occ.addLine(point_id8, point_id5)

        line_id9 = gmsh.model.occ.addLine(point_id1, point_id5)
        line_id10 = gmsh.model.occ.addLine(point_id2, point_id6)
        line_id11 = gmsh.model.occ.addLine(point_id3, point_id7)
        line_id12 = gmsh.model.occ.addLine(point_id4, point_id8)

        curve_loop_id = gmsh.model.occ.addCurveLoop([line_id1, line_id2, line_id3, line_id4])
        curve_loop_id2 = gmsh.model.occ.addCurveLoop([line_id5, line_id6, line_id7, line_id8])
        curve_loop_id3 = gmsh.model.occ.addCurveLoop([line_id1, line_id9, line_id5, line_id10])
        curve_loop_id4 = gmsh.model.occ.addCurveLoop([line_id2, line_id10, line_id6, line_id11])
        curve_loop_id5 = gmsh.model.occ.addCurveLoop([line_id3, line_id11, line_id7, line_id12])
        curve_loop_id6 = gmsh.model.occ.addCurveLoop([line_id4, line_id12, line_id8, line_id9])

        surface_id_1: int = gmsh.model.occ.addPlaneSurface([curve_loop_id])
        surface_id_2: int = gmsh.model.occ.addPlaneSurface([curve_loop_id2])
        surface_id_3: int = gmsh.model.occ.addPlaneSurface([curve_loop_id3])
        surface_id_4: int = gmsh.model.occ.addPlaneSurface([curve_loop_id4])
        surface_id_5: int = gmsh.model.occ.addPlaneSurface([curve_loop_id5])
        surface_id_6: int = gmsh.model.occ.addPlaneSurface([curve_loop_id6])

        shell_id = gmsh.model.occ.addSurfaceLoop([surface_id_1, surface_id_2, surface_id_3,
                                                  surface_id_4, surface_id_5, surface_id_6])

        volume_id: int = gmsh.model.occ.addVolume([shell_id])

        gmsh.model.addPhysicalGroup(3, [volume_id], tag=-1, name="volume")

        # create intersection surface
        # create surface lines
        line_id13 = gmsh.model.occ.addLine(point_id1, point_id8)
        line_id14 = gmsh.model.occ.addLine(point_id2, point_id7)

        # create surface
        curve_loop_id7 = gmsh.model.occ.addCurveLoop([line_id1, line_id14, line_id7,line_id13])
        surface_id_7: int = gmsh.model.occ.addPlaneSurface([curve_loop_id7])
        gmsh.model.addPhysicalGroup(2, [surface_id_7], name="new_surface")

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()

        # extract geo data
        gmsh_io.extract_geo_data()
        filled_geo_data = gmsh_io.geo_data

        expected_geo_data = {'points': {1: [0.0, 0.0, 0.0], 2: [1.0, 0.0, 0.0], 3: [1.0, 1.0, 0.0], 4: [0.0, 1.0, 0.0],
                                        5: [0.0, 0.0, 1.0], 6: [1.0, 0.0, 1.0], 7: [1.0, 1.0, 1.0], 8: [0.0, 1.0, 1.0]},
                             'lines': {1: [1, 2], 2: [2, 3], 3: [3, 4], 4: [4, 1],
                                       5: [5, 6], 6: [6, 7], 7: [7, 8], 8: [8, 5],
                                       9: [1, 5], 10: [2, 6], 11: [3, 7], 12: [4, 8],
                                       13: [1, 8], 14: [2, 7]},
                             'surfaces': {1: [1, 2, 3, 4], 2: [5, 6, 7, 8], 3: [1, 10, -5, -9], 5: [3, 12, -7, -11],
                                          7: [1, 14, 7, -13], 8: [-14, 2, 11], 9: [-12, 4, 13], 10: [-10, 14, -6],
                                          11: [-13, 9, -8]},
                             'volumes': {1: [-1, 7, 8, 5, 9], 2: [3, -7, 10, 2, 11]},
                             'physical_groups': {'new_surface': {'geometry_ids': [7], 'id': 2, 'ndim': 2},
                                                 'volume': {'geometry_ids': [1, 2], 'id': 1, 'ndim': 3}}}

        # check if geo data is as expected
        TestUtils.assert_dictionary_almost_equal(filled_geo_data, expected_geo_data)

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()

        # extract geo data
        gmsh_io.extract_geo_data()
        filled_geo_data = gmsh_io.geo_data

        # check if geo data hasn't changed after re-synchronizing
        TestUtils.assert_dictionary_almost_equal(filled_geo_data, expected_geo_data)

    def test_reset_gmsh(self):
        """
        Checks whether gmsh is reset after calling reset_gmsh.
        """

        # initialize gmsh
        gmsh_io = GmshIO()
        gmsh.initialize()

        # create two points and add to physical group
        point_id1 = gmsh.model.occ.addPoint(0, 0, 0)
        point_id2 = gmsh.model.occ.addPoint(1, 0, 0)

        gmsh.model.addPhysicalGroup(0, [point_id1, point_id2], name="test")

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()

        # reset gmsh
        gmsh_io.reset_gmsh_instance()

        # extract geo data
        gmsh_io.extract_geo_data()
        empty_geo_data = gmsh_io.geo_data

        # check if geo data is empty after resetting
        expected_empty_geo_data = {'lines': {}, 'physical_groups': {}, 'points': {}, 'surfaces': {}, 'volumes': {}}

        TestUtils.assert_dictionary_almost_equal(empty_geo_data, expected_empty_geo_data)

    def test_clear_geo_data(self):
        """
        Checks whether geo data is cleared after calling clear_geo_data.
        """

        # initialize gmsh
        gmsh_io = GmshIO()
        gmsh.initialize()

        # create two points and add to physical group
        point_id1 = gmsh.model.occ.addPoint(0, 0, 0)
        point_id2 = gmsh.model.occ.addPoint(1, 0, 0)

        gmsh.model.addPhysicalGroup(0, [point_id1, point_id2], name="test")

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()

        # extract geo data
        gmsh_io.extract_geo_data()

        # clear geo data
        gmsh_io.clear_geo_data()
        empty_geo_data = gmsh_io.geo_data

        # check if geo data is empty after resetting
        expected_empty_geo_data = {'lines': {}, 'physical_groups': {}, 'points': {}, 'surfaces': {}, 'volumes': {}}
        TestUtils.assert_dictionary_almost_equal(empty_geo_data, expected_empty_geo_data)

        # check if geo data is present after re-extracting
        gmsh_io.extract_geo_data()
        filled_geo_data = gmsh_io.geo_data

        # check if geo data is filled after synchronizing
        expected_filled_geo_data = {'points': {1: [0., 0., 0.], 2: [1., 0., 0.]},
                                    'lines': {},
                                    'surfaces': {},
                                    'volumes': {},
                                    'physical_groups': {'test': {'geometry_ids': [1, 2], 'id': 1, 'ndim': 0}}}
        TestUtils.assert_dictionary_almost_equal(filled_geo_data, expected_filled_geo_data)

    def test_clear_mesh_data(self):
        """
        Checks whether mesh data is cleared after calling clear_mesh_data.
        """

        # initialize gmsh
        gmsh_io = GmshIO()
        gmsh.initialize()

        # create a line and add to physical group
        point_id1 = gmsh.model.occ.addPoint(0, 0, 0)
        point_id2 = gmsh.model.occ.addPoint(1, 0, 0)

        line = gmsh_io.create_line([point_id1, point_id2])

        gmsh.model.addPhysicalGroup(1, [line], name="test")

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()
        gmsh_io.extract_geo_data()

        # generate mesh
        gmsh_io.generate_mesh(1, element_size=0.5)

        expected_filled_mesh_data = {
            "ndim": 1,
            'elements': {'LINE_2N': {1: [1, 3],
                                     2: [3, 2]},
                         'POINT_1N': {3: [1],
                                      4: [2]}},
            'nodes': {1: [0., 0., 0.],
                      2: [1., 0., 0.],
                      3: [0.5, 0., 0.]},
            'physical_groups': {'test': {"ndim": 1, 'element_ids': [1, 2], "node_ids": [1, 2, 3], "element_type": "LINE_2N"}}}

        # check if mesh data is filled after generating mesh
        TestUtils.assert_dictionary_almost_equal(gmsh_io.mesh_data, expected_filled_mesh_data)

        # clear mesh data
        gmsh_io.clear_mesh_data()

        # check if mesh data is empty after clearing
        expected_empty_mesh_data = {}
        TestUtils.assert_dictionary_almost_equal(gmsh_io.mesh_data, expected_empty_mesh_data)

        # regenerate mesh
        gmsh_io.generate_mesh(1, element_size=0.5)
        TestUtils.assert_dictionary_almost_equal(gmsh_io.mesh_data, expected_filled_mesh_data)

    def test_make_geometry_0D(self):
        """
        Checks whether a 0D geometry is created correctly.
        """

        # define point coordinates
        point_coordinates = [(0, 0, 0), (1, 0, 0), (2, 1, 0)]

        # initialize gmsh
        gmsh_io = GmshIO()
        gmsh.initialize()

        # create multiple points and add to physical group
        gmsh_io.make_geometry_0d(point_coordinates, "point_group")

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()
        gmsh_io.extract_geo_data()

        # check if geo data is filled after synchronizing
        expected_filled_geo_data = {'points': {1: [0., 0., 0.],
                                               2: [1., 0., 0.],
                                               3: [2., 1., 0.]},
                                    'lines': {},
                                    'surfaces': {},
                                    'volumes': {},
                                    'physical_groups': {'point_group': {'geometry_ids': [1, 2, 3], 'id': 1, 'ndim': 0}}}

        TestUtils.assert_dictionary_almost_equal(gmsh_io.geo_data, expected_filled_geo_data)

    def test_make_geometry_1D(self):
        """
        Checks whether a 1D geometry is created correctly.
        """

        # define point coordinates
        point_coordinates = [(0, 0, 0), (1, 0, 0), (2, 1, 0)]

        # initialize gmsh
        gmsh_io = GmshIO()
        gmsh.initialize()

        # create multiple points and add to physical group
        gmsh_io.make_geometry_1d(point_coordinates, "line_group")

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()
        gmsh_io.extract_geo_data()

        # check if geo data is filled after synchronizing
        expected_filled_geo_data = {'points': {1: [0., 0., 0.],
                                               2: [1., 0., 0.],
                                               3: [2., 1., 0.]},
                                    'lines': {1: [1, 2], 2: [2, 3]},
                                    'surfaces': {},
                                    'volumes': {},
                                    'physical_groups': {'line_group': {'geometry_ids': [1, 2], 'id': 1, 'ndim': 1}}}

        TestUtils.assert_dictionary_almost_equal(gmsh_io.geo_data, expected_filled_geo_data)

    def test_make_geometry_2D(self):
        """
        Checks whether a 2D geometry is created correctly.
        """

        # define point coordinates
        point_coordinates = [(0, 0, 0), (1, 0, 0), (2, 1, 0)]

        # initialize gmsh
        gmsh_io = GmshIO()
        gmsh.initialize()

        # create multiple points and add to physical group
        gmsh_io.make_geometry_2d(point_coordinates, "surface_group")

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()
        gmsh_io.extract_geo_data()

        # check if geo data is filled after synchronizing
        expected_filled_geo_data = {'points': {1: [0., 0., 0.],
                                               2: [1., 0., 0.],
                                               3: [2., 1., 0.]},
                                    'lines': {1: [1, 2], 2: [2, 3], 3: [3, 1]},
                                    'surfaces': {1: [1, 2, 3]},
                                    'volumes': {},
                                    'physical_groups': {'surface_group': {'geometry_ids': [1], 'id': 1, 'ndim': 2}}}

        TestUtils.assert_dictionary_almost_equal(gmsh_io.geo_data, expected_filled_geo_data)

    def test_make_geometry_3D_by_extrusion(self):
        """
        Checks whether a 3D geometry is created correctly by extruding a 2D geometry.
        """

        # define point coordinates
        point_coordinates = [(0, 0, 0), (1, 0, 0), (2, 1, 0)]

        # initialize gmsh
        gmsh_io = GmshIO()
        gmsh.initialize()

        extrusion_length = [0, 0, 1]
        # create multiple points and add to physical group
        gmsh_io.make_geometry_3d_by_extrusion(point_coordinates, extrusion_length, "volume_group")

        # synchronize gmsh
        gmsh_io.synchronize_gmsh()
        gmsh_io.extract_geo_data()

        # check if geo data is filled after synchronizing
        expected_filled_geo_data = {'points': {1: [0., 0., 0.],
                                               2: [1., 0., 0.],
                                               3: [2., 1., 0.],
                                               4: [0., 0., 1.],
                                               5: [1., 0., 1.],
                                               6: [2., 1., 1.]},
                                    'lines': {1: [1, 2], 2: [2, 3], 3: [3, 1], 4: [1, 4], 5: [2, 5], 6: [4, 5],
                                              7: [3, 6], 8: [5, 6], 9: [6, 4]},
                                    'surfaces': {1: [1, 2, 3], 2: [4, 6, -5, -1], 3: [5, 8, -7, -2], 4: [7, 9, -4, -3],
                                                 5: [6, 8, 9]},
                                    'volumes': {1: [-2, -3, -4, -1, 5]},
                                    'physical_groups': {'volume_group': {'geometry_ids': [1], 'id': 1, 'ndim': 3}}}

        TestUtils.assert_dictionary_almost_equal(gmsh_io.geo_data, expected_filled_geo_data)

    def test_mesh_file_not_found(self):
        """
        Checks whether an error is raised if the mesh file is not found.

        """
        gmsh_io = GmshIO()
        pytest.raises(FileNotFoundError, gmsh_io.read_gmsh_msh, "not_existing_file.msh")

    def test_geo_file_not_found(self):
        """
        Checks whether an error is raised if the geo file is not found.

        """
        gmsh_io = GmshIO()
        pytest.raises(FileNotFoundError, gmsh_io.read_gmsh_geo, "not_existing_file.geo")

    def test_add_physical_group(self):
        """
        Checks whether physical groups are added correctly. This test checks whether the physical group is added
        correctly for a 0D geometry, 1D geometry, 2D geometry and 3D geometry.

        """

        # initialize gmsh and create a 3D geometry
        gmsh_io = GmshIO()
        gmsh.initialize()
        gmsh_io.make_geometry_3d_by_extrusion([(0, 0, 0), (1, 0, 0), (2, 1, 0)],[0,0,1] , "volume_group")

        # add physical group for all dimensions
        gmsh_io.add_physical_group("new_volume_group", 3, [1])
        gmsh_io.add_physical_group("new_surface_group", 2, [1])
        gmsh_io.add_physical_group("new_line_group", 1, [1])
        gmsh_io.add_physical_group("new_point_group", 0, [1])

        geo_data = gmsh_io.geo_data

        # set expected data
        expected_group_data = {"new_volume_group": {"geometry_ids": [1], "id": 2, "ndim": 3},
                               "new_surface_group": {"geometry_ids": [1], "id": 3, "ndim": 2},
                               "new_line_group": {"geometry_ids": [1], "id": 4, "ndim": 1},
                               "new_point_group": {"geometry_ids": [1], "id": 5, "ndim": 0},
                               "volume_group": {"geometry_ids": [1], "id": 1, "ndim": 3}}

        # check if physical groups are added correctly
        TestUtils.assert_dictionary_almost_equal(geo_data["physical_groups"], expected_group_data)

    def test_add_point_at_surface_point(self):
        """
        Checks whether a point is added at the surface point correctly. Both the physical groups of the surface and the
        point should be maintained.
        """

        gmsh_io = GmshIO()

        # create surface input
        input_surface = {'surface': {"coordinates": [(0, 0, 0), (3, 0, 0), (3, 1, 0), (0, 1, 0)],
                                     "ndim": 2}}

        # create point input
        input_point = {'point': {"coordinates": [(0, 0, 0)],
                                 "ndim": 0}}

        # generate point after surface
        gmsh_io.generate_geometry(input_surface, "")
        gmsh_io.generate_geometry(input_point, "")

        output_group_names = list(gmsh_io.geo_data["physical_groups"].keys())
        expected_group_names = ['surface', 'point']

        # check if all groups are added
        for group_name in expected_group_names:
            assert group_name in output_group_names

        assert gmsh_io.geo_data["physical_groups"]["surface"]["geometry_ids"] == [1]
        assert gmsh_io.geo_data["physical_groups"]["point"]["geometry_ids"] == [1]

    def test_add_multiple_points_at_surface_points(self):
        """
        Checks whether multiple points are added at the surface points correctly. Both the physical groups of the
        surface and the points should be maintained.
        """

        gmsh_io = GmshIO()

        # create surface input
        input_surface = {'surface': {"coordinates": [(0, 0, 0), (3, 0, 0), (3, 1, 0), (0, 1, 0)],
                                     "ndim": 2}}

        # create point input
        input_point = {'point': {"coordinates": [(3, 1, 0), (0, 1, 0)],
                                 "ndim": 0}}

        # generate point after surface
        gmsh_io.generate_geometry(input_surface, "")
        gmsh_io.generate_geometry(input_point, "")

        output_group_names = list(gmsh_io.geo_data["physical_groups"].keys())
        expected_group_names = ['surface', 'point']

        # check if all groups are added
        for group_name in expected_group_names:
            assert group_name in output_group_names

        assert gmsh_io.geo_data["physical_groups"]["surface"]["geometry_ids"] == [1]
        assert gmsh_io.geo_data["physical_groups"]["point"]["geometry_ids"] == [3, 4]

    def test_add_multiple_point_groups_to_surface_point(self):
        """
        Checks whether multiple points each in a different group added at the surface same points correctly.
        All the physical groups of the surface and the points should be maintained.

        """

        gmsh_io = GmshIO()

        # create surface input
        input_surface = {'surface': {"coordinates": [(0, 0, 0), (3, 0, 0), (3, 1, 0), (0, 1, 0)],
                                     "ndim": 2}}

        # create point input
        input_point = {'point1': {"coordinates": [(3, 1, 0)],
                                  "ndim": 0}}

        # create point input
        input_point2 = {'point2': {"coordinates": [(3, 1, 0)],
                                   "ndim": 0}}

        # generate point after surface
        gmsh_io.generate_geometry(input_surface, "")
        gmsh_io.generate_geometry(input_point, "")
        gmsh_io.generate_geometry(input_point2, "")

        output_group_names = list(gmsh_io.geo_data["physical_groups"].keys())
        expected_group_names = ['surface', 'point1','point2']

        # check if all groups are added
        for group_name in expected_group_names:
            assert group_name in output_group_names

        assert gmsh_io.geo_data["physical_groups"]["surface"]["geometry_ids"] == [1]
        assert gmsh_io.geo_data["physical_groups"]["point1"]["geometry_ids"] == [3]
        assert gmsh_io.geo_data["physical_groups"]["point2"]["geometry_ids"] == [3]

    def test_validate_layer_parameters(self, capfd: CaptureFixture[str]):
        """
        Checks whether the layer parameters are validated correctly. And errors should be raised if the layer parameters
        are not valid.

        Args:
            - capfd: Pytest fixture to capture stdout and stderr.

        """

        gmsh_io = GmshIO()

        # valid layer parameters, no error should be raised
        layer_parameters = {"layer_1": {"coordinates": [[0,0,0], [1,1,1]],
                                        "element_size": 1,
                                        "ndim": 1,
                                        "extrusion_length": 1}}
        gmsh_io.validate_layer_parameters(layer_parameters)

        # missing coordinates, error should be raised
        layer_parameters = {"layer_1": {"element_size": 1,
                                        "ndim": 1,
                                        "extrusion_length": 1}}

        with pytest.raises(ValueError, match=r"Layer layer_1 must contain the key 'coordinates'"):
            gmsh_io.validate_layer_parameters(layer_parameters)

        # missing ndim, error should be raised
        layer_parameters = {"layer_1": {"coordinates": [[0,0,0], [1,1,1]],
                                        "element_size": 1,
                                        "extrusion_length": 1}}

        with pytest.raises(ValueError, match=r"Layer layer_1 must contain the key 'ndim'"):
            gmsh_io.validate_layer_parameters(layer_parameters)

        # missing extrusion length in 3D, error should be raised
        layer_parameters = {"layer_1": {"coordinates": [[0,0,0], [1,1,1]],
                                        "element_size": 1,
                                        "ndim": 3}}

        with pytest.raises(ValueError, match=r"Layer layer_1 must contain the key 'extrusion_length', "
                                             r"which is needed for 3D geometries"):
            gmsh_io.validate_layer_parameters(layer_parameters)

        # missing extrusion length in non-3D, should not raise an error
        layer_parameters = {"layer_1": {"coordinates": [[0,0,0], [1,1,1]],
                                        "element_size": 1,
                                        "ndim": 2}}

        gmsh_io.validate_layer_parameters(layer_parameters)

        # missing element size, warning should be printed
        layer_parameters = {"layer_1": {"coordinates": [[0,0,0], [1,1,1]],
                                        "ndim": 1,
                                        "extrusion_length": 1}}

        gmsh_io.validate_layer_parameters(layer_parameters)
        console_output, _ = capfd.readouterr()

        assert ("Warning: Layer layer_1 does not contain the key 'element_size'. "
                "The element size will be determined by gmsh.") == console_output.strip()

        # non supported ndim value, error should be raised
        layer_parameters = {"layer_1": {"coordinates": [[0, 0, 0], [1, 1, 1]],
                                        "element_size": 1,
                                        "extrusion_length": 1,
                                        "ndim": 4}}

        with pytest.raises(ValueError, match=f"ndim must be 0, 1, 2 or 3. ndim=4"):
            gmsh_io.validate_layer_parameters(layer_parameters)

    def test_write_mesh(self):
        """
        Checks whether the mesh is written correctly to a file.
        """

        # create line geometry
        layer_parameters = {"line": {"coordinates": [[0, 0, 0], [1, 1, 1]],
                                     "element_size": 1,
                                     "ndim": 1}}

        gmsh_io = GmshIO()

        gmsh_io.generate_geometry(layer_parameters, "")

        gmsh_io.generate_extract_mesh(2, "test_mesh_line.msh", ".", True, False)

        # check if file is created
        assert Path("test_mesh_line.msh").is_file()

        # open generated mesh file
        with open("test_mesh_line.msh", "r") as file:
            generate_mesh = file.readlines()

        # open expected mesh file
        with open("tests/test_data/expected_mesh_line.msh", "r") as file:
            expected_mesh = file.readlines()

        # check if generated mesh file is equal to expected mesh file
        for generated_line, expected_line in zip(generate_mesh, expected_mesh):
            temp = re.findall(r'\d+', generated_line)
            generated_numbers = list(map(int, temp))

            temp = re.findall(r'\d+', expected_line)
            expected_numbers = list(map(int, temp))

            # assert string if line does not contain numbers, else assert numbers
            if len(expected_numbers) == 0:
                assert generated_line == expected_line
            else:
                np.testing.assert_array_almost_equal(generated_numbers, expected_numbers)

        # remove test file
        Path("test_mesh_line.msh").unlink()

    def test_two_point_groups_with_same_name(self):
        """
        Checks whether the points are added correctly to the existing physical group.
        """

        gmsh_io = GmshIO()

        # create first point input
        input_first_point = {'point': {"coordinates": [(0, 0, 0), (1,0,0), (0,1,0)],
                                       "ndim": 0}}

        # create second point input, note that one coordinate already exists
        input_second_point = {'point': {"coordinates": [(3, 0, 0), (2,0,0), (0,1,0)],
                                        "ndim": 0}}

        # generate points separately
        gmsh_io.generate_geometry(input_first_point, "")
        gmsh_io.generate_geometry(input_second_point, "")

        geo_data = gmsh_io.geo_data

        # set expected data
        expected_group_data = {"point": {"geometry_ids": [1, 2, 3, 4, 5], "id": 1, "ndim": 0}}

        TestUtils.assert_dictionary_almost_equal(geo_data["physical_groups"], expected_group_data)

    def test_two_line_groups_with_same_name(self):
        """
        Checks whether the lines are added correctly to the existing physical group.
        """

        gmsh_io = GmshIO()

        # create first line input
        input_first_line = {'line': {"coordinates": [(0, 0, 0), (3, 0, 0)],
                                     "ndim": 1}}

        # create second line input
        input_second_line = {'line': {"coordinates": [(0, 0, 1), (3, 0, 1)],
                                      "ndim": 1}}

        # generate lines separately
        gmsh_io.generate_geometry(input_first_line, "")
        gmsh_io.generate_geometry(input_second_line, "")

        geo_data = gmsh_io.geo_data

        # set expected data
        expected_group_data = {"line": {"geometry_ids": [1, 2], "id": 1, "ndim": 1}}

        TestUtils.assert_dictionary_almost_equal(geo_data["physical_groups"], expected_group_data)

    def test_two_surface_groups_with_same_name(self):
        """
        Checks whether the surfaces are added correctly to the existing physical group.
        """

        gmsh_io = GmshIO()

        # create first surface input
        input_first_surface = {'surface': {"coordinates": [(0, 0, 0), (3, 0, 0), (3, 1, 0), (0, 1, 0)],
                                           "ndim": 2}}

        # create second surface input
        input_second_surface = {'surface': {"coordinates": [(0, 0, 1), (3, 0, 1), (3, 1, 1), (0, 1, 1)],
                                            "ndim": 2}}

        # generate surfaces separately
        gmsh_io.generate_geometry(input_first_surface, "")
        gmsh_io.generate_geometry(input_second_surface, "")

        geo_data = gmsh_io.geo_data

        # set expected data
        expected_group_data = {"surface": {"geometry_ids": [1, 2], "id": 1, "ndim": 2}}

        TestUtils.assert_dictionary_almost_equal(geo_data["physical_groups"], expected_group_data)

    def test_two_volume_groups_with_same_name(self):
        """
        Checks whether the volumes are added correctly to the existing physical group.
        """

        gmsh_io = GmshIO()

        # create first volume input
        input_first_volume = {'volume': {"coordinates": [(0, 0, 0), (3, 0, 0), (3, 1, 0), (0, 1, 0)],
                                            "ndim": 3,
                                            "extrusion_length": [0, 0, 1]}}

        # create second volume input
        input_second_volume = {'volume': {"coordinates": [(0, 0, 2), (3, 0, 2), (3, 1, 2), (0, 1, 2)],
                                            "ndim": 3,
                                            "extrusion_length": [0, 0, 1]}}

        # generate volumes separately
        gmsh_io.generate_geometry(input_first_volume, "")
        gmsh_io.generate_geometry(input_second_volume, "")

        geo_data = gmsh_io.geo_data

        # set expected data
        expected_group_data = {"volume": {"geometry_ids": [1, 2], "id": 1, "ndim": 3}}

        TestUtils.assert_dictionary_almost_equal(geo_data["physical_groups"], expected_group_data)

    def test_two_different_dimension_groups_with_same_name(self):
        """
        Checks whether an exception is raised that all items in a physical group must have the same dimension.
        """

        gmsh_io = GmshIO()

        # create first line input
        input_first_line = {'group_1': {"coordinates": [(0, 0, 0), (3, 0, 0)],
                                        "ndim": 1}}

        # create first point input
        input_first_point = {'group_1': {"coordinates": [(0, 0, 0)],
                                         "ndim": 0}}

        # generate groups separately, exception is raised
        gmsh_io.generate_geometry(input_first_line, "")
        with pytest.raises(ValueError, match=r"Cannot add geometry ids to physical group group_1 with dimension 0 as "
                                             r"the physical group already exists with dimension 1."):
            gmsh_io.generate_geometry(input_first_point, "")











