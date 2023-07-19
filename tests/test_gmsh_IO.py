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

        expected_mesh_data = {'nodes': {1: [0., 0., 0.],
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
        expected_mesh_data = {'nodes': {1: [0., 0., 0.],
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
                              'physical_groups': {'group_1': {"node_ids": [1, 2, 3, 4, 5],
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
        expected_groups_in_mesh_data = {'group_1': {"node_ids": [1, 2, 3, 4, 7],
                                                    "element_ids": [5, 6, 7, 8],
                                                    "element_type": "TRIANGLE_3N"},
                                        'group_2': {"node_ids": [3, 4, 5, 6, 8],
                                                    "element_ids": [9, 10, 11, 12],
                                                    "element_type": "TRIANGLE_3N"},
                                        "combined_group": {"node_ids": [1, 2, 3, 4, 5, 6, 7, 8],
                                                           "element_ids": [5, 6, 7, 8, 9, 10, 11, 12],
                                                           "element_type": "TRIANGLE_3N"},
                                        "line_group": {"node_ids": [1, 2, 3, 6],
                                                       "element_ids": [3, 4],
                                                       "element_type": "LINE_2N"},
                                        "point_group": {"node_ids": [1, 2],
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
            'elements': {'LINE_2N': {1: [1, 3],
                                     2: [3, 2]},
                         'POINT_1N': {3: [1],
                                      4: [2]}},
            'nodes': {1: [0., 0., 0.],
                      2: [1., 0., 0.],
                      3: [0.5, 0., 0.]},
            'physical_groups': {'test': {'element_ids': [1, 2], "node_ids": [1, 2, 3], "element_type": "LINE_2N"}}}

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