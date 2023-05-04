from gmsh_utils.gmsh_IO import GmshIO
import numpy as np
import numpy.testing as npt
class TestGmshIO:
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

    def test_read_gmsh_geo_3D(self):
        """
        Checks whether a gmsh .geo file is read correctly.
        """
        geo_file = r"test_data/column_3D_tetra4.geo"

        gmsh_io = GmshIO()
        gmsh_io.read_gmsh_geo(geo_file)

        geo_data = gmsh_io.geo_data

        # check if the number of points is correct
        assert len(geo_data["points"]) == 12

        # check if the coordinates of the points are correct
        point_coordinates = np.array(list(geo_data["points"].values()))
        min_x, max_x = np.min(point_coordinates[:, 0]), np.max(point_coordinates[:, 0])
        min_y, max_y = np.min(point_coordinates[:, 1]), np.max(point_coordinates[:, 1])
        min_z, max_z = np.min(point_coordinates[:, 2]), np.max(point_coordinates[:, 2])

        npt.assert_allclose([min_x, max_x], [0, 0.5])
        npt.assert_allclose([min_y, max_y], [0, 2])
        npt.assert_allclose([min_z, max_z], [-0.5, 0.0])

        # check if the number of lines is correct
        assert len(geo_data["lines"]) == 20

        # check if the number of surfaces is correct
        assert len(geo_data["surfaces"]) == 11

        # check if the number of volumes is correct
        assert len(geo_data["volumes"]) == 2

        # check if the number of physical groups is correct
        assert len(geo_data["physical_groups"]) == 2

        # check if geometry ids are correct
        group_1_id = geo_data["physical_groups"]["group_1"]["geometry_id"]
        group_2_id = geo_data["physical_groups"]["group_2"]["geometry_id"]

        assert group_1_id == 1
        assert group_2_id == 1

        # check if group dimension is correct
        assert geo_data["physical_groups"]["group_1"]["ndim"] == 3

        # get surface ids of group 1, any direction is fine
        group_1_surface_ids = np.abs(geo_data["volumes"][group_1_id])

        # get line ids of group 1, any direction is fine
        group_1_line_ids = np.abs([geo_data["surfaces"][id] for id in group_1_surface_ids])

        # check if number of lines of the volume in group 1 is correct
        assert group_1_line_ids.shape == (6, 4)





