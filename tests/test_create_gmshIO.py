import pytest
from gmsh_utils.gmsh_IO import GmshIO
import numpy as np
import pytest

def test_generate_mesh_2d():
    """
    Checks whether mesh data generated for 2D geometries is not empty.

    """
    # define the default mesh size
    default_mesh_size = 1
    # define the points of the surface as a list of tuples
    input_points_list = [[(0, 0, 0), (3, 0, 0), (3, 1, 0), (0, 1, 0)],
                         [(3, 0, 0), (5, 0, 0), (5, 1, 0), (4, 1.5, 0), (3, 1, 0)],
                         [(0, 1, 0), (2, 1, 0), (2, 3, 0), (0, 3, 0)],
                         [(2, 1, 0), (3, 1, 0), (4, 1.5, 0), (5, 1, 0), (5, 3, 0), (2, 3, 0)],
                         [(0, 3, 0), (2.5, 3, 0), (2, 4, 0), (0, 4, 0)],
                         [(0.8, 4, 0), (1.2, 4, 0), (1.2, 4.1, 0), (0.8, 4.1, 0)]]
    # define the name labels for the surfaces
    name_label_list = ["First Soil Layer", "FSL", "Second Soil Layer", "SSL", "Soil Ballast", "Line Track"]

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
    mesh_output_dir = "./"

    gmsh_io = GmshIO()

    gmsh_io.generate_gmsh_mesh(input_points_list, extrusion_length, default_mesh_size, dims,
                               name_label_list, mesh_output_name, mesh_output_dir,
                               save_file, open_gmsh_gui)

    mesh_data = gmsh_io.mesh_data

    assert mesh_data["nodes"]["coordinates"].size > 0  # check if node_coords is not empty
    assert mesh_data["nodes"]["ids"].size > 0  # check if node_tags is not empty
    assert list(mesh_data["elements"].keys()) == ["LINE_2N", "TRIANGLE_3N",
                                                  "POINT_1N"]  # check if correct elements are present

    # check each element type contains ids and nodes
    for value in mesh_data["elements"].values():
        assert value["element_ids"].size > 0
        assert value["connectivities"].size > 0

#
def test_generate_mesh_3d():
    """
    Checks whether mesh data generated for 3D geometries is not empty.

    """

    # define the default mesh size
    default_mesh_size = 1
    # define the points of the surface as a list of tuples
    input_points_list = [[(0, 0, 0), (3, 0, 0), (3, 1, 0), (0, 1, 0)],
                         [(3, 0, 0), (5, 0, 0), (5, 1, 0), (4, 1.5, 0), (3, 1, 0)],
                         [(0, 1, 0), (2, 1, 0), (2, 3, 0), (0, 3, 0)],
                         [(2, 1, 0), (3, 1, 0), (4, 1.5, 0), (5, 1, 0), (5, 3, 0), (2, 3, 0)],
                         [(0, 3, 0), (2.5, 3, 0), (2, 4, 0), (0, 4, 0)],
                         [(0.8, 4, 0), (1.2, 4, 0), (1.2, 4.1, 0), (0.8, 4.1, 0)]]
    # define the name labels for the surfaces
    name_label_list = ["First Soil Layer", "FSL", "Second Soil Layer", "SSL", "Soil Ballast", "Line Track"]

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
    mesh_output_dir = "./"

    gmsh_io = GmshIO()

    gmsh_io.generate_gmsh_mesh(input_points_list, extrusion_length, default_mesh_size, dims,
                               name_label_list, mesh_output_name, mesh_output_dir,
                               save_file, open_gmsh_gui)

    mesh_data = gmsh_io.mesh_data


    assert mesh_data["nodes"]["coordinates"].size > 0  # check if node_coords is not empty
    assert mesh_data["nodes"]["ids"].size > 0  # check if node_tags is not empty
    assert list(mesh_data["elements"].keys()) == ["LINE_2N", "TRIANGLE_3N", "TETRAHEDRON_4N",
                                                  "POINT_1N"]  # check if correct elements are present

    # check each element type contains ids and nodes
    for value in mesh_data["elements"].values():
        assert value["element_ids"].size > 0
        assert value["connectivities"].size > 0
