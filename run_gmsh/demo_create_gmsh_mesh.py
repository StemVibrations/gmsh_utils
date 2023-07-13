# todo group input into geometry, settings, ...
from gmsh_utils.gmsh_IO import GmshIO
import warnings

# define a default mesh size, if set to -1 mesh size is logically chosen by Gmsh itself based on the geometry
default_mesh_size: float = -1
# define the name labels of the layers and points coordinates of the surface in order /
# (regardless of clockwise or anticlockwise) and mesh sizes for each layer as a dictionary. If arbitrary_mesh_size is
# false, the default_mesh_size is assigned to all layers.
input_dict = {'First left Soil Layer': {"mesh_size": default_mesh_size,
                                        "coordinates": [(0, 0, 0), (3, 0, 0), (5, 1.5, 0), (2, 1, 0), (0, 1, 0)]},
              'Second right Soil Layer': {"mesh_size": default_mesh_size,
                                          "coordinates": [(3, 0, 0), (5, 0, 0), (5, 1.5, 0)]},
              'Third top Soil Layer': {"mesh_size": 0.1,
                                       "coordinates": [(0, 1, 0), (0, 3, 0), (2, 3, 0), (2, 1, 0)]}}

# input_dict = {'First left Soil Layer': {"mesh_size": default_mesh_size,
#                                         "coordinates": [(0, 0, 0), (3, 0, 0), (5, 1.5, 0), (2, 1, 0), (0, 1, 0)]},
#               'Second right Soil Layer': (0.1, [(3, 0, 0), (5, 0, 0), (5, 1.5, 0)]),
#               'Third top Soil Layer': (0.2, [(0, 1, 0), (0, 3, 0), (2, 3, 0), (2, 1, 0)])}
# define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
dims = 3
# if 3D, input depth of geometry to be extruded from 2D surface
extrusion_length = [0, 0, 3]
# if "True", saves mesh data to separate mdpa files; otherwise "False"
save_file = True
# if "True", opens gmsh interface; otherwise "False"
open_gmsh_gui = True
# set a name for mesh output file
mesh_output_name = "geometry"
# set output directory
mesh_output_dir = "./"


# Extract the points, mesh sizes and name labels from the dictionary
input_points_list = []
mesh_size_list = []
name_label_list = []
for value in input_dict.values():
    mesh_size_list.append(value["mesh_size"])
    input_points_list.append(value["coordinates"])

# for value in input_dict.values():
#     input_points_list.append(list(value[1]))  # Extract the points
#     if value[0] <= 0:  # Check if mesh size is negative
#         warnings.warn('Warning Message: Mesh size cannot be negative!')
#         mesh_size = 1
#         warnings.warn('Warning Message: Mesh size changed to 1!')
#     else:
#         mesh_size = value[0] if value[0] != default_mesh_size else default_mesh_size
#     mesh_size_list.append(mesh_size)  # Extract the mesh sizes

# Directly access the dictionary keys
keys = input_dict.keys()
# Print the keys
for key in keys:
    name_label_list.append(key)  # Extract the name label


gmsh_io = GmshIO()

gmsh_io.generate_geometry(input_points_list, extrusion_length, dims, mesh_output_name,
                          name_label_list, default_mesh_size)

# for group in groups:
#
#     gmsh_io.set_mesh_size_of_group(group, 5)
#     gmsh_io.set_mesh_size_of_group("soil_1", 5)

gmsh_io.generate_extract_mesh(dims, mesh_output_name, mesh_output_dir, mesh_size_list, save_file,
                              open_gmsh_gui)
mesh_data = gmsh_io.mesh_data
