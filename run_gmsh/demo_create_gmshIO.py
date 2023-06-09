# todo group input into geometry, settings, ...
from gmsh_utils.gmsh_IO import GmshIO

# define the default mesh size
default_mesh_size = 0.3
# define the points of the surface and mesh sizes as a dictionary
input_dict = {'First Soil Layer': (5, [(0, 0, 0), (3, 0, 0), (5, 1.5, 0), (2, 1, 0), (0, 1, 0)]),
              'Second left Soil Layer': (1, [(0, 1, 0), (2, 1, 0), (2, 3, 0), (0, 3, 0)])}

input_points_list = []
mesh_size_list = []
name_label_list = []
number_of_layers = len(input_dict)
for value in input_dict.values():
    mesh_size_list.append(value[0])  # Extract the mesh size
    input_points_list.append(value[1])
# Directly access the dictionary keys
keys = input_dict.keys()
# Print the keys
for key in keys:
    name_label_list.append(key)  # Extract the name label

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


gmsh_io = GmshIO()

gmsh_io.generate_gmsh_mesh(input_points_list, extrusion_length, default_mesh_size, dims, name_label_list,
                           mesh_output_name, mesh_output_dir, save_file, open_gmsh_gui)

mesh_data = gmsh_io.mesh_data
