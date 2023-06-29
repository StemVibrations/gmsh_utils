# todo group input into geometry, settings, ...
from gmsh_utils.gmsh_IO import GmshIO

# define the default mesh size, if -1, the mesh size is logically chosen by Gmsh itself based on the geometry
default_mesh_size: float = -1
# define the points of the surface and mesh sizes as a dictionary
input_dict = {'Second left Soil Layer': (default_mesh_size, [(3, 0, 0), (5, 0, 0), (5, 1.5, 0)]),
              'First Soil Layer': (default_mesh_size, [(0, 0, 0), (3, 0, 0), (5, 1.5, 0), (2, 1, 0), (0, 1, 0)]),
              'Third left Soil Layer': (default_mesh_size, [(0, 1, 0), (0, 3, 0), (2, 3, 0), (2, 1, 0)])}

input_points_list = []
mesh_size_list = []
name_label_list = []
number_of_layers = len(input_dict)
for value in input_dict.values():
    input_points_list.append(list(value[1]))
# Directly access the dictionary keys
keys = input_dict.keys()
# Print the keys
for key in keys:
    name_label_list.append(key)  # Extract the name label

# define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
dims = 2
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

gmsh_io.generate_geometry(input_points_list, extrusion_length, dims,
                          mesh_output_name, name_label_list, default_mesh_size)
gmsh_io.generate_extract_mesh(dims, mesh_output_name, mesh_output_dir, save_file, open_gmsh_gui)
mesh_data = gmsh_io.mesh_data
