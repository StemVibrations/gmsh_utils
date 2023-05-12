# todo group input into geometry, settings, ...
from gmsh_utils.gmsh_IO import GmshIO

number_of_layers = 6
default_mesh_size = 1
# define the points of the surface as a list of tuples
input_points_list = [[(0, 0, 0), (3, 0, 0), (3, 1, 0),  (0, 1, 0)],
                     [(3, 0, 0), (5, 0, 0), (5, 1, 0), (4, 1.5, 0), (3, 1, 0)],
                     [(0, 1, 0), (2, 1, 0),  (2, 3, 0), (0, 3, 0)],
                     [(2, 1, 0), (3, 1, 0), (4, 1.5, 0), (5, 1, 0),  (5, 3, 0), (2, 3, 0)],
                     [(0, 3, 0), (2.5, 3, 0), (2, 4, 0), (0, 4, 0)],
                     [(0.8, 4, 0), (1.2, 4, 0), (1.2, 4.1, 0), (0.8, 4.1, 0)]]
name_label_list = ["First Soil Layer", "FSL", "Second Soil Layer", "SSL", "Soil Ballast", "Line Track"]
# # define the element size
# element_size = 2
# define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
dims = 3
# if 3D, input depth of geometry to be extruded from 2D surface
extrusion_length = [0, 0, 3]
# set a name label for the surface

# if "True", saves mesh data to separate mdpa files; otherwise "False"
save_file = False
# if "True", opens gmsh interface; otherwise "False"
open_gmsh_gui = True
# set a name for mesh output file
mesh_output_name = "geometry"
# set output directory
mesh_output_dir = "./"


gmsh_io = GmshIO()

gmsh_io.generate_gmsh_mesh(input_points_list, number_of_layers, extrusion_length, default_mesh_size, dims,
                           name_label_list, mesh_output_name, mesh_output_dir,
                           save_file, open_gmsh_gui)

mesh_data = gmsh_io.mesh_data
