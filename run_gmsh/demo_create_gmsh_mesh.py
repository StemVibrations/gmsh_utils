# todo group input into geometry, settings, ...
from gmsh_utils.gmsh_IO import GmshIO

# define a global mesh size, if set to -1 mesh size is logically chosen by Gmsh itself based on the geometry
global_mesh_size: float = -1
# define the name labels of the layers and points coordinates of the surface in anticlockwise order as a dictionary.
# The global_mesh_size is assigned to all layers if a specific mesh size is not given in the function
# set_mesh_size_of_group below
input_dict = {'soil_1': {"mesh_size": global_mesh_size,
                         "coordinates": [(0, 0, 0), (3, 0, 0), (5, 1.5, 0), (2, 1, 0), (0, 1, 0)]},
              'soil_2': {"mesh_size": global_mesh_size,
                         "coordinates": [(3, 0, 0), (5, 0, 0), (5, 1.5, 0)]},
              'soil_3': {"mesh_size": global_mesh_size,
                         "coordinates": [(0, 1, 0), (0, 3, 0), (2, 3, 0), (2, 1, 0)]}}

# define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
dims = 2
# if 3D, input depth of geometry to be extruded from 2D surface
extrusion_length = [0, 0, 3]

# define the name labels of the layers and points coordinates of the surface in order /
# (regardless of clockwise or anticlockwise) and mesh sizes for each layer as a dictionary
input_dict = {'First left Soil Layer': {"element_size": default_mesh_size,
                                        "coordinates": [(0, 0, 0), (3, 0, 0), (5, 1.5, 0), (2, 1, 0), (0, 1, 0)],
                                        "ndim": dims,
                                        "extrusion_length": extrusion_length},
              'Second right Soil Layer': {"element_size": default_mesh_size,
                                          "coordinates": [(3, 0, 0), (5, 0, 0), (5, 1.5, 0)],
                                          "ndim": dims,
                                          "extrusion_length": extrusion_length},
              'Third top Soil Layer': {"element_size": default_mesh_size,
                                       "coordinates": [(0, 1, 0), (0, 3, 0), (2, 3, 0), (2, 1, 0)],
                                       "ndim": dims,
                                       "extrusion_length": extrusion_length}}

# if "True", saves mesh data to separate mdpa files; otherwise "False"
save_file = True
# if "True", opens gmsh interface; otherwise "False"
open_gmsh_gui = True
# set a name for mesh output file
mesh_output_name = "geometry"
# set output directory
mesh_output_dir = "./"


gmsh_io = GmshIO()
coordinates_list, name_label_list = gmsh_io.input_dict_to_list(input_dict)
gmsh_io.generate_geometry(coordinates_list, name_label_list, extrusion_length, dims, mesh_output_name, global_mesh_size)

gmsh_io.generate_geometry(input_dict, mesh_output_name)
# set mesh size of a group by defining the name label of the group and the desired mesh size
gmsh_io.set_mesh_size_of_group("soil_1", 0.1)

gmsh_io.generate_extract_mesh(dims, mesh_output_name, mesh_output_dir, save_file, open_gmsh_gui)

mesh_data = gmsh_io.mesh_data
