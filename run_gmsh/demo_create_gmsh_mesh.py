# todo group input into geometry, settings, ...
from gmsh_utils.gmsh_IO import GmshIO

# define the default mesh size, if -1, the mesh size is logically chosen by Gmsh itself based on the geometry
default_mesh_size: float = -1

# define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
dims = 3

# if 3D, input depth of geometry to be extruded from 2D surface
extrusion_length = [0, 0, 3]

# define the name labels of the layers; points coordinates of the surface in clockwise or anti-clockwise order;
# and mesh sizes for each layer as a dictionary
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

gmsh_io.generate_geometry(input_dict, mesh_output_name)
gmsh_io.generate_extract_mesh(dims, mesh_output_name, mesh_output_dir, save_file, open_gmsh_gui)
mesh_data = gmsh_io.mesh_data
