import gmsh
import embankment_functions as efn


inputs = efn.InputGenerator(number_of_layers=7,
                            mesh_size_list=[2, 0.1, 0.5, 5, 0.5, 0.1, 0.01],
                            input_points_list=[[(0, 0, 0), (3, 0, 0), (3, 1, 0),  (0, 1, 0)],
                                               [(3, 0, 0), (5, 0, 0), (5, 1, 0), (4, 1.5, 0), (3, 1, 0)],
                                               [(0, 1, 0), (2, 1, 0),  (2, 3, 0), (0, 3, 0)],
                                               [(2, 1, 0), (3, 1, 0), (4, 1.5, 0), (5, 1, 0),  (5, 3, 0), (2, 3, 0)],
                                               [(0, 3, 0), (2.5, 3, 0), (2, 4, 0), (0, 4, 0)],
                                               [(2.5, 3, 0), (4, 3, 0), (4, 4, 0), (2, 4, 0)],
                                               [(0.8, 4, 0), (1.2, 4, 0), (1.2, 4.1, 0), (0.8, 4.1, 0)]],
                            name_label_list=["First Soil Layer", "FSL", "Second Soil Layer", "SSL", "Soil Ballast",
                                             "Soil Ballast2", "Line Track"],
                            eps=1e-3,
                            default_mesh_size=0.5, dims=3, depth=2, save_file=True, gmsh_interface=True,
                            mesh_output_name="Embankment")

layers = inputs.prepare_inputs()  # List of all points list


if __name__ == '__main__':
    gmsh.initialize()
    gmsh.model.add(inputs.mesh_output_name)

    efn.make_geometry(layers, inputs.default_mesh_size, inputs.name_label_list, inputs.dims, inputs.depth)

    efn.fragment(inputs.dims, inputs.number_of_layers)

    gmsh.model.occ.synchronize()

    efn.set_mesh_size(layers, inputs.number_of_layers, inputs.mesh_size_list, inputs.eps, inputs.depth)

    gmsh.model.mesh.generate(inputs.dims)
    file_extension = ".msh"
    mesh_output_file = str(inputs.mesh_output_name) + str(inputs.dims) + "D" + file_extension
    gmsh.write(mesh_output_file)
    if inputs.gmsh_interface:
        gmsh.fltk.run()
    gmsh.finalize()
