import gmsh
import numpy as np


def init():
    """
    gets user input
    :return: input_points, depth, mesh_size, dims, save_file, name_label, mesh_output_name, gmsh_interface
    """
    # define the points of the surface as a list of tuples
    input_points = [(0, 0, 0), (1, 0, 0), (1, 3, 0), (0, 3, 0), (-1, 1.5, 0)]
    # define the mesh size
    mesh_size = 2
    # define geometry dimension; input "3" for 3D to extrude the 2D surface, input "2" for 2D
    dims = 3
    # if 3D, input depth of geometry to be extruded from 2D surface
    depth = 2
    # set a name label for the surface
    name_label = "Soil Layer"
    # if "True", saves mesh data to separate mdpa files; otherwise "False"
    save_file = True
    # if "True", opens gmsh interface; otherwise "False"
    gmsh_interface = True
    # set a name for mesh output file
    mesh_output_name = "geometry"

    return input_points, depth, mesh_size, dims, save_file, name_label, mesh_output_name, gmsh_interface


def create_point(coordinates, mesh_size):
    """
    creates points in gmsh
    :param coordinates: gets points coordinates in order from user
    :param mesh_size: gets mesh size from user
    :return: -
    """
    x = coordinates[0]
    y = coordinates[1]
    z = coordinates[2]
    lc = mesh_size
    gmsh.model.geo.addPoint(x, y, z, lc)


def create_line(point_ids):
    """
    creates lines in gmsh
    :param point_ids: gets point tags in order
    :return: -
    """
    point1 = point_ids[0]
    point2 = point_ids[1]
    gmsh.model.geo.addLine(point1, point2)


def create_surface(line_ids, name_label):
    """
    creates curve and then surface in gmsh by using line tags
    :param line_ids: gets line tags in order
    :param name_label: surface name label from user input
    :return: returns the surface tag and surface name label
    """
    gmsh.model.geo.addCurveLoop(line_ids, 1)
    surfaces = gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.setPhysicalName(2, surfaces, name_label)
    return surfaces


def create_volume(surface_id, depth):
    """
    creates volume by extruding 2D surface
    :param surface_id: surface tag
    :param depth: depth of 3D geometry
    :return: -
    """
    gmsh.model.geo.extrude([(2, surface_id)], 0, 0, depth)


def make_geometry_2D(point_coordinates, point_pairs, mesh_size, name_label):
    """
    takes point_pairs and puts their tags as the beginning and end of line in gmsh to create line
    then creates surface
    :param point_coordinates: point coordinates
    :param point_pairs: has been generated by storing point tags of two consecutive points in an array
    :param mesh_size: mesh size
    :param name_label: surface name label from user input
    :return: creates 2D geometries
    """
    for point in point_coordinates:
        coordinate = [point[0], point[1], point[2]]
        create_point(coordinate, mesh_size)

    line_lists = []

    for i in range(len(point_pairs)):
        line = [point_pairs[i][0], point_pairs[i][1]]
        line_lists.append(i + 1)
        create_line(line)

    surfaces = create_surface(line_lists, name_label)
    return surfaces


def make_geometry_3D(point_coordinates, point_pairs, mesh_size, depth, name_label):
    """
    creates 3D geometries by extruding the 2D surface
    :param point_coordinates: geometry points coordinates
    :param point_pairs: points paired for lines
    :param mesh_size: mesh size
    :param depth: depth of 3D geometry
    :param name_label: surface name label from user input
    :return: -
    """
    surfaces = make_geometry_2D(point_coordinates, point_pairs, mesh_size, name_label)
    create_volume(surfaces, depth)


def extract_mesh_data(dims):
    """
    gets gmsh output data
    :param dims: geometry dimension (2=2D or 3=3D)
    :return: geometry and mesh data: node tags, node coordinates, element types, element tags 0D, 1D, 2D, 3D
    """
    node_tags, node_coords, node_params = gmsh.model.mesh.getNodes()  # nodes, elements
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements()

    # get number of nodes per element type
    node_shape = [get_num_nodes_from_elem_type(elem_type) for elem_type in elem_types]
    num_elem = sum(len(i) for i in elem_tags)
    print(" - Mesh has " + str(len(node_tags)) + " nodes and " + str(num_elem) +
          " elements")

    num_nodes = len(node_tags)
    coord = np.reshape(node_coords, (num_nodes, 3))

    num_1d_elements = len(elem_tags[0])
    node_tag_1D = np.reshape(elem_node_tags[0], (num_1d_elements, node_shape[0]))

    num_2d_elements = len(elem_tags[1])
    node_tag_2D = np.reshape(elem_node_tags[1], (num_2d_elements, node_shape[1]))

    if dims == 3:
        num_3d_elements = len(elem_tags[2])
        node_tag_3D = np.reshape(elem_node_tags[2], (num_3d_elements, node_shape[2]))
        return coord, node_tags, elem_types, elem_tags, node_tag_1D, node_tag_2D, node_tag_3D

    if dims == 2:
        return coord, node_tags, elem_types, elem_tags, node_tag_1D, node_tag_2D


def get_num_nodes_from_elem_type(elem_type):
    """
    gets number of nodes from element types
    :param elem_type: int that defines type of elements
    :return: number of nodes needed for a type of element
    """

    # 2 node line
    if elem_type == 1:
        return 2  # number of nodes needed for 2-node line
    # 3 node triangle
    if elem_type == 2:
        return 3  # number of nodes needed for 3-node triangle
    # 4 node quadrangle
    if elem_type == 3:
        return 4  # number of nodes needed for 4-node quadrangle
    # 4 node tetrahedron
    if elem_type == 4:
        return 4  # number of nodes needed for 4-node tetrahedron
        # 1 node
    if elem_type == 15:
        return 1  # number of nodes needed for 1-node


def generate_gmsh_mesh(point_coordinates, depth, mesh_size, dims, save_file, name_label, mesh_output_name,
                       gmsh_interface):
    """
    crates point pairs by storing point tags of two consecutive points in an array then
    generates mesh for geometries in gmsh
    :param point_coordinates: user input points of the surface as a list of tuples
    :param depth: depth of 3D geometry
    :param mesh_size: mesh size
    :param dims: geometry dimension (2=2D or 3=3D)
    :param save_file: if True saves mesh data to mdpa files
    :param name_label: surface name label from user input
    :param mesh_output_name: name of mesh output file
    :param gmsh_interface: user indicates whether to open gmsh interface
    :return: saves mesh data and opens gmsh interface
    """

    point_pairs = []
    for i in range(len(point_coordinates) - 1):
        # puts two consecutive points tags as the beginning and end of line in an array
        point_pair = [i + 1, i + 2]
        point_pairs.append(point_pair)
    # make a pair that connects last point to first point
    point_pairs.append([len(point_coordinates), 1])

    gmsh.initialize()
    gmsh.model.add(mesh_output_name)

    volumes = []
    if dims == 3:
        make_geometry_3D(point_coordinates, point_pairs, mesh_size, depth, name_label)
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(dims)
        # extract mesh data
        node_coords, node_tags, elem_types, elem_tags, node_tag_1D, node_tag_2D, node_tag_3D = extract_mesh_data(dims)
        nodes, lines, surfaces = get_data_for_kratos(node_coords, node_tags, elem_tags, node_tag_1D, node_tag_2D)
        volumes = np.concatenate((elem_tags[2][:, None], np.array(node_tag_3D)), axis=1)

    elif dims == 2:
        make_geometry_2D(point_coordinates, point_pairs, mesh_size, name_label)
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(dims)
        # extract mesh data
        node_coords, node_tags, elem_types, elem_tags, node_tag_1D, node_tag_2D = extract_mesh_data(dims)
        nodes, lines, surfaces = get_data_for_kratos(node_coords, node_tags, elem_tags, node_tag_1D, node_tag_2D)

    if save_file:
        save_mesh_data(nodes, lines, surfaces, volumes)

    # writes mesh file output in .msh format
    file_extension = ".msh"
    mesh_output_file = mesh_output_name + file_extension
    gmsh.write(mesh_output_file)

    # opens GMsh interface
    if gmsh_interface:
        gmsh.fltk.run()
    gmsh.finalize()

    if dims == 2:
        return node_coords, node_tags, elem_types, elem_tags, node_tag_1D, node_tag_2D
    if dims == 3:
        return node_coords, node_tags, elem_types, elem_tags, node_tag_1D, node_tag_2D, node_tag_3D


def get_data_for_kratos(node_coords, node_tags, elem_tags, node_tag_1D, node_tag_2D):
    """
    gets mesh data for Kratos
    :param node_coords: node coordinates
    :param node_tags: node tags
    :param elem_tags: all element tags in an array separated by element type
    :param node_tag_1D: node tags of start and end of line
    :param node_tag_2D: node tags of surface
    :return: node tag followed by node coordinates and element tag followed by node tags in an array
    """
    nodes = np.concatenate((node_tags[:, None], np.array(node_coords)), axis=1)
    lines = np.concatenate((elem_tags[0][:, None], np.array(node_tag_1D)), axis=1)
    surfaces = np.concatenate((elem_tags[1][:, None], np.array(node_tag_2D)), axis=1)

    return nodes, lines, surfaces


def save_mesh_data(nodes, lines, surfaces, volumes):
    """
    saves mesh data to mdpa file
    :param nodes: node tag followed by node coordinates in an array
    :param lines: line tag followed by node tags of start and end of line in an array
    :param surfaces: surface tag followed by node tags of surface in an array
    :param volumes: volume tag followed by node tags of volume in an array
    :return: -
    """
    np.savetxt('0.nodes.mdpa', nodes, fmt=['%.f', '%.10f', '%.10f', '%.10f'], delimiter=' ')
    np.savetxt('1.lines.mdpa', lines, delimiter=' ')
    np.savetxt('2.surfaces.mdpa', surfaces, delimiter=' ')
    np.savetxt('3.volumes.mdpa', volumes, delimiter=' ')


if __name__ == '__main__':
    input_points, depth, mesh_size, dims, save_file, name_label, mesh_output_name, gmsh_interface = init()
    generate_gmsh_mesh(input_points, depth, mesh_size, dims, save_file, name_label, mesh_output_name, gmsh_interface)
