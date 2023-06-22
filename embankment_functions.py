import gmsh
import numpy as np

class InputGenerator:
    def __init__(self, number_of_layers, mesh_size_list, input_points_list, name_label_list, eps,
                 default_mesh_size, dims, depth, save_file, gmsh_interface,
                 mesh_output_name):
        self.number_of_layers = number_of_layers
        self.mesh_size_list = mesh_size_list
        self.input_points_list = input_points_list
        self.name_label_list = name_label_list
        self.eps = eps
        self.default_mesh_size = default_mesh_size
        self.dims = dims
        self.depth = depth
        self.save_file = save_file
        self.gmsh_interface = gmsh_interface
        self.mesh_output_name = mesh_output_name

    def prepare_inputs(self):
        layer_list = []
        for i in range(self.number_of_layers):
            layer = self.input_points_list[i]
            layer_list.append(layer)
        return layer_list


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
    gmsh.model.occ.addPoint(x, y, z, lc)


def create_line(point_ids):
    """
    creates lines in gmsh
    :param point_ids: gets point tags in order
    :return: -
    """
    point1 = point_ids[0]
    point2 = point_ids[1]
    gmsh.model.occ.addLine(point1, point2)


def create_surface(line_ids, name_label):
    """
    creates curve and then surface in gmsh by using line tags
    :param line_ids: gets line tags in order
    :param name_label: surface name label from user input
    :return: returns the surface tag and surface name label
    """
    curve_loop = gmsh.model.occ.addCurveLoop(line_ids)
    surfaces = gmsh.model.occ.addPlaneSurface([curve_loop])
    gmsh.model.setPhysicalName(2, surfaces, name_label)
    return surfaces


def create_volume(surface_id, depth, name_label, volume_index):
    """
        extrude the surface to create a volume
        :param surface_id: surface tag
        :param depth: depth of 3D geometry
        :param name_label: surface name label from user input
        :param volume_index: volume index
        :return: -
        """
    gmsh.model.occ.extrude([(2, surface_id)], 0, 0, depth)
    gmsh.model.setPhysicalName(2, volume_index+1, name_label)


def generate_point_pairs(coord_list):
    """
    Generates pairs of point IDs which form a line
    :param coord_list: geometry points coordinates
    :return: pairs of point ids which create a line
    """
    list_point_pairs = []
    counter = 0
    first_point_tag = 0
    for i in range(len(coord_list)):  # number of layers
        point_pairs = []
        for j in range(len(coord_list[i]) - 1):
            if j == 0:
                first_point_tag = counter+i+1  # saving the first point tag in order to return from last point to it
            # puts two consecutive points tags as the beginning and end of line in an array
            point_pair = [counter + i + 1, counter + i + 2]
            point_pairs.append(point_pair)
            counter += 1
        # make a pair that connects last point to first point
        point_pairs.append([counter+i+1, first_point_tag])
        list_point_pairs.append(point_pairs)

    return list_point_pairs


def make_points(point_coordinates, default_mesh_size):
    for point in point_coordinates:
        coordinate = [point[0], point[1], point[2]]
        create_point(coordinate, default_mesh_size)


def make_lines(point_pairs):
    counter = 0
    list_lines = []
    for i in range(len(point_pairs)):
        lines = []
        for j in range(len(point_pairs[i])):
            line = [point_pairs[i][j][0], point_pairs[i][j][1]]
            lines.append(counter + 1)
            create_line(line)
            counter += 1
        list_lines.append(lines)
    return list_lines


def make_surfaces(line_list, name_label):
    surfaces = []
    for i in range(len(line_list)):
        surfaces.append(create_surface(line_list[i], name_label[i]))

    return surfaces


def make_volume(surface_id, depth, name_label):
    volumes = []
    for i in range(len(surface_id)):
        volumes.append(create_volume(surface_id[i], depth, name_label[i], i))


def make_geometry(point_lists, default_mesh_size, name_label_list, dims, depth):
    for i in range(len(point_lists)):
        make_points(point_lists[i], default_mesh_size)

    pair_lists = generate_point_pairs(point_lists)
    line_lists = make_lines(pair_lists)
    surface_tags = make_surfaces(line_lists, name_label_list)
    if dims == 2:
        return surface_tags
    if dims == 3:
        make_volume(surface_tags, depth, name_label_list)


def fragment(dims, number_of_layers):
    """
    fragments two surfaces to remove the double line and connect the two surfaces
    :param dims: dimension
    :param number_of_layers: number of layers
    :return:
    """
    # gmsh.model.occ.synchronize()
    # entities = gmsh.model.get_entities(dims)
    # print(entities)
    # gmsh.model.occ.fragment([(dims, 1)], [(dims, 2)])
    # gmsh.model.occ.synchronize()
    # gmsh.model.occ.fragment([(dims, 1)], [(dims, 3)])
    # gmsh.model.occ.synchronize()
    # common_surface = gmsh.model.occ.intersect([(3, 1), (3, 2)])
    # gmsh.model.occ.remove([(2, common_surface)])
    gmsh.model.occ.removeAllDuplicates()

    # for i in range(len(entities)):
    #     for j in range(i, len(entities)):
    #         if i == j: continue
    #         print("i=",i,"\nj=",j,"\nentities[i][1]=", entities[i][1], "\nentities[j][1]=", entities[j][1])
    #         gmsh.model.occ.fragment([(dims, entities[i][1])], [(dims, entities[j][1])])

    # for i in range(number_of_layers):
    #     for j in range(i, number_of_layers):
    #         if i == j: continue
    #         print("i=",i,"\nj=",j,"\nentities[i][1]=", i+1, "\nentities[j][1]=", j+1)
    #         gmsh.model.occ.fragment([(dims, i+1)], [(dims, j+1)])

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


def set_mesh_size(layer_list, number_of_layers, mesh_size_lists, eps, depth):
    """
    Sets the mesh size by considering the left bottom point and the right top point and drawing a box within epsilon of
    these points to select all the points in that box to assign the related mesh size
    :param layer_list: list of layers' points
    :param number_of_layers: number of layers
    :param mesh_size_lists: mesh sizes for each layer
    :param eps: epsilon defined for the box to be big enough to cover the whole layer
    :param depth: depth of geometry to be extruded
    :return: -
    """
    entities_list = []
    for i in range(number_of_layers):
        min_x = min(layer_list[i], key=lambda p: p[0])
        max_x = max(layer_list[i], key=lambda p: p[0])
        min_y = min(layer_list[i], key=lambda p: p[1])
        max_y = max(layer_list[i], key=lambda p: p[1])
        min_z = min(layer_list[i], key=lambda p: p[2])
        max_z = max(layer_list[i], key=lambda p: p[2])

        entities_list.append(gmsh.model.getEntitiesInBoundingBox(min_x[0] - eps, min_y[1] - eps, min_z[2] - eps,
                                                                 max_x[0] + eps,
                                                                 max_y[1] + eps,
                                                                 max_z[2] + depth + eps, 0))

    for i in range(len(entities_list)):
        gmsh.model.mesh.setSize(entities_list[i], mesh_size_lists[i])
