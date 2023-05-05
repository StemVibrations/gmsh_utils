import gmsh
import numpy as np

class input_generator:
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


def create_volume(surface_id, depth, name_label, i):
    """
        extrude the surface to create a volume
        :param surface_id: surface tag
        :param depth: depth of 3D geometry
        :return: -
        """
    volumes = gmsh.model.occ.extrude([(2, surface_id)], 0, 0, depth)
    gmsh.model.setPhysicalName(2, i+1, name_label)


def generate_point_pairs(coord_list):
    """
    Generates pairs of point IDs which form a line
    :param coord_list: geometry points coordinates
    :return: pairs of point ids which create a line
    """
    list_point_pairs = []
    counter = 0
    for i in range(len(coord_list)): # number of layers
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


def make_surfaces(line_list, name_label, dims):
    surfaces = []
    for i in range(len(line_list)):
        surfaces.append(create_surface(line_list[i], name_label[i]))

    return surfaces

def make_volume(surface_id, depth, name_label):
    volumes =[]
    for i in range(len(surface_id)):
        volumes.append(create_volume(surface_id[i], depth, name_label[i], i))
    # vol_tag = gmsh.model.addPhysicalGroup(3, [volume_id])
    return volumes

def make_geometry(point_lists, default_mesh_size, name_label_list, dims, depth):
    for i in range (len(point_lists)):
        make_points(point_lists[i], default_mesh_size)

    pair_lists = generate_point_pairs(point_lists)
    line_lists = make_lines(pair_lists)
    # print(line_lists)
    # input()
    surface_tags  = make_surfaces(line_lists, name_label_list, dims)
    # print(surface_tags)
    # input()
    if dims == 2:
        return surface_tags
    if dims == 3:
        volume_tags = make_volume(surface_tags, depth, name_label_list)
        # print (volume_tags)
        # input()
        return volume_tags


def fragment(tags, dims, number_of_layers):
    """
    fragments two surfaces to remove the dubble line and connect the two surfaces
    :param dims: dimension
    :param number_of_layers: number of layers
    :return:
    """
    for i in range(number_of_layers - 1):
        # gmsh.model.occ.fragment([(3, tags[i])], [(3, tags[i + 1])])
        gmsh.model.occ.fragment([(dims,i+1)], [(dims,i+2)])

def set_mesh_size(layer_list, number_of_layers, mesh_size_lists, eps, depth):
    """
    Sets the mesh size by considering the left bottom point and the right top point and drawing a box within epsilon of
    these points to select all the points in that box to assign the related mesh size
    :param layer_list: list of layers' points
    :param number_of_layers: number of layers
    :param mesh_size_lists: mesh sizes for each layer
    :param eps: epsilon defined for the box to be big enough to cover the whole layer
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