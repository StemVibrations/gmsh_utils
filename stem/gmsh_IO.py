import gmsh
import numpy as np

from enum import Enum
import re

#todo Put this file in its own package, e.g. GmshUtils

class ElementType(Enum):
    """
    Enum of the element types as present in Gmsh, where the enum value corresponds to the element type number in gmsh

    """

    LINE_2N = 1
    TRIANGLE_3N = 2
    QUADRANGLE_4N = 3
    TETRAHEDRON_4N = 4
    HEXAHEDRON_8N = 5
    PRISM_6N = 6
    PYRAMID_5N = 7
    LINE_3N = 8
    TRIANGLE_6N = 9
    QUADRANGLE_9N = 10
    TETRAHEDRON_10N = 11
    HEXAHEDRON_27N = 12
    PRISM_18N = 13
    PYRAMID_14N = 14
    POINT_1N = 15
    QUADRANGLE_8N = 16
    HEXAHEDRON_20N = 17
    PRISM_15N = 18
    PYRAMID_13N = 19
    #todo complete list


class GmshIO:

    def __init__(self):
        pass

    def create_point(self, coordinates, mesh_size):
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


    def create_line(self, point_ids):
        """
        creates lines in gmsh
        :param point_ids: gets point tags in order
        :return: -
        """
        point1 = point_ids[0]
        point2 = point_ids[1]
        gmsh.model.geo.addLine(point1, point2)


    def create_surface(self, line_ids, name_label):
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

    def create_volume_by_extruding_surface(self, surface_id, depth):
        """
        creates volume by extruding 2D surface
        :param surface_id: surface tag
        :param depth: depth of 3D geometry
        :return: -
        """
        surface_ndim = 2
        gmsh.model.geo.extrude([(surface_ndim, surface_id)], 0, 0, depth)

    def make_geometry_2D(self, point_coordinates, point_pairs, mesh_size, name_label):
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
            self.create_point(coordinate, mesh_size)

        line_lists = []

        for i in range(len(point_pairs)):
            line = [point_pairs[i][0], point_pairs[i][1]]
            line_lists.append(i + 1)
            self.create_line(line)

        surfaces = self.create_surface(line_lists, name_label)
        return surfaces

    def make_geometry_3D(self, point_coordinates, point_pairs, mesh_size, depth, name_label):
        """
        creates 3D geometries by extruding the 2D surface
        :param point_coordinates: geometry points coordinates
        :param point_pairs: points paired for lines
        :param mesh_size: mesh size
        :param depth: depth of 3D geometry
        :param name_label: surface name label from user input
        :return: -
        """
        surfaces = self.make_geometry_2D(point_coordinates, point_pairs, mesh_size, name_label)
        self.create_volume_by_extruding_surface(surfaces, depth)

    def generate_point_pairs(self, point_coordinates):
        """
        Generates pairs of point IDs which form a line

        :param point_coordinates: geometry points coordinates
        :return: pairs of point ids which create a line
        """
        point_pairs = []
        for i in range(len(point_coordinates) - 1):
            # puts two consecutive points tags as the beginning and end of line in an array
            point_pair = [i + 1, i + 2]
            point_pairs.append(point_pair)
        # make a pair that connects last point to first point
        point_pairs.append([len(point_coordinates), 1])

        return point_pairs

    def get_num_nodes_from_elem_type(self, elem_type):
        """
        gets number of nodes from element types
        :param elem_type: int that defines type of elements
        :return: number of nodes needed for a type of element
        """

        # get name from element type enum
        element_name = ElementType(elem_type).name

        # get number of nodes from the enum name
        num_nodes = int(re.findall(r'\d+', element_name)[0])

        return num_nodes

    def generate_gmsh_mesh(self, point_coordinates, depth, mesh_size, dims, name_label, mesh_output_name,
                           save_file=False, open_gmsh_interface=False):
        """
        crates point pairs by storing point tags of two consecutive points in an array then
        generates mesh for geometries in gmsh
        :param point_coordinates: user input points of the surface as a list of tuples
        :param depth: depth of 3D geometry
        :param mesh_size: mesh size
        :param dims: geometry dimension (2=2D or 3=3D)
        :param name_label: surface name label from user input
        :param mesh_output_name: name of mesh output file
        :param save_file: if True saves mesh data to gmsh msh file
        :param open_gmsh_interface: user indicates whether to open gmsh interface
        :return: saves mesh data and opens gmsh interface
        """

        point_pairs = self.generate_point_pairs(point_coordinates)

        gmsh.initialize()
        gmsh.model.add(mesh_output_name)

        volumes = []
        if dims == 3:
            self.make_geometry_3D(point_coordinates, point_pairs, mesh_size, depth, name_label)
            # # extract mesh data
            # node_coords, node_tags, elem_types, elem_tags, node_tag_1D, node_tag_2D, node_tag_3D = extract_mesh_data(
            #     dims)
            # nodes, lines, surfaces = get_data_for_kratos(node_coords, node_tags, elem_tags, node_tag_1D, node_tag_2D)
            # volumes = np.concatenate((elem_tags[2][:, None], np.array(node_tag_3D)), axis=1)

        elif dims == 2:
            self.make_geometry_2D(point_coordinates, point_pairs, mesh_size, name_label)

        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(dims)

        mesh_data = self.extract_mesh_data(gmsh.model.mesh)

        if save_file:
            # writes mesh file output in .msh format
            file_extension = ".msh"
            mesh_output_file = mesh_output_name + file_extension
            gmsh.write(mesh_output_file)

        # opens GMsh interface
        if open_gmsh_interface:
            gmsh.fltk.run()

        gmsh.finalize()

        return mesh_data

    def extract_element_data(self,elem_type, elem_tags, elem_node_tags):
        """
        Gets gmsh data belonging to a single element

        :param elem_type: gmsh id for element type
        :param elem_tags:  gmsh id for element number
        :param elem_node_tags: gmsh node ids which belong the the element

        :return: dictionary of element data
        """

        element_name = ElementType(elem_type).name
        n_nodes_per_element = self.get_num_nodes_from_elem_type(elem_type)
        num_elements = len(elem_tags)
        element_node_ids = np.reshape(elem_node_tags, (num_elements, n_nodes_per_element))

        return {element_name: {"element_ids": elem_tags,
                               "element_nodes": element_node_ids}}


    def extract_mesh_data(self, gmsh_mesh):
        """
        gets gmsh output data

        :param gmsh_mesh: the mesh as generated by gmsh
        :return: dictionary which contains nodal and elemental information
        """

        mesh_data = {"nodes": {},
                     "elements": {}}

        # get nodal information
        node_tags, node_coords, node_params = gmsh_mesh.getNodes()  # nodes, elements

        # reshape nodal coordinate array to [num nodes, 3]
        num_nodes = len(node_tags)
        node_coordinates = np.reshape(node_coords, (num_nodes, 3))

        mesh_data["nodes"]["coordinates"] = node_coordinates
        mesh_data["nodes"]["ids"] = node_tags

        # get all elemental information
        elem_types, elem_tags, elem_node_tags = gmsh_mesh.getElements()

        # todo, this is unhandy for the future and the connection to kratos, handier would be to group elements by physical group
        for elem_type, elem_tag, elem_node_tag in zip(elem_types, elem_tags, elem_node_tags):
            element_dict = self.extract_element_data(elem_type, elem_tag, elem_node_tag)
            mesh_data["elements"].update(element_dict)

        return mesh_data


if __name__ == '__main__':
    pass