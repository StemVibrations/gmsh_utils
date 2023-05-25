from typing import Dict, List, Union, Type
from enum import Enum
import re
import gmsh
import numpy as np
import numpy.typing as npt


# todo Put this file in its own package, e.g. GmshUtils

class ElementType(Enum):
    """
    Enum of the element types as present in Gmsh, where the enum value corresponds to the element type number in gmsh

    """

    LINE_2N = 1
    TRIANGLE_3N = 2
    QUADRANGLE_4N = 3
    TETRAHEDRON_4N = 4
    HEXAHEDRON_8N = 5
    LINE_3N = 8
    TRIANGLE_6N = 9
    TETRAHEDRON_10N = 11
    POINT_1N = 15
    QUADRANGLE_8N = 16
    HEXAHEDRON_20N = 17


class GmshIO:
    """
    Class for reading and writing mesh data to and from Gmsh

    Attributes
    ----------
    mesh_data : Dict
        Dictionary containing the mesh data, i.e. nodal ids and coordinates; and elemental ids, connectivity's
        and element types.

    """

    def __init__(self):
        self.__mesh_data = {}

    @property
    def mesh_data(self) -> Dict[str, object]:
        """
        Returns the mesh data dictionary

        Returns:
            Dict: Dictionary containing the mesh data, i.e. nodal ids and coordinates; and elemental ids, connectivity's
            and element types.
        """

        return self.__mesh_data

    # def prepare_inputs(self, input_points_list):
    #     """
    #     Prepares the input points list for the mesh generation.
    #
    #     Args:
    #         input_points_list (List[List[float]]): The list of input points provided by user input.
    #
    #     Returns:
    #         List[List[float]]: The list of input points for the mesh generation.
    #     """
    #     layer_list = []
    #     number_of_layers = len(input_points_list)
    #     for i in range(number_of_layers):
    #         layer = input_points_list[i]
    #         layer_list.append(layer)
    #     return layer_list

    def create_point(self, coordinates: Union[List[float], npt.NDArray[np.float64]], element_size: float) -> None:
        """
        Creates points in gmsh.

        Args:
            coordinates (Union[List[float], npt.NDArray[float]]): An Iterable of point x,y,z coordinates.
            mesh_size (float): The element size.

        Returns
        -------
        None
        """
        x = coordinates[0]
        y = coordinates[1]
        z = coordinates[2]
        point_id = gmsh.model.occ.addPoint(x, y, z, element_size)
        return point_id

    def create_line(self, point_ids: Union[List[int], npt.NDArray[np.int_]]) -> None:
        """
        Creates lines in gmsh.

        Args:
            point_ids (Union[List[int], npt.NDArray[int]]): A list of point tags in order.

        Returns:
            None
        """

        point1 = point_ids[0]
        point2 = point_ids[1]
        line_id = gmsh.model.occ.addLine(point1, point2)
        return line_id

    def create_surface(self, line_ids: Union[List[int], npt.NDArray[np.int_]], name_label: str) -> int:
        """
        Creates curve and then surface in gmsh by using line tags.

        Args:
            line_ids (Union[List[int], npt.NDArray[int]]): A list of line tags in order.
            name_label (str): The surface name label provided by user input.

        Returns:
            int: surface id
        """

        curve_loop_id = gmsh.model.occ.addCurveLoop(line_ids)
        surface_id: int = gmsh.model.occ.addPlaneSurface([curve_loop_id])
        surface_ndim = 2
        gmsh.model.setPhysicalName(surface_ndim, surface_id, name_label)
        return surface_id

    def create_volume_by_extruding_surface(self, surface_id: int,
                                           extrusion_length: Union[List[float], npt.NDArray[np.float64]],
                                           name_label: str,
                                           volume_ids: int) -> None:
        """
        Creates volume by extruding a 2D surface

        Args:
            surface_id (int): The surface tag.
            extrusion_length (Union[List[float], npt.NDArray[float]]): The extrusion length in x, y and z direction.
            name_label (str): The volume name label provided by user input
            volume_ids (int): The volume tag.
        Returns:
            None
        """

        surface_ndim = 2
        gmsh.model.occ.extrude([(surface_ndim, surface_id)], extrusion_length[0], extrusion_length[1],
                               extrusion_length[2])
        gmsh.model.setPhysicalName(surface_ndim, volume_ids+1, name_label)

    def generate_point_pairs(self, number_of_layers,
                             number_of_points_in_layers: Union[List[List[float]], npt.NDArray[np.float64]], point_ids)\
                             -> List[List[List[int]]]:
        """
        Generates pairs of point IDs which form a line

        Args:

        Returns:
            List[List[int]]: A list of pairs of point IDs which create a line.
        """
        list_point_pairs = []
        first_point_tag = 0
        for i in range(number_of_layers):  # number of layers
            point_pairs = []
            for j in range(number_of_points_in_layers[i] - 1): # number of points in each layer (-1 for turning back)
                if j == 0:
                    # saving the first point tag in order to return from last point to it
                    first_point_tag = point_ids[i][j]
                # puts two consecutive points tags as the beginning and end of line in an array
                point_pair = [point_ids[i][j], point_ids[i][j+1]]
                point_pairs.append(point_pair)
            # make a pair that connects last point to first point
            last_point_tag = len(point_ids[i])-1
            point_pairs.append([point_ids[i][last_point_tag], first_point_tag])
            list_point_pairs.append(point_pairs)

        return list_point_pairs

    def make_points(self, point_coordinates: Union[List[List[float]], npt.NDArray[np.float64]],
                    default_mesh_size:  float) -> List[None]:
        """Makes points with point tags by getting coordinates.

        Args:
            point_coordinates (Union[List[List[float]], npt.NDArray[np.float64]]): An Iterable of point x,y,z coordinates.
            default_mesh_size (float): The element size.

        Returns:
            List[None]: A list of point tags.
            """
        list_point_ids = []
        for point in point_coordinates:
            point_id = self.create_point(point, default_mesh_size)
            list_point_ids.append(point_id)
        return list_point_ids

    def make_lines(self, point_pairs: Union[List[List[int]], npt.NDArray[np.int_]]):
        """Makes lines with line tags by getting point pairs.

        Args:
            point_pairs (Union[List[List[int]], npt.NDArray[np.int_]]): A list of pairs of point tags which create a line.

        Returns:
            List[List[int]]: A list of line tags.
            """
        list_lines = []
        for i in range(len(point_pairs)):
            lines = []
            for j in range(len(point_pairs[i])):
                line = [point_pairs[i][j][0], point_pairs[i][j][1]]
                line_id = self.create_line(line)
                lines.append(line_id)
            list_lines.append(lines)
        return list_lines

    def make_surfaces(self, line_list: Union[List[List[float]], npt.NDArray[np.float64]], name_label: List[str]):
        """Makes surfaces with surface tags by getting line tags.

        Args:
            line_list (Union[List[List[float]], npt.NDArray[np.float64]]): A list of line tags in order.
            name_label (List[str]): A list of surface name labels provided by user input.

        Returns:
            List[int]: A list of surface tags.
            """
        surfaces = []
        for i in range(len(line_list)):
            surfaces.append(self.create_surface(line_list[i], name_label[i]))

        return surfaces

    def make_volume(self, surface_id: [List[int], npt.NDArray[np.int_]],
                    extrusion_length: Union[List[float], npt.NDArray[np.float64]], name_label: List[str]):
        volumes = []
        for volume in range(len(surface_id)):
            volumes.append(self.create_volume_by_extruding_surface(surface_id[volume], extrusion_length,
                                                                   name_label[volume], volume))

    def make_geometry_2d(self, point_coordinates: Union[List[List[float]], npt.NDArray[np.float64]],
                         default_mesh_size: float, name_label_list: List[str]) -> List[int]:
        """
        Takes point_pairs and puts their tags as the beginning and end of line in gmsh to create line,
        then creates surface to make 2D geometry.

        Args:
            point_coordinates (Union[List[float], npt.NDArray[np.float64]]): A list of point coordinates.
            default_mesh_size (float): The default mesh size provided by user.
            name_label_list (List[str]): A list of surface name labels provided by user input.

        Returns:
            int: Surface id.
        """
        number_of_layers = len(point_coordinates)
        n_points_in_layers = [] #number of points in each layer
        point_ids_of_layers = [] #list of point ids of coordinates of all layers
        for i in range(number_of_layers):
            n_points_in_layers.append(len(point_coordinates[i]))
            list_point_id = self.make_points(point_coordinates[i], default_mesh_size)
            point_ids_of_layers.append(list_point_id)

        pair_lists = self.generate_point_pairs(number_of_layers, n_points_in_layers, point_ids_of_layers)
        line_lists = self.make_lines(pair_lists)
        surface_ids = self.make_surfaces(line_lists, name_label_list)

        return surface_ids

    def make_geometry_3d(self, point_coordinates: Union[List[List[float]], npt.NDArray[np.float64]],
                         default_mesh_size: float, name_label_list: List[str],
                         extrusion_length: Union[List[float], npt.NDArray[np.float64]], ) -> None:
        """
        Creates 3D geometries by extruding the 2D surface

        Args:
            point_coordinates (Union[List[float], npt.NDArray[float]]): Geometry points coordinates.
                points in an array.
            default_mesh_size (float): The default mesh size provided by user.
            name_label_list (List[str]): A list of labels provided by user input.
            extrusion_length (Union[List[float], npt.NDArray[float]]): The extrusion length in x, y and z direction.

        Returns:
            None
        """

        surfaces = self.make_geometry_2d(point_coordinates, default_mesh_size, name_label_list)
        self.make_volume(surfaces, extrusion_length, name_label_list)

    def remove_duplicates(self):
        """
        Removes duplicate entities from the geometry.

        """
        gmsh.model.occ.removeAllDuplicates()

    @staticmethod
    def get_num_nodes_from_elem_type(elem_type: int) -> int:
        """
        Gets number of nodes from element types

        Args:
            elem_type (int): An integer that defines the type of element.

        Returns:
            int: The number of nodes needed for a type of element.
        """

        # get name from element type enum
        element_name = ElementType(elem_type).name

        # get number of nodes from the enum name
        num_nodes = int(re.findall(r'\d+', element_name)[0])

        return num_nodes

    def generate_gmsh_mesh(self, point_coordinates: Union[List[List[float]], npt.NDArray[np.float64]],
                           extrusion_length: Union[List[float], npt.NDArray[np.float64]],
                           mesh_size: float, dims: int, name_label: List[str], mesh_name: str,
                           mesh_output_dir: str, save_file: bool = False, open_gmsh_gui: bool = False) -> None:
        """
        Creates point pairs by storing point tags of two consecutive points in an array,
        then generates mesh for geometries in gmsh.

        Args:
            point_coordinates (Union[List[List[float]], npt.NDArray[np.float64]]): User input points of the surface as
                a list or ndarray.
            extrusion_length (Union[List[float], npt.NDArray[float]]): The depth of 3D geometry.
            mesh_size (float): The mesh size provided by user.
            dims (int): The dimension of geometry (2=2D or 3=3D).
            name_label (str): The surface name label provided by user input.
            mesh_name (str): Name of gmsh model and mesh output file.
            mesh_output_dir (str): Output directory of mesh file.
            save_file (bool, optional): If True, saves mesh data to gmsh msh file. (default is False)
            open_gmsh_gui (bool, optional): User indicates whether to open gmsh interface (default is False)

        Returns:
            None
        """

        # todo add check for clockwise or anticlockwise

        gmsh.initialize()
        gmsh.model.add(mesh_name)

        # layer_list = self.prepare_inputs(point_coordinates)
        if dims == 3:
            self.make_geometry_3d(point_coordinates, mesh_size, name_label, extrusion_length)

        elif dims == 2:
            self.make_geometry_2d(point_coordinates, mesh_size, name_label)

        self.remove_duplicates()
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(dims)

        # extracts mesh data from gmsh
        self.__mesh_data = self.extract_mesh_data(gmsh.model.mesh)

        if save_file:
            # writes mesh file output in .msh format
            file_extension = ".msh"
            mesh_output_file = mesh_output_dir + mesh_name + file_extension
            gmsh.write(mesh_output_file)

        # opens Gmsh interface
        if open_gmsh_gui:
            gmsh.fltk.run()

        gmsh.finalize()

    def extract_element_data(self, elem_type: int, elem_tags: List[int], element_connectivities: List[int]) -> \
            Dict[str, object]:
        """
        Extracts element data from gmsh mesh

        Args:
            elem_type (int): Element type.
            elem_tags (List[int]): Element tags.
            element_connectivities (List[int]): Element node tags.

        Returns:
            dict: Dictionary which contains element data.
        """

        element_name = ElementType(elem_type).name

        n_nodes_per_element = self.get_num_nodes_from_elem_type(elem_type)
        num_elements = len(elem_tags)
        connectivities = np.reshape(element_connectivities, (num_elements, n_nodes_per_element))
        # print("element_name=",element_name)
        # print("n_nodes_per_element=",n_nodes_per_element)
        # print("num_elements=",num_elements)
        # print("elem_tags=",elem_tags)
        # print("connectivities=",connectivities)
        return {element_name: {"element_ids": elem_tags,
                               "connectivities": connectivities}}

    def extract_mesh_data(self, gmsh_mesh: Type[gmsh.model.mesh]):
        """
        Gets gmsh output data

        Args:
            gmsh_mesh (gmsh.model.mesh): The mesh as generated by gmsh.

        Returns:
            dict: Dictionary which contains nodal and elemental information.
        """

        mesh_data: Dict[str, Dict[str, object]] = {"nodes": {},
                                                   "elements": {}}

        # get nodal information
        node_tags, node_coords, node_params = gmsh_mesh.getNodes()  # nodes, elements

        # reshape nodal coordinate array to [num nodes, 3]
        num_nodes = len(node_tags)
        node_coordinates = np.reshape(node_coords, (num_nodes, 3))

        mesh_data["nodes"]["coordinates"] = node_coordinates
        mesh_data["nodes"]["ids"] = node_tags
        # print("node_coordinates=", node_coordinates)
        # print("node_tags=", node_tags)
        # get all elemental information
        elem_types, elem_tags, elem_node_tags = gmsh_mesh.getElements()

        # todo, this is unhandy for the future and the connection to kratos,
        #  handier would be to group elements by physical group
        for elem_type, elem_tag, elem_node_tag in zip(elem_types, elem_tags, elem_node_tags):
            element_dict = self.extract_element_data(elem_type, elem_tag, elem_node_tag)
            mesh_data["elements"].update(element_dict)

        return mesh_data
