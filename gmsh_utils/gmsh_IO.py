import pathlib
from typing import Dict, List, Union, Type, Any
from enum import Enum
import re

import gmsh
import numpy as np
import numpy.typing as npt


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
    geo_data : Dict
        Dictionary containing the geometry data, the geometry data contains: points, lines, surfaces, volumes
        and the physical groups.


    """

    def __init__(self):
        self.__mesh_data = {}
        self.__geo_data = {}

    @property
    def mesh_data(self) -> Dict[str, Dict[str, Any]]:
        """
        Returns the mesh data dictionary

        Returns:
            Dict: Dictionary containing the mesh data, i.e. nodal ids and coordinates; and elemental ids, connectivity's
            and element types.
        """

        return self.__mesh_data

    @mesh_data.setter
    def mesh_data(self, mesh_data: Dict[str, Dict[str, Any]]) -> None:
        """
        Sets the mesh data dictionary. For now, an exception is raised if this method is called, this is because the
        mesh data can only be set by internal methods

        Args:
            mesh_data (Dict): Dictionary containing the mesh data, i.e. nodal ids and coordinates; and elemental ids,
            connectivity's and element types.

        Raises:
            Exception: Mesh data can only be set by internal methods.

        Returns:
            None

        """
        raise Exception("Mesh data can only be set by internal methods.")

    @property
    def geo_data(self) -> Dict[str, Dict[str, Any]]:
        """
        Returns the geometry data dictionary

        Returns:
            Dict: Dictionary containing the geometry data, the geometry data contains: points, lines, surfaces, volumes
            and the physical groups.
        """

        return self.__geo_data

    @geo_data.setter
    def geo_data(self, geo_data: Dict[str, Dict[str, Any]]) -> None:
        """
        Sets the geometry data dictionary. For now, an exception is raised if this method is called, this is because the
        geometry data can only be set by internal method.

        Args:
            geo_data (Dict): Dictionary containing the geometry data, the geometry data contains: points, lines,
            surfaces, volumes and the physical groups.

        Raises:
            Exception: Geometry data can only be set by internal methods.

        Returns:
            None
        """

        raise Exception("Geometry data can only be set by internal methods.")

    def create_point(self, coordinates: Union[List[float], npt.NDArray[np.float64]], element_size: float) -> int:
        """
        Creates lines in gmsh.

        Args:
            coordinates (Union[List[float], npt.NDArray[np.float64]]): A list of point tags in order.
            element_size (float): The element size provided by user input.

        Returns:
            int: point tag
        """

        x = coordinates[0]
        y = coordinates[1]
        z = coordinates[2]
        point_id: int = gmsh.model.occ.addPoint(x, y, z, element_size)
        return point_id

    def create_line(self, point_ids: Union[List[int], npt.NDArray[np.int_]]) -> int:
        """
        Creates lines in gmsh.

        Args:
            point_ids (Union[List[int], npt.NDArray[int]]): A list of point tags in order.

        Returns:
            int: line tag
        """

        point1 = point_ids[0]
        point2 = point_ids[1]
        line_id: int = gmsh.model.occ.addLine(point1, point2)
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
                                           volume_ids: int) -> int:
        """
        Creates volume by extruding a 2D surface

        Args:
            surface_id (int): The surface tag.
            extrusion_length (Union[List[float], npt.NDArray[float]]): The extrusion length in x, y and z direction.
            name_label (str): The volume name label provided by user input
            volume_ids (int): The volume tag.
        Returns:
            int: volume tag
        """

        surface_ndim = 2
        volume_tag: int = gmsh.model.occ.extrude([(surface_ndim, surface_id)], extrusion_length[0], extrusion_length[1],
                               extrusion_length[2])
        gmsh.model.setPhysicalName(surface_ndim, volume_ids+1, name_label)
        return volume_tag

    def generate_point_pairs(self, point_ids) -> List[List[int]]:
        """
        Generates pairs of point IDs which form a line

        Args:

        Returns:
            List[List[int]]: A list of pairs of point IDs which create a line.
        """

        first_point_tag = 0
        point_pairs = []
        for point_id in range(len(point_ids) - 1):  # number of points in each layer (-1 for turning back)
            if point_id == 0:
                # saving the first point tag in order to return from last point to it
                first_point_tag = point_ids[point_id]
            # puts two consecutive points tags as the beginning and end of line in an array
            point_pair = [point_ids[point_id], point_ids[point_id+1]]
            point_pairs.append(point_pair)
        # make a pair that connects last point to first point
        last_point_tag = point_ids[len(point_ids)-1]
        point_pairs.append([last_point_tag, first_point_tag])

        return point_pairs

    def make_points(self, point_coordinates: Union[List[List[float]], npt.NDArray[np.float64]],
                    default_mesh_size:  float) -> List[int]:
        """
        Makes points with point tags by getting coordinates.

        Args:
            point_coordinates (Union[List[List[float]], npt.NDArray[np.float64]]): An Iterable of point x,y,z
                                                                                   coordinates.
            default_mesh_size (float): The element size.

        Returns:
            List[None]: A list of point tags.
        """

        list_point_ids = []
        for point in point_coordinates:
            point = [point[0], point[1], point[2]]
            point_id = self.create_point(point, default_mesh_size)
            list_point_ids.append(point_id)
        return list_point_ids

    def make_lines(self, point_pairs: Union[List[List[int]], npt.NDArray[np.int_]])\
            -> Union[List[int], npt.NDArray[np.int_]]:
        """
        Makes lines with line tags by getting point pairs.

        Args:
            point_pairs (Union[List[List[int]], npt.NDArray[np.int_]]): A list of pairs of point tags which create a
            line.

        Returns:
            List[List[int]]: A list of line tags.
        """

        list_lines = []
        for point_pair in range(len(point_pairs)):
            line = [point_pairs[point_pair][0], point_pairs[point_pair][1]]
            line_id = self.create_line(line)
            list_lines.append(line_id)
        return list_lines

    def make_surfaces(self, line_list: Union[List[int], npt.NDArray[np.int_]], name_label: str)\
            -> Union[List[int], npt.NDArray[np.int_]]:
        """
        Makes surfaces with surface tags by getting line tags.

        Args:
            line_list (Union[List[int], npt.NDArray[np.int_]]): A list of line tags in order.
            name_label (str): surface name labels provided by user input.

        Returns:
            List[int]: A list of surface tags.
        """

        surfaces = []
        surfaces.append(self.create_surface(line_list, name_label))

        return surfaces

    def make_volume(self, surface_id: Union[List[int], npt.NDArray[np.int_]],
                    extrusion_length: Union[List[float], npt.NDArray[np.float64]], name_label: str) -> None:
        """
        Makes volume with volume tags by getting surface tags.

        Args:
            surface_id (Union[List[int], npt.NDArray[np.int_]]): A list of surface tags.
            extrusion_length (Union[List[float], npt.NDArray[np.float64]]): The extrusion length in x, y and z direction
            name_label (str): A list of volume name labels provided by user input.

        Returns:
            None
        """

        volumes = []
        for volume in range(len(surface_id)):
            volume_tag = self.create_volume_by_extruding_surface(surface_id[volume], extrusion_length,
                                                                   name_label[volume], volume)
            volumes.append(volume_tag)

    def make_geometry_2d(self, point_coordinates: Union[List[List[float]], npt.NDArray[np.float64]],
                         default_mesh_size: float, name_label_list: str) -> Union[List[int], npt.NDArray[np.int_]]:
        """
        Takes point_pairs and puts their tags as the beginning and end of line in gmsh to create line,
        then creates surface to make 2D geometry.

        Args:
            point_coordinates (Union[List[List[float]], npt.NDArray[np.float64]]): A list of point coordinates.
            default_mesh_size (float): The default mesh size provided by user.
            name_label_list (List[str]): A list of surface name labels provided by user input.

        Returns:
            List[int]: List of Surface ids
        """

        list_point_id = self.make_points(point_coordinates, default_mesh_size)
        pair_lists = self.generate_point_pairs(list_point_id)
        line_lists = self.make_lines(pair_lists)
        surface_ids = self.make_surfaces(line_lists, name_label_list)

        return surface_ids

    def make_geometry_3d(self, point_coordinates: Union[List[List[float]], npt.NDArray[np.float64]],
                         default_mesh_size: float, name_label_list: str,
                         extrusion_length: Union[List[float], npt.NDArray[np.float64]], ) -> None:
        """
        Creates 3D geometries by extruding the 2D surface

        Args:
            point_coordinates (Union[List[List[float]], npt.NDArray[float]]): Geometry points coordinates.
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
                           mesh_size: float, dims: int, name_label: List[str],
                           mesh_name: str, mesh_output_dir: str, save_file: bool = False,
                           open_gmsh_gui: bool = False) -> None:
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
        self.generate_point_pairs(point_coordinates)

        gmsh.initialize()
        gmsh.model.add(mesh_name)

        for layer in range(len(point_coordinates)):
            # layer_list = self.prepare_inputs(point_coordinates)
            if dims == 3:
                self.make_geometry_3d(point_coordinates[layer], mesh_size, name_label[layer], extrusion_length)

            elif dims == 2:
                self.make_geometry_2d(point_coordinates[layer], mesh_size, name_label[layer])

        self.remove_duplicates()
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(dims)

        self.extract_mesh_data(gmsh.model.mesh)

        if save_file:
            # writes mesh file output in .msh format
            file_extension = ".msh"
            mesh_output_file = mesh_output_dir + mesh_name + file_extension
            gmsh.write(mesh_output_file)

        # opens Gmsh interface
        if open_gmsh_gui:
            gmsh.fltk.run()

        gmsh.finalize()

    def extract_node_data(self, node_tags: npt.NDArray[np.int_],
                          node_coordinates: npt.NDArray[np.float64]) \
            -> Dict[str, Union[npt.NDArray[np.int_], npt.NDArray[np.float64]]]:
        """
        Gets gmsh data belonging to nodal data

        Args:
            node_tags (npt.NDArray[np.int_]): gmsh node ids
            node_coordinates (npt.NDArray[float]) : gmsh node coordinates

        Returns:
            Dict[str, Union[npt.NDArray[np.int_], npt.NDArray[np.float64]]]: A dictionary containing node ids and
            coordinates

        """

        # reshape nodal coordinate array to [num nodes, 3]
        num_nodes = len(node_tags)
        node_coordinates = np.reshape(node_coordinates, (num_nodes, 3))

        return {"coordinates": node_coordinates,
                "ids": node_tags}

    def extract_elements_data(self, elem_types: npt.NDArray[np.int_], elem_tags: List[npt.NDArray[np.int_]],
                              elem_node_tags: List[npt.NDArray[np.int_]]) -> Dict[str, Dict[str, npt.NDArray[np.int_]]]:
        """
        Extracts element data from gmsh mesh

        Args:
            elem_types (npt.NDArray[np.int_]): Element types.
            elem_tags (List[npt.NDArray[np.int_]]): Element tags.
            elem_node_tags (List[npt.NDArray[np.int_]]): Element node tags.

        Returns:
            Dict (Dict[str, Dict[str, npt.NDArray[np.int_]]]): Dictionary which contains element data.

        """

        # initialize empty dictionary
        elements_data: Dict[str, Dict[str, npt.NDArray[np.int_]]] = {}

        # fill dictionary with element data
        for elem_type, elem_tag, elem_node_tag in zip(elem_types, elem_tags, elem_node_tags):
            element_dict = self.extract_element_data(elem_type, elem_tag, elem_node_tag)
            elements_data.update(element_dict)

        return elements_data

    def extract_element_data(self, elem_type: int, elem_tags: npt.NDArray[np.int_],
                             element_connectivities: npt.NDArray[np.int_]) -> \
            Dict[str, Dict[str, npt.NDArray[np.int_]]]:
        """
        Extracts element data from gmsh mesh
        Gets gmsh data belonging to a single element type

        Args:
            elem_type (int): Element type.
            elem_tags (npt.NDArray[np.int_]): Element ids.
            element_connectivities (npt.NDArray[np.int_]): Element node tags.

        Returns:
            dict: Dictionary which contains element data.
        """

        element_name = ElementType(elem_type).name
        n_nodes_per_element = self.get_num_nodes_from_elem_type(elem_type)
        num_elements = len(elem_tags)
        connectivities = np.reshape(element_connectivities, (num_elements, n_nodes_per_element))

        return {element_name: {"element_ids": elem_tags,
                               "connectivities": connectivities}}

    def extract_mesh_data(self, gmsh_mesh: Type[gmsh.model.mesh]):
        """
        Gets gmsh mesh data and stores it in a dictionary

        Args:
            gmsh_mesh (gmsh.model.mesh): The mesh as generated by gmsh.

        """

        mesh_data: Dict[str, Dict[str, Any]] = {"nodes": {},
                                                "elements": {}}

        # get nodal information
        node_tags, node_coords, node_params = gmsh_mesh.getNodes()  # nodes
        nodes_dict = self.extract_node_data(node_tags, node_coords)
        mesh_data["nodes"].update(nodes_dict)

        # get all elemental information
        elem_types, elem_tags, elem_node_tags = gmsh_mesh.getElements()

        # todo, this is unhandy for the future and the connection to kratos, handier would be to group
        #  elements by physical group
        mesh_data["elements"] = self.extract_elements_data(elem_types, elem_tags, elem_node_tags)

        self.__mesh_data = mesh_data

    def read_gmsh_msh(self, filename: str):
        """
        Reads a Gmsh .msh file and stores the data in a dictionary

        Args:
            filename (str): name of the Gmsh .msh file

        """
        gmsh.initialize()
        gmsh.open(filename)

        self.extract_mesh_data(gmsh.model.mesh)

        gmsh.finalize()

    def get_nodes_in_group(self, group_name: str) -> Dict[str, Union[npt.NDArray[np.int_], npt.NDArray[np.float64]]]:
        """
        Gets all nodes which are part of a certain group

        Args:
            group_name (str): Name of the requested group.

        Returns:
             Dict[str, Union[npt.NDArray[np.int_], npt.NDArray[np.float64]]]: Dictionary which contains nodal data.

        """

        groups = gmsh.model.getPhysicalGroups()

        for group in groups:
            name = gmsh.model.getPhysicalName(group[0], group[1])
            if name == group_name:
                # gets nodes per group
                nodes = gmsh.model.mesh.get_nodes_for_physical_group(group[0], group[1])

                nodes_data = self.extract_node_data(nodes[0], nodes[1])

                return nodes_data

        return {}

    def get_elements_in_group(self, group_name: str) -> Dict[str, Dict[str, npt.NDArray[np.int_]]]:
        """
        Gets all elements which are part of a certain group

        Args:
            group_name (str): Name of the requested group.

        Returns:
            Dict[str, npt.NDArray[np.int_]]: Dictionary which contains element data.

        """

        # Get group dimensions and ids
        groups = gmsh.model.getPhysicalGroups()

        for group in groups:

            # get name of the group
            name = gmsh.model.getPhysicalName(group[0], group[1])

            # if the requested group name is equal to the group, retrieve element data
            if name == group_name:
                # gets elements per group
                entity = gmsh.model.getEntitiesForPhysicalGroup(group[0], group[1])[0]
                elements = gmsh.model.mesh.getElements(entity)

                # extract element data
                element_data = self.extract_elements_data(elements[0], elements[1], elements[2])

                return element_data

        return {}

    def get_boundary_data(self, entity_ndim: int, entity_id: int) -> List[int]:
        """
        Gets lower entities of a certain entity, i.e. get surfaces from a volume; lines from a surface;
        points from a line.

        Args:
            entity_ndim (int): Dimension of the entity.
            entity_id (int): Id of the entity.

        Returns:
            List[int]: List of lower entities.

        """

        # get boundary entities of current entity
        lower_entities = gmsh.model.getBoundary([(entity_ndim, entity_id)], combined=False)

        # get ids of lower entities
        lower_entity_ids = [entity[1] for entity in lower_entities]

        return lower_entity_ids

    def extract_geo_data(self):
        """
        Extracts geometry data from gmsh model

        """

        # get all entities
        entities = gmsh.model.get_entities()

        geo_data: Dict[str, Dict[str, Any]] = {"points":   {},
                                               "lines":    {},
                                               "surfaces": {},
                                               "volumes":  {},
                                               "physical_groups": {}}

        # loop over all entities
        for entity in entities:

            # get dimension and id of entity
            entity_ndim, entity_id = entity[0], entity[1]

            # get point data
            if entity_ndim == 0:
                geo_data["points"][entity_id] = gmsh.model.get_value(entity_ndim, entity_id, [])
            # get line data
            if entity_ndim == 1:
                geo_data["lines"][entity_id] = self.get_boundary_data(entity_ndim, entity_id)
            # get surface data
            if entity_ndim == 2:
                geo_data["surfaces"][entity_id] = self.get_boundary_data(entity_ndim, entity_id)
            # get volume data
            if entity_ndim == 3:
                geo_data["volumes"][entity_id] = self.get_boundary_data(entity_ndim, entity_id)

        # Get group dimensions and ids
        groups = gmsh.model.getPhysicalGroups()

        # loop over all physical groups
        for group in groups:
            # get name of the group
            name = gmsh.model.getPhysicalName(group[0], group[1])

            # gets entity per group
            entity = gmsh.model.getEntitiesForPhysicalGroup(group[0], group[1])[0]

            # add group to dictionary
            geo_data["physical_groups"][name] = {"ndim": group[0],
                                                 "id": group[1],
                                                 "geometry_id": entity}

        self.__geo_data = geo_data

    def read_gmsh_geo(self, filename: str):
        """
        Reads a Gmsh .geo file and extracts the geometry data.

        Args:
            filename (str): Name of the Gmsh .geo file.

        Raises:
            FileNotFoundError: If the file does not exist.


        """

        if pathlib.Path(filename).exists():
            gmsh.initialize()
            gmsh.open(filename)

            self.extract_geo_data()

            gmsh.finalize()
        else:
            raise FileNotFoundError(f"File {filename} does not exist!")

    def generate_geo_from_geo_data(self):
        """
        Generates a Gmsh geometry data from a geometry data dictionary.

        The curve loops which form the surface are reoriented such that the surfaces are valid.

        """

        # initialize gmsh
        gmsh.initialize()

        # add points to the geometry
        for k, v in self.__geo_data["points"].items():
            gmsh.model.occ.addPoint(v[0], v[1], v[2], tag=k)

        # add lines to the geometry
        for k, v in self.__geo_data["lines"].items():
            print(self.__geo_data["lines"].items())
            gmsh.model.occ.addLine(v[0], v[1], tag=k)

        # add surfaces to the geometry
        for k, v in self.__geo_data["surfaces"].items():
            curve_key = gmsh.model.occ.addCurveLoop(v)
            gmsh.model.occ.addPlaneSurface([curve_key], tag=k)

        # add volumes to the geometry
        for k, v in self.__geo_data["volumes"].items():
            gmsh.model.occ.addSurfaceLoop(np.abs(v), tag=k)
            gmsh.model.occ.add_volume([k], tag=k)

        # add physical groups to the geometry
        for k, v in self.__geo_data["physical_groups"].items():
            gmsh.model.addPhysicalGroup(v["ndim"], [v["geometry_id"]], tag=v["id"], name=k)

        # synchronize the geometry for generating the mesh
        gmsh.model.occ.synchronize()

        # synchronize the geo geometry such that physical groups are added, important is that this is done after
        # synchronizing the occ geometry :-D
        gmsh.model.geo.synchronize()

    def generate_mesh(self, ndim: int, element_size: float = 0.0, order: int = 1):
        """
        Generates a mesh from the geometry data.

        Args:
            ndim (int): Dimension of the mesh.
            element_size (float, optional): Element size. Defaults to 0.0.
            order (int, optional): Order of the mesh. Defaults to 1.

        """

        # sets gmsh geometry from a geometry data dictionary
        self.generate_geo_from_geo_data()

        if element_size > 0.0:
            gmsh.model.mesh.setSize(gmsh.model.getEntities(), element_size)

        # set mesh order
        gmsh.model.mesh.setOrder(order)

        # generate mesh
        gmsh.model.mesh.generate(ndim)

        # parses gmsh mesh data into a mesh data dictionary
        self.extract_mesh_data(gmsh.model.mesh)

        # finalize gmsh
        gmsh.finalize()
