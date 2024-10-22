import pathlib
from typing import Dict, List, Union, Type, Any, Sequence, Tuple
from enum import Enum
import re

import gmsh
import numpy as np

from gmsh_utils.math_utils import MathUtils
# force gmsh io to not use numpy for gmsh
gmsh.use_numpy = False


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

    Attributes:
        - mesh_data (dict): Dictionary containing the mesh data, i.e. nodal ids and coordinates; and elemental ids, \
            connectivity's and element types.
        - geo_data (dict): Dictionary containing the geometry data, the geometry data contains: points, lines, \
            surfaces, volumes and the physical groups.

    """

    def __init__(self):
        """
        Constructor of GmshIO class

        """
        self.__mesh_data = {}
        self.__geo_data: Dict[str, Any] = {"points": {},
                                           "lines": {},
                                           "surfaces": {},
                                           "volumes": {},
                                           "physical_groups": {}}

    @property
    def mesh_data(self) -> Dict[str, Any]:
        """
        Returns the mesh data dictionary

        Returns:
            - Dict[str, Any]: Dictionary containing the mesh data, i.e. nodal ids and coordinates; and elemental ids, \
                connectivity's and element types.
        """

        return self.__mesh_data

    @mesh_data.setter
    def mesh_data(self, mesh_data: Dict[str, Any]):
        """
        Sets the mesh data dictionary. For now, an exception is raised if this method is called, this is because the
        mesh data can only be set by internal methods

        Args:
            - mesh_data (Dict[str, Any]): Dictionary containing the mesh data, i.e. nodal ids and coordinates; and \
                elemental ids, connectivity's and element types.

        Raises:
            - Exception: Mesh data can only be set by internal methods.
        """

        raise Exception("Mesh data can only be set by internal methods.")

    @property
    def geo_data(self) -> Dict[str, Dict[Any, Any]]:
        """
        Returns the geometry data dictionary

        Returns:
            - Dict[str, Dict[Any, Any]]: Dictionary containing the geometry data, the geometry data contains: points,
            lines, surfaces, volumes and the physical groups.
        """

        return self.__geo_data

    @geo_data.setter
    def geo_data(self, geo_data: Dict[str, Dict[Any, Any]]):
        """
        Sets the geometry data dictionary. For now, an exception is raised if this method is called, this is because the
        geometry data can only be set by internal method.

        Args:
            - geo_data (Dict[str, Dict[Any, Any]]): Dictionary containing the geometry data, the geometry data contains:
            points, lines, surfaces, volumes and the physical groups.

        Raises:
            - Exception: Geometry data can only be set by internal methods.
        """

        raise Exception("Geometry data can only be set by internal methods.")

    def create_point(self, coordinates: Sequence[float], mesh_size: float = -1) -> int:
        """
        Creates points in gmsh.

        Args:
            - coordinates (Sequence[float]): A list of point tags in order.
            - mesh_size (float): The element size provided by user input (default is -1, \
                which let gmsh choose the size).

        Returns:
            - int: point tag
        """

        x = coordinates[0]
        y = coordinates[1]
        z = coordinates[2]

        # define typing point id
        point_id: int

        if [x, y, z] in self.__geo_data["points"].values():
            point_id = next((point_tag for point_tag, point_coordinates in self.__geo_data["points"].items() if
                             point_coordinates == [x, y, z]))
        else:
            point_id = gmsh.model.occ.addPoint(x, y, z, mesh_size)

        return point_id

    def create_line(self, point_ids: Sequence[int]) -> int:
        """
        Creates lines in gmsh.

        Args:
            - point_ids (Sequence[int]): A list of point tags in order.

        Returns:
            - int: line tag
        """

        point1 = point_ids[0]
        point2 = point_ids[1]
        line_id: int = gmsh.model.occ.addLine(point1, point2)
        return line_id

    def create_surface(self, line_ids: Sequence[int]) -> int:
        """
        Creates curve and then surface in gmsh by using line tags.

        Args:
            - line_ids (Sequence[int]): A list of line tags in order.

        Returns:
            - int: surface id
        """

        curve_loop_id = gmsh.model.occ.addCurveLoop(line_ids)
        surface_id: int = gmsh.model.occ.addPlaneSurface([curve_loop_id])

        return surface_id

    def create_volume_by_extruding_surface(self, surface_id: int, extrusion_length: Sequence[float]) -> int:
        """
        Creates volume by extruding a 2D surface

        Args:
            - surface_id (int): The surface tag.
            - extrusion_length (Sequence[float]): The extrusion length in x, y and z direction.
        Returns:
            - int: volume tag
        """

        surface_dim = 2
        volume_dim = 3
        new_dim_tags = gmsh.model.occ.extrude([(surface_dim, surface_id)], extrusion_length[0], extrusion_length[1],
                                              extrusion_length[2])
        # gets the first volume tag from the list of new dimension tags
        volume_tag: int = next((dim_tag[1] for dim_tag in new_dim_tags if dim_tag[0] == volume_dim))

        return volume_tag

    def __generate_point_pairs_for_closed_loop(self, point_ids: List[int]) -> List[List[int]]:
        """
        Generates list of consecutive pairs of point tags which is needed to form the lines and closed surfaces.

        Args:
            - point_ids (List[int]): A list of two ordered point tags

        Returns:
            - List[List[int]]: A list of pairs of point tags which is needed to create a line
        """

        first_point_tag = 0
        point_pairs = []
        for i in range(len(point_ids) - 1):  # number of points in each layer (-1 for turning back)
            if i == 0:
                # saving the first point tag in order to return from last point to it
                first_point_tag = point_ids[i]
            # puts two consecutive points tags as the beginning and end of line in an array
            point_pair = [point_ids[i], point_ids[i + 1]]
            point_pairs.append(point_pair)
        # make a pair that connects last point to first point
        last_point_tag = point_ids[len(point_ids) - 1]
        point_pairs.append([last_point_tag, first_point_tag])

        return point_pairs

    def make_points(self, point_coordinates: Sequence[Sequence[float]], element_size: float = -1) -> List[int]:
        """
        Makes points with point tags from coordinates.

        Args:
           - point_coordinates (Sequence[Sequence[float]]): An Iterable of point x,y,z coordinates.
           - element_size (float): The element size (default is -1, which let gmsh choose the size).

        Returns:
            - List[int]: A list of point tags.
        """

        list_point_ids = [self.create_point(point, element_size) for point in point_coordinates]
        return list_point_ids

    def make_lines(self, point_pairs: Sequence[Sequence[int]]) -> List[int]:
        """
        Makes lines with line tags by getting point pairs.

        Args:
            - point_pairs (Sequence[Sequence[int]]): A sequence of pairs of point tags which create a line.

        Returns:
            - List[int]: A list of line tags.
        """

        line_ids = [self.create_line(point_pair) for point_pair in point_pairs]
        return line_ids

    def make_geometry_0d(self, point_coordinates: Sequence[Sequence[float]], name_label: str = "",
                         element_size: float = -1) -> List[int]:
        """
        Makes 0D geometry by creating points in gmsh.

        Args:
            - point_coordinates (Sequence[Sequence[float]]): A sequence of point x,y,z coordinates.
            - name_label (str): The name of the physical group containing the points.
            - element_size (float): The element size (default is -1, which let gmsh choose the size).

        Returns:
            - List[int]: A list of point tags.
        """

        point_ids = self.make_points(point_coordinates, element_size)

        # only add physical group if name label is not empty
        if name_label != "":
            point_ndim = 0
            self.__add_or_append_to_physical_group(name_label, point_ndim, point_ids)

        return point_ids

    def make_geometry_1d(self, point_coordinates: Sequence[Sequence[float]],
                         name_label: str = "", element_size: float = -1) -> List[int]:
        """
        Makes 1D geometry by creating points and lines in gmsh.

        Args:
            - point_coordinates (Sequence[Sequence[float]]): A sequence of point x,y,z coordinates.
            - name_label (str): The name of the physical group containing the lines.
            - element_size (float): The element size (default is -1, which let gmsh choose the size).

        Returns:
            - List[int]: A list of line tags.
        """

        # create point ids
        point_ids = self.make_points(point_coordinates, element_size=element_size)

        # create lines
        line_ids = [self.create_line([point_ids[i], point_ids[i + 1]]) for i in range(len(point_ids) - 1)]

        # only add physical group if name label is not empty
        if name_label != "":
            line_ndim = 1
            self.__add_or_append_to_physical_group(name_label, line_ndim, line_ids)

        return line_ids

    def __generate_closed_line_loop(self, point_coordinates: Sequence[Sequence[float]],
                                    element_size: float = -1) -> List[int]:
        """
        Generates a closed line loop from a list of point coordinates.

        Args:
            - point_coordinates (Sequence[Sequence[float]]): A sequence of point x,y,z coordinates.
            - element_size (float): The element size.

        Returns:
            - List[int]: A list of line tags.
        """

        point_ids = self.make_points(point_coordinates, element_size=element_size)
        pair_list = self.__generate_point_pairs_for_closed_loop(point_ids)
        line_ids = self.make_lines(pair_list)
        return line_ids

    def make_geometry_2d(self, point_coordinates: Sequence[Sequence[float]],
                         name_label: str = "", element_size: float = -1) -> int:
        """
        Takes point_pairs and puts their tags as the beginning and end of line in gmsh to create line,
        then creates surface to make 2D geometry.

        Args:
            - point_coordinates (Sequence[Sequence[float]]): A list of point coordinates.
            - name_label (str): A name label provided for the volume by user input.
            - element_size (float): The default mesh size provided by user.

        Returns:
            - int: Surface id
        """

        lined_ids = self.__generate_closed_line_loop(point_coordinates, element_size=element_size)
        surface_id = self.create_surface(lined_ids)

        # only add physical group if name label is not empty
        if name_label != "":
            surface_ndim = 2
            self.__add_or_append_to_physical_group(name_label, surface_ndim, [surface_id])

        return surface_id

    def make_geometry_3d_by_extrusion(self, point_coordinates: Sequence[Sequence[float]],
                                      extrusion_length: Sequence[float], name_label: str = "",
                                      element_size: float = -1) -> int:
        """
        Creates 3D geometries by extruding the 2D surface

        Args:
            - point_coordinates (Sequence[Sequence[float]]): Geometry points coordinates.
            - extrusion_length (Sequence[float]): The extrusion length in x, y and z direction.
            - name_label (str): A name label provided for the volume by user input.
            - element_size (float): The default mesh size provided by user.

        Returns:
            - int: Volume id
        """

        # create surface without a name label, which is used for extrusion
        surface_id = self.make_geometry_2d(point_coordinates, element_size=element_size)
        volume_id = self.create_volume_by_extruding_surface(surface_id, extrusion_length)

        if name_label != "":
            volume_dim = 3
            self.__add_or_append_to_physical_group(name_label, volume_dim, [volume_id])

        return volume_id

    def remove_duplicate_entities_from_geometry(self):
        """
        Removes duplicate entities from the geometry.

        """

        gmsh.model.occ.removeAllDuplicates()

    @staticmethod
    def get_num_nodes_from_elem_type(elem_type: int) -> int:
        """
        Gets number of nodes from element types

        Args:
            - elem_type (int): An integer that defines the type of element.

        Returns:
            - int: The number of nodes needed for a type of element.
        """

        # get name from element type enum
        element_name = ElementType(elem_type).name

        # get number of nodes from the enum name
        num_nodes = int(re.findall(r"\d+", element_name)[0])

        return num_nodes

    @staticmethod
    def validate_layer_parameters(layer_parameters: Dict[str, Any]):
        """
        Validates the layer parameters

        Args:
            - layer_parameters (Dict[str, Any]): A dictionary containing the layer information.

        """

        for key, layer in layer_parameters.items():
            if "coordinates" not in layer:
                raise ValueError(f"Layer {key} must contain the key 'coordinates'")
            if "ndim" not in layer:
                raise ValueError(f"Layer {key} must contain the key 'ndim'")

            if layer["ndim"] == 3 and "extrusion_length" not in layer:
                raise ValueError(f"Layer {key} must contain the key 'extrusion_length', which is needed "
                                 f"for 3D geometries.")

            if "element_size" not in layer:
                layer["element_size"] = -1
                print(f"Warning: Layer {key} does not contain the key 'element_size'. The element size will be "
                      "determined by gmsh.")

            if layer["ndim"] not in [0, 1, 2, 3]:
                raise ValueError(f"ndim must be 0, 1, 2 or 3. ndim={layer['ndim']}")

    def generate_geometry(self, layer_parameters: Dict[str, Any], model_name: str):
        """
        Generates the geometry. Gmsh is initialized if it is not already initialized. Then for each layer in the
        layer_settings dictionary, the geometry is created. The layer_settings dictionary must contain the following
        keys:
            - name: The name of the layer.\

        Per item in the layer_settings dictionary, the following keys are required:
            - coordinates: A list of point coordinates which make up the layer.\
            - ndim: The number of dimensions of the layer geometry.\
            - in 3D, extrusion_length: The extrusion length in x, y and z direction.\
            - Optional[element_size]: The element size. If not provided, the element size is determined by gmsh.\

        Args:
            - layer_parameters (Dict[str, Any]): A dictionary containing the layer information.
            - model_name (str): Name of gmsh model and mesh output file.
        """

        # validate the layer_dictionary
        self.validate_layer_parameters(layer_parameters)

        if not gmsh.isInitialized():
            gmsh.initialize()
            gmsh.model.add(model_name)

        for layer_name, layer in layer_parameters.items():
            ndim = layer["ndim"]

            if ndim == 0:
                self.make_geometry_0d(layer["coordinates"], layer_name, layer["element_size"])
            elif ndim == 1:
                self.make_geometry_1d(layer["coordinates"], layer_name, layer["element_size"])
            elif ndim == 2:
                self.make_geometry_2d(layer["coordinates"], layer_name, layer["element_size"])
            elif ndim == 3:
                self.make_geometry_3d_by_extrusion(layer["coordinates"], layer["extrusion_length"], layer_name,
                                                   layer["element_size"])

        self.remove_duplicate_entities_from_geometry()
        self.synchronize_gmsh()

        self.extract_geo_data()

        # add element size to geo data
        for layer_name, layer in layer_parameters.items():
            if "element_size" in layer and layer["element_size"] > 0:
                self.__geo_data["physical_groups"][layer_name]["element_size"] = layer["element_size"]

    def set_mesh_size_of_group(self, group_name: str, mesh_size: float) -> None:
        """
        Sets the customized mesh size of a specific group

        Args:
             - group_name (str): The name of the group provided by user input
             - mesh_size (float): The customized mesh size provided by user for a specific group.

        Returns:
            - None
        """

        if mesh_size <= 0:
            raise ValueError(f"The mesh size of {group_name} is smaller than or equal to zero. "
                             f"Please provide a mesh size larger than zero.")
        group = self.geo_data["physical_groups"][group_name]

        # add data to physical group
        group["element_size"] = mesh_size

        geometry_ids = group["geometry_ids"]
        for geometry_id in geometry_ids:
            ndim = group["ndim"]
            if ndim == 0:
                gmsh.model.mesh.setSize([(ndim, geometry_id)], mesh_size)
            else:
                entities_list = gmsh.model.getBoundary([(ndim, geometry_id)], recursive=True)
                gmsh.model.mesh.setSize(entities_list, mesh_size)

    def create_node_data_dict(self, node_tags: List[int],
                              node_coordinates: List[float]) -> Dict[int, Sequence[float]]:
        """
        Creates a dictionary containing node ids as keys and coordinates as values

        Args:
           - node_tags (List[int]): gmsh node ids
           - node_coordinates (List[float]) : gmsh node coordinates

        Returns:
           - Dict[int, Sequence[float]]: A dictionary containing node ids and
            coordinates

        """

        # reshape nodal coordinate array to [num nodes, 3]
        num_nodes = len(node_tags)
        ndim = 3

        # create dictionary of nodal data with node ids as keys and coordinates as values
        nodal_data: Dict[int, Any] = {}
        for i in range(num_nodes):
            nodal_data[node_tags[i]] = node_coordinates[i * ndim: (i + 1) * ndim]

        return nodal_data

    def extract_all_elements_data(self, elem_types: List[int], elem_tags: List[List[int]],
                                  elem_node_tags: List[List[int]]) -> Dict[str, Any]:
        """
        Extracts element data for all the elements in the gmsh mesh. Element data is defined as:
            - element type
            - element ids
            - element node connectivities

        Args:
            - elem_types (List[int]): Element types.
            - elem_tags (List[List[int]]): Element tags.
            - elem_node_tags (List[List[int]]): Element node tags.

        Returns:
            - Dict (Dict[str, Any]): Dictionary which contains element data.

        """

        # initialize empty dictionary
        elements_data: Dict[str, Any] = {}

        # fill dictionary with element data
        for elem_type, elem_tag, elem_node_tag in zip(elem_types, elem_tags, elem_node_tags):
            element_dict: Dict[str, Any] = self.extract_element_data(elem_type, elem_tag, elem_node_tag)
            elements_data.update(element_dict)

        return elements_data

    def extract_element_data(self, elem_type: int, elem_tags: List[int],
                             element_connectivities: List[int]) -> \
            Dict[str, Any]:
        """
        Extracts element data from gmsh mesh belonging to a single element type. Element data is defined as:
            - element type
            - element ids
            - element node connectivities

        Args:
            - elem_type (int): Element type.
            - elem_tags (List[int]): Element ids.
            - element_connectivities (List[int]): Element node tags.

        Returns:
            - Dict: Dictionary which contains element data.
        """

        element_name = ElementType(elem_type).name
        n_nodes_per_element = self.get_num_nodes_from_elem_type(elem_type)
        num_elements = len(elem_tags)

        # add element data to dictionary with key = element id and value = element connectivity
        element_data: Dict[int, Any] = {}
        for i in range(num_elements):
            element_data[elem_tags[i]] = element_connectivities[i * n_nodes_per_element: (i + 1) * n_nodes_per_element]

        return {element_name: element_data}

    def extract_mesh_data(self):
        """
        Gets gmsh mesh data and stores it in a dictionary. The dictionary contains nodal data, elemental data and
        physical group data. Each physical group contains the node ids and element ids, which are part of the group

        """

        mesh_data: Dict[str, Any] = {"ndim": gmsh.model.getDimension(),
                                     "nodes": {},
                                     "elements": {},
                                     "physical_groups": {}}

        # get nodal information
        node_tags, node_coords, node_params = gmsh.model.mesh.getNodes()
        mesh_data["nodes"] = self.create_node_data_dict(node_tags, node_coords)

        # get all elemental information
        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements()
        mesh_data["elements"] = self.extract_all_elements_data(elem_types, elem_tags, elem_node_tags)

        # get all physical group information
        physical_groups = gmsh.model.getPhysicalGroups()

        # loop over all physical groups
        for group_dim, group_id in physical_groups:
            # get name of the group
            name = gmsh.model.getPhysicalName(group_dim, group_id)

            # get the node ids belonging to the group
            node_ids = gmsh.model.mesh.get_nodes_for_physical_group(group_dim, group_id)[0]

            # gets elements per group
            entities = gmsh.model.getEntitiesForPhysicalGroup(group_dim, group_id)

            element_ids = []
            element_type = None
            for entity in entities:

                # gets element ids belonging to physical group
                element_ids.extend(gmsh.model.mesh.getElements(dim=group_dim, tag=entity)[1][0])

                # gets element type of elements in physical group, note that gmsh makes sure that all elements in a
                # physical group are of the same type
                element_type = gmsh.model.mesh.getElements(dim=group_dim, tag=entity)[0][0]

            # store group information in dictionary
            mesh_data["physical_groups"][name] = {"ndim": group_dim,
                                                  "node_ids": node_ids,
                                                  "element_ids": element_ids}

            # store element type in dictionary if it exists
            if element_type is not None:
                mesh_data["physical_groups"][name]["element_type"] = ElementType(element_type).name

        self.__mesh_data = mesh_data

    def read_gmsh_msh(self, filename: str):
        """
        Reads a Gmsh .msh file and stores the data in a dictionary

        Args:
            - filename (str): name of the Gmsh .msh file

        """

        # check if file exists
        if not pathlib.Path(filename).exists():
            raise FileNotFoundError(f"File {filename} not found")

        self.reset_gmsh_instance()

        gmsh.open(filename)

        self.extract_mesh_data()

        self.finalize_gmsh()

    def get_boundary_data(self, entity_ndim: int, entity_id: int) -> List[int]:
        """
        Gets lower entities of a certain entity, i.e. get surfaces from a volume; lines from a surface;
        points from a line.

        Args:
            - entity_ndim (int): Dimension of the entity.
            - entity_id (int): Id of the entity.

        Returns:
            - List[int]: List of lower entities.

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
        entities = gmsh.model.occ.getEntities()

        # if no occ entities are found, it means that geo_data is read from a geo file, in this case entities are stored
        # directly on the gmsh model
        if len(entities) == 0:
            entities = gmsh.model.get_entities()


        geo_data: Dict[str, Dict[str, Any]] = {"points": {},
                                               "lines": {},
                                               "surfaces": {},
                                               "volumes": {},
                                               "physical_groups": {}}

        #todo, it is possible that duplicated points exist

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
            entities = gmsh.model.getEntitiesForPhysicalGroup(group[0], group[1])

            # add group to dictionary
            geo_data["physical_groups"][name] = {"ndim": group[0],
                                                 "id": group[1],
                                                 "geometry_ids": entities}

            # add previously defined element size to physical group
            if name in self.geo_data["physical_groups"]:
                if "element_size" in self.geo_data["physical_groups"][name]:
                    geo_data["physical_groups"][name]["element_size"] = (
                        self.geo_data)["physical_groups"][name]["element_size"]

        self.__geo_data = geo_data

    def read_gmsh_geo(self, filename: str):
        """
        Reads a Gmsh .geo file and extracts the geometry data.

        Args:
            - filename (str): Name of the Gmsh .geo file.

        Raises:
            FileNotFoundError: If the file does not exist.

        """

        if pathlib.Path(filename).exists():
            self.reset_gmsh_instance()

            gmsh.open(filename)

            self.extract_geo_data()

            self.finalize_gmsh()
        else:
            raise FileNotFoundError(f"File {filename} does not exist!")

    def generate_geo_from_geo_data(self):
        """
        Generates a Gmsh geometry data from a geometry data dictionary.

        The curve loops which form the surface are reoriented such that the surfaces are valid.

        """

        # reset gmsh
        self.reset_gmsh_instance()

        # add points to the geometry
        for k, v in self.__geo_data["points"].items():
            gmsh.model.occ.addPoint(v[0], v[1], v[2], tag=k)

        # add lines to the geometry
        for k, v in self.__geo_data["lines"].items():
            gmsh.model.occ.addLine(v[0], v[1], tag=k)

        # add surfaces to the geometry
        for k, v in self.__geo_data["surfaces"].items():
            curve_key = gmsh.model.occ.addCurveLoop(v)
            gmsh.model.occ.addPlaneSurface([curve_key], tag=k)

        # add volumes to the geometry
        for k, v in self.__geo_data["volumes"].items():
            gmsh.model.occ.addSurfaceLoop([abs(surface_id) for surface_id in v], tag=k)
            gmsh.model.occ.add_volume([k], tag=k)

        # add physical groups to the geometry
        for k, v in self.__geo_data["physical_groups"].items():
            gmsh.model.addPhysicalGroup(v["ndim"], v["geometry_ids"], tag=v["id"], name=k)

        self.synchronize_gmsh()

    def __add_or_append_to_physical_group(self, name: str, ndim: int, geometry_ids: Sequence[int]):
        """
        Adds or appends geometry ids to a physical group. If the physical group does not exist, it is created.
        If the physical group already exists, the dimension of the physical group is checked. If the dimension of the
        physical group is different from the dimension of the geometry ids, an error is raised. Otherwise, the
        geometry ids are appended to the existing physical group.

        Args:
            - name (str): Name of the physical group.
            - ndim (int): Dimension of the physical group.
            - geometry_ids (Sequence[int]): Sequence of geometry ids belonging to the physical group.

        Raises:
            - ValueError: If the dimension of the existing physical group is different from the dimension of the
            new physical group in case the names of the physical groups are the same.
        """

        if "physical_groups" in self.__geo_data.keys() and name in self.__geo_data["physical_groups"]:
            ndim_existing_group = self.__geo_data["physical_groups"][name]["ndim"]
            if ndim != ndim_existing_group:
                raise ValueError(f"Cannot add geometry ids to physical group {name} with dimension {ndim} as the "
                                 f"physical group already exists with dimension "
                                 f"{ndim_existing_group}.")

            existing_geometry_ids = self.__geo_data["physical_groups"][name]["geometry_ids"]
            new_geometry_ids = existing_geometry_ids + list(geometry_ids)

            # remove existing physical group
            gmsh.model.removePhysicalGroups([(ndim_existing_group, self.__geo_data["physical_groups"][name]["id"])])

            # re-add new physical group
            gmsh.model.addPhysicalGroup(ndim_existing_group, new_geometry_ids,
                                        tag=self.__geo_data["physical_groups"][name]["id"], name=name)

        else:
            # add physical group to the geometry
            gmsh.model.addPhysicalGroup(ndim, geometry_ids, name=name)

    def add_physical_group(self, name: str, ndim: int, geometry_ids: Sequence[int]):
        """
        Adds a physical group to the existing geometry.

        Args:
            - name (str): Name of the physical group.
            - ndim (int): Dimension of the physical group.
            - geometry_ids (Sequence[int]): Sequence of geometry ids belonging to the physical group.
        """

        # add physical group to the geometry
        self.__add_or_append_to_physical_group(name, ndim, geometry_ids)

        # synchronize the geometry
        self.synchronize_gmsh()

        # extract the geometry data
        self.extract_geo_data()

    def get_surface_ids_at_plane(self, plane_vertices: Sequence[Sequence[float]]) -> List[int]:
        """
        Gets surface ids at a plane defined by three vertices.

        Args:
            - plane_vertices (Sequence[Sequence[float]]): A list of three vertices defining the plane.

        Returns:
            - List[int]: A list of surface ids at the plane.
        """

        surface_ids_at_plane = []

        # get normal of the plane
        plane_vertices_array  = np.array(plane_vertices, dtype=float)
        normal_plane = MathUtils.calculate_normal_polygon(plane_vertices_array)

        # check all surfaces
        for surface_id in self.geo_data["surfaces"].keys():

            # get point ids of the surface
            surface_point_dim_ids = gmsh.model.getBoundary([(2, surface_id)], combined=False, recursive=True)

            # get surface coordinates
            surface_coordinates = np.array([self.geo_data["points"][point_id] for _, point_id in surface_point_dim_ids])

            # check if all surface coordinates are on the plane, if so add surface id to list
            if all(MathUtils.is_point_on_plane(point, plane_vertices_array[0], normal_plane)
                   for point in surface_coordinates):
                surface_ids_at_plane.append(surface_id)

        return surface_ids_at_plane

    def get_surface_ids_at_polygon(self, polygon_vertices: Sequence[Sequence[float]]) -> List[int]:
        """
        Gets all surface ids at a convex or concave polygon

        Args:
            - polygon_vertices (Sequence[Sequence[float]]): A list of all vertices defining the polygon.

        Returns:
            - List[int]: A list of surface ids at the plane.
        """

        # check if all surface_vertices lie on the same plane
        surface_ids_at_plane = []

        for surface_id in self.geo_data["surfaces"].keys():

            # get point ids of the surface
            surface_point_dim_ids = gmsh.model.getBoundary([(2, surface_id)], combined=False, recursive=True)

            # get surface coordinates
            surface_coordinates = [self.geo_data["points"][point_id] for _, point_id in surface_point_dim_ids]

            # check if all surface coordinates are in the polygon, if so add surface id to list
            if all(MathUtils.is_point_in_polygon(point, polygon_vertices) for point in surface_coordinates):
                surface_ids_at_plane.append(surface_id)

        return surface_ids_at_plane

    def __sort_groups_by_element_size(self) -> Tuple[List[Dict[str, Any]], List[str]]:
        """
        Sorts physical groups by element size in descending order. Physical groups without element size information
        are sorted by their original order and are placed at the beginning of the list.

        Returns:
            - Tuple[List[Dict[str, Any]], List[str]]: Tuple of sorted groups and sorted group names.
        """

        # get all physical groups from geo data
        group_names, groups = self.geo_data["physical_groups"].keys(), self.geo_data["physical_groups"].values()

        # Filter groups with element_size information
        groups_with_info = [(name, group) for name, group in zip(group_names, groups) if "element_size" in group]

        if len(list(groups_with_info)) > 0:

            # Sort groups by element_size in descending order
            sorted_groups_with_info: List[Tuple[str, Dict[Any, Any]]] = sorted(groups_with_info,
                                                                               key=lambda x:
                                                                               x[1].get("element_size", -1),
                                                                               reverse=True)

            sorted_group_names_with_info, sorted_groups_with_info = zip(*sorted_groups_with_info)

            # Find groups without element_size information
            groups_without_info = [(name, group) for name, group in zip(group_names, groups) if
                                   "element_size" not in group]

            # Combine groups with and without info, preserving their original order
            sorted_groups = [group for _, group in groups_without_info] + list(sorted_groups_with_info)
            sorted_group_names = [name for name, _ in groups_without_info] + list(sorted_group_names_with_info)

        else:
            sorted_groups, sorted_group_names = list(groups), list(group_names)

        return sorted_groups, sorted_group_names

    def generate_mesh(self, ndim: int, element_size: float = -1, order: int = 1, save_file: bool = False,
                      mesh_name: str = "mesh_file", mesh_output_dir: str = "./", open_gmsh_gui: bool = False):
        """
        Generates a mesh from the geometry data.

        Args:
            - ndim (int): Dimension of the mesh.
            - element_size (float): Element size. Defaults to -1, which lets gmsh choose the element size.
            - order (int, optional): Order of the mesh. Defaults to 1.
            - save_file (bool): If True, saves mesh data to gmsh msh file. Defaults to False.
            - mesh_name (str): Name of gmsh model and mesh output file. Defaults to `mesh_file`.
            - mesh_output_dir (str): Output directory of mesh file. Defaults to working directory.
            - open_gmsh_gui (bool): User indicates whether to open gmsh interface. Defaults to False.

        """

        # sets gmsh geometry from a geometry data dictionary
        self.generate_geo_from_geo_data()

        # sort physical groups by element size in descending order
        sorted_groups, sorted_group_names = self.__sort_groups_by_element_size()

        # set mesh size per group, if element size is provided in group, else use the global element size
        for group_name, group in zip(sorted_group_names, sorted_groups):
            if "element_size" in group:
                self.set_mesh_size_of_group(group_name, group["element_size"])
            elif element_size > 0.0:
                self.set_mesh_size_of_group(group_name, element_size)

        # set mesh order
        gmsh.model.mesh.setOrder(order)

        # generate mesh
        gmsh.model.mesh.generate(ndim)

        # parses gmsh mesh data into a mesh data dictionary
        self.extract_mesh_data()

        if save_file:
            # writes mesh file output in .msh format

            # create directory if it does not exist
            pathlib.Path(mesh_output_dir).mkdir(parents=True, exist_ok=True)
            mesh_output_file = (pathlib.Path(mesh_output_dir) / mesh_name).with_suffix(".msh")
            gmsh.write(str(mesh_output_file))

        # opens Gmsh interface
        if open_gmsh_gui:
            gmsh.fltk.run()

        # finalize gmsh
        self.finalize_gmsh()

    @staticmethod
    def finalize_gmsh():
        """
        Finalizes gmsh.

        """

        if gmsh.isInitialized():
            gmsh.finalize()

    def __synchronize_intersection(self):
        """
        Synchronizes the intersection of entities. In the intersection of entities, the entities are split into
        smaller entities. Both of the new entities are then manually added to the existing physical group, as gmsh
        does not do this automatically.

        """

        # in 3D geometries, firstly intersect surfaces with lines
        if gmsh.model.getDimension() == 3:

            # get all entities
            occ_entities = gmsh.model.occ.get_entities()

            # get all lines and surfaces
            lines = [entity for entity in occ_entities if entity[0] == 1]
            surfaces = [entity for entity in occ_entities if entity[0] == 2]

            # intersect all lines with all surfaces
            new_entities, new_entities_map = gmsh.model.occ.fragment(surfaces, lines, removeTool=True,
                                                                     removeObject=True)

            # create the new entities map in the correct order, including all points, lines, surfaces and volumes
            new_surfaces = new_entities_map[:len(surfaces)]
            new_lines = new_entities_map[len(surfaces):]

            new_occ_entities = gmsh.model.occ.get_entities()
            points = [[entity] for entity in new_occ_entities if entity[0] == 0]
            volumes = [[entity] for entity in occ_entities if entity[0] == 3]
            new_entities_map = points + new_lines + new_surfaces + volumes

            # in the original entities list, exchange original points with new points, other entities remain the same
            unique_original_entities = ([entity for entity in new_occ_entities if entity[0] == 0] +
                                        [entity for entity in occ_entities if entity[0] != 0])

            # re-add physical groups on split entities
            self.__readd_physical_group_on_split_entities(unique_original_entities, new_entities_map)

        # get all entities
        occ_entities = gmsh.model.occ.get_entities()

        # intersect all entities with each other
        new_entities, new_entities_map = gmsh.model.occ.fragment(occ_entities, occ_entities, removeTool=True,
                                                                 removeObject=True)

        # the new entities map duplicates the entities, therefore only the first half of the entities map is used
        filtered_entities_map = new_entities_map[: len(occ_entities)]

        # re-add physical groups on split entities
        if len(filtered_entities_map) > 0:
            self.__readd_physical_group_on_split_entities(occ_entities, filtered_entities_map)

    def __readd_physical_group_on_split_entities(self, original_entities: List[Tuple[int, int]],
                                                 split_entities: List[List[Tuple[int, int]]]):
        """
        Re-adds physical groups on split entities. The physical groups are removed and then re-added with the new
        entities.

        Args:
            - original_entities (List[Tuple[int, int]]): List of original entities.
            - split_entities (List[List[Tuple[int, int]]]): List of the new split entities.
        """

        # get all physical groups
        physical_groups = gmsh.model.getPhysicalGroups()

        # loop over the filtered physical groups
        for group_dim, group_id in physical_groups:
            # get name of the group
            name = gmsh.model.getPhysicalName(group_dim, group_id)

            # gets elements per group
            entities_group = gmsh.model.getEntitiesForPhysicalGroup(group_dim, group_id)

            # get indices within the entities array of the entities which belong to the group
            indices = [original_entities.index((group_dim, entity_id)) for entity_id in entities_group]

            # get new entities which belong to the group
            new_entities_group = [split_entities[index] for index in indices]

            # get all new geometry ids belonging to the group
            new_geom_ids = [dimtag[1] for new_entity in new_entities_group for dimtag in new_entity]

            # remove existing physical group
            gmsh.model.removePhysicalGroups([(group_dim, group_id)])

            # re-add new physical group
            gmsh.model.addPhysicalGroup(group_dim, new_geom_ids, tag=group_id, name=name)

        self.__synchronize_geometry()

    @staticmethod
    def __synchronize_geometry():
        """
        Synchronizes the geometry.
        """
        # synchronize the geometry for generating the mesh
        gmsh.model.occ.synchronize()

        # synchronize the geo geometry such that physical groups are added, important is that this is done after
        # synchronizing the occ geometry :-D
        gmsh.model.geo.synchronize()

    def synchronize_gmsh(self):
        """
        Synchronizes the gmsh geometry and takes care of the intersections.

        """
        if gmsh.isInitialized():
            # synchronize the geometry for generating the mesh
            self.__synchronize_geometry()

            # synchronize intersections of entities
            self.__synchronize_intersection()

    @staticmethod
    def reset_gmsh_instance():
        """
        Resets gmsh object. Finalizes gmsh if the gmsh object is initialized and initializes gmsh.

        """
        if gmsh.isInitialized():
            gmsh.finalize()
        gmsh.initialize()

    def clear_geo_data(self):
        """
        Clears the geometry data.

        """

        self.__geo_data = {"points": {},
                           "lines": {},
                           "surfaces": {},
                           "volumes": {},
                           "physical_groups": {}}

    def clear_mesh_data(self):
        """
        Clears the mesh data.

        """

        self.__mesh_data = {}

    @staticmethod
    def set_verbosity_level(verbosity_level: int):
        """
        Sets the verbosity level of gmsh.
        Level of information printed on the terminal and the message console (0: silent except for fatal errors,
        1: errors, 2: warnings, 3: direct, 4: information, 5: status, 99: debug) Default value: 5

        """

        if verbosity_level not in [0, 1, 2, 3, 4, 5, 99]:
            raise ValueError(f"Verbosity level must be 0, 1, 2, 3, 4, 5 or 99. Verbosity level is {verbosity_level}")

        gmsh.option.setNumber("General.Verbosity", verbosity_level)
