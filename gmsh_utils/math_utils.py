from typing import Sequence

import numpy as np
import numpy.typing as npt

from gmsh_utils.globals import PRECISION

class MathUtils:
    """
    Class with mathematical utilities.
    """

    @staticmethod
    def calculate_rotation_matrix_polygon(polygon_vertices: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """
        Calculate the rotation matrix to align the polygon with the x-y plane. Using the Rodrigues' rotation formula.

        Args:
            - polygon_vertices (npt.NDArray[np.float64]): Vertices of the polygon.

        Returns:
            - npt.NDArray[np.float64]: Rotation matrix.
        """

        polygon_normal = MathUtils.calculate_normal_polygon(polygon_vertices)

        # Define the target normal (Z-axis in this case)
        target_normal = np.array([0, 0, 1])

        # Calculate the rotation vector
        rotation_vector = np.cross(polygon_normal, target_normal)

        # Calculate sin and cos of the angle between the polygon normal and the target normal
        sin_theta = np.linalg.norm(rotation_vector)
        cos_theta = np.dot(polygon_normal, target_normal)

        # Normalize the rotation vector
        if not np.isclose(sin_theta, 0.0):
            rotation_vector = rotation_vector / sin_theta

        # Skew-symmetric cross-product matrix of rotation_vector
        v_cross = np.array([[0, -rotation_vector[2], rotation_vector[1]],
                            [rotation_vector[2], 0, -rotation_vector[0]],
                            [-rotation_vector[1], rotation_vector[0], 0]])

        # Rotation matrix using the Rodrigues' rotation formula
        rotation_matrix: npt.NDArray[np.float64] = (np.eye(3) + sin_theta * v_cross +
                                                    (1 - cos_theta) * np.dot(v_cross, v_cross))

        return rotation_matrix

    @staticmethod
    def is_point_in_polygon(point: Sequence[float], polygon_vertices: Sequence[Sequence[float]]) \
            -> bool:
        """
        Check if a point is inside a concave polygon. With the ray tracing algorithm.

        Args:
            - point (Sequence[float]): Point to check.
            - polygon_vertices (Sequence[Sequence[float]]): Vertices of the polygon.

        Returns:
            - bool: True if the point is inside the polygon, False otherwise.
        """

        # convert the input to numpy arrays
        polygon_vertices_array = np.array(polygon_vertices, dtype=float)
        point_array = np.array(point)

        # calculate the normal of the polygon
        polygon_normal = MathUtils.calculate_normal_polygon(polygon_vertices_array)

        # check if the point is on the plane of the polygon
        if not MathUtils.is_point_on_plane(point_array, polygon_vertices_array[0], polygon_normal):
            return False

        # reduce the polygon to 2D.
        # If the polygon is in the x-y plane, keep the first two coordinates
        if np.isclose(polygon_vertices_array[:, 2], polygon_vertices_array[0, 2]).all():
            projected_polygon = polygon_vertices_array[:, :2]
            eta, xi = point_array[:2]
        # if the polygon is in the x-z plane, keep the first and third coordinates
        elif np.isclose(polygon_vertices_array[:, 1], polygon_vertices_array[0, 1]).all():
            projected_polygon = polygon_vertices_array[:, [0, 2]]
            eta, xi = point_array[[0, 2]]
        # if the polygon is in the y-z plane, keep the second and third coordinates
        elif np.isclose(polygon_vertices_array[:, 0], polygon_vertices_array[0, 0]).all():
            projected_polygon = polygon_vertices_array[:, [1, 2]]
            eta, xi = point_array[[1, 2]]
        # if the polygon is not aligned with any of the planes, rotate the polygon to align it with the x-y plane
        else:
            # Calculate the rotation matrix to align the polygon with the x-y plane
            rotation_matrix = MathUtils.calculate_rotation_matrix_polygon(polygon_vertices_array)

            # Rotate the polygon vertices
            projected_polygon = np.dot(polygon_vertices_array, rotation_matrix.T)[:,:2]

            # Rotate the point
            point_array_rot = point_array.dot(rotation_matrix.T)

            # Unpack 2D coordinates
            eta, xi = point_array_rot[[0, 1]]

        # Shift projected_polygon by one point to the left (circular shift)
        shifted_polygon = np.roll(projected_polygon, -1, axis=0)

        # Unpack coordinates
        p1eta, p1xi = projected_polygon[:, 0], projected_polygon[:, 1]
        p2eta, p2xi = shifted_polygon[:, 0], shifted_polygon[:, 1]

        # Check if the point is within the xi-bounds of the polygon edges
        between_xi_bounds = (((p1xi < xi + PRECISION) & (xi < p2xi + PRECISION)) |
                             ((p2xi < xi + PRECISION) & (xi < p1xi + PRECISION)))

        # Compute the intersection points of the polygon edges with the horizontal line at y
        eta_intersections = (p2eta - p1eta) * (xi - p1xi) / (p2xi - p1xi) + p1eta

        # check if point is on the edge
        on_edge = np.isclose(eta, eta_intersections) & between_xi_bounds
        if on_edge.any():
            return True

        # Check if the point is to the left of the intersection points
        left_of_intersection = between_xi_bounds & (eta < eta_intersections+PRECISION)

        # point is inside if an odd number of edges are to the right of the point
        is_point_inside: bool = np.any(left_of_intersection) & np.sum(left_of_intersection) % 2 == 1
        return is_point_inside

    @staticmethod
    def is_point_on_plane(point: npt.NDArray[np.float64],
                          plane_point: npt.NDArray[np.float64],
                          plane_normal: npt.NDArray[np.float64]) -> bool:
        """
        Check if a point is on a plane.

        Args:
            - point (npt.NDArray[np.float64]): Point to check.
            - plane_point (npt.NDArray[np.float64]): Point on the plane.
            - plane_normal (npt.NDArray[np.float64]): Vector normal of the plane.

        Returns:
            - bool: True if the point is on the plane, False


        """

        # check if the point is on the plane
        is_point_on_plane: bool = np.isclose(np.dot(point - plane_point, plane_normal), 0.0)
        return is_point_on_plane

    @staticmethod
    def calculate_normal_polygon(polygon_vertices: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """
        Calculate the normal of a polygon defined by a minimum three vertices.

        Args:
            - polygon_vertices (npt.NDArray[np.float64]): Vertices of the polygon.

        Raises:
            - ValueError: If all the vertices are collinear.

        Returns:
            - npt.NDArray[np.float64]: Normal of the polygon.
        """

        v1 = polygon_vertices[1] - polygon_vertices[0]

        # make sure that the vertices are not collinear
        for i in range(2, len(polygon_vertices)):
            v2 = polygon_vertices[i] - polygon_vertices[0]
            normal: npt.NDArray[np.float64] = np.cross(v1, v2)

            # return the normal if it is not zero
            if not np.allclose(normal, 0):
                return normal / np.linalg.norm(normal)

        raise ValueError("All polygon vertices are collinear.")
