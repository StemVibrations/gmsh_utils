import numpy as np
import numpy.typing as npt

class MathUtils:

    @staticmethod
    def is_point_in_polygon(point, polygon_vertices, eps=1e-8):
        """
        Check if a point is inside a concave polygon. With the ray tracing algorithm.
        """

        polygon_vertices_array = np.array(polygon_vertices, dtype=float)
        point_array = np.array(point)

        polygon_normal = MathUtils.calculate_normal_plane(polygon_vertices_array)

        # check if the point is on the plane of the polygon
        if not MathUtils.is_point_on_plane(point_array, polygon_vertices_array[0], polygon_normal):
            return False

        # reduce the polygon to 2D. If the polygon is in the x-y plane, keep the first two coordinates
        # otherwise, keep the first and third coordinates
        if np.isclose(polygon_vertices_array[:, 2], polygon_vertices_array[0, 2]).all():
            projected_polygon = polygon_vertices_array[:, :2]
            eta, xi = point_array[:2]
        else:
            projected_polygon = polygon_vertices_array[:, [0, 2]]
            eta, xi = point_array[[0, 2]]

        # Shift projected_polygon by one point to the left (circular shift)
        shifted_polygon = np.roll(projected_polygon, -1, axis=0)

        # Unpack coordinates
        p1eta, p1xi = projected_polygon[:, 0], projected_polygon[:, 1]
        p2eta, p2xi = shifted_polygon[:, 0], shifted_polygon[:, 1]

        # Check if the point is within the y-bounds of the polygon edges
        inside = ((p1xi < xi + eps) & (xi < p2xi + eps)) | ((p2xi < xi + eps) & (xi < p1xi + eps))

        # Compute the intersection points of the polygon edges with the horizontal line at y
        eta_intersections = (p2eta - p1eta) * (xi - p1xi) / (p2xi - p1xi) + p1eta

        # check if point is on the edge
        on_edge = np.isclose(eta, eta_intersections) & inside
        if on_edge.any():
            return True

        # Check if the point is to the left of the intersection points
        inside = inside & (eta < eta_intersections+eps)

        # point is inside if an odd number of edges are to the right of the point
        return np.any(inside) & np.sum(inside) % 2 == 1

    @staticmethod
    def is_point_on_plane(point, plane_point, plane_normal):
        """
        Check if a point is on a plane.
        """

        # check if the point is on the plane
        return np.isclose(np.dot(point - plane_point, plane_normal), 0.0)

    @staticmethod
    def calculate_normal_plane(plane_vertices: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """
        Calculate the normal of a plane defined by three vertices.
        """

        v1 = plane_vertices[1] - plane_vertices[0]

        # make sure that the vertices are not collinear
        for i in range(2, len(plane_vertices)):
            v2 = plane_vertices[i] - plane_vertices[0]
            normal = np.cross(v1, v2)

            # return the normal if it is not zero
            if not np.allclose(normal, 0):
                return normal / np.linalg.norm(normal)

        raise ValueError("All plane vertices are collinear.")
