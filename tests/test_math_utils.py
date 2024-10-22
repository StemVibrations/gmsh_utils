import numpy as np
import pytest

from gmsh_utils.math_utils import MathUtils

class TestMathUtils:
    """
    Test the MathUtils class.
    """

    def test_calculate_normal_2d_plane(self):
        """
        Test the calculation of the normal vector of a 2D plane.

        """

        plane_vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        normal = MathUtils.calculate_normal_polygon(plane_vertices)
        expected_normal = np.array([0, 0, 1])
        assert np.allclose(normal, expected_normal)

        # vertices with collinear points
        plane_vertices = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [0,1, 0]])
        normal = MathUtils.calculate_normal_polygon(plane_vertices)
        expected_normal = np.array([0, 0, 1])
        assert np.allclose(normal, expected_normal)

        # vertices with only collinear points
        plane_vertices = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]])
        with pytest.raises(ValueError, match="All polygon vertices are collinear."):
            MathUtils.calculate_normal_polygon(plane_vertices)



    def test_calculate_normal_3d_plane(self):
        """
        Test the calculation of the normal vector of a 3D plane.
        """

        plane_vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 1]])
        normal = MathUtils.calculate_normal_polygon(plane_vertices)
        expected_normal = np.array([0, -1/np.sqrt(2), 1/np.sqrt(2)])
        assert np.allclose(normal, expected_normal)

    def test_is_point_on_2d_plane(self):
        """
        Test if a point is on a 2D plane.
        """
        point = np.array([14, 99, 0])
        plane_point = np.array([1, 4, 0])
        plane_normal = np.array([0, 0, 1])
        assert MathUtils.is_point_on_plane(point, plane_point, plane_normal)

        point_not_on_plane = np.array([14, 99, -99])
        assert not MathUtils.is_point_on_plane(point_not_on_plane, plane_point, plane_normal)


    def test_is_point_on_3d_plane(self):
        """
        Test if a point is on a 3D plane.
        """
        point = np.array([2, -4, 2])
        plane_point = np.array([0, 0, 0])
        plane_normal = np.array([1, 1, 1]) / np.sqrt(3)
        assert MathUtils.is_point_on_plane(point, plane_point, plane_normal)

        point_not_on_plane = np.array([2, -4, 3])
        assert not MathUtils.is_point_on_plane(point_not_on_plane, plane_point, plane_normal)

    def test_is_point_in_convex_polygon_xy_plane(self):
        """
        Test if a point is inside a convex polygon in the x-y plane.
        """
        polygon = [(0, 0, 0), (0, 1, 0), (1, 2, 0), (2, 1, 0), (2, 0, 0)]

        # point inside the polygon
        point = (0.5, 0.5, 0)
        assert MathUtils.is_point_in_polygon(point, polygon)

        # point on left edge of the polygon
        point_left_edge = (0, 0.5, 0)
        assert MathUtils.is_point_in_polygon(point_left_edge, polygon)

        # point on right edge of the polygon
        point_right_edge = (2, 0.5, 0)
        assert MathUtils.is_point_in_polygon(point_right_edge, polygon)

        # point on corner of the polygon
        point_corner = (1, 2, 0)
        assert MathUtils.is_point_in_polygon(point_corner, polygon)

        point_outside = (1.6, 1.6, 0)
        assert not MathUtils.is_point_in_polygon(point_outside, polygon)

    def test_is_point_in_convex_polygon_xz_plane(self):
        """
        Test if a point is inside a convex polygon in the x-z plane.
        """
        polygon = [(0, 0, 0), (0, 0, 1), (1, 0, 2), (2, 0, 1), (2, 0, 0)]

        # point inside the polygon
        point = (0.5, 0, 0.5)
        assert MathUtils.is_point_in_polygon(point, polygon)

        # point outside the polygon
        point_outside = (0.5, 0.5, 0.0)
        assert not MathUtils.is_point_in_polygon(point_outside, polygon)

    def test_is_point_in_concave_polygon(self):
        """
        Test if a point is inside a concave polygon.
        """

        polygon = [(0, 0, 0), (0, 1, 0), (0.5, 0.5, 0),(1,1,0), (1, 0, 0)]

        # point inside the polygon
        point = (0.75, 0.5, 0)
        assert MathUtils.is_point_in_polygon(point, polygon)

        # point in concavity of the polygon
        point_concave = (0.5, 0.75, 0)
        assert not MathUtils.is_point_in_polygon(point_concave, polygon)

        # point on the left edge of the polygon
        point_left_edge = (0, 0.9, 0)
        assert MathUtils.is_point_in_polygon(point_left_edge, polygon)

        # point on the right edge of the polygon
        point_right_edge = (1, 0.2, 0)
        assert MathUtils.is_point_in_polygon(point_right_edge, polygon)

        # point on corner of the polygon
        point_corner = (1, 1, 0)
        assert MathUtils.is_point_in_polygon(point_corner, polygon)

    def test_is_point_in_polygon_inclined_plane(self):
        """
        Test if a point is inside a polygon in an inclined plane.
        """
        polygon = [(0, 0, 0), (0, 1, 1), (1, 1, 1), (1, 0, 0)]

        # point inside the polygon
        point = (0.5, 0.5, 0.5)
        assert MathUtils.is_point_in_polygon(point, polygon)

        # point on edge the polygon
        point_edge = (1, 0.5, 0.5)
        assert MathUtils.is_point_in_polygon(point_edge, polygon)

        # point outside the polygon
        point_outside = (1.4, 0.5, 0.5)
        assert not MathUtils.is_point_in_polygon(point_outside, polygon)

    def test_is_point_in_polygon_on_yz_plane(self):
        """
        Test if a point is on a polygon plane.
        """
        polygon = [(0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0)]

        # point inside the polygon
        point = (0, 0.5, 0.5)
        assert MathUtils.is_point_in_polygon(point, polygon)

        # point on yz plane outside of the the polygon
        point = (0.0, 2.0, 1)
        assert not MathUtils.is_point_in_polygon(point, polygon)

    def test_calculate_rotation_matrix_polygon_with_xy_plane(self):
        """
        Test the calculation of the rotation matrix from a xy plane polygon to the x-y plane.
        """
        polygon_vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        rotation_matrix = MathUtils.calculate_rotation_matrix_polygon(polygon_vertices)
        expected_matrix = np.eye(3)
        assert np.allclose(rotation_matrix, expected_matrix)

    def test_calculate_rotation_matrix_polygon_with_yz_plane(self):
        """
        Test the calculation of the rotation matrix from a yz plane polygon to the x-y plane.
        """
        polygon_vertices = np.array([[0, 0, 0], [0, 1, 0], [0, 0, 1]])
        rotation_matrix = MathUtils.calculate_rotation_matrix_polygon(polygon_vertices)
        expected_matrix = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
        assert np.allclose(rotation_matrix, expected_matrix)

    def test_calculate_rotation_matrix_polygon_with_inclined_plane(self):
        """
        Test the calculation of the rotation matrix from an inclined plane polygon to the x-y plane.
        """
        polygon_vertices = np.array([[0, 0, 0], [0, 1, 1], [1, 1, 1]])
        rotation_matrix = MathUtils.calculate_rotation_matrix_polygon(polygon_vertices)
        expected_matrix = np.array([[1, 0, 0], [0, -1/np.sqrt(2), -1/np.sqrt(2)], [0, 1/np.sqrt(2), -1/np.sqrt(2)]])
        assert np.allclose(rotation_matrix, expected_matrix)

    def test_calculate_rotation_matrix_polygon_with_collinear_vertices(self):
        """
        Test the calculation of the rotation matrix from a polygon with collinear vertices.
        """
        polygon_vertices = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
        with pytest.raises(ValueError, match="All polygon vertices are collinear."):
            MathUtils.calculate_rotation_matrix_polygon(polygon_vertices)
