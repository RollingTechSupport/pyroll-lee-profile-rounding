"""Shapely-based geometric helpers for the three-roll bulging model."""
import numpy as np
from shapely import Point, LineString
from shapely.ops import nearest_points


def polar_point(radius: float, angle: float) -> Point:
    return Point(radius * np.cos(angle), radius * np.sin(angle))


def circle_radius_at_angle(center_radius: float, center_angle: float, radius: float, angle: float) -> float:
    """Distance from the origin to the nearer crossing of the ray at ``angle`` with the circle of ``radius`` centered at ``polar_point(center_radius, center_angle)``, or ``inf`` if the ray misses it."""
    center_x = center_radius * np.cos(center_angle)
    center_y = center_radius * np.sin(center_angle)
    b = -2 * (center_x * np.cos(angle) + center_y * np.sin(angle))
    c = center_x ** 2 + center_y ** 2 - radius ** 2
    discriminant = b ** 2 - 4 * c

    if discriminant < 0:
        return np.inf

    root = np.sqrt(discriminant)
    crossing_distances = [t for t in ((-b + root) / 2, (-b - root) / 2) if t > 0]
    return min(crossing_distances) if crossing_distances else np.inf


def boundary_radius_at_angle(boundary, angle: float, fallback_radius: float) -> float:
    """Distance from the origin to the nearest crossing of ``boundary`` with the ray at ``angle``, or ``fallback_radius`` if it never crosses."""
    ray = LineString([(0, 0), (fallback_radius * np.cos(angle), fallback_radius * np.sin(angle))])
    crossings = ray.intersection(boundary)

    if crossings.is_empty:
        return fallback_radius

    crossing_points = list(crossings.geoms) if hasattr(crossings, "geoms") else [crossings]
    return min(np.hypot(p.x, p.y) for p in crossing_points if hasattr(p, "x"))


def nearest_point_in(geometry, near: Point) -> Point:
    return nearest_points(geometry, near)[0]


def point_on_line_at_radius(line: LineString, radius: float, near: Point) -> Point:
    """Point where ``line`` crosses the circle of ``radius`` around the origin, nearest ``near``; falls back to the closest point on ``line`` if it never reaches that radius."""
    circle = Point(0, 0).buffer(radius, quad_segs=256).boundary
    crossings = line.intersection(circle)

    if crossings.is_empty:
        return nearest_points(line, near)[0]

    return nearest_point_in(crossings, near)


def two_nearest_lines(lines, point: Point):
    """Indices of the two entries of ``lines`` closest to ``point``."""
    distances = [line.distance(point) for line in lines]
    return sorted(range(len(lines)), key=lambda i: distances[i])[:2]
