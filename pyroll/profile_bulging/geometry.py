"""Shapely-based geometric helpers shared by the two- and three-roll bulging models."""
import numpy as np
from shapely import Point, Polygon, LineString
from shapely.ops import nearest_points


def polar_point(radius: float, angle: float) -> Point:
    """A point at the given ``radius`` and ``angle`` (rad) from the origin."""
    return Point(radius * np.cos(angle), radius * np.sin(angle))


def circle_radius_at_angle(center_radius: float, center_angle: float, radius: float, angle: float) -> float:
    """The distance from the origin to the point where the ray at ``angle`` crosses the circle
    of the given ``radius`` centered at ``polar_point(center_radius, center_angle)`` - the
    closer of the (up to two) crossings on the positive ray, or ``inf`` if the ray misses the
    circle entirely. Used to evaluate a bulge circle's contribution to the free surface at a
    given angle without constructing (and repeatedly re-discretizing) its polygon boundary."""
    cx = center_radius * np.cos(center_angle)
    cy = center_radius * np.sin(center_angle)
    b = -2 * (cx * np.cos(angle) + cy * np.sin(angle))
    c = cx ** 2 + cy ** 2 - radius ** 2
    discriminant = b ** 2 - 4 * c

    if discriminant < 0:
        return np.inf

    root = np.sqrt(discriminant)
    candidates = [t for t in ((-b + root) / 2, (-b - root) / 2) if t > 0]
    return min(candidates) if candidates else np.inf


def boundary_radius_at_angle(boundary, angle: float, far: float) -> float:
    """The distance from the origin to the nearest point where ``boundary`` (a LineString or
    MultiLineString, e.g. a polygon's exterior) crosses the ray at ``angle`` - or ``far`` if it
    never does."""
    ray = LineString([(0, 0), (far * np.cos(angle), far * np.sin(angle))])
    intersection = ray.intersection(boundary)

    if intersection.is_empty:
        return far

    points = list(intersection.geoms) if hasattr(intersection, "geoms") else [intersection]
    return min(np.hypot(p.x, p.y) for p in points if hasattr(p, "x"))


def nearest_point_in(geometry, near: Point) -> Point:
    """The point of ``geometry`` (of any type - Point, MultiPoint, LineString, a boundary
    intersection, ...) closest to ``near``."""
    return nearest_points(geometry, near)[0]


def point_on_line_at_radius(line: LineString, radius: float, near: Point) -> Point:
    """The point where ``line`` crosses the circle of the given ``radius`` around the origin,
    choosing the crossing closest to ``near`` where there is more than one. Falls back to the
    point on ``line`` nearest to ``near`` if the line never reaches that radius, e.g. because
    the roll face itself ends before it would cross the circle."""
    circle = Point(0, 0).buffer(radius, quad_segs=256).boundary
    intersection = line.intersection(circle)

    if intersection.is_empty:
        return nearest_points(line, near)[0]

    return nearest_point_in(intersection, near)


def two_nearest_lines(lines, point: Point):
    """Indices of the two entries of ``lines`` (a sequence of geometries) closest to ``point``."""
    distances = [line.distance(point) for line in lines]
    return sorted(range(len(lines)), key=lambda i: distances[i])[:2]
