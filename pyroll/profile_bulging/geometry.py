"""Shapely-based geometric helpers shared by the two- and three-roll bulging models."""
import numpy as np
from shapely import Point, LineString
from shapely.ops import nearest_points


def polar_point(radius: float, angle: float) -> Point:
    """A point at the given ``radius`` and ``angle`` (rad) from the origin."""
    return Point(radius * np.cos(angle), radius * np.sin(angle))


def point_on_line_at_radius(line: LineString, radius: float, near: Point) -> Point:
    """The point where ``line`` crosses the circle of the given ``radius`` around the origin,
    choosing the crossing closest to ``near`` where there is more than one. Falls back to the
    point on ``line`` nearest to ``near`` if the line never reaches that radius, e.g. because
    the roll face itself ends before it would cross the circle."""
    circle = Point(0, 0).buffer(radius, quad_segs=256).boundary
    intersection = line.intersection(circle)

    if intersection.is_empty:
        return nearest_points(line, near)[0]

    candidates = list(intersection.geoms) if hasattr(intersection, "geoms") else [intersection]
    return min(candidates, key=lambda candidate: candidate.distance(near))


def two_nearest_lines(lines, point: Point):
    """Indices of the two entries of ``lines`` (a sequence of geometries) closest to ``point``."""
    distances = [line.distance(point) for line in lines]
    return sorted(range(len(lines)), key=lambda i: distances[i])[:2]
