"""Free-surface (bulging/rounding) model for three-roll passes, after Min et al. (2003) and Byon et al. (2017)."""
import logging
import math

import numpy as np
from shapely import Polygon
from pyroll.core import Hook, Unit, Profile as BaseProfile, ThreeRollPass
from pyroll.core.roll_pass.hookimpls.helpers import out_cross_section3

from . import geometry

log = logging.getLogger(__name__)


class BulgeModelNotApplicable(Exception):
    """Raised when the predicted width makes a corner's free-surface arc inconsistent with the roll's own contour."""


ThreeRollPass.OutProfile.bulge_radius = Hook[float]()
"""Radius of curvature of the free surface, or ``None`` for a straight-chord pass or a fallback."""

ThreeRollPass.OutProfile.bulge_eccentricity = Hook[float]()
"""Offset of the free surface's center of curvature from the profile's center, or ``None`` for a fallback."""

CORNER_ANGLES = np.deg2rad([90, -30, 210])
"""Angular positions of the three corners of a 3-fold symmetric profile."""

MIN_ECCENTRICITY_COEFFICIENT = 2.40
"""Min et al.'s ``e / B0 = 2.40 * FH / FO`` coefficient, fitted to FE results of a flat-roll three-roll mill."""

BYON_ECCENTRICITY_COEFFICIENT = 3.133
"""Byon et al.'s ``Pc / B0 = 3.133 * FH / F0`` coefficient, fitted to FE results of a Kocks mill."""

CIRCULAR_CROSS_SECTION_RESOLUTION = 360
CROSSING_SEARCH_RESOLUTION = 121


def nearest_corner_angle(angle: float) -> float:
    return min(CORNER_ANGLES, key=lambda corner_angle: abs((angle - corner_angle + np.pi) % (2 * np.pi) - np.pi))


class ThreeRollBulgingModel(Unit):
    def __init__(self, roll_pass: ThreeRollPass):
        self.roll_pass = roll_pass
        super().__init__(label=f"Three-Roll Bulging Model for {self.roll_pass}")

    @property
    def is_flat_roll(self) -> bool:
        return "flat" in self.roll_pass.classifiers

    @property
    def in_is_round(self) -> bool:
        return "round" in self.roll_pass.in_profile.classifiers

    @property
    def out_is_round(self) -> bool:
        return "round" in self.roll_pass.classifiers

    @property
    def displaced_area_fraction(self) -> float:
        """``FH / F0``: fraction of the incoming cross-section that must be displaced by rolling."""
        displaced_area = self.roll_pass.displaced_cross_section.area
        incoming_area = self.roll_pass.in_profile.cross_section.area
        return displaced_area / incoming_area

    def eccentricity(self, out_profile: BaseProfile) -> float:
        in_radius = self.roll_pass.in_profile.width / 2

        if self.in_is_round:
            coefficient = MIN_ECCENTRICITY_COEFFICIENT if self.is_flat_roll else BYON_ECCENTRICITY_COEFFICIENT
            return coefficient * self.displaced_area_fraction * in_radius

        if self.is_flat_roll:
            return None

        if self.out_is_round:
            out_radius = out_profile.width / 2
            inscribed_diameter = self.roll_pass.inscribed_circle_diameter
            # Byon Eqs. (11)-(12) solved for Pc: Rs = DS/2 + Pc and (B1 sin60)^2 + (B1 cos60 + Pc)^2 = Rs^2.
            return ((inscribed_diameter / 2) ** 2 - out_radius ** 2) / (out_radius - inscribed_diameter)

        return 0.0

    def bulge_radius(self, out_profile: BaseProfile, eccentricity: float) -> float:
        if eccentricity is None:
            return None

        if not self.is_flat_roll and self.out_is_round:
            inscribed_diameter = self.roll_pass.inscribed_circle_diameter
            radius = inscribed_diameter / 2 + eccentricity
        else:
            radius = out_profile.width / 2 - eccentricity

        if radius <= 0:
            raise BulgeModelNotApplicable(
                f"{self.roll_pass}: free-surface radius of curvature {radius * 1e3:.2f} mm is not positive."
            )

        return radius

    def corner_crossing_offset(self, eccentricity: float, bulge_radius: float, roll_contour_boundary, fallback_radius: float) -> float:
        """Angular offset from a corner where its bulge circle first grows back to the roll's own contour, or ``None`` if that never happens within the corner's own half-sector."""
        half_sector = np.pi / len(CORNER_ANGLES)
        corner_angle = CORNER_ANGLES[0]

        for offset in np.linspace(0, half_sector, CROSSING_SEARCH_RESOLUTION):
            angle = corner_angle + offset
            circle_radius = geometry.circle_radius_at_angle(eccentricity, corner_angle, bulge_radius, angle)
            contour_radius = geometry.boundary_radius_at_angle(roll_contour_boundary, angle, fallback_radius)
            if circle_radius >= contour_radius:
                return offset

        return None

    def circular_cross_section(self, out_profile: BaseProfile, eccentricity: float, bulge_radius: float) -> Polygon:
        roll_contour = out_cross_section3(self.roll_pass, math.inf)
        roll_contour_boundary = roll_contour.boundary
        fallback_radius = 10 * out_profile.width

        crossing_offset = self.corner_crossing_offset(eccentricity, bulge_radius, roll_contour_boundary, fallback_radius)
        if crossing_offset is None:
            raise BulgeModelNotApplicable(
                f"{self.roll_pass}: free-surface arc (eccentricity={eccentricity * 1e3:.2f} mm, "
                f"radius={bulge_radius * 1e3:.2f} mm) never reaches the roll's own contour within a corner's half-sector."
            )

        angles = np.linspace(-np.pi, np.pi, CIRCULAR_CROSS_SECTION_RESOLUTION, endpoint=False)
        points = []
        for angle in angles:
            corner_angle = nearest_corner_angle(angle)
            offset_from_corner = abs((angle - corner_angle + np.pi) % (2 * np.pi) - np.pi)
            if offset_from_corner <= crossing_offset:
                radius = geometry.circle_radius_at_angle(eccentricity, corner_angle, bulge_radius, angle)
            else:
                radius = geometry.boundary_radius_at_angle(roll_contour_boundary, angle, fallback_radius)
            points.append((radius * np.cos(angle), radius * np.sin(angle)))

        return Polygon(points)

    def linear_chamfer_cross_section(self, out_profile: BaseProfile) -> Polygon:
        """Min's straight-chord approximation: each corner is replaced by chords from W1 (where roll contact ends) to the spread tip B1."""
        roll_faces = list(self.roll_pass.contour_lines.geoms)
        spread_tip_radius = out_profile.width / 2

        corner_points = [geometry.polar_point(spread_tip_radius, angle) for angle in CORNER_ANGLES]
        corner_faces = [geometry.two_nearest_lines(roll_faces, point) for point in corner_points]

        vertices = []
        for corner_index in range(len(CORNER_ANGLES)):
            spread_tip = corner_points[corner_index]
            incoming_faces = set(corner_faces[corner_index - 1])
            outgoing_faces = set(corner_faces[(corner_index + 1) % len(CORNER_ANGLES)])

            face_in = next(i for i in corner_faces[corner_index] if i in incoming_faces)
            face_out = next(i for i in corner_faces[corner_index] if i in outgoing_faces)

            chord_start = geometry.point_on_line_at_radius(roll_faces[face_in], spread_tip_radius, spread_tip)
            chord_end = geometry.point_on_line_at_radius(roll_faces[face_out], spread_tip_radius, spread_tip)

            vertices.extend([chord_start, spread_tip, chord_end])

        return Polygon([vertex.coords[0] for vertex in vertices])

    def cross_section(self, out_profile: BaseProfile, eccentricity: float, bulge_radius: float) -> Polygon:
        if bulge_radius is None:
            return self.linear_chamfer_cross_section(out_profile)
        return self.circular_cross_section(out_profile, eccentricity, bulge_radius)

    def solve(self, in_profile: BaseProfile) -> BaseProfile:
        try:
            eccentricity = self.eccentricity(in_profile)
            bulge_radius = self.bulge_radius(in_profile, eccentricity)
            cross_section = self.cross_section(in_profile, eccentricity, bulge_radius)
        except BulgeModelNotApplicable as e:
            log.warning("%s Falling back to the un-bulged cross-section from pyroll-core.", e)
            in_profile.bulge_eccentricity = None
            in_profile.bulge_radius = None
            return in_profile

        in_profile.bulge_eccentricity = eccentricity
        in_profile.bulge_radius = bulge_radius
        in_profile.cross_section = cross_section
        return in_profile


ThreeRollPass.post_processors.append(ThreeRollBulgingModel)
