"""
Free-surface (bulging / rounding) model for three-roll passes.

Implements the two analytical models this plugin is built around:

- Min, Kwon, Lee, Woo and Im (2003), for three-roll mills with **flat** rolls, distinguishing
  the round -> curved-hexagonal pass (circular free-surface arc) from the curved-hexagonal ->
  hexagonal and hexagonal -> hexagonal passes (straight free-surface chord).
- Byon, Kim, Kim and Lee (2017), for Kocks mills with a **curved (grooved)** roll surface,
  distinguishing round stock -> Koval, Koval -> Koval and Koval -> round passes, all modeled
  by a circular free-surface arc.

Both models describe the free (out-of-contact) surface of the workpiece given the maximum
width ``B1`` of the profile after rolling; how ``B1`` itself is predicted is the concern of a
separate spread model (e.g. ``pyroll-wusatowski-spreading``), not of this plugin. Since that is
a separate, independently-calibrated model, the width it predicts is not guaranteed to be
consistent with the eccentricity this plugin computes from it (e.g. a spread formula not
calibrated for three-roll mills can under-predict the width for a given reduction); see
:py:exc:`BulgeModelNotApplicable`.
"""
import logging
import math

import numpy as np
from shapely import Polygon
from pyroll.core import Hook, Unit, Profile as BaseProfile, ThreeRollPass
from pyroll.core.roll_pass.hookimpls.helpers import out_cross_section3

from . import geometry

log = logging.getLogger(__name__)


class BulgeModelNotApplicable(Exception):
    """Raised when the eccentricity/radius of a corner's free-surface arc are inconsistent
    with the predicted width, so that the arc would have to be trusted well outside the
    (roughly one corner's own 60 degree half-sector) range the source papers fitted and
    validated it for - producing an unphysical, concave "necked" cross-section rather than a
    plausible rolled shape. See :py:meth:`ThreeRollBulgingModel._corner_crossing_offset`."""

ThreeRollPass.OutProfile.bulge_radius = Hook[float]()
"""Radius of curvature of the free (out-of-contact) surface, or ``None`` where the model
approximates the free surface by a straight chord instead of a circular arc
(Min's curved-hexagonal -> hexagonal and hexagonal -> hexagonal passes)."""

ThreeRollPass.OutProfile.bulge_eccentricity = Hook[float]()
"""Distance of the free-surface arc's center of curvature from the profile's center
(``Pc`` in Byon's notation, ``e`` in Min's)."""

CORNER_ANGLES = np.deg2rad([90, -30, 210])
"""Angular positions of the three corners of a 3-fold symmetric profile, matching the
orientation of :py:attr:`ThreeRollPass.contour_lines`."""

MIN_ECCENTRICITY_COEFFICIENT = 2.40
"""Coefficient of Min et al.'s eccentricity formula ``e / B0 = 2.40 * F_H / F_O``, fitted to
FE results of a three-roll mill with flat rolls."""

BYON_ECCENTRICITY_COEFFICIENT = 3.133
"""Coefficient of Byon et al.'s eccentricity formula ``Pc / B0 = 3.133 * F_H / F_0``, fitted to
FE results of a Kocks mill (curved-groove three-roll mill)."""

CIRCULAR_CROSS_SECTION_RESOLUTION = 360
"""Number of angular samples used to trace the bulged free surface in :py:meth:`_circular_cross_section`."""

CROSSING_SEARCH_RESOLUTION = 121
"""Number of angular samples used to search for the corner/natural-contour crossing in
:py:meth:`ThreeRollBulgingModel._corner_crossing_offset`."""


def _nearest_corner_angle(angle: float) -> float:
    return min(CORNER_ANGLES, key=lambda corner: abs((angle - corner + np.pi) % (2 * np.pi) - np.pi))


class ThreeRollBulgingModel(Unit):
    """Computes the free-surface (bulge/rounding) shape of profiles rolled in a
    :py:class:`ThreeRollPass`."""

    def __init__(self, roll_pass: ThreeRollPass):
        self.roll_pass = roll_pass
        super().__init__(label=f"Three-Roll Bulging Model for {self.roll_pass}")

    @property
    def is_flat_roll(self) -> bool:
        """True for Min's flat-roll three-roll mill, false for Byon's Kocks mill
        (curved/grooved roll surface)."""
        return "flat" in self.roll_pass.classifiers

    @property
    def in_is_round(self) -> bool:
        return "round" in self.roll_pass.in_profile.classifiers

    @property
    def out_is_round(self) -> bool:
        return "round" in self.roll_pass.classifiers

    @property
    def area_fraction(self) -> float:
        """``F_H / F_0``: fraction of the incoming cross-section that does not fit within the
        new pass's usable cross-section and must therefore spread."""
        f_h = self.roll_pass.displaced_cross_section.area
        f_0 = self.roll_pass.in_profile.cross_section.area
        return f_h / f_0

    def eccentricity(self, profile: BaseProfile) -> float:
        """The free-surface arc's eccentricity ``Pc`` (Byon) / ``e`` (Min).
        ``None`` for Min's straight-chord passes, where no circular arc is used."""
        b0 = self.roll_pass.in_profile.width / 2

        if self.in_is_round:
            coefficient = MIN_ECCENTRICITY_COEFFICIENT if self.is_flat_roll else BYON_ECCENTRICITY_COEFFICIENT
            return coefficient * self.area_fraction * b0

        if self.is_flat_roll:
            return None

        if self.out_is_round:
            # Byon Eqs. (11)-(12): solve Rs = DS/2 + Pc and
            # (B1 sin60)^2 + (B1 cos60 + Pc)^2 = Rs^2 simultaneously for Pc.
            b1 = profile.width / 2
            ds = self.roll_pass.inscribed_circle_diameter
            return ((ds / 2) ** 2 - b1 ** 2) / (b1 - ds)

        # Byon Koval -> Koval pass (Eq. 9): the surface profile is a circle through the
        # section's own center, i.e. zero eccentricity.
        return 0.0

    def bulge_radius(self, profile: BaseProfile, eccentricity: float) -> float:
        """The free-surface arc's radius of curvature ``Rs`` (Byon) / ``R1`` (Min).
        ``None`` for Min's straight-chord passes."""
        if eccentricity is None:
            return None

        if not self.is_flat_roll and self.out_is_round:
            ds = self.roll_pass.inscribed_circle_diameter
            radius = ds / 2 + eccentricity
        else:
            radius = profile.width / 2 - eccentricity

        if radius <= 0:
            raise BulgeModelNotApplicable(
                f"{self.roll_pass}: the free-surface arc's radius of curvature is not positive "
                f"({radius * 1e3:.2f} mm) for eccentricity {eccentricity * 1e3:.2f} mm; the "
                f"predicted width is inconsistent with this model."
            )

        return radius

    def _corner_crossing_offset(self, eccentricity: float, bulge_radius: float, boundary, far: float) -> float:
        """The angular distance (rad) from a corner at which its bulge circle first becomes
        as large as the roll's own natural (un-bulged) contour, searching outward up to one
        corner's own 60 degree half-sector. By construction (3-fold symmetry) this is the same
        for every corner and on both sides of each corner, so it only needs to be found once,
        for the canonical corner (:py:data:`CORNER_ANGLES`\\ [0]).

        Beyond this offset the free surface is the natural contour (full roll contact); within
        it, the corner's own circle. Without this cutoff, a corner's circle - fitted and valid
        only close to its own tip - can still be numerically smaller than the natural contour
        far away near the *next* corner if its eccentricity is large enough, incorrectly
        pulling the boundary inward there and producing an unphysical, concave "necked" shape
        instead of the plausible, convexity-preserving cross-section the two source papers
        actually depict. Returns ``None`` if no such crossing is found within the half-sector,
        meaning the circle stays smaller than the natural contour all the way to the *next*
        corner's own territory - i.e. the eccentricity/radius (and thus the predicted width
        they were computed from) are inconsistent with this model, see
        :py:exc:`BulgeModelNotApplicable`."""
        half_sector = np.pi / len(CORNER_ANGLES)
        corner = CORNER_ANGLES[0]

        for offset in np.linspace(0, half_sector, CROSSING_SEARCH_RESOLUTION):
            angle = corner + offset
            circle_r = geometry.circle_radius_at_angle(eccentricity, corner, bulge_radius, angle)
            natural_r = geometry.boundary_radius_at_angle(boundary, angle, far)
            if circle_r >= natural_r:
                return offset

        return None

    def _circular_cross_section(self, profile: BaseProfile, eccentricity: float, bulge_radius: float) -> Polygon:
        """Traces the bulged free surface directly, angle by angle: within
        :py:meth:`_corner_crossing_offset` of its nearest corner, the boundary follows that
        corner's bulge circle; beyond it, the roll's own natural (un-bulged) contour (full roll
        contact resumes). Built by direct sampling rather than polygon boolean ops
        (union/intersection of the three, generally overlapping, bulge circles) because the
        latter is prone to spurious reentrant cusps wherever two corners' circles, or a circle
        and the natural contour, cross without being tangent there."""
        rp = self.roll_pass
        max_cross_section = out_cross_section3(rp, math.inf)
        boundary = max_cross_section.boundary
        far = 10 * profile.width

        crossing_offset = self._corner_crossing_offset(eccentricity, bulge_radius, boundary, far)
        if crossing_offset is None:
            raise BulgeModelNotApplicable(
                f"{self.roll_pass}: the free-surface arc (eccentricity="
                f"{eccentricity * 1e3:.2f} mm, radius={bulge_radius * 1e3:.2f} mm) "
                f"never becomes as large as the roll's natural contour within a 60 degree "
                f"half-sector of its corner; the predicted width is inconsistent with this model."
            )

        angles = np.linspace(-np.pi, np.pi, CIRCULAR_CROSS_SECTION_RESOLUTION, endpoint=False)
        points = []
        for angle in angles:
            corner = _nearest_corner_angle(angle)
            offset = abs((angle - corner + np.pi) % (2 * np.pi) - np.pi)
            if offset <= crossing_offset:
                radius = geometry.circle_radius_at_angle(eccentricity, corner, bulge_radius, angle)
            else:
                radius = geometry.boundary_radius_at_angle(boundary, angle, far)
            points.append((radius * np.cos(angle), radius * np.sin(angle)))

        return Polygon(points)

    def _linear_chamfer_cross_section(self, profile: BaseProfile) -> Polygon:
        """Min's straight-chord free-surface approximation (curved-hexagonal -> hexagonal and
        hexagonal -> hexagonal passes): each sharp corner of the plain flat-roll contour is
        replaced by the two straight chords from the points ``W1`` (where roll contact ends, on
        each of the two roll faces flanking the corner) to the spread tip ``B1``, with ``W1``
        taken as the point on the respective roll face nearest to ``B1``. Between corners, the
        boundary follows the (straight) roll face exactly, so the whole free surface is built
        directly as a polygon through the ``W1``/``B1`` points, without needing the unbounded
        sharp-cornered contour at all."""
        rp = self.roll_pass
        b1_radius = profile.width / 2
        faces = list(rp.contour_lines.geoms)

        corner_points = [geometry.polar_point(b1_radius, angle) for angle in CORNER_ANGLES]
        corner_faces = [geometry.two_nearest_lines(faces, point) for point in corner_points]

        vertices = []
        for k in range(len(CORNER_ANGLES)):
            b1 = corner_points[k]
            incoming_faces = set(corner_faces[k - 1])
            outgoing_faces = set(corner_faces[(k + 1) % len(CORNER_ANGLES)])

            face_in = next(i for i in corner_faces[k] if i in incoming_faces)
            face_out = next(i for i in corner_faces[k] if i in outgoing_faces)

            w1_in = geometry.point_on_line_at_radius(faces[face_in], b1_radius, b1)
            w1_out = geometry.point_on_line_at_radius(faces[face_out], b1_radius, b1)

            vertices.extend([w1_in, b1, w1_out])

        return Polygon([v.coords[0] for v in vertices])

    def cross_section(self, profile: BaseProfile, eccentricity: float, bulge_radius: float) -> Polygon:
        if bulge_radius is None:
            return self._linear_chamfer_cross_section(profile)
        return self._circular_cross_section(profile, eccentricity, bulge_radius)

    def solve(self, in_profile: BaseProfile) -> BaseProfile:
        try:
            eccentricity = self.eccentricity(in_profile)
            bulge_radius = self.bulge_radius(in_profile, eccentricity)
            cross_section = self.cross_section(in_profile, eccentricity, bulge_radius)
        except BulgeModelNotApplicable as e:
            log.warning(
                "%s Falling back to the un-bulged cross-section from pyroll-core for this profile.", e
            )
            in_profile.bulge_eccentricity = None
            in_profile.bulge_radius = None
            return in_profile

        in_profile.bulge_eccentricity = eccentricity
        in_profile.bulge_radius = bulge_radius
        in_profile.cross_section = cross_section
        return in_profile


ThreeRollPass.post_processors.append(ThreeRollBulgingModel)
