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
separate spread model (e.g. ``pyroll-wusatowski-spreading``), not of this plugin.
"""
import math

import numpy as np
from shapely import Polygon, unary_union
from pyroll.core import Hook, Unit, Profile as BaseProfile, ThreeRollPass
from pyroll.core.roll_pass.hookimpls.helpers import out_cross_section3

from . import geometry

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

CORE_WIDTH_FRACTION = 0.9
"""Fraction of the profile width up to which the plain (non-bulged) contour is trusted to
already coincide with the bulged shape, used to fill the area between the three corner arcs."""


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
            return ds / 2 + eccentricity

        return profile.width / 2 - eccentricity

    def _circular_cross_section(self, profile: BaseProfile) -> Polygon:
        rp = self.roll_pass
        max_cross_section = out_cross_section3(rp, math.inf)
        core = out_cross_section3(rp, profile.width * CORE_WIDTH_FRACTION)

        bulges = [
            max_cross_section.intersection(
                geometry.polar_point(profile.bulge_eccentricity, angle).buffer(profile.bulge_radius)
            )
            for angle in CORNER_ANGLES
        ]

        return unary_union([core, *bulges])

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

    def cross_section(self, profile: BaseProfile) -> Polygon:
        if profile.bulge_radius is None:
            return self._linear_chamfer_cross_section(profile)
        return self._circular_cross_section(profile)

    def solve(self, in_profile: BaseProfile) -> BaseProfile:
        eccentricity = self.eccentricity(in_profile)
        in_profile.bulge_eccentricity = eccentricity
        in_profile.bulge_radius = self.bulge_radius(in_profile, eccentricity)
        in_profile.cross_section = self.cross_section(in_profile)
        return in_profile


ThreeRollPass.post_processors.append(ThreeRollBulgingModel)
