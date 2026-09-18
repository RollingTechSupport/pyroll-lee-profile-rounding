"""
Free-surface (bulging) model for two-roll passes.

Implements the round-oval-round model of Lee and Choi (2000) / Lee and Goldhahn (2001) /
Lee (2002), and the square-diamond-square / oval-square construction attributed to Schmidt.
These models are not the focus of this plugin's three-roll restructuring; they are kept
functionally unchanged from the previous implementation, only reorganized into their own module.
"""
import math
import logging

import numpy as np
from shapely import Point, intersection, unary_union
from pyroll.core import Hook, Unit, Profile as BaseProfile, TwoRollPass
from pyroll.core.roll_pass.hookimpls.helpers import out_cross_section

TwoRollPass.OutProfile.bulge_radius = Hook[float]()


class TwoRollBulgingModel(Unit):
    """Computes the free-surface (bulge) shape of profiles rolled in a two-roll pass."""

    def __init__(self, roll_pass: TwoRollPass):
        self.roll_pass = roll_pass
        super().__init__(label=f"Two-Roll Bulging Model for {self.roll_pass}")

    def bulge_radius_round_oval_lee(self, profile: BaseProfile):
        weight = (self.roll_pass.roll.groove.usable_width - profile.width) / (
            self.roll_pass.roll.groove.usable_width - self.roll_pass.in_profile.width
        )
        usable_radius = (
            self.roll_pass.roll.groove.r2 * self.roll_pass.height
            - (self.roll_pass.roll.groove.usable_width ** 2 + self.roll_pass.height ** 2) / 4
        ) / (2 * self.roll_pass.roll.groove.r2 - self.roll_pass.roll.groove.usable_width)
        return self.roll_pass.in_profile.equivalent_radius * weight + usable_radius * (1 - weight)

    def bulge_radius_oval_round_lee(self, profile: BaseProfile):
        weight = (self.roll_pass.roll.groove.usable_width - profile.width) / (
            self.roll_pass.roll.groove.usable_width - self.roll_pass.in_profile.width
        )
        oval_radius = self.roll_pass.prev_of(TwoRollPass).roll.groove.r2

        if self.roll_pass.height == 2 * self.roll_pass.roll.groove.r2:
            usable_radius = 2 * self.roll_pass.roll.groove.r2
        else:
            usable_radius = self.roll_pass.roll.groove.r2 + (self.roll_pass.height - 2 * self.roll_pass.roll.groove.r2)

        return oval_radius * weight + usable_radius * (1 - weight)

    def bulge_radius_schmidt(self, profile: BaseProfile):
        return (
            self.roll_pass.height / 2 ** 2
            + profile.width / 2 ** 2
            - 2 * self.roll_pass.roll.groove.r2 * self.roll_pass.height / 2
        ) / (2 * (profile.width / 2 - self.roll_pass.roll.groove.r2))

    def bulge_radius(self, profile: BaseProfile):
        in_classifiers = self.roll_pass.in_profile.classifiers
        classifiers = self.roll_pass.classifiers

        if "round" in in_classifiers and "oval" in classifiers:
            return self.bulge_radius_round_oval_lee(profile)
        elif "oval" in in_classifiers and "round" in classifiers:
            return self.bulge_radius_oval_round_lee(profile)
        elif "square" in in_classifiers and "diamond" in classifiers:
            return self.bulge_radius_schmidt(profile)
        elif "diamond" in in_classifiers and "square" in classifiers:
            return self.bulge_radius_schmidt(profile)
        elif "oval" in in_classifiers and "square" in classifiers:
            return self.bulge_radius_schmidt(profile)
        elif "square" in in_classifiers and "oval" in classifiers:
            return self.bulge_radius_schmidt(profile)

    def cross_section_round_oval_round(self, profile: BaseProfile):
        circle_center = profile.width / 2 - profile.bulge_radius
        right_circle = Point(circle_center, 0).buffer(profile.bulge_radius)
        left_circle = Point(-circle_center, 0).buffer(profile.bulge_radius)
        max_cross_section = out_cross_section(self.roll_pass, math.inf)
        intersection_points = max_cross_section.boundary.intersection(right_circle.boundary)

        if intersection_points.is_empty:
            logging.getLogger(__name__).info("No intersection point found. Continuing without bulging.")
            return None

        elif (profile.bulge_radius * 2) > (abs(max_cross_section.bounds[0]) + max_cross_section.bounds[2]):
            circle_intersection = intersection(left_circle, right_circle)
            return intersection(circle_intersection, max_cross_section)

        else:
            intersection_points = list(intersection_points.geoms)
            first_intersection_point = min(intersection_points, key=lambda point: abs(point.y))
            cross_section_till_intersection = out_cross_section(
                self.roll_pass, abs(first_intersection_point.x) * 2
            )
            left_side_cross_section = intersection(max_cross_section, left_circle)
            right_side_cross_section = intersection(max_cross_section, right_circle)
            return unary_union([left_side_cross_section, cross_section_till_intersection, right_side_cross_section])

    def cross_section_square_diamond_square(self, profile: BaseProfile):
        separation_point_angle = np.arcsin(
            (profile.width / 2 - profile.bulge_radius) / (self.roll_pass.roll.groove.r2 - profile.bulge_radius)
        )
        separation_point_z_coordinate = self.roll_pass.roll.groove.r2 * np.sin(separation_point_angle)

        left_bulge = Point(-profile.width / 2 + profile.bulge_radius, 0).buffer(profile.bulge_radius)
        right_bulge = Point(profile.width / 2 - profile.bulge_radius, 0).buffer(profile.bulge_radius)
        intersection_cross_section = out_cross_section(self.roll_pass, 2 * np.abs(separation_point_z_coordinate))

        if (2 * profile.bulge_radius) < self.roll_pass.height:
            return unary_union([left_bulge, intersection_cross_section, right_bulge])
        else:
            helper_cs = out_cross_section(self.roll_pass, profile.width)
            left_bulge_with_intersection = left_bulge.intersection(helper_cs)
            right_bulge_with_intersection = right_bulge.intersection(helper_cs)
            return unary_union([left_bulge_with_intersection, right_bulge_with_intersection])

    def cross_section_square_oval_square(self, profile: BaseProfile):
        separation_point_angle = np.arcsin(
            (profile.width / 2 - profile.bulge_radius) / (self.roll_pass.roll.groove.r2 - profile.bulge_radius)
        )
        separation_point_z_coordinate = self.roll_pass.roll.groove.r2 * np.sin(separation_point_angle)

        left_bulge = Point(-profile.width / 2 + profile.bulge_radius, 0).buffer(profile.bulge_radius)
        right_bulge = Point(profile.width / 2 - profile.bulge_radius, 0).buffer(profile.bulge_radius)
        intersection_cross_section = out_cross_section(self.roll_pass, 2 * np.abs(separation_point_z_coordinate))

        if (2 * profile.bulge_radius) < self.roll_pass.height:
            return unary_union([left_bulge, intersection_cross_section, right_bulge])
        else:
            helper_cs = out_cross_section(self.roll_pass, profile.width)
            left_bulge_with_intersection = left_bulge.intersection(helper_cs)
            right_bulge_with_intersection = right_bulge.intersection(helper_cs)
            return left_bulge_with_intersection.intersection(right_bulge_with_intersection)

    def cross_section(self, profile: BaseProfile):
        in_classifiers = self.roll_pass.in_profile.classifiers
        classifiers = self.roll_pass.classifiers

        if "square" in in_classifiers and "diamond" in classifiers:
            return self.cross_section_square_diamond_square(profile)
        elif "diamond" in in_classifiers and "square" in classifiers:
            return self.cross_section_square_diamond_square(profile)
        elif "round" in in_classifiers and "oval" in classifiers:
            return self.cross_section_round_oval_round(profile)
        elif "oval" in in_classifiers and "round" in classifiers:
            return self.cross_section_round_oval_round(profile)
        elif "square" in in_classifiers and "oval" in classifiers:
            return self.cross_section_square_oval_square(profile)
        elif "oval" in in_classifiers and "square" in classifiers:
            return self.cross_section_square_oval_square(profile)
        else:
            return profile.cross_section

    def solve(self, in_profile: BaseProfile) -> BaseProfile:
        in_profile.bulge_radius = self.bulge_radius(profile=in_profile)
        in_profile.cross_section = self.cross_section(profile=in_profile)
        return in_profile


TwoRollPass.post_processors.append(TwoRollBulgingModel)
