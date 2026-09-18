"""
Validates the flat-roll three-roll model of Min et al. (2003) against the digitized shapes of
Fig. 2 / Fig. 13 of that paper (4-pass round -> curved-hexagonal -> hexagonal -> hexagonal
sequence). The predicted out-profile width is pinned to the digitized shape's own spread tip
radius for each pass (see ``tests/visual.py``), so that only the free-surface (bulge/rounding)
model is under test, independent of spread prediction.

Produces both side-by-side and overlaid comparison plots under ``tests/output/`` for human
visual review; the overlay makes it easy to see that the small change in slope where the
bulge circle hands off to the plain contour (see ``three_roll_pass.py``) is not visible at the
scale of an actual rolled profile.
"""
from shapely import Polygon

from pyroll.core import FlatGroove, Profile, Roll, ThreeRollPass

import pyroll.profile_bulging  # noqa: F401  (registers the bulging post-processors)

from data.min2003 import PASS_1_BOUNDARY_MM, PASS_2_BOUNDARY_MM, PROCESS_PARAMS
from visual import max_radius_mm, mirror_right_half, pin_width, plot_comparison, plot_overlay, rotate_points


def digitized_polygon_area_mm2(right_half_boundary_mm):
    """Area of the full (mirrored) polygon described by a right-half boundary trace, via the
    shoelace formula."""
    return Polygon(mirror_right_half(right_half_boundary_mm)).area


def solve_pass_1():
    in_profile = Profile.round(
        diameter=PROCESS_PARAMS["in_diameter"],
        temperature=PROCESS_PARAMS["temperature_celsius"] + 273.15,
        strain=0,
        material=["C45", "steel"],
        flow_stress=100e6,
        length=1,
        density=7.5e3,
        specific_heat_capacity=690,
        thermal_conductivity=23,
    )

    roll_pass = ThreeRollPass(
        label="Min Pass 1 (round -> curved hexagonal)",
        roll=Roll(
            groove=FlatGroove(r1=5e-3, usable_width=45e-3, pad_angle=30),
            nominal_radius=PROCESS_PARAMS["roll_radius"],
            rotational_frequency=PROCESS_PARAMS["roll_speed_rpm"] / 60,
        ),
        inscribed_circle_diameter=29e-3,
    )

    target_width = 2 * max_radius_mm(PASS_1_BOUNDARY_MM) * 1e-3
    out_profile = pin_width(roll_pass, in_profile, target_width)

    assert out_profile.cross_section.is_valid
    assert "3fold" in roll_pass.classifiers

    reference_area = digitized_polygon_area_mm2(PASS_1_BOUNDARY_MM)
    predicted_area = out_profile.cross_section.area * 1e6
    relative_error = abs(predicted_area - reference_area) / reference_area
    assert relative_error < 0.20, (
        f"Predicted area {predicted_area:.1f} mm^2 deviates "
        f"{relative_error:.1%} from the digitized Fig. 2(a)/13(a) area {reference_area:.1f} mm^2"
    )

    plot_comparison(
        "min_pass_1_round_to_curved_hexagonal.png",
        "Min et al. (2003) Fig. 2(a)/13(a): round -> curved-hexagonal",
        out_profile,
        PASS_1_BOUNDARY_MM,
        reference_label="Digitized Fig. 2(a)/13(a) (right half)",
    )
    plot_overlay(
        "min_pass_1_round_to_curved_hexagonal_overlay.png",
        "Min et al. (2003) Fig. 2(a)/13(a): round -> curved-hexagonal (overlay)",
        out_profile,
        rotate_points(mirror_right_half(PASS_1_BOUNDARY_MM), 180),
        reference_label="Digitized Fig. 2(a)/13(a) (rotated 180 deg to this plugin's own corner convention)",
    )

    return out_profile


def test_min_pass_1_round_to_curved_hexagonal():
    solve_pass_1()


def test_min_pass_2_curved_hexagonal_to_hexagonal():
    out_profile_1 = solve_pass_1()

    roll_pass = ThreeRollPass(
        label="Min Pass 2 (curved hexagonal -> hexagonal)",
        orientation="AntiY",
        roll=Roll(
            groove=FlatGroove(r1=5e-3, usable_width=34e-3, pad_angle=30),
            nominal_radius=PROCESS_PARAMS["roll_radius"],
            rotational_frequency=PROCESS_PARAMS["roll_speed_rpm"] / 60,
        ),
        inscribed_circle_diameter=22e-3,
    )

    target_width = 2 * max_radius_mm(PASS_2_BOUNDARY_MM) * 1e-3
    out_profile = pin_width(roll_pass, out_profile_1, target_width)

    assert out_profile.cross_section.is_valid
    assert out_profile.bulge_radius is None, "the curved-hexagonal -> hexagonal pass uses Min's linear (straight-chord) approximation, not a circular bulge"

    reference_area = digitized_polygon_area_mm2(PASS_2_BOUNDARY_MM)
    predicted_area = out_profile.cross_section.area * 1e6
    relative_error = abs(predicted_area - reference_area) / reference_area
    assert relative_error < 0.20, (
        f"Predicted area {predicted_area:.1f} mm^2 deviates "
        f"{relative_error:.1%} from the digitized Fig. 2(b)/13(b) area {reference_area:.1f} mm^2"
    )

    plot_comparison(
        "min_pass_2_curved_hexagonal_to_hexagonal.png",
        "Min et al. (2003) Fig. 2(b)/13(b): curved-hexagonal -> hexagonal",
        out_profile,
        PASS_2_BOUNDARY_MM,
        reference_label="Digitized Fig. 2(b)/13(b) (right half)",
    )
    plot_overlay(
        "min_pass_2_curved_hexagonal_to_hexagonal_overlay.png",
        "Min et al. (2003) Fig. 2(b)/13(b): curved-hexagonal -> hexagonal (overlay)",
        out_profile,
        rotate_points(mirror_right_half(PASS_2_BOUNDARY_MM), 180),
        reference_label="Digitized Fig. 2(b)/13(b) (rotated 180 deg to this plugin's own corner convention)",
    )
