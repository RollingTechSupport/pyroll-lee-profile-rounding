"""
Validates the curved-groove (Kocks mill) three-roll model of Byon et al. (2017) against the
digitized surface profiles of Fig. 9 / Fig. 10 of that paper, using the same real 7-pass mill
schedule (Seah Chang Won Special Steel Corp.) the paper itself validates against: round stock
-> Koval (1st pass) -> Koval (2nd, 3rd passes) -> ... -> round (7th pass). Passes 1-3 use the
exact groove geometry of the existing ``test_solve_3rp_round_oval_oval`` test, which already
matches Byon's own Table/Fig. 5 values. The predicted out-profile width is pinned to each
digitized shape's own spread tip radius (see ``tests/visual.py``), so that only the
free-surface (bulge/rounding) model is under test, independent of spread prediction.

Byon's Figs. 8-10 plot one sixth-symmetric sector of the shape in the paper's own local axis
orientation, which does not necessarily align with this plugin's corner-angle convention;
the comparison plots are therefore drawn side by side rather than overlaid, for a qualitative
human visual review, with light quantitative sanity checks (validity, area ordering,
non-negative eccentricity/radius) supplementing them.

Produces comparison plots under ``tests/output/`` for human visual review.
"""
from pyroll.core import CircularOvalGroove, Profile, Roll, RoundGroove, ThreeRollPass

import pyroll.profile_bulging  # noqa: F401  (registers the bulging post-processors)

from data.byon2017 import (
    MILL_SCHEDULE,
    PASS_1_BOUNDARY_MM,
    PASS_2_BOUNDARY_MM,
    PASS_3_BOUNDARY_MM,
    PASS_7_BOUNDARY_MM,
)
from visual import max_radius_mm, pin_width, plot_comparison


def in_profile():
    return Profile.round(
        diameter=MILL_SCHEDULE["in_diameter"],
        temperature=1000 + 273.15,
        strain=0,
        material=[MILL_SCHEDULE["material"], "steel"],
        flow_stress=100e6,
        length=1,
        density=7.5e3,
        specific_heat_capacity=690,
        thermal_conductivity=23,
    )


def solve_pass_1():
    roll_pass = ThreeRollPass(
        label="Byon Pass 1 (round -> Koval)",
        roll=Roll(
            groove=CircularOvalGroove(usable_width=67.12e-3, r1=0.1e-3, r2=129e-3 / 2, pad_angle=30),
            nominal_radius=195e-3 / 2,
            rotational_frequency=130 / 60,
        ),
        inscribed_circle_diameter=59.9e-3,
    )
    target_width = 2 * max_radius_mm(PASS_1_BOUNDARY_MM) * 1e-3
    return roll_pass, pin_width(roll_pass, in_profile(), target_width)


def solve_pass_2():
    _, out_profile_1 = solve_pass_1()
    roll_pass = ThreeRollPass(
        label="Byon Pass 2 (Koval -> Koval)",
        orientation="AntiY",
        roll=Roll(
            groove=CircularOvalGroove(usable_width=62.37e-3, r1=0.1e-3, r2=129e-3 / 2, pad_angle=30),
            nominal_radius=197e-3 / 2,
            rotational_frequency=100 / 60,
        ),
        inscribed_circle_diameter=54.4e-3,
    )
    target_width = 2 * max_radius_mm(PASS_2_BOUNDARY_MM) * 1e-3
    return roll_pass, pin_width(roll_pass, out_profile_1, target_width)


def solve_pass_3():
    _, out_profile_2 = solve_pass_2()
    roll_pass = ThreeRollPass(
        label="Byon Pass 3 (Koval -> Koval)",
        roll=Roll(
            groove=CircularOvalGroove(usable_width=56.25e-3, r1=0.1e-3, r2=129e-3 / 2, pad_angle=30),
            nominal_radius=201e-3 / 2,
            rotational_frequency=100 / 60,
        ),
        inscribed_circle_diameter=47.7e-3,
    )
    target_width = 2 * max_radius_mm(PASS_3_BOUNDARY_MM) * 1e-3
    return roll_pass, pin_width(roll_pass, out_profile_2, target_width)


def test_byon_pass_1_round_to_koval():
    roll_pass, out_profile = solve_pass_1()

    assert out_profile.cross_section.is_valid
    assert out_profile.bulge_eccentricity >= 0
    assert out_profile.bulge_radius > 0
    assert out_profile.cross_section.area < in_profile().cross_section.area

    plot_comparison(
        "byon_pass_1_round_to_koval.png",
        "Byon et al. (2017) Fig. 9(a): round -> Koval",
        out_profile,
        PASS_1_BOUNDARY_MM,
        reference_label="Digitized Fig. 9(a) (1/6 sector)",
    )


def test_byon_pass_2_koval_to_koval():
    roll_pass, out_profile = solve_pass_2()

    assert out_profile.cross_section.is_valid
    assert out_profile.bulge_eccentricity == 0, "Koval -> Koval passes use a circle through the section's own center"
    assert out_profile.bulge_radius == out_profile.width / 2

    plot_comparison(
        "byon_pass_2_koval_to_koval.png",
        "Byon et al. (2017) Fig. 9(b): Koval -> Koval",
        out_profile,
        PASS_2_BOUNDARY_MM,
        reference_label="Digitized Fig. 9(b) (1/6 sector)",
    )


def test_byon_pass_3_koval_to_koval():
    roll_pass, out_profile = solve_pass_3()

    assert out_profile.cross_section.is_valid
    assert out_profile.bulge_eccentricity == 0
    assert out_profile.bulge_radius == out_profile.width / 2

    plot_comparison(
        "byon_pass_3_koval_to_koval.png",
        "Byon et al. (2017) Fig. 10(a): Koval -> Koval",
        out_profile,
        PASS_3_BOUNDARY_MM,
        reference_label="Digitized Fig. 10(a) (1/6 sector)",
    )


def test_byon_pass_7_koval_to_round():
    _, out_profile_3 = solve_pass_3()

    roll_pass = ThreeRollPass(
        label="Byon Pass 7 (Koval -> round)",
        roll=Roll(
            groove=RoundGroove(r1=2e-3, r2=20e-3, depth=8e-3, pad_angle=30),
            nominal_radius=100e-3,
            rotational_frequency=100 / 60,
        ),
        inscribed_circle_diameter=MILL_SCHEDULE["out_diameter"],
    )
    target_width = 2 * max_radius_mm(PASS_7_BOUNDARY_MM) * 1e-3
    out_profile = pin_width(roll_pass, out_profile_3, target_width)

    assert out_profile.cross_section.is_valid
    # Byon Eqs. (11)-(12): Rs = DS/2 + Pc, solved together with the geometric relation between
    # B1, Pc and Rs; check that the returned (Pc, Rs) pair actually satisfies both.
    ds = roll_pass.inscribed_circle_diameter
    b1 = out_profile.width / 2
    pc = out_profile.bulge_eccentricity
    rs = out_profile.bulge_radius
    assert abs(rs - (ds / 2 + pc)) < 1e-9
    assert abs((b1 * (3 ** 0.5 / 2)) ** 2 + (b1 * 0.5 + pc) ** 2 - rs ** 2) < 1e-9 * rs ** 2

    plot_comparison(
        "byon_pass_7_koval_to_round.png",
        "Byon et al. (2017) Fig. 10(e): Koval -> round",
        out_profile,
        PASS_7_BOUNDARY_MM,
        reference_label="Digitized Fig. 10(e) (1/6 sector)",
    )
