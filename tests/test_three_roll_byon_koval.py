"""Validates the Kocks-mill model of Byon et al. (2017) against digitized Fig. 9/10 shapes, using the paper's own 7-pass mill schedule."""
from pyroll.core import CircularOvalGroove, Profile, Roll, RoundGroove, ThreeRollPass

import pyroll.profile_bulging

from data.byon2017 import (
    MILL_SCHEDULE,
    PASS_1_BOUNDARY_MM,
    PASS_2_BOUNDARY_MM,
    PASS_3_BOUNDARY_MM,
    PASS_7_BOUNDARY_MM,
)
from visual import max_radius_mm, pin_width, plot_comparison, plot_overlay, reconstruct_threefold_sector


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
    """This large first-pass reduction pushes Byon's eccentricity formula outside its valid range, so the model must fall back."""
    roll_pass, out_profile = solve_pass_1()

    assert out_profile.cross_section.is_valid
    assert out_profile.bulge_eccentricity is None
    assert out_profile.bulge_radius is None
    assert out_profile.cross_section.area < in_profile().cross_section.area

    plot_comparison(
        "byon_pass_1_round_to_koval.png",
        "Byon et al. (2017) Fig. 9(a): round -> Koval (falls back to the core model, see docstring)",
        out_profile,
        PASS_1_BOUNDARY_MM,
        reference_label="Digitized Fig. 9(a) (1/6 sector)",
    )
    plot_overlay(
        "byon_pass_1_round_to_koval_overlay.png",
        "Byon et al. (2017) Fig. 9(a): round -> Koval (fallback shape vs. digitized bulge)",
        out_profile,
        reconstruct_threefold_sector(PASS_1_BOUNDARY_MM),
        reference_label="Digitized Fig. 9(a) (true bulge, for reference)",
    )


def test_byon_pass_2_koval_to_koval():
    roll_pass, out_profile = solve_pass_2()

    assert out_profile.cross_section.is_valid
    assert out_profile.bulge_eccentricity == 0
    assert out_profile.bulge_radius == out_profile.width / 2

    plot_comparison(
        "byon_pass_2_koval_to_koval.png",
        "Byon et al. (2017) Fig. 9(b): Koval -> Koval",
        out_profile,
        PASS_2_BOUNDARY_MM,
        reference_label="Digitized Fig. 9(b) (1/6 sector)",
    )
    plot_overlay(
        "byon_pass_2_koval_to_koval_overlay.png",
        "Byon et al. (2017) Fig. 9(b): Koval -> Koval (overlay)",
        out_profile,
        reconstruct_threefold_sector(PASS_2_BOUNDARY_MM),
        reference_label="Digitized Fig. 9(b)",
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
    plot_overlay(
        "byon_pass_3_koval_to_koval_overlay.png",
        "Byon et al. (2017) Fig. 10(a): Koval -> Koval (overlay)",
        out_profile,
        reconstruct_threefold_sector(PASS_3_BOUNDARY_MM),
        reference_label="Digitized Fig. 10(a)",
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
    # Byon Eqs. (11)-(12): Rs = DS/2 + Pc together with the geometric relation between B1, Pc, Rs.
    inscribed_diameter = roll_pass.inscribed_circle_diameter
    out_radius = out_profile.width / 2
    eccentricity = out_profile.bulge_eccentricity
    bulge_radius = out_profile.bulge_radius
    assert abs(bulge_radius - (inscribed_diameter / 2 + eccentricity)) < 1e-9
    assert abs((out_radius * (3 ** 0.5 / 2)) ** 2 + (out_radius * 0.5 + eccentricity) ** 2 - bulge_radius ** 2) < 1e-9 * bulge_radius ** 2

    plot_comparison(
        "byon_pass_7_koval_to_round.png",
        "Byon et al. (2017) Fig. 10(e): Koval -> round",
        out_profile,
        PASS_7_BOUNDARY_MM,
        reference_label="Digitized Fig. 10(e) (1/6 sector)",
    )
    plot_overlay(
        "byon_pass_7_koval_to_round_overlay.png",
        "Byon et al. (2017) Fig. 10(e): Koval -> round (overlay)",
        out_profile,
        reconstruct_threefold_sector(PASS_7_BOUNDARY_MM),
        reference_label="Digitized Fig. 10(e)",
    )
