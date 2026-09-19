"""End-to-end validation: lets PyRolL predict the out-profile width via Wusatowski + Lendl instead of pinning it."""
from pyroll.core import CircularOvalGroove, FlatGroove, Profile, Roll, ThreeRollPass

import pyroll.wusatowski_spreading
import pyroll.lendl_equivalent_method
import pyroll.profile_bulging

from data.min2003 import PROCESS_PARAMS
from data.byon2017 import MILL_SCHEDULE
from visual import plot_comparison


def test_min_sequence_with_predicted_spread():
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

    stand_1 = ThreeRollPass(
        label="Stand I",
        roll=Roll(
            groove=FlatGroove(r1=5e-3, usable_width=45e-3, pad_angle=30),
            nominal_radius=PROCESS_PARAMS["roll_radius"],
            rotational_frequency=PROCESS_PARAMS["roll_speed_rpm"] / 60,
        ),
        inscribed_circle_diameter=29e-3,
    )
    stand_2 = ThreeRollPass(
        label="Stand II",
        roll=Roll(
            groove=FlatGroove(r1=5e-3, usable_width=34e-3, pad_angle=30),
            nominal_radius=PROCESS_PARAMS["roll_radius"],
            rotational_frequency=PROCESS_PARAMS["roll_speed_rpm"] / 60,
        ),
        inscribed_circle_diameter=22e-3,
    )

    profile_1 = stand_1.solve(in_profile)
    profile_2 = stand_2.solve(profile_1)

    for profile, label in [(profile_1, "stand_1"), (profile_2, "stand_2")]:
        assert profile.cross_section.is_valid, f"{label}: invalid cross-section"
        assert profile.width > 0

    assert profile_2.cross_section.area < profile_1.cross_section.area, "material must lose area from pass to pass"

    plot_comparison(
        "min_sequence_predicted_spread.png",
        "Min-style 2-pass sequence with predicted (Wusatowski + Lendl) spread",
        profile_2,
        [(x * 1000, y * 1000) for x, y in profile_1.cross_section.exterior.coords],
        reference_label="Stand I out-profile (for scale reference)",
    )


def test_byon_sequence_with_predicted_spread():
    in_profile = Profile.round(
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

    stand_1 = ThreeRollPass(
        label="Oval I",
        roll=Roll(
            groove=CircularOvalGroove(usable_width=67.12e-3, r1=0.1e-3, r2=129e-3 / 2, pad_angle=30),
            nominal_radius=195e-3 / 2,
            rotational_frequency=130 / 60,
        ),
        inscribed_circle_diameter=59.9e-3,
    )
    stand_2 = ThreeRollPass(
        label="Oval II",
        orientation="AntiY",
        roll=Roll(
            groove=CircularOvalGroove(usable_width=62.37e-3, r1=0.1e-3, r2=129e-3 / 2, pad_angle=30),
            nominal_radius=197e-3 / 2,
            rotational_frequency=100 / 60,
        ),
        inscribed_circle_diameter=54.4e-3,
    )

    profile_1 = stand_1.solve(in_profile)
    profile_2 = stand_2.solve(profile_1)

    for profile, label in [(profile_1, "stand_1"), (profile_2, "stand_2")]:
        assert profile.cross_section.is_valid, f"{label}: invalid cross-section"
        # None means the model fell back (see BulgeModelNotApplicable); stand_1 hits this here.
        assert profile.bulge_radius is None or profile.bulge_radius > 0, f"{label}: invalid bulge_radius"

    assert profile_2.cross_section.area < profile_1.cross_section.area

    plot_comparison(
        "byon_sequence_predicted_spread.png",
        "Byon-style 2-pass Kocks-mill sequence with predicted (Wusatowski + Lendl) spread",
        profile_2,
        [(x * 1000, y * 1000) for x, y in profile_1.cross_section.exterior.coords],
        reference_label="Oval I out-profile (for scale reference)",
    )
