"""
End-to-end validation: unlike the other three-roll tests, which pin the out-profile width to
a value read off the source papers' own figures (to validate the free-surface/bulge model in
isolation from spread prediction, see ``tests/visual.py``), these tests let PyRolL predict the
width itself, using Wusatowski's generic spread equation
(``pyroll.wusatowski_spreading``) fed by Lendl's equivalent rectangle method
(``pyroll.lendl_equivalent_method``, which supplies the ``equivalent_height``/
``equivalent_width`` a generic spread formula needs, including for the 3-fold symmetric
geometry of a ``ThreeRollPass``, where PyRolL core does not define them on its own).

This is what a real user of the plugin does: solve a pass with a spread-prediction plugin
loaded and this plugin's bulging model applying on top, with no manual intervention. Since
Wusatowski's equation is a generic flat-rolling formula rather than Min's or Byon's own
three-roll-specific spread equation, the predicted width is not expected to reproduce the
source papers' values exactly (see ``tests/test_three_roll_min_flat.py`` and
``tests/test_three_roll_byon_koval.py`` for that, tighter, comparison) - the point here is to
exercise and visualize the full, realistic model chain.
"""
from pyroll.core import CircularOvalGroove, FlatGroove, Profile, Roll, ThreeRollPass

import pyroll.wusatowski_spreading  # noqa: F401  (generic spread prediction)
import pyroll.lendl_equivalent_method  # noqa: F401  (equivalent_height/width for the spread model)
import pyroll.profile_bulging  # noqa: F401  (registers the bulging post-processors)

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
        # Either the bulge model applied (bulge_radius > 0), or it correctly recognized the
        # predicted width/eccentricity combination as outside its valid range and fell back to
        # the plain pyroll-core cross-section instead of an unphysical shape (bulge_radius is
        # None then) - see BulgeModelNotApplicable and test_three_roll_byon_koval.py's
        # test_byon_pass_1_round_to_koval, which hits exactly this for stand 1's own pass type.
        assert profile.bulge_radius is None or profile.bulge_radius > 0, f"{label}: invalid bulge_radius"

    assert profile_2.cross_section.area < profile_1.cross_section.area

    plot_comparison(
        "byon_sequence_predicted_spread.png",
        "Byon-style 2-pass Kocks-mill sequence with predicted (Wusatowski + Lendl) spread",
        profile_2,
        [(x * 1000, y * 1000) for x, y in profile_1.cross_section.exterior.coords],
        reference_label="Oval I out-profile (for scale reference)",
    )
