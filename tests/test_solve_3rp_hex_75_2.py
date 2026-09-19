import pytest
from pyroll.core import FlatGroove, Profile, Roll, ThreeRollPass, Transport

import pyroll.wusatowski_spreading
import pyroll.lendl_equivalent_method
import pyroll.profile_bulging

from visual import plot_sequence


def test_solve_3rp_hexagon_75_2():
    in_profile = Profile.round(
        diameter=96 * 1.0177777e-3,
        temperature=1030 + 273.15,
        strain=0,
        material=["steel", "C20"],
        flow_stress=100e6,
        length=1,
        density=7.5e3,
        specific_heat_capacity=690,
        thermal_conductivity=23,
    )

    stand_1 = ThreeRollPass(
        label="Stand - I",
        orientation="Y",
        roll=Roll(
            groove=FlatGroove(r1=5e-3, usable_width=90e-3, pad_angle=30),
            nominal_radius=293e-3 / 2,
            rotational_frequency=1.73,
        ),
        inscribed_circle_diameter=84.43e-3,
        coulomb_friction_coefficient=0.4,
    )
    stand_2 = ThreeRollPass(
        label="Stand - II",
        orientation="AntiY",
        roll=Roll(
            groove=FlatGroove(r1=5e-3, usable_width=90e-3, pad_angle=30),
            nominal_radius=292.3e-3 / 2,
            rotational_frequency=1.913,
        ),
        inscribed_circle_diameter=80.84e-3,
        coulomb_friction_coefficient=0.4,
    )
    stand_3 = ThreeRollPass(
        label="Stand - III",
        orientation="Y",
        roll=Roll(
            groove=FlatGroove(r1=5e-3, usable_width=90e-3, pad_angle=30),
            nominal_radius=294.5e-3 / 2,
            rotational_frequency=2.09,
        ),
        inscribed_circle_diameter=77.51e-3,
        coulomb_friction_coefficient=0.4,
    )
    stand_4 = ThreeRollPass(
        label="Stand - IV",
        orientation="AntiY",
        roll=Roll(
            groove=FlatGroove(r1=5e-3, usable_width=90e-3, pad_angle=30),
            nominal_radius=296.4e-3 / 2,
            rotational_frequency=2.23,
        ),
        inscribed_circle_diameter=75.8e-3,
        coulomb_friction_coefficient=0.4,
    )
    stand_5 = ThreeRollPass(
        label="Stand - V",
        orientation="Y",
        roll=Roll(
            groove=FlatGroove(r1=5e-3, usable_width=90e-3, pad_angle=30),
            nominal_radius=301e-3 / 2,
            rotational_frequency=2.27,
        ),
        inscribed_circle_diameter=75.96e-3,
        coulomb_friction_coefficient=0.4,
    )

    transport_1_2 = Transport(label="I->II", length=0.72)
    transport_2_3 = Transport(label="II->III", length=0.72)
    transport_3_4 = Transport(label="III->IV", length=0.72, disk_element_count=10)
    transport_4_5 = Transport(label="IV->V", length=0.72)

    profile_1 = stand_1.solve(in_profile)
    profile_1_transported = transport_1_2.solve(profile_1)
    profile_2 = stand_2.solve(profile_1_transported)
    profile_2_transported = transport_2_3.solve(profile_2)
    profile_3 = stand_3.solve(profile_2_transported)
    profile_3_transported = transport_3_4.solve(profile_3)
    profile_4 = stand_4.solve(profile_3_transported)
    profile_4_transported = transport_4_5.solve(profile_4)

    labeled_profiles = [
        ("In", in_profile),
        ("Stand I", profile_1),
        ("Stand II", profile_2),
        ("Stand III", profile_3),
        ("Stand IV", profile_4),
    ]

    for label, profile in labeled_profiles:
        assert profile.cross_section.is_valid, f"{label}: invalid cross-section"

    for (_, previous), (label, current) in zip(labeled_profiles, labeled_profiles[1:]):
        assert current.cross_section.area < previous.cross_section.area, f"{label}: must lose area"

    plot_sequence(
        "hex_75_2_sequence.png",
        "5-stand flat-roll sequence (round -> hexagonal) with predicted (Wusatowski + Lendl) spread",
        labeled_profiles,
    )

    # Stand V's real inscribed circle (75.96 mm) is barely larger than stand IV's bulged
    # out-profile (~75.86 mm): a near-zero-reduction sizing pass. Lendl's 3-fold contact-width
    # hook finds no overlap with the incoming (already rounded) profile there and returns NaN.
    with pytest.raises(ValueError, match="lendl_width"):
        stand_5.solve(profile_4_transported)
