"""
Direct numeric validation of the eccentricity formulas against the FE calibration data
digitized from the source papers (Min et al. 2003 Fig. 12, Byon et al. 2017 Fig. 7), and a
closed-form self-consistency check of the Koval -> round eccentricity formula, independent of
any PyRolL geometry construction.
"""
import numpy as np

from data.min2003 import ECCENTRICITY_COEFFICIENT as MIN_COEFFICIENT
from data.min2003 import ECCENTRICITY_FIG12
from data.byon2017 import ECCENTRICITY_COEFFICIENT as BYON_COEFFICIENT
from data.byon2017 import ECCENTRICITY_FIG7


def test_min_eccentricity_coefficient_matches_fig12():
    area_fractions, relative_eccentricities = np.transpose(ECCENTRICITY_FIG12)
    predicted = MIN_COEFFICIENT * area_fractions
    relative_error = np.abs(predicted - relative_eccentricities) / relative_eccentricities
    assert np.all(relative_error < 0.15), (
        f"e/B0 = {MIN_COEFFICIENT} * FH/FO deviates by up to "
        f"{relative_error.max():.1%} from the digitized Fig. 12 data"
    )


def test_byon_eccentricity_coefficient_matches_fig7():
    area_fractions, relative_eccentricities = np.transpose(ECCENTRICITY_FIG7)
    predicted = BYON_COEFFICIENT * area_fractions
    relative_error = np.abs(predicted - relative_eccentricities) / relative_eccentricities
    assert np.all(relative_error < 0.15), (
        f"Pc/B0 = {BYON_COEFFICIENT} * FH/F0 deviates by up to "
        f"{relative_error.max():.1%} from the digitized Fig. 7 data"
    )


def test_koval_to_round_eccentricity_boundary_condition():
    """When the spread B1 exactly matches the roll's own inscribed radius DS/2 (no additional
    free bulging beyond the round contact surface), Byon's Eqs. (11)-(12) must reduce to
    Pc = 0 and Rs = DS/2 exactly (a perfectly round profile, fully on the roll contact)."""
    ds = 0.040

    b1 = ds / 2
    pc = ((ds / 2) ** 2 - b1 ** 2) / (b1 - ds)
    rs = ds / 2 + pc

    assert abs(pc) < 1e-12
    assert abs(rs - ds / 2) < 1e-12


def test_koval_to_round_eccentricity_satisfies_geometric_relation():
    """For a range of B1 < DS, the closed-form Pc solving Byon's Eqs. (11)-(12) must satisfy
    both equations simultaneously (see ThreeRollBulgingModel.eccentricity)."""
    ds = 0.040

    for b1 in np.linspace(0.030, 0.039, 10):
        pc = ((ds / 2) ** 2 - b1 ** 2) / (b1 - ds)
        rs = ds / 2 + pc

        lhs = (b1 * np.sin(np.deg2rad(60))) ** 2 + (b1 * np.cos(np.deg2rad(60)) + pc) ** 2
        assert abs(lhs - rs ** 2) < 1e-9 * rs ** 2
