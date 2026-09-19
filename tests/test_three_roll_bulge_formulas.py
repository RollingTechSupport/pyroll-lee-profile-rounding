"""Numeric validation of the eccentricity formulas against digitized FE calibration data, independent of any PyRolL geometry."""
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
    """A spread tip radius exactly at the roll's own inscribed radius must give zero eccentricity."""
    inscribed_diameter = 0.040

    out_radius = inscribed_diameter / 2
    eccentricity = ((inscribed_diameter / 2) ** 2 - out_radius ** 2) / (out_radius - inscribed_diameter)
    bulge_radius = inscribed_diameter / 2 + eccentricity

    assert abs(eccentricity) < 1e-12
    assert abs(bulge_radius - inscribed_diameter / 2) < 1e-12


def test_koval_to_round_eccentricity_satisfies_geometric_relation():
    """The closed-form eccentricity must satisfy Byon's Eqs. (11)-(12) simultaneously."""
    inscribed_diameter = 0.040

    for out_radius in np.linspace(0.030, 0.039, 10):
        eccentricity = ((inscribed_diameter / 2) ** 2 - out_radius ** 2) / (out_radius - inscribed_diameter)
        bulge_radius = inscribed_diameter / 2 + eccentricity

        lhs = (out_radius * np.sin(np.deg2rad(60))) ** 2 + (out_radius * np.cos(np.deg2rad(60)) + eccentricity) ** 2
        assert abs(lhs - bulge_radius ** 2) < 1e-9 * bulge_radius ** 2
