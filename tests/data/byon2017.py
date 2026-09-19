"""Digitized data from Byon, Kim, Kim and Lee (2017), J. Mech. Sci. Technol. 31(7), 3489-3497."""

# Fig. 5/Sec. 5.1: 7-pass Kocks mill, 71 mm round stock into 6 Koval passes then round at 35.4 mm.
MILL_SCHEDULE = dict(
    in_diameter=71e-3,
    out_diameter=35.4e-3,
    material="STS304",
)

# Table 1: fitted alpha/beta of B1 = alpha*B0*{(FH/F0)*(Ldm/B0) - beta*(1 - DS/(2*R0))} + B0.
SPREAD_COEFFICIENTS = {
    "round_to_koval": dict(alpha=7.118, beta=0.145),
    "koval_to_koval": dict(alpha=9.963, beta=0.145),
    "koval_to_round": dict(alpha=9.568, beta=0.145),
}

# Fig. 7: gamma calibration, Pc/B0 vs. FH/F0, with fitted line Pc/B0 = 3.133 * FH/F0.
ECCENTRICITY_FIG7 = [
    (0.08, 0.265),
    (0.12, 0.370),
    (0.15, 0.460),
    (0.19, 0.605),
    (0.23, 0.730),
]
ECCENTRICITY_COEFFICIENT = 3.133

# Fig. 9(a), pass 1 (round -> Koval, 14.7% reduction): one sixth-symmetric sector (mm).
PASS_1_BOUNDARY_MM = [
    (0, 30.0),
    (3.5, 29.6),
    (7, 29.2),
    (11, 28.5),
    (14.5, 27.5),
    (18, 26.3),
    (21.5, 24.8),
    (24.5, 22.5),
    (27, 20.5),
    (29, 18.7),
    (31, 18.0),
]

# Fig. 9(b), pass 2 (Koval -> Koval, 19.4% reduction).
PASS_2_BOUNDARY_MM = [
    (0, 31.0),
    (3, 30.7),
    (6.5, 30.0),
    (9.5, 28.6),
    (12, 26.5),
    (15, 24.3),
    (16, 22.0),
    (18, 19.3),
    (20, 17.0),
    (22, 14.5),
    (23.5, 13.8),
]

# Fig. 10(a), pass 3 (Koval -> Koval, 20.4% reduction).
PASS_3_BOUNDARY_MM = [
    (0, 24.0),
    (3, 23.9),
    (6, 23.6),
    (9, 23.1),
    (12, 22.3),
    (15, 21.2),
    (18, 19.8),
    (20.5, 17.8),
    (22.5, 15.8),
    (24.5, 14.2),
]

# Fig. 10(e), pass 7 (Koval -> round, 11.2% reduction, final pass).
PASS_7_BOUNDARY_MM = [
    (0, 17.6),
    (2.5, 17.5),
    (5, 17.2),
    (7.5, 16.7),
    (10, 15.8),
    (12.5, 14.4),
    (14, 12.5),
    (15, 10.5),
    (15.5, 8.7),
]
