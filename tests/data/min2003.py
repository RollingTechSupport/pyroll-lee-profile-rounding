"""Digitized data from Min, Kwon, Lee, Woo and Im (2003), J. Mater. Process. Technol. 140, 471-477."""

# Table 1: rolling conditions of the 4-pass experiment (Figs. 2 and 13).
PROCESS_PARAMS = dict(
    roll_radius=185e-3,
    roll_speed_rpm=248,
    temperature_celsius=1050,
    in_diameter=35.8e-3,
)

# Fig. 2(a)/13(a), pass 1 (round -> curved-hexagonal): right-half boundary trace (mm).
PASS_1_BOUNDARY_MM = [
    (0, 14.5),
    (6, 14.5),
    (12, 14.3),
    (13.5, 13.5),
    (16, 10.0),
    (17.5, 5.0),
    (17.0, 0.0),
    (15.0, -6.0),
    (11.0, -12.0),
    (6.0, -16.5),
    (0.0, -19.0),
]

# Fig. 2(b)/13(b), pass 2 (curved-hexagonal -> hexagonal): right-half boundary trace (mm).
PASS_2_BOUNDARY_MM = [
    (0, 13.3),
    (9, 13.3),
    (14.5, 4.0),
    (14.5, -4.0),
    (9, -13.3),
    (0, -16.2),
]

# Fig. 9: fitted alpha/beta of Delta B / B0 = alpha * (Ldm/B0) * (FH/FO) + beta, by pass type.
SPREAD_COEFFICIENTS = {
    "round_to_curved_hexagonal": dict(alpha=8.22, beta=1.03),
    "curved_hexagonal_to_hexagonal": dict(alpha=5.76, beta=1.00),
    "hexagonal_to_hexagonal": dict(alpha=7.17, beta=0.97),
}

# Fig. 12: e/B0 vs. FH/FO FE data points (round -> curved-hexagonal), fitted line e/B0 = 2.40 * FH/FO.
ECCENTRICITY_FIG12 = [
    (0.095, 0.250),
    (0.115, 0.275),
    (0.130, 0.320),
    (0.140, 0.355),
    (0.145, 0.360),
    (0.160, 0.385),
    (0.163, 0.395),
    (0.165, 0.400),
    (0.168, 0.405),
    (0.185, 0.445),
    (0.190, 0.460),
    (0.205, 0.490),
    (0.215, 0.505),
    (0.230, 0.540),
]
ECCENTRICITY_COEFFICIENT = 2.40
