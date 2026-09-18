"""
Digitized data from Min, Kwon, Lee, Woo and Im (2003), "Analytical model for prediction of
deformed shape in three-roll rolling process", J. Mater. Process. Technol. 140, 471-477.

Coordinates were read off the printed figures (300 dpi page renders) by eye against the
printed axis gridlines; see the extraction notes for per-figure uncertainty estimates
(typically +/-0.5-1 mm for the mm-scale shape plots, +/-0.01-0.03 for the dimensionless
spread/eccentricity plots). Boundary traces are given as one representative half/quarter of
the true 3-fold symmetric shape, in the paper's own local axis orientation (not aligned with
this plugin's own corner-angle convention), in millimeters.
"""

# Table 1: rolling conditions of the 4-pass experiment (Figs. 2 and 13).
PROCESS_PARAMS = dict(
    roll_radius=185e-3,
    roll_speed_rpm=248,
    temperature_celsius=1050,
    in_diameter=35.8e-3,
)

# Fig. 2(a) / Fig. 13(a): 1st pass, round -> curved-hexagonal. Boundary trace of the right half
# of the shape (mm), starting at the top of the flat cap and ending at the bottom tip.
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

# Fig. 2(b) / Fig. 13(b): 2nd pass, curved-hexagonal -> hexagonal. Boundary trace of the right
# half of the shape (mm): flat top, slightly bulged sides/bottom.
PASS_2_BOUNDARY_MM = [
    (0, 13.3),
    (9, 13.3),
    (14.5, 4.0),
    (14.5, -4.0),
    (9, -13.3),
    (0, -16.2),
]

# Fig. 9: least-squares fit constants of the spread formula
# Delta B / B0 = alpha * (Ldm/B0) * (FH/FO) + beta, by pass-type.
SPREAD_COEFFICIENTS = {
    "round_to_curved_hexagonal": dict(alpha=8.22, beta=1.03),
    "curved_hexagonal_to_hexagonal": dict(alpha=5.76, beta=1.00),
    "hexagonal_to_hexagonal": dict(alpha=7.17, beta=0.97),
}

# Fig. 12: eccentricity e/B0 vs. area fraction FH/FO (round -> curved-hexagonal pass), read off
# the plotted FE data points, with the fitted line e/B0 = 2.40 * FH/FO.
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
