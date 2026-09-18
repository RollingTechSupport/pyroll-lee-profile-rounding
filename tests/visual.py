"""Shared helpers for the three-roll validation tests: pinning a profile's width to a known
value (decoupling the bulge/profile-shape model under test from spread prediction, which is a
separate concern - see the plugin docs), and rendering side-by-side comparison plots of the
digitized paper data against PyRolL's prediction for human visual review.
"""
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)


def pin_width(roll_pass, in_profile, width):
    """Solves ``roll_pass`` with its out-profile width fixed to ``width``, bypassing spread
    prediction so that only the free-surface (bulge) model is under test."""
    roll_pass.out_profile = roll_pass.OutProfile(roll_pass, in_profile)
    roll_pass.out_profile.width = width
    return roll_pass.solve(in_profile)


def max_radius_mm(points_mm):
    """Distance from the origin of the farther of the two endpoints of a digitized boundary
    trace, in millimeters, i.e. the spread tip radius B1. The endpoints (rather than the
    global maximum over all points) are used because a traced arc runs from a roll-contact
    reference point to the free-surface spread tip (or vice versa), and its extremes - not
    necessarily an interior point - represent that geometry; interior points can locally
    exceed the endpoints due to digitization noise without being the true spread tip."""
    return max(np.hypot(*points_mm[0]), np.hypot(*points_mm[-1]))


def rotate_points(points_mm, degrees):
    """Rotates a list of (x, y) points by ``degrees`` about the origin - used to align a
    digitized reference shape's own rotational phase with this plugin's corner-angle
    convention (``CORNER_ANGLES`` in ``three_roll_pass.py``) for an overlay plot; the source
    papers are not drawn in a consistent orientation relative to it."""
    theta = np.radians(degrees)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    return [(x * cos_t - y * sin_t, x * sin_t + y * cos_t) for x, y in points_mm]


def mirror_right_half(points_mm):
    """Reconstructs a full, closed boundary from a right-half trace (as digitized for Min et
    al.'s Fig. 2/13, mirror-symmetric about the y axis / z=0), by appending its mirror image."""
    left_half = [(-x, y) for x, y in reversed(points_mm)]
    return list(points_mm) + left_half


def reconstruct_threefold_sector(sector_points_mm):
    """Reconstructs a full, closed boundary from one sixth-symmetric sector (as digitized for
    Byon et al.'s Figs. 8-10: one sector, from a corner at 90 degrees to the mirror line
    towards the next corner at 30 degrees, given corner-to-valley) by mirroring it about the
    corner's own axis (the y axis) to get one full 120 degree "lobe" through that corner, then
    rotating that lobe by 0/+-120 degrees to complete the other two corners, matching this
    plugin's own corner-angle convention (see ``CORNER_ANGLES`` in ``three_roll_pass.py``)."""
    corner_to_valley = np.array(sector_points_mm, dtype=float)
    valley_to_corner_mirrored = corner_to_valley[::-1] * np.array([-1, 1])
    lobe = np.vstack([valley_to_corner_mirrored, corner_to_valley[1:]])

    full = []
    for angle_deg in (0, -120, -240):
        theta = np.radians(angle_deg)
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        rotation = np.array([[cos_t, -sin_t], [sin_t, cos_t]])
        full.extend((lobe @ rotation.T).tolist())

    return full


def plot_overlay(filename, title, pyroll_profile, reference_points_mm, reference_label="Digitized from paper"):
    """Saves a plot overlaying PyRolL's predicted cross-section boundary directly on top of a
    (already reconstructed to a full, closed boundary - see ``mirror_right_half`` /
    ``reconstruct_threefold_sector``) digitized reference boundary from the source paper, for
    human visual review of how closely, and how smoothly, the prediction matches it."""
    fig, ax = plt.subplots(figsize=(5, 5))

    reference = np.array(reference_points_mm)
    ax.plot(reference[:, 0], reference[:, 1], "-", color="tab:orange", lw=3, alpha=0.6, label=reference_label)

    xs, ys = pyroll_profile.cross_section.exterior.xy
    ax.plot(np.array(xs) * 1000, np.array(ys) * 1000, "-", color="tab:blue", lw=1.2, label="PyRolL (this plugin)")

    ax.set_aspect("equal")
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.set_xlabel("z (mm)")
    ax.set_ylabel("y (mm)")
    ax.legend(loc="lower center", fontsize=8)
    ax.set_title(title, fontsize=10)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / filename, dpi=130)
    plt.close(fig)


def plot_comparison(filename, title, pyroll_profile, reference_points_mm, reference_label="Digitized from paper"):
    """Saves a side-by-side comparison plot: PyRolL's predicted cross-section boundary vs.
    a digitized reference boundary trace from the source paper, for human visual review."""
    fig, (ax_ref, ax_pyroll) = plt.subplots(1, 2, figsize=(9, 4.5))

    reference = np.array(reference_points_mm)
    ax_ref.plot(reference[:, 0], reference[:, 1], "o-", color="tab:orange")
    ax_ref.set_aspect("equal")
    ax_ref.axhline(0, color="gray", lw=0.5)
    ax_ref.axvline(0, color="gray", lw=0.5)
    ax_ref.set_title(reference_label, fontsize=9)
    ax_ref.set_xlabel("z (mm)")
    ax_ref.set_ylabel("y (mm)")

    xs, ys = pyroll_profile.cross_section.exterior.xy
    ax_pyroll.plot(np.array(xs) * 1000, np.array(ys) * 1000, "o-", ms=2, color="tab:blue")
    ax_pyroll.set_aspect("equal")
    ax_pyroll.axhline(0, color="gray", lw=0.5)
    ax_pyroll.axvline(0, color="gray", lw=0.5)
    ax_pyroll.set_title("PyRolL (this plugin)", fontsize=9)
    ax_pyroll.set_xlabel("z (mm)")

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / filename, dpi=130)
    plt.close(fig)
