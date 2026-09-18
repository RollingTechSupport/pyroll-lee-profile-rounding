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
