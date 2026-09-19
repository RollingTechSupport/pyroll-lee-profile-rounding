"""Helpers for pinning a profile's width and rendering comparison plots for human review."""
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTPUT_DIR = Path(__file__).parent / "output"
OUTPUT_DIR.mkdir(exist_ok=True)


def pin_width(roll_pass, in_profile, width):
    """Solves ``roll_pass`` with its out-profile width fixed, bypassing spread prediction."""
    roll_pass.out_profile = roll_pass.OutProfile(roll_pass, in_profile)
    roll_pass.out_profile.width = width
    return roll_pass.solve(in_profile)


def max_radius_mm(points_mm):
    """Distance from the origin of the farther endpoint of a digitized boundary trace."""
    return max(np.hypot(*points_mm[0]), np.hypot(*points_mm[-1]))


def rotate_points(points_mm, degrees):
    """Rotates a list of (x, y) points by ``degrees`` about the origin."""
    theta = np.radians(degrees)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    return [(x * cos_t - y * sin_t, x * sin_t + y * cos_t) for x, y in points_mm]


def mirror_right_half(points_mm):
    """Reconstructs a full, closed boundary from a right-half trace by appending its mirror."""
    left_half = [(-x, y) for x, y in reversed(points_mm)]
    return list(points_mm) + left_half


def reconstruct_threefold_sector(sector_points_mm):
    """Reconstructs a full, closed boundary from one sixth-symmetric sector by mirroring and threefold rotation."""
    corner_to_valley = np.array(sector_points_mm, dtype=float)
    valley_to_corner_mirrored = corner_to_valley[::-1] * np.array([-1, 1])
    lobe = np.vstack([valley_to_corner_mirrored, corner_to_valley[1:]])

    full_boundary = []
    for angle_deg in (0, -120, -240):
        theta = np.radians(angle_deg)
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        rotation = np.array([[cos_t, -sin_t], [sin_t, cos_t]])
        full_boundary.extend((lobe @ rotation.T).tolist())

    return full_boundary


def plot_overlay(filename, title, pyroll_profile, reference_points_mm, reference_label="Digitized from paper"):
    """Saves a plot overlaying the digitized reference boundary directly on the PyRolL prediction."""
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


def plot_sequence(filename, title, labeled_profiles):
    """Saves a row of cross-section plots, one per ``(label, profile)`` pair, for visual review of a pass sequence."""
    fig, axes = plt.subplots(1, len(labeled_profiles), figsize=(3.2 * len(labeled_profiles), 3.6))

    for ax, (label, profile) in zip(axes, labeled_profiles):
        xs, ys = profile.cross_section.exterior.xy
        ax.plot(np.array(xs) * 1000, np.array(ys) * 1000, "-", color="tab:blue")
        ax.set_aspect("equal")
        ax.axhline(0, color="gray", lw=0.5)
        ax.axvline(0, color="gray", lw=0.5)
        ax.set_title(f"{label}\nwidth={profile.width * 1000:.1f} mm, area={profile.cross_section.area * 1e6:.0f} mm2", fontsize=8)
        ax.set_xlabel("z (mm)")

    axes[0].set_ylabel("y (mm)")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / filename, dpi=130)
    plt.close(fig)


def plot_comparison(filename, title, pyroll_profile, reference_points_mm, reference_label="Digitized from paper"):
    """Saves a side-by-side comparison plot: digitized reference boundary vs. PyRolL prediction."""
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
