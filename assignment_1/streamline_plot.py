#!/usr/bin/env python3
"""AE 683 (Fluid Dynamics) 2026-1 — Assignment 1

Combined potential flow: sink + irrotational vortex (bathtub vortex).

Generates three figures:
  1. High-resolution streamlines coloured by speed.
  2. Orthogonal streamline/equipotential flow net.
  3. Four-panel analysis: streamlines, velocity directions, radial speed
     profiles, and an idealised free-surface depression.

Model (m > 0, k > 0):
  u_r     = -m/r                 (sink, radially inward)
  u_theta =  k/r                 (counter-clockwise swirl)
  u       = -(m*x + k*y)/r^2
  v       =  (k*x - m*y)/r^2
  phi     = -m*ln(r/r_ref) + k*theta
  psi     = -k*ln(r/r_ref) - m*theta

The ideal point sink/vortex is singular at r=0. The code masks a small core
(r < r_min) in the 2-D velocity plots rather than drawing a false zero-speed
point there. The free-surface plot separately uses a finite core and clips the
ideal Bernoulli estimate for display.
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


STUDENT_NAME = "Aryan Meena"
ROLL_NUMBER = "23B0074"
DEFAULT_M = 1.0
DEFAULT_K = 1.5
DEFAULT_DOMAIN_LIMIT = 3.5
DEFAULT_GRID_SIZE = 320
DEFAULT_CORE_RADIUS = 0.08
SPEED_LIMIT = 6.0  # Common colour scale limit for the two speed-coloured plots
R_REF = 1.0        # Reference length in metres inside the logarithms


def compute_flow_field(X, Y, m=DEFAULT_M, k=DEFAULT_K, r_min=DEFAULT_CORE_RADIUS):
    """Return velocity components, speed, potential, streamfunction, r and theta.

    Values inside r_min are masked as NaN because the point-flow model is not
    physically/numerically defined at its singular core.
    """
    if m <= 0 or k <= 0:
        raise ValueError("For this sink + vortex model, m and k must both be positive.")
    if r_min <= 0:
        raise ValueError("r_min must be positive.")

    R = np.hypot(X, Y)
    TH = np.arctan2(Y, X)
    core_mask = R < r_min
    R_safe = np.maximum(R, r_min)
    R2_safe = R_safe**2

    U = -(m * X + k * Y) / R2_safe
    V = (k * X - m * Y) / R2_safe
    speed = np.hypot(U, V)

    # The logarithms use a reference length so their arguments are dimensionless.
    phi = -m * np.log(R_safe / R_REF) + k * TH
    psi = -k * np.log(R_safe / R_REF) - m * TH

    # Do not fabricate a finite velocity at the mathematical singularity.
    for field in (U, V, speed, phi, psi):
        field[core_mask] = np.nan

    return U, V, speed, phi, psi, R, TH


def _save_figure(fig, output_path):
    """Save a figure at print-friendly resolution and close it."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {output_path}")


def _add_domain_axis_format(ax, L, label_size=11):
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [m]", fontsize=label_size, fontweight="bold")
    ax.set_ylabel("y [m]", fontsize=label_size, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.4)


def _plot_analytical_streamlines(ax, m, k, L, r_min, count=8):
    """Overlay analytical logarithmic spirals as a light visual cross-check.

    For each seed angle theta0 at R_outer:
        theta(r) = theta0 + (k/m)*ln(R_outer/r)
    Points outside the plotted square are replaced by NaN, which breaks the
    curve cleanly rather than joining separated in-domain segments.
    """
    R_outer = L * np.sqrt(2.0)
    r = np.geomspace(r_min * 1.05, R_outer, 1200)
    for i, theta0 in enumerate(np.linspace(0, 2 * np.pi, count, endpoint=False)):
        theta = theta0 + (k / m) * np.log(R_outer / r)
        xs = r * np.cos(theta)
        ys = r * np.sin(theta)
        outside = (np.abs(xs) > L) | (np.abs(ys) > L)
        xs = np.where(outside, np.nan, xs)
        ys = np.where(outside, np.nan, ys)
        ax.plot(
            xs, ys,
            color="cyan", linestyle="--", linewidth=0.85, alpha=0.75,
            label="Analytical streamlines" if i == 0 else None,
            zorder=3,
        )


def plot_streamlines(
    m=DEFAULT_M, k=DEFAULT_K, L=DEFAULT_DOMAIN_LIMIT,
    N=DEFAULT_GRID_SIZE, r_min=DEFAULT_CORE_RADIUS,
    output_path="bathtub_vortex_streamlines.png",
):
    """High-resolution streamline plot with a consistent speed colour scale."""
    x = np.linspace(-L, L, N)
    y = np.linspace(-L, L, N)
    X, Y = np.meshgrid(x, y)
    U, V, speed, _, _, _, _ = compute_flow_field(X, Y, m=m, k=k, r_min=r_min)
    norm = Normalize(vmin=0.0, vmax=SPEED_LIMIT, clip=True)
    speed_colours = np.ma.masked_invalid(np.clip(speed, 0.0, SPEED_LIMIT))

    fig, ax = plt.subplots(figsize=(9.5, 8.5), constrained_layout=True)
    ax.contourf(
        X, Y, speed_colours, levels=np.linspace(0, SPEED_LIMIT, 50),
        cmap="viridis", norm=norm, alpha=0.32, extend="max",
    )
    stream = ax.streamplot(
        x, y, np.ma.masked_invalid(U), np.ma.masked_invalid(V),
        color=speed_colours, cmap="inferno", norm=norm,
        density=1.75, linewidth=1.15, arrowsize=1.1, arrowstyle="->",
        zorder=2,
    )
    _plot_analytical_streamlines(ax, m, k, L, r_min)
    ax.scatter(
        [0], [0], color="red", s=75, zorder=5,
        edgecolors="black", linewidths=1.0,
        label="Sink + vortex singularity",
    )
    _add_domain_axis_format(ax, L, label_size=12)
    ax.set_title(
        "AE 683 Assignment 1: Sink + Vortex (Bathtub Vortex)\n"
        r"$\psi=-k\ln(r/r_{ref})-m\theta$; "
        r"$r=C\exp[-(m/k)\theta]$"
        + f"\n{STUDENT_NAME} | Roll No: {ROLL_NUMBER}",
        fontsize=12, fontweight="bold", pad=12,
    )
    cbar = fig.colorbar(stream.lines, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label(r"Speed $|\mathbf{V}|$ [m/s]", fontsize=11, fontweight="bold")
    ax.legend(loc="upper right", framealpha=0.95, fontsize=8.5)
    _save_figure(fig, output_path)


def _add_inward_arrow(ax, xs, ys, valid):
    """Add one clearly visible arrow pointing toward smaller radius."""
    valid_indices = np.flatnonzero(valid)
    if len(valid_indices) < 12:
        return
    # Choose a contiguous portion of the curve in the plotted domain.
    mid = len(valid_indices) // 2
    idx_outer = valid_indices[mid + 5] if mid + 5 < len(valid_indices) else valid_indices[-1]
    idx_inner = valid_indices[mid - 5] if mid >= 5 else valid_indices[0]
    if idx_outer == idx_inner:
        return
    # r_vals increases with index, so the arrow points from outer to inner.
    ax.annotate(
        "", xy=(xs[idx_inner], ys[idx_inner]),
        xytext=(xs[idx_outer], ys[idx_outer]),
        arrowprops=dict(arrowstyle="->", color="#0055d4", lw=1.5),
        zorder=4,
    )


def plot_flow_net(
    m=DEFAULT_M, k=DEFAULT_K, L=DEFAULT_DOMAIN_LIMIT,
    r_min=DEFAULT_CORE_RADIUS,
    output_path="bathtub_vortex_flow_net.png",
):
    """Plot analytical streamlines and equipotential lines (an orthogonal net)."""
    fig, ax = plt.subplots(figsize=(9.5, 8.5), constrained_layout=True)
    R_outer = L * np.sqrt(2.0)
    r_vals = np.geomspace(r_min * 1.1, R_outer, 900)
    num_lines = 16

    # psi = constant: logarithmic spirals, moving inward as theta increases.
    for i, theta0 in enumerate(np.linspace(0, 2 * np.pi, num_lines, endpoint=False)):
        theta_stream = theta0 + (k / m) * np.log(R_outer / r_vals)
        xs = r_vals * np.cos(theta_stream)
        ys = r_vals * np.sin(theta_stream)
        outside = (np.abs(xs) > L) | (np.abs(ys) > L)
        valid = ~outside
        ax.plot(
            np.where(valid, xs, np.nan), np.where(valid, ys, np.nan),
            color="#0055d4", linewidth=1.35, alpha=0.88,
            label=r"Streamlines ($\psi=\mathrm{const}$)" if i == 0 else None,
        )
        if i % 2 == 0:
            _add_inward_arrow(ax, xs, ys, valid)

    # phi = constant: another family of logarithmic spirals, orthogonal to psi.
    for i, theta0 in enumerate(np.linspace(0, 2 * np.pi, num_lines, endpoint=False)):
        theta_phi = theta0 - (m / k) * np.log(R_outer / r_vals)
        xp = r_vals * np.cos(theta_phi)
        yp = r_vals * np.sin(theta_phi)
        outside = (np.abs(xp) > L) | (np.abs(yp) > L)
        ax.plot(
            np.where(outside, np.nan, xp), np.where(outside, np.nan, yp),
            color="#d9381e", linestyle="--", linewidth=1.2, alpha=0.85,
            label=r"Equipotential lines ($\phi=\mathrm{const}$)" if i == 0 else None,
        )

    ax.scatter(
        [0], [0], color="black", s=80, zorder=6,
        label="Sink + vortex singularity",
    )
    _add_domain_axis_format(ax, L, label_size=12)
    ax.set_title(
        "Orthogonal Flow Net of a Bathtub Vortex\n"
        r"Streamlines ($\psi=\mathrm{const}$) and equipotentials "
        r"($\phi=\mathrm{const}$)" +
        f"\n{STUDENT_NAME} | Roll No: {ROLL_NUMBER}",
        fontsize=12, fontweight="bold", pad=12,
    )
    ax.legend(loc="upper right", framealpha=0.95, fontsize=9)
    _save_figure(fig, output_path)


def plot_comprehensive_analysis(
    m=DEFAULT_M, k=DEFAULT_K, L=3.0, r_min=DEFAULT_CORE_RADIUS,
    output_path="bathtub_vortex_analysis.png",
):
    """Create a four-panel analysis figure."""
    fig = plt.figure(figsize=(15, 12.5), constrained_layout=True)
    norm = Normalize(vmin=0.0, vmax=SPEED_LIMIT, clip=True)

    # (a) Streamlines and speed distribution.
    ax1 = fig.add_subplot(2, 2, 1)
    x = np.linspace(-L, L, 260)
    y = np.linspace(-L, L, 260)
    X, Y = np.meshgrid(x, y)
    U, V, speed, _, _, _, _ = compute_flow_field(X, Y, m=m, k=k, r_min=r_min)
    speed_colours = np.ma.masked_invalid(np.clip(speed, 0.0, SPEED_LIMIT))
    ax1.contourf(
        X, Y, speed_colours, levels=np.linspace(0, SPEED_LIMIT, 40),
        cmap="viridis", norm=norm, alpha=0.32, extend="max",
    )
    stream1 = ax1.streamplot(
        x, y, np.ma.masked_invalid(U), np.ma.masked_invalid(V),
        color=speed_colours, cmap="inferno", norm=norm,
        density=1.45, linewidth=1.05, arrowsize=1.0,
    )
    ax1.scatter([0], [0], color="red", s=55, zorder=5)
    _add_domain_axis_format(ax1, L, label_size=9)
    ax1.set_title("(a) Streamlines and Speed Distribution", fontweight="bold")
    fig.colorbar(stream1.lines, ax=ax1, fraction=0.046, pad=0.03, label=r"$|\mathbf{V}|$ [m/s]")

    # (b) Velocity directions; arrows are unit-length direction indicators,
    # while their colour represents the local speed.
    ax2 = fig.add_subplot(2, 2, 2)
    xq = np.linspace(-L, L, 23)
    yq = np.linspace(-L, L, 23)
    Xq, Yq = np.meshgrid(xq, yq)
    Uq, Vq, speed_q, _, _, _, _ = compute_flow_field(Xq, Yq, m=m, k=k, r_min=r_min)
    valid_q = np.isfinite(speed_q) & (speed_q > 0)
    Uq_norm = np.divide(Uq, speed_q, out=np.zeros_like(Uq), where=valid_q)
    Vq_norm = np.divide(Vq, speed_q, out=np.zeros_like(Vq), where=valid_q)
    q = ax2.quiver(
        Xq, Yq, Uq_norm, Vq_norm, np.ma.masked_invalid(np.clip(speed_q, 0, SPEED_LIMIT)),
        cmap="turbo", norm=norm, pivot="mid", scale=25, width=0.0045,
        headwidth=3.8,
    )
    ax2.scatter([0], [0], color="red", s=55, zorder=5)
    _add_domain_axis_format(ax2, L, label_size=9)
    ax2.set_title("(b) Velocity Direction (colour = speed)", fontweight="bold")
    fig.colorbar(q, ax=ax2, fraction=0.046, pad=0.03, label=r"Speed $|\mathbf{V}|$ [m/s]")

    # (c) Analytical radial velocity profiles.
    ax3 = fig.add_subplot(2, 2, 3)
    r_profile = np.linspace(max(r_min * 1.1, 0.12), L, 300)
    u_r_mag = m / r_profile
    u_theta_mag = k / r_profile
    total_speed = np.sqrt(m**2 + k**2) / r_profile
    ax3.plot(r_profile, u_r_mag, label=r"$|u_r|=m/r$ (inward radial)", linewidth=2.2)
    ax3.plot(r_profile, u_theta_mag, label=r"$|u_\theta|=k/r$ (swirl)", linewidth=2.2)
    ax3.plot(r_profile, total_speed, "--", label=r"$|\mathbf{V}|=\sqrt{m^2+k^2}/r$", linewidth=2.2)
    ax3.set_xlim(0, L)
    ax3.set_ylim(bottom=0)
    ax3.set_xlabel("Radial distance r [m]", fontweight="bold")
    ax3.set_ylabel("Speed [m/s]", fontweight="bold")
    ax3.set_title("(c) Velocity Decay with Radial Distance", fontweight="bold")
    ax3.grid(True, linestyle=":", alpha=0.55)
    ax3.legend(framealpha=0.95, fontsize=9)

    # (d) Idealised free-surface depression from Bernoulli. A finite core and
    # display clipping prevent an unbounded, unreadable point at r=0.
    ax4 = fig.add_subplot(2, 2, 4, projection="3d")
    x3 = np.linspace(-L, L, 100)
    y3 = np.linspace(-L, L, 100)
    X3, Y3 = np.meshgrid(x3, y3)
    R3 = np.hypot(X3, Y3)
    surface_core_radius = max(0.3, r_min)
    R3_safe = np.maximum(R3, surface_core_radius)
    g = 9.81
    Z_surface = -(m**2 + k**2) / (2 * g * R3_safe**2)
    Z_surface = np.clip(Z_surface, -1.0, 0.0)
    surf = ax4.plot_surface(
        X3, Y3, Z_surface, cmap="YlGnBu_r", edgecolor="none",
        alpha=0.94, antialiased=True,
    )
    ax4.contour(X3, Y3, Z_surface, zdir="z", offset=-1.0, cmap="YlGnBu_r", linewidths=0.7)
    ax4.set_title("(d) Idealised Free-Surface Depression", fontweight="bold", pad=10)
    ax4.set_xlabel("x [m]", labelpad=2)
    ax4.set_ylabel("y [m]", labelpad=2)
    ax4.set_zlabel("Surface elevation Δz [m]", labelpad=3)
    ax4.set_zlim(-1.0, 0.05)
    ax4.view_init(elev=27, azim=-52)
    fig.colorbar(surf, ax=ax4, shrink=0.62, pad=0.08, label="Surface elevation [m]")

    fig.suptitle(
        "AE 683 (Fluid Dynamics) 2026-1 | Assignment 1: Combined Potential Flow\n"
        f"Sink + Vortex (Bathtub Vortex) | {STUDENT_NAME} ({ROLL_NUMBER})\n"
        f"Parameters: m = {m:g}, k = {k:g}; core cutoff r_min = {r_min:g} m",
        fontsize=13, fontweight="bold",
    )
    _save_figure(fig, output_path)


def main():
    parser = argparse.ArgumentParser(
        description="Generate Assignment 1 figures for a sink + vortex potential flow."
    )
    parser.add_argument("--sink-strength", "-m", type=float, default=DEFAULT_M,
                        help="Positive sink parameter m (default: 1.0).")
    parser.add_argument("--vortex-strength", "-k", type=float, default=DEFAULT_K,
                        help="Positive vortex parameter k (default: 1.5).")
    parser.add_argument("--domain-limit", "-L", type=float, default=DEFAULT_DOMAIN_LIMIT,
                        help="Half-width of the 2-D plotting domain in metres (default: 3.5).")
    parser.add_argument("--grid-size", "-N", type=int, default=DEFAULT_GRID_SIZE,
                        help="Grid points per direction for the high-resolution plot (default: 320).")
    parser.add_argument("--core-radius", type=float, default=DEFAULT_CORE_RADIUS,
                        help="Radius masked around the point singularity in 2-D plots (default: 0.08 m).")
    parser.add_argument("--output-dir", "-o", type=Path, default=Path("."),
                        help="Directory for generated PNGs (default: current directory).")
    args = parser.parse_args()

    if args.sink_strength <= 0 or args.vortex_strength <= 0:
        parser.error("--sink-strength and --vortex-strength must both be positive.")
    if args.domain_limit <= 0:
        parser.error("--domain-limit must be positive.")
    if args.grid_size < 30:
        parser.error("--grid-size must be at least 30.")
    if args.core_radius <= 0:
        parser.error("--core-radius must be positive.")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Generating figures for {STUDENT_NAME} ({ROLL_NUMBER})")
    print(f"Parameters: m = {args.sink_strength}, k = {args.vortex_strength}")
    print(f"2-D core cutoff: r_min = {args.core_radius} m")

    plot_streamlines(
        m=args.sink_strength, k=args.vortex_strength, L=args.domain_limit,
        N=args.grid_size, r_min=args.core_radius,
        output_path=args.output_dir / "bathtub_vortex_streamlines.png",
    )
    plot_flow_net(
        m=args.sink_strength, k=args.vortex_strength, L=args.domain_limit,
        r_min=args.core_radius,
        output_path=args.output_dir / "bathtub_vortex_flow_net.png",
    )
    plot_comprehensive_analysis(
        m=args.sink_strength, k=args.vortex_strength, L=min(args.domain_limit, 3.0),
        r_min=args.core_radius,
        output_path=args.output_dir / "bathtub_vortex_analysis.png",
    )
    print("All three figures generated successfully.")


if __name__ == "__main__":
    main()
