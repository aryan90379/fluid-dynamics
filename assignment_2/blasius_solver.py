#!/usr/bin/env python3
"""
AE 683 (Fluid Dynamics) 2026-1 - Assignment 2
Numerical solution of the Blasius boundary-layer equation

Student Name: Aryan Meena
Roll Number: 23B0074

Equation:
    f''' + 0.5 * f * f'' = 0

Boundary conditions:
    f(0) = 0, f'(0) = 0, f'(infinity) = 1

The program saves fig1.png, fig2.png and a combined figure.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar


def blasius_ode(eta, y):
    """Return the first-order form of the Blasius equation."""
    f, fp, fpp = y
    fppp = -0.5 * f * fpp
    return [fp, fpp, fppp]


def solve_blasius(eta_max=10.0, n_points=1000):
    """Solve the equation using shooting and the RK45 integrator."""

    def shooting_error(slope):
        sol = solve_ivp(
            blasius_ode,
            (0.0, eta_max),
            [0.0, 0.0, slope],
            rtol=1e-10,
            atol=1e-12,
        )
        if not sol.success:
            raise RuntimeError("The Blasius integration did not converge.")
        return sol.y[1, -1] - 1.0

    # Find f''(0) so that f'(eta_max) is close to 1.
    root = root_scalar(
        shooting_error,
        bracket=(0.1, 0.5),
        method="brentq",
        xtol=1e-12,
    )
    if not root.converged:
        raise RuntimeError("Could not find the initial slope f''(0).")

    eta = np.linspace(0.0, eta_max, n_points)
    sol = solve_ivp(
        blasius_ode,
        (0.0, eta_max),
        [0.0, 0.0, root.root],
        t_eval=eta,
        rtol=1e-10,
        atol=1e-12,
    )
    if not sol.success:
        raise RuntimeError("Could not calculate the final Blasius solution.")

    f, fp, fpp = sol.y

    # Dimensionless normal velocity from the similarity solution.
    v_norm = 0.5 * (eta * fp - f)

    # Boundary-layer thickness coefficients.
    delta_star = np.trapezoid(1.0 - fp, eta)
    theta = np.trapezoid(fp * (1.0 - fp), eta)
    delta_99 = np.interp(0.99, fp, eta)
    shape_factor = delta_star / theta

    results = {
        "fpp0": root.root,
        "delta_star": delta_star,
        "theta": theta,
        "delta_99": delta_99,
        "shape_factor": shape_factor,
        "v_inf": v_norm[-1],
    }
    return eta, f, fp, v_norm, results


def mark_thicknesses(ax, results):
    """Add the momentum and displacement thickness markers to Fig. 1."""
    for value, label in (
        (results["theta"], r"$\theta$"),
        (results["delta_star"], r"$\delta^*$"),
    ):
        ax.annotate(
            label,
            xy=(value, 0.0),
            xycoords="data",
            xytext=(value, -0.12),
            textcoords=("data", "axes fraction"),
            ha="center",
            va="top",
            fontsize=10,
            annotation_clip=False,
            arrowprops={"arrowstyle": "->", "color": "black", "lw": 1.0},
        )


def format_axes(ax, x_max, y_label, x_ticks):
    ax.set_xlim(0, x_max)
    ax.set_ylim(0, 1.05)
    ax.set_xticks(x_ticks)
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_xlabel(r"$\eta = y\sqrt{\frac{U}{\nu x}}$", fontsize=11, labelpad=8)
    ax.set_ylabel(y_label, fontsize=12, rotation=0, labelpad=15)
    ax.tick_params(direction="in", top=True, right=True)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)


def plot_fig1(eta, fp, results, output_path):
    fig, ax = plt.subplots(figsize=(6.5, 5.0), dpi=300)
    mask = eta <= 7.0
    ax.plot(eta[mask], fp[mask], color="black", linewidth=2.0)
    format_axes(ax, 7.0, r"$u/U = f'(\eta)$", np.arange(0, 8, 1))
    ax.set_title("Fig. 1: Blasius Velocity Profile", fontsize=11, pad=12)
    mark_thicknesses(ax, results)
    fig.subplots_adjust(bottom=0.23)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")


def plot_fig2(eta, v_norm, results, output_path):
    fig, ax = plt.subplots(figsize=(6.5, 5.0), dpi=300)
    mask = eta <= 6.0
    ax.plot(eta[mask], v_norm[mask], color="black", linewidth=2.0)
    ax.axhline(results["v_inf"], color="black", linestyle="--", linewidth=1.1)
    ax.text(0.15, results["v_inf"] + 0.025, f'{results["v_inf"]:.2f}', fontsize=10)
    format_axes(ax, 6.0, r"$v\sqrt{\frac{x}{\nu U}}$", [0, 3, 6])
    ax.set_yticks([0.0, 0.5, round(results["v_inf"], 2), 1.0])
    ax.set_title("Fig. 2: Dimensionless Normal Velocity", fontsize=11, pad=12)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")


def plot_combined(eta, fp, v_norm, results, output_path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2), dpi=300)

    mask1 = eta <= 7.0
    ax1.plot(eta[mask1], fp[mask1], color="black", linewidth=2.0)
    format_axes(ax1, 7.0, r"$u/U = f'(\eta)$", np.arange(0, 8, 1))
    ax1.set_title("Fig. 1: Longitudinal Velocity", fontsize=10)
    mark_thicknesses(ax1, results)
    ax1.grid(True, linestyle=":", alpha=0.35)

    mask2 = eta <= 6.0
    ax2.plot(eta[mask2], v_norm[mask2], color="black", linewidth=2.0)
    ax2.axhline(results["v_inf"], color="black", linestyle="--", linewidth=1.1)
    ax2.text(0.15, results["v_inf"] + 0.025, f'{results["v_inf"]:.2f}', fontsize=9)
    format_axes(ax2, 6.0, r"$v\sqrt{\frac{x}{\nu U}}$", [0, 3, 6])
    ax2.set_yticks([0.0, 0.5, round(results["v_inf"], 2), 1.0])
    ax2.set_title("Fig. 2: Transverse Velocity", fontsize=10)
    ax2.grid(True, linestyle=":", alpha=0.35)

    fig.suptitle(
        "AE 683 Assignment 2: Blasius Equation\n"
        r"$f''' + \frac{1}{2}ff'' = 0$  |  Aryan Meena (23B0074)",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Solve the Blasius boundary-layer equation.")
    parser.add_argument("--eta-max", type=float, default=10.0,
                        help="Maximum similarity coordinate (default: 10)")
    parser.add_argument("--n-points", type=int, default=1000,
                        help="Number of points in the solution (default: 1000)")
    parser.add_argument("--output-dir", "-o", default=".",
                        help="Folder for the output figures (default: current folder)")
    args = parser.parse_args()

    if args.eta_max < 7.0:
        parser.error("--eta-max must be at least 7 to cover the required plots.")
    if args.n_points < 100:
        parser.error("--n-points must be at least 100.")

    os.makedirs(args.output_dir, exist_ok=True)
    eta, f, fp, v_norm, results = solve_blasius(args.eta_max, args.n_points)

    print("AE 683 Assignment 2 - Blasius Boundary Layer")
    print("Student: Aryan Meena | Roll No: 23B0074")
    print(f"f''(0)                 = {results['fpp0']:.8f} (expected about 0.33205734)")
    print(f"Displacement thickness = {results['delta_star']:.4f} sqrt(nu*x/U)")
    print(f"Momentum thickness     = {results['theta']:.4f} sqrt(nu*x/U)")
    print(f"99% thickness          = {results['delta_99']:.4f} sqrt(nu*x/U)")
    print(f"Shape factor           = {results['shape_factor']:.4f}")
    print(f"Normal velocity limit  = {results['v_inf']:.4f} (expected about 0.8604)")

    fig1_path = os.path.join(args.output_dir, "fig1.png")
    fig2_path = os.path.join(args.output_dir, "fig2.png")
    combined_path = os.path.join(args.output_dir, "blasius_fig1_fig2_combined.png")

    plot_fig1(eta, fp, results, fig1_path)
    plot_fig2(eta, v_norm, results, fig2_path)
    plot_combined(eta, fp, v_norm, results, combined_path)
    print("Done. All figures were generated.")


if __name__ == "__main__":
    main()
