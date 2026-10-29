#!/usr/bin/env python3

"""
AE 683 - Fluid Dynamics
Assignment 4: Numerical Solutions for Compressible Aerodynamics

Part (a): Area-Mach number relation
Part (b): Oblique shock wave relation

Name: Aryan Meena
Roll No.: 23B0074
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar


# Part (a): Area-Mach number relation

def area_mach(M, gamma=1.4):
    """Calculate the area ratio A/A* for a given Mach number."""
    term = 2 / (gamma + 1) + ((gamma - 1) / (gamma + 1)) * M**2
    power = (gamma + 1) / (2 * (gamma - 1))

    return (1 / M) * term**power


def find_mach(area_ratio, branch="subsonic", gamma=1.4):
    """Find the Mach number for a given area ratio."""
    if area_ratio < 1:
        raise ValueError("Area ratio must be greater than or equal to 1.")

    if np.isclose(area_ratio, 1.0):
        return 1.0

    if branch == "subsonic":
        limits = [1e-6, 1.0]
    else:
        limits = [1.0, 15.0]

    result = root_scalar(
        lambda M: area_mach(M, gamma) - area_ratio,
        bracket=limits,
        method="brentq"
    )

    return result.root


def plot_area_mach(gamma=1.4, filename="fig4.png"):
    """Plot the area-Mach number relation."""
    fig, ax = plt.subplots(figsize=(6.5, 6), dpi=300)

    # Subsonic branch
    M1 = np.linspace(0.023, 1.0, 400)
    A1 = area_mach(M1, gamma)

    # Supersonic branch
    M2 = np.linspace(1.0, 5.0, 600)
    A2 = area_mach(M2, gamma)

    ax.plot(A1[A1 <= 30], M1[A1 <= 30], "r", linewidth=2.2)
    ax.plot(A2[A2 <= 30], M2[A2 <= 30], "r", linewidth=2.2)

    ax.set_xlim(0, 30)
    ax.set_ylim(0, 5)

    ax.set_xticks([0, 10, 20, 30])
    ax.set_yticks([0, 1, 2, 3, 4, 5])

    ax.set_xlabel(r"$\frac{A}{A^*}$", fontsize=15, labelpad=8)
    ax.set_ylabel(r"$M$", fontsize=15, rotation=0, labelpad=10)
    ax.yaxis.set_label_coords(-0.08, 0.5)

    ax.set_title("Fig. 4: Area-Mach Number Relation", fontsize=12)

    ax.tick_params(direction="in", top=True, right=True)
    ax.grid(True, linestyle=":", alpha=0.55)

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

    print(f"Saved {filename}")


# Part (b): Oblique shock relation

def oblique_shock(M1, gamma=1.4, n_points=800):
    """Calculate the flow deflection angle and shock angle."""
    if np.isinf(M1):
        beta = np.linspace(1e-4, np.pi / 2, n_points)
        tan_theta = np.sin(2 * beta) / (gamma + np.cos(2 * beta))

    else:
        mu = np.arcsin(1 / M1)
        beta = np.linspace(mu + 1e-6, np.pi / 2, n_points)

        numerator = M1**2 * np.sin(beta)**2 - 1
        denominator = M1**2 * (gamma + np.cos(2 * beta)) + 2

        tan_theta = 2 / np.tan(beta) * numerator / denominator

    theta = np.degrees(np.arctan(np.maximum(tan_theta, 0)))
    beta = np.degrees(beta)

    return theta, beta


def plot_oblique_shock(gamma=1.4, filename="fig5.png"):
    """Plot the theta-beta-M diagram."""
    fig, ax = plt.subplots(figsize=(8, 7), dpi=300)

    mach_numbers = [1.05, 1.1, 1.2, 1.4, 1.6, 2, 3, 4, 5, 10, np.inf]

    labels = {
        1.2: (4, 56, r"$M_1 = 1.2$"),
        1.4: (9, 55, r"$M_1 = 1.4$"),
        1.6: (14, 56, r"$M_1 = 1.6$"),
        2.0: (17.5, 33, r"$M_1 = 2.0$"),
        3.0: (17.5, 23.5, r"$M_1 = 3$"),
        4.0: (17.5, 18, r"$M_1 = 4$"),
        5.0: (17.5, 14, r"$M_1 = 5$"),
        10.0: (14, 8.5, r"$M_1 = 10$"),
        np.inf: (14, 3.5, r"$M_1 = \infty$")
    }

    for M1 in mach_numbers:
        theta, beta = oblique_shock(M1, gamma)

        ax.plot(theta, beta, color="black", linewidth=1.2)

        if M1 in labels:
            x, y, label = labels[M1]

            ax.text(
                x, y, label,
                fontsize=7.5,
                rotation=-20,
                ha="center",
                va="center",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85)
            )

    ax.set_xlim(0, 45.6)
    ax.set_ylim(0, 90)

    ax.set_xticks(np.arange(0, 46, 5))
    ax.set_yticks(np.arange(10, 91, 10))

    ax.set_xlabel("Deflection angle $\\theta$, deg", fontsize=11)
    ax.set_ylabel("Wave angle $\\beta$, deg", fontsize=11)

    ax.set_title(
        r"Fig. 5: Oblique Shock Wave Solution ($\theta-\beta-M$ Diagram)",
        fontsize=12
    )

    ax.tick_params(direction="in", top=True, right=True)
    ax.grid(True, linestyle=":", alpha=0.55)

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

    print(f"Saved {filename}")


def main():
    parser = argparse.ArgumentParser(
        description="AE 683 Assignment 4"
    )

    parser.add_argument(
        "--gamma",
        type=float,
        default=1.4,
        help="Ratio of specific heats"
    )

    parser.add_argument(
        "--output-dir",
        "-o",
        default=".",
        help="Directory for saving the plots"
    )

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    gamma = args.gamma

    print("\nAE 683 - Assignment 4")
    print("Student: Aryan Meena | Roll No.: 23B0074")
    print(f"Ratio of specific heats: {gamma}\n")

    # Check the area-Mach relation
    print("Part (a): Area-Mach Number Relation")
    print(f"{'A/A*':>8} | {'Subsonic M':>12} | {'Supersonic M':>14}")
    print("-" * 42)

    area_ratios = [1, 1.5, 2, 5, 10, 25]

    for area in area_ratios:
        M_sub = find_mach(area, "subsonic", gamma)
        M_sup = find_mach(area, "supersonic", gamma)

        print(f"{area:8.2f} | {M_sub:12.5f} | {M_sup:14.5f}")

    plot_area_mach(
        gamma,
        os.path.join(args.output_dir, "fig4.png")
    )

    # Check the oblique shock relation
    print("\nPart (b): Oblique Shock Relation")

    print(
        f"{'M1':>6} | {'Mach angle':>12} | "
        f"{'Max theta':>12} | {'Beta at max':>12}"
    )
    print("-" * 50)

    mach_numbers = [1.2, 1.4, 1.6, 2, 3, 4, 5, 10, np.inf]

    for M1 in mach_numbers:
        theta, beta = oblique_shock(M1, gamma, 1000)

        i = np.argmax(theta)

        if np.isinf(M1):
            mu = 0.0
            label = "inf"
        else:
            mu = np.degrees(np.arcsin(1 / M1))
            label = f"{M1:.2f}"

        print(
            f"{label:>6} | {mu:12.2f} | "
            f"{theta[i]:12.2f} | {beta[i]:12.2f}"
        )

    plot_oblique_shock(
        gamma,
        os.path.join(args.output_dir, "fig5.png")
    )

    print("\nAssignment completed. Both figures have been saved.")


if __name__ == "__main__":
    main()