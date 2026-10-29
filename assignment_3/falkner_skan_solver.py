#!/usr/bin/env python3
"""AE 683 - Assignment 3: Falkner-Skan boundary-layer equation."""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp


def falkner_skan_equations(eta, y, n):
    """Return the three first-order equations for f, f' and f''."""
    f = y[0]
    fp = y[1]
    fpp = y[2]

    # f''' = -((n + 1)/2) f f'' + n (f')^2 - n
    fppp = -0.5 * (n + 1.0) * f * fpp + n * (fp**2 - 1.0)
    return np.vstack((fp, fpp, fppp))


def boundary_conditions(ya, yb):
    """f(0) = 0, f'(0) = 0 and f'(eta_max) = 1."""
    return np.array([ya[0], ya[1], yb[1] - 1.0])


def solve_for_n(n, previous_solution, eta_max, tol=1e-7):
    """Solve for one value of n, using the previous result as a starting guess."""
    if previous_solution is None:
        eta = np.linspace(0.0, eta_max, 200)
        y_guess = np.zeros((3, eta.size))
        # A smooth first guess that approaches f' = 1 away from the wall.
        y_guess[0] = eta - 1.0 + np.exp(-eta)
        y_guess[1] = 1.0 - np.exp(-eta)
        y_guess[2] = np.exp(-eta)
    else:
        eta = previous_solution.x
        y_guess = previous_solution.y

    solution = solve_bvp(
        lambda e, y: falkner_skan_equations(e, y, n),
        boundary_conditions,
        eta,
        y_guess,
        tol=tol,
        max_nodes=20000,
    )

    if not solution.success:
        raise RuntimeError(f"Solver did not converge for n={n:g}: {solution.message}")
    return solution


def solve_all_values(n_values, eta_max):
    """Solve the requested values, using continuation between nearby values."""
    solutions = {}
    base = solve_for_n(0.0, None, eta_max)
    solutions[0.0] = base

    # Solve positive values in increasing order, starting from the Blasius case.
    current_n = 0.0
    current_solution = base
    for target_n in sorted(n for n in n_values if n > 0.0):
        steps = max(1, int(np.ceil((target_n - current_n) / 0.25)))
        for n in np.linspace(current_n, target_n, steps + 1)[1:]:
            current_solution = solve_for_n(float(n), current_solution, eta_max)
        solutions[target_n] = current_solution
        current_n = target_n

    # Continue to negative values in small steps; this helps near separation.
    current_n = 0.0
    current_solution = base
    for target_n in sorted((n for n in n_values if n < 0.0), reverse=True):
        steps = max(1, int(np.ceil(abs(target_n - current_n) / 0.002)))
        for n in np.linspace(current_n, target_n, steps + 1)[1:]:
            current_solution = solve_for_n(float(n), current_solution, eta_max, tol=1e-6)
        solutions[target_n] = current_solution
        current_n = target_n

    return solutions


def plot_fig3(solutions, output_path):
    """Plot the velocity profiles for the requested values of n."""
    fig, ax = plt.subplots(figsize=(7.5, 6.0), dpi=300)

    # Put the larger positive n curves first, like the assignment figure.
    n_values = [4.0, 1.0, 0.5, 0.0, -0.0654, -0.0904]
    labels = {
        4.0: "n = 4",
        1.0: "1",
        0.5: r"$\frac{1}{2}$",
        0.0: "0",
        -0.0654: "-0.0654",
        -0.0904: "-0.0904",
    }

    for n in n_values:
        sol = solutions[n]
        eta = sol.x
        fp = sol.y[1]
        x = np.sqrt((n + 1.0) / 2.0) * eta
        visible = x <= 4.8
        ax.plot(x[visible], fp[visible], color="black", linewidth=1.8)

    # Put each label near its curve so they do not overlap.
    label_positions = {
        4.0: ((0.80, 0.94), 0.70),
        1.0: ((1.15, 0.87), 0.88),
        0.5: ((1.27, 0.68), 1.12),
        0.0: ((1.50, 0.49), 1.35),
        -0.0654: ((2.10, 0.42), 2.08),
        -0.0904: ((2.72, 0.23), 2.62),
    }
    for n in n_values:
        (text_x, text_y), point_x = label_positions[n]
        sol = solutions[n]
        x_curve = np.sqrt((n + 1.0) / 2.0) * sol.x
        y_curve = sol.y[1]
        point_y = np.interp(point_x, x_curve, y_curve)
        ax.annotate(
            labels[n],
            xy=(point_x, point_y),
            xytext=(text_x, text_y),
            fontsize=10,
            fontweight="bold",
            arrowprops={"arrowstyle": "->", "color": "black", "lw": 0.9},
        )

    ax.set_xlim(0, 4.8)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_xlabel(r"$\sqrt{\frac{n+1}{2}}\,\eta$", fontsize=13)
    ax.set_ylabel(r"$u/U = f'(\eta)$", fontsize=13, rotation=0, labelpad=25)
    ax.set_title("Fig. 3: Falkner-Skan Velocity Profiles", fontsize=12)
    ax.tick_params(direction="in", top=True, right=True)
    ax.grid(True, linestyle=":", alpha=0.45)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Solve the Falkner-Skan equation.")
    parser.add_argument("--eta-max", type=float, default=10.0,
                        help="Upper limit for eta (default: 10)")
    parser.add_argument("--output-dir", "-o", default=".",
                        help="Folder for the output plot (default: current folder)")
    args = parser.parse_args()

    if args.eta_max <= 0:
        parser.error("--eta-max must be positive")

    os.makedirs(args.output_dir, exist_ok=True)
    n_values = [4.0, 1.0, 0.5, 0.0, -0.0654, -0.0904]
    solutions = solve_all_values(n_values, args.eta_max)

    print("AE 683 - Assignment 3: Falkner-Skan equation")
    print("{:>9} | {:>10} | {:>14}".format("n", "f''(0)", "Solver status"))
    print("-" * 42)
    for n in n_values:
        sol = solutions[n]
        print(f"{n:9.4f} | {sol.y[2, 0]:10.6f} | {'converged':>14}")

    output_path = os.path.join(args.output_dir, "fig3.png")
    plot_fig3(solutions, output_path)
    print("Assignment 3 complete.")


if __name__ == "__main__":
    main()
