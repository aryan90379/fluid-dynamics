# AE 683: Fluid Dynamics (Autumn 2026-1)
**Department of Aerospace Engineering, IIT Bombay**

**Student:** Aryan Meena  
**Roll Number:** 23B0074  
**Email:** 23b0074@iitb.ac.in  

---

## Repository Overview

This repository contains numerical solvers, theoretical formulation, and output visualizations for all four computational assignments in **AE 683 (Fluid Dynamics)**.

Each directory is completely self-contained with its dedicated Python solver, dependencies, and high-resolution output figures.

```
.
├── assignment_1/             # Sink + Vortex Potential Flow (Bathtub Vortex)
│   ├── streamline_plot.py    # Potential flow streamline & flow net solver
│   ├── bathtub_vortex_streamlines.png
│   ├── bathtub_vortex_flow_net.png
│   ├── bathtub_vortex_analysis.png
│   ├── README.md
│   └── requirements.txt
├── assignment_2/             # Blasius Boundary Layer Equation
│   ├── blasius_solver.py     # RK4 + Secant shooting method solver
│   ├── fig1.png              # u/U vs eta profile
│   ├── fig2.png              # Normal velocity v*sqrt(x/(nu*U)) vs eta
│   ├── blasius_fig1_fig2_combined.png
│   ├── README.md
│   └── requirements.txt
├── assignment_3/             # Falkner-Skan Boundary Layer Family
│   ├── falkner_skan_solver.py # Multi-pressure-gradient similarity solver
│   ├── fig3.png              # Velocity profiles across accelerated & adverse n
│   ├── README.md
│   └── requirements.txt
└── assignment_4/             # Compressible Aerodynamics
    ├── compressible_flow_solver.py # Area-Mach & theta-beta-M relations
    ├── fig4.png              # Area ratio A/A* vs Mach number M
    ├── fig5.png              # Oblique shock theta-beta-M chart
    ├── README.md
    └── requirements.txt
```

---

## Assignment Summaries & Direct Code Links

### [Assignment 1: Streamline Plotting of Combined Potential Flow](assignment_1/)
- **Assigned Flow:** Sink + Vortex (*Bathtub Vortex*)
- **Code:** [`assignment_1/streamline_plot.py`](assignment_1/streamline_plot.py)
- **Key Outputs:** 2D Streamlines with velocity vector field, orthogonal flow net ($\psi$ and $\phi$ contours), and radial velocity/pressure coefficient profiles.

### [Assignment 2: Numerical Solution of Blasius Equation](assignment_2/)
- **Problem:** $f f'' + 2 f''' = 0$ with $f(0) = 0$, $f'(0) = 0$, $f'(\infty) = 1$.
- **Code:** [`assignment_2/blasius_solver.py`](assignment_2/blasius_solver.py)
- **Key Outputs:** Fig. 1 ($u/U$ vs $\eta$), Fig. 2 ($v \sqrt{x/(\nu U)}$ vs $\eta$), and boundary layer thicknesses ($\delta^* = 1.7208$, $\theta = 0.6641$, $\delta_{99} = 4.91$, $H = 2.591$).

### [Assignment 3: Numerical Solution of Falkner-Skan Equation](assignment_3/)
- **Problem:** $f''' + \frac{n+1}{2} f f'' - n (f')^2 + n = 0$ with $f(0) = 0$, $f'(0) = 0$, $f'(\infty) = 1$.
- **Code:** [`assignment_3/falkner_skan_solver.py`](assignment_3/falkner_skan_solver.py)
- **Key Outputs:** Fig. 3 velocity profiles across favorable ($n=4.0, 1.0, 0.5$), zero ($n=0.0$), and adverse ($n=-0.0654, -0.0904$) pressure gradients down to separation ($f''(0) \approx 0$).

### [Assignment 4: Compressible Aerodynamics Solutions](assignment_4/)
- **Problem:**
  - Part (a): Area - Mach Number relation $\frac{A}{A^*} = \frac{1}{M}\left[\frac{2}{\gamma+1} + \frac{\gamma-1}{\gamma+1}M^2\right]^{\frac{\gamma+1}{2(\gamma-1)}}$
  - Part (b): Oblique shock wave $\theta$-$\beta$-$M$ relations $\tan\theta = 2\cot\beta\left[\frac{M_1^2\sin^2\beta - 1}{M_1^2(\gamma + \cos 2\beta) + 2}\right]$
- **Code:** [`assignment_4/compressible_flow_solver.py`](assignment_4/compressible_flow_solver.py)
- **Key Outputs:** Fig. 4 (Area-Mach curve with subsonic/supersonic branches and sonic throat at $M=1$) and Fig. 5 ($\theta$-$\beta$-$M$ oblique shock chart with detachment limit and sonic boundary).

---

## Quick Execution

To run any assignment solver and reproduce all figures:

```bash
# Assignment 1
python3 assignment_1/streamline_plot.py

# Assignment 2
python3 assignment_2/blasius_solver.py

# Assignment 3
python3 assignment_3/falkner_skan_solver.py

# Assignment 4
python3 assignment_4/compressible_flow_solver.py
```
