# AE 683 (Fluid Dynamics) 2026-1
## Assignment 2: Numerical Solution of the Blasius Boundary Layer Equation

**Student Name:** Aryan Meena  
**Roll Number:** 23B0074  

---

## 1. Problem Formulation

The laminar boundary layer over a semi-infinite flat plate aligned with a uniform freestream $U$ in steady, incompressible, 2D flow without pressure gradient ($dp/dx = 0$) is governed by the Prandtl boundary layer equations:

$$\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0$$
$$u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} = \nu \frac{\partial^2 u}{\partial y^2}$$

Boundary conditions:
- At the wall ($y = 0$): $u(x, 0) = 0, \quad v(x, 0) = 0$ (no-slip and impermeability)
- In the freestream ($y \to \infty$): $u(x, \infty) = U$

### 1.1 Similarity Transformation
Introducing the similarity variable $\eta$ and streamfunction $\psi(x, y)$:
$$\eta = y \sqrt{\frac{U}{\nu x}}$$
$$\psi(x, y) = \sqrt{\nu x U} \, f(\eta)$$

Velocity components:
$$u = \frac{\partial \psi}{\partial y} = U f'(\eta)$$
$$v = -\frac{\partial \psi}{\partial x} = \frac{1}{2} \sqrt{\frac{\nu U}{x}} \left( \eta f'(\eta) - f(\eta) \right)$$

Substituting into the momentum equation yields the celebrated **Blasius Equation**:
$$f f'' + 2 f''' = 0 \quad \Longleftrightarrow \quad f''' + \frac{1}{2} f f'' = 0$$

Boundary conditions transform to:
$$f(0) = 0, \quad f'(0) = 0, \quad \lim_{\eta \to \infty} f'(\eta) = 1$$

---

## 2. Numerical Technique: Shooting Method with Runge-Kutta 4th Order

The Blasius equation is a nonlinear third-order two-point Boundary Value Problem (BVP). We convert it into a system of three first-order ODEs:
$$y_1 = f, \quad y_2 = f', \quad y_3 = f''$$

$$\frac{d}{d\eta} \begin{bmatrix} y_1 \\ y_2 \\ y_3 \end{bmatrix} = \begin{bmatrix} y_2 \\ y_3 \\ -\frac{1}{2} y_1 y_3 \end{bmatrix}$$

Initial conditions at $\eta = 0$:
$$y_1(0) = 0, \quad y_2(0) = 0, \quad y_3(0) = s$$
where $s = f''(0)$ is the unknown initial curvature (wall shear stress parameter).

### 2.1 Shooting Algorithm
1. Choose an asymptotic boundary layer infinity $\eta_\infty = 10.0$.
2. For an initial guess $s$, integrate the initial value problem from $\eta = 0$ to $\eta = \eta_\infty$ using 4th-order Runge-Kutta (RK4 / DOP853).
3. Evaluate the shooting residual:
   $$g(s) = y_2(\eta_\infty; s) - 1.0 = 0$$
4. Update $s$ using Brent's method / Newton iteration until $|g(s)| < 10^{-10}$.

### 2.2 Computed Key Results
| Parameter | Definition / Formula | Numerical Value | Literature Value |
| :--- | :--- | :--- | :--- |
| **Initial Curvature** | $f''(0)$ | **`0.33205734`** | `0.33205734` |
| **Displacement Thickness** | $\frac{\delta^*}{x}\sqrt{Re_x} = \int_0^\infty (1 - f') d\eta$ | **`1.7208`** | `1.7208` |
| **Momentum Thickness** | $\frac{\theta}{x}\sqrt{Re_x} = \int_0^\infty f'(1 - f') d\eta = 2 f''(0)$ | **`0.6641`** | `0.6641` |
| **99% Boundary Layer Thickness** | $\delta_{99} / \sqrt{\frac{\nu x}{U}}$ (where $u/U = 0.99$) | **`4.9149`** | `~4.91` |
| **Shape Factor** | $H = \delta^* / \theta$ | **`2.5911`** | `2.591` |
| **Skin Friction Coefficient** | $C_f = \frac{2 \tau_w}{\rho U^2} = \frac{2 f''(0)}{\sqrt{Re_x}}$ | $\frac{0.6641}{\sqrt{Re_x}}$ | $\frac{0.664}{\sqrt{Re_x}}$ |
| **Normal Velocity Asymptote** | $\lim_{\eta \to \infty} v \sqrt{\frac{x}{\nu U}} = \frac{1}{2}(\eta f' - f)_\infty$ | **`0.8604`** | `0.86` |

---

## 3. Output Plots

The code produces two primary figures matching the assignment specifications:

1. **`fig1.png` (Fig. 1 in PDF):**
   - Longitudinal velocity profile $\frac{u}{U} = f'(\eta)$ vs $\eta = y \sqrt{\frac{U}{\nu x}}$.
   - $x$-axis spanning $0 \le \eta \le 7$.
   - Explicit markers for momentum thickness $\theta \approx 0.664$ and displacement thickness $\delta^* \approx 1.721$.

2. **`fig2.png` (Fig. 2 in PDF):**
   - Transverse (normal) velocity component $v \sqrt{\frac{x}{\nu U}} = \frac{1}{2}(\eta f' - f)$ vs $\eta = y \sqrt{\frac{U}{\nu x}}$.
   - $x$-axis spanning $0 \le \eta \le 6$.
   - Horizontal asymptotic dashed line at $v \sqrt{\frac{x}{\nu U}} = 0.86$ showing the displacement-induced outward flow from the boundary layer.

3. **`blasius_fig1_fig2_combined.png`:**
   - Side-by-side display of Fig. 1 and Fig. 2 for comprehensive report presentation.

---

## 4. File Structure

```text
assignment_2/
├── blasius_solver.py                 # Numerical solver (Shooting method + RK4)
├── requirements.txt                  # Dependencies
├── fig1.png                          # Output Fig. 1 (u/U vs eta)
├── fig2.png                          # Output Fig. 2 (v*sqrt(x/nuU) vs eta)
├── blasius_fig1_fig2_combined.png    # Combined Fig. 1 & Fig. 2
└── README.md                         # Detailed derivation and report
```

---

## 5. How to Run

Install requirements:
```bash
pip install -r requirements.txt
```

Run the solver:
```bash
python3 blasius_solver.py
```
