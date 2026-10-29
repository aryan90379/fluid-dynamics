# AE 683 (Fluid Dynamics) 2026-1
## Assignment 3: Numerical Solution of the Falkner-Skan Equation

**Student Name:** Aryan Meena  
**Roll Number:** 23B0074  

---

## 1. Problem Formulation

The **Falkner-Skan boundary layer family** describes self-similar laminar flows over wedges where the potential external velocity follows a power law:
$$U(x) = C x^n$$

From Bernoulli's equation, the pressure gradient in the external stream is:
$$-\frac{1}{\rho}\frac{dp}{dx} = U \frac{dU}{dx} = n C^2 x^{2n - 1}$$
- If $n > 0$: $\frac{dp}{dx} < 0$ (**favorable pressure gradient**, accelerated flow).
- If $n = 0$: $\frac{dp}{dx} = 0$ (**zero pressure gradient**, Blasius flat-plate flow).
- If $n < 0$: $\frac{dp}{dx} > 0$ (**adverse pressure gradient**, decelerated flow prone to boundary layer separation).

### 1.1 Similarity Transformation
Using the similarity variable:
$$\eta = y \sqrt{\frac{U(x)}{\nu x}} = y \sqrt{\frac{C x^n}{\nu x}}$$
$$\psi(x, y) = \sqrt{\nu x U(x)} \, f(\eta)$$

The momentum boundary layer equation reduces to the **Falkner-Skan Equation**:
$$f''' + \frac{n + 1}{2} f f'' - n (f')^2 + n = 0$$

Equivalently written as:
$$f''' + \frac{n + 1}{2} f f'' + n \left(1 - (f')^2\right) = 0$$

Boundary conditions:
$$f(0) = 0, \quad f'(0) = 0, \quad \lim_{\eta \to \infty} f'(\eta) = 1$$

---

## 2. Relation to Hartree Parameter $\beta$

By defining the scaled similarity variable:
$$\xi = \left\{\frac{n + 1}{2}\right\}^{1/2} \eta, \qquad F(\xi) = \left\{\frac{n + 1}{2}\right\}^{1/2} f(\eta)$$

The equation transforms to the standard **Hartree form**:
$$F''' + F F'' + \beta \left(1 - (F')^2\right) = 0$$
where the Hartree pressure gradient parameter is:
$$\beta = \frac{2n}{n + 1} \iff n = \frac{\beta}{2 - \beta}$$

The wedge angle corresponding to $\beta$ is $\theta_{\text{wedge}} = \beta \pi$.

### 2.1 Investigated Cases
| $n$ | Hartree $\beta = \frac{2n}{n+1}$ | Wedge Angle $\beta\pi$ | Flow Description | $f''(0)$ |
| :---: | :---: | :---: | :--- | :---: |
| **`4.0`** | `1.6000` | $1.600\pi$ | Strongly accelerating flow | **`2.40572`** |
| **`1.0`** | `1.0000` | $\pi$ | 2D Hiemenz stagnation-point flow | **`1.23259`** |
| **`0.5`** | `0.6667` | $\frac{2}{3}\pi$ | Moderately accelerating wedge flow | **`0.89972`** |
| **`0.0`** | `0.0000` | $0$ | Blasius flat-plate boundary layer | **`0.33206`** |
| **`-0.0654`** | `-0.1400` | $-0.140\pi$ | Decelerating flow with adverse pressure gradient | **`0.16395`** |
| **`-0.0904`** | `-0.1988` | $-0.199\pi$ | **Incipient Separation Profile** ($f''(0) \approx 0$) | **`0.00498`** |

---

## 3. Numerical Implementation

### 3.1 Numerical Continuation
Solving Falkner-Skan for adverse pressure gradients ($n < 0$) is notoriously ill-conditioned for standard shooting due to solution non-uniqueness (separated vs unseparated branches). 

To ensure convergence:
1. The third-order ODE is formulated as a first-order system of 3 equations:
   $$y_1' = y_2, \quad y_2' = y_3, \quad y_3' = -\frac{n+1}{2} y_1 y_3 + n (y_2^2 - 1)$$
2. For $n \ge 0$, solved directly with `scipy.integrate.solve_bvp`.
3. For $n < 0$, a **numerical parameter continuation** scheme is employed, marching in small steps $\Delta n$ from $n = 0$ down to $n = -0.0904$, using the previous converged solution as the initial guess for the next.

---

## 4. Output Plot

- **`fig3.png` (Fig. 3 in PDF):**
  - Plots normalized velocity $\frac{u}{U} = f'(\eta)$ on the vertical axis ($0 \le u/U \le 1.0$).
  - Plots scaled similarity coordinate $\left\{\frac{1}{2}(n+1)\right\}^{1/2}\eta$ on the horizontal axis ($0 \le \xi \le 4.8$).
  - Displays all 6 profiles with annotations:
    - $n = 4$
    - $n = 1$
    - $n = 1/2$
    - $n = 0$
    - $n = -0.0654$
    - $n = -0.0904$
  - Highlights how adverse pressure gradients cause the profile to become S-shaped (inflectional) with zero wall slope at separation ($f''(0) \to 0$).

---

## 5. File Structure

```text
assignment_3/
├── falkner_skan_solver.py            # BVP solver with numerical continuation
├── requirements.txt                  # Dependencies (numpy, matplotlib, scipy)
├── fig3.png                          # Output Fig. 3 matching Assignment 3
└── README.md                         # Mathematical derivation and documentation
```

---

## 6. How to Run

Install dependencies:
```bash
pip install -r requirements.txt
```

Run solver:
```bash
python3 falkner_skan_solver.py
```
