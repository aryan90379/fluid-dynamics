# AE 683 (Fluid Dynamics) 2026-1
## Assignment 4: Compressible Aerodynamics Numerical Solutions

**Student Name:** Aryan Meena  
**Roll Number:** 23B0074  

---

## 1. Part (a): Area - Mach Number Relation (Fig. 4)

### 1.1 Governing Equation
For 1D steady, isentropic flow of a calorically perfect gas through a variable-area duct (de Laval nozzle), conservation of mass, momentum, and energy yields the **Area - Mach number relation**:

$$\frac{A}{A^*} = \frac{1}{M} \left[ \frac{2}{\gamma + 1} + \frac{\gamma - 1}{\gamma + 1} M^2 \right]^{\frac{\gamma + 1}{2(\gamma - 1)}}$$

For air ($\gamma = 1.4$):
$$\frac{\gamma + 1}{2(\gamma - 1)} = \frac{2.4}{2(0.4)} = 3.0$$
$$\frac{A}{A^*} = \frac{1}{M} \left( \frac{5 + M^2}{6} \right)^3$$

### 1.2 Physical Properties
1. **Throat Condition ($A = A^*$):** At $M = 1$, the area ratio reaches its minimum $\frac{A}{A^*} = 1$. The curve exhibits a vertical tangent $\frac{d(A/A^*)}{dM} = 0$.
2. **Dual Solutions:** For every area ratio $A/A^* > 1$, there exist exactly two physical Mach numbers:
   - **Subsonic solution ($0 < M < 1$):** Converging section accelerates flow; diverging section decelerates flow (diffuser).
   - **Supersonic solution ($M > 1$):** Diverging section accelerates flow (supersonic nozzle).
3. **Asymptotic Behavior:**
   - As $M \to 0$: $\frac{A}{A^*} \approx \left(\frac{5}{6}\right)^3 \frac{1}{M} \to \infty$.
   - As $M = 5.0$: $\frac{A}{A^*} = \frac{1}{5} \left(\frac{5 + 25}{6}\right)^3 = \frac{1}{5} (5)^3 = 25.0$.

### 1.3 Numerical Values (Air, $\gamma = 1.4$)
| $A / A^*$ | Subsonic Root $M$ | Supersonic Root $M$ |
| :---: | :---: | :---: |
| **`1.00`** | `1.00000` | `1.00000` |
| **`1.50`** | `0.43026` | `1.85412` |
| **`2.00`** | `0.30590` | `2.19720` |
| **`5.00`** | `0.11669` | `3.17478` |
| **`10.00`**| `0.05799` | `3.92255` |
| **`25.00`**| `0.02316` | **`5.00000`** |

---

## 2. Part (b): Oblique Shock Wave Solution ($\theta - \beta - M$ Relation, Fig. 5)

### 2.1 Governing Equation
Across an oblique shock wave with deflection angle $\theta$, shock wave angle $\beta$, and freestream Mach number $M_1$:

$$\tan\theta = 2 \cot\beta \left[ \frac{M_1^2 \sin^2\beta - 1}{M_1^2 (\gamma + \cos 2\beta) + 2} \right]$$

### 2.2 Limiting Behaviors
1. **Infinitesimal Wave (Mach Wave):**
   When the shock strength approaches zero, $M_{n1} = M_1 \sin\beta \to 1$:
   $$\beta \to \mu = \arcsin\left(\frac{1}{M_1}\right), \qquad \theta \to 0$$
2. **Normal Shock Wave:**
   When $\beta = 90^\circ$ ($\frac{\pi}{2}$ rad), $\cot\beta = 0$:
   $$\theta = 0$$
3. **Infinite Mach Number Limit ($M_1 \to \infty$):**
   Dividing numerator and denominator by $M_1^2$:
   $$\tan\theta = \frac{\sin 2\beta}{\gamma + \cos 2\beta}$$
   The theoretical maximum deflection angle for $\gamma = 1.4$ is:
   $$\theta_{\max,\infty} = \arctan\left(\frac{1}{\sqrt{\gamma^2 - 1}}\right) = \arctan\left(\frac{1}{\sqrt{0.96}}\right) \approx 45.58^\circ \quad \text{at } \beta \approx 67.79^\circ$$

### 2.3 Physical Interpretation of the $\theta - \beta - M$ Diagram
- For any given deflection $\theta < \theta_{\max}$, there exist two possible oblique shock angles:
  - **Weak Shock Solution (Lower $\beta$):** Smaller wave angle, downstream flow usually remains supersonic ($M_2 > 1$, except near $\theta_{\max}$). This is the solution naturally observed in unconfined external supersonic flows.
  - **Strong Shock Solution (Upper $\beta$):** Larger wave angle approaching $90^\circ$, downstream flow is always subsonic ($M_2 < 1$). Occurs under high downstream back-pressure or nozzle blockage.
- **Detached Shock Regime:** If the body wedge angle $\theta > \theta_{\max}(M_1)$, no attached oblique shock solution exists. The shock detaches and stands upstream of the body as a curved bow shock.

### 2.4 Computed Maximum Deflection Angles ($\gamma = 1.4$)
| $M_1$ | Mach Angle $\mu$ | $\theta_{\max}$ | $\beta$ at $\theta_{\max}$ |
| :---: | :---: | :---: | :---: |
| **`1.2`** | $56.44^\circ$ | **$3.94^\circ$** | $71.96^\circ$ |
| **`1.4`** | $45.58^\circ$ | **$9.43^\circ$** | $67.73^\circ$ |
| **`1.6`** | $38.68^\circ$ | **$14.65^\circ$** | $65.81^\circ$ |
| **`2.0`** | $30.00^\circ$ | **$22.97^\circ$** | $64.65^\circ$ |
| **`3.0`** | $19.47^\circ$ | **$34.07^\circ$** | $65.22^\circ$ |
| **`4.0`** | $14.48^\circ$ | **$38.77^\circ$** | $66.04^\circ$ |
| **`5.0`** | $11.54^\circ$ | **$41.12^\circ$** | $66.59^\circ$ |
| **`10.0`**| $5.74^\circ$ | **$44.43^\circ$** | $67.48^\circ$ |
| **`$\infty$`** | $0.00^\circ$ | **$45.58^\circ$** | $67.75^\circ$ |

---

## 3. Output Plots

1. **`fig4.png` (Fig. 4 in PDF):**
   - Area Mach number curve in red spanning $A/A^* \in [0, 30]$ and $M \in [0, 5]$.
   - Highlights the sonic throat at $(A/A^* = 1, M = 1)$, subsonic lower branch, and supersonic upper branch ending at $(25, 5.0)$.

2. **`fig5.png` (Fig. 5 in PDF):**
   - $\theta - \beta - M$ shock polar diagram with deflection angle $\theta \in [0^\circ, 45.6^\circ]$ and wave angle $\beta \in [0^\circ, 90^\circ]$.
   - Curves for $M_1 = 1.05, 1.1, 1.2, 1.4, 1.6, 2.0, 3, 4, 5, 10, \infty$ with labeled contours.

---

## 4. File Structure

```text
assignment_4/
├── compressible_flow_solver.py       # Numerical solver for Area-Mach and Oblique Shock
├── requirements.txt                  # Dependencies (numpy, matplotlib, scipy)
├── fig4.png                          # Output Fig. 4 (Area - Mach relation)
├── fig5.png                          # Output Fig. 5 (theta - beta - M relation)
└── README.md                         # Mathematical derivations and documentation
```

---

## 5. How to Run

Install dependencies:
```bash
pip install -r requirements.txt
```

Run solver:
```bash
python3 compressible_flow_solver.py
```
