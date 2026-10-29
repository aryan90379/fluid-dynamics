# AE 683 (Fluid Dynamics) 2026-1
## Assignment 1: Streamline Plotting of Combined Potential Flow

**Student Name:** Aryan Meena  
**Roll Number:** 23B0074  
**Assigned Flow:** Sink + Vortex (*A bathtub vortex*)  

---

## 1. Overview & Physical Context

In classical potential flow theory (ideal, incompressible, irrotational, inviscid flow), linear superposition of elementary elementary potential flows allows modeling complex fluid phenomena.

The **Sink + Vortex** flow combination models the classic **"Bathtub Vortex"** (or spiral vortex). In this flow:
- A **sink** at the origin draws fluid radially inward (modeling the drain).
- A **free (irrotational) vortex** at the origin induces a circulating tangential velocity (modeling background angular momentum or swirl).

The combined action creates inward spiraling streamlines that form **logarithmic (equiangular) spirals**, perfectly matching the swirling funnel observed when water drains from a bathtub or container.

---

## 2. Mathematical Formulation

### 2.1 Elementary Potential Flows in Polar Coordinates $(r, \theta)$

#### A. Sink Flow (Strength $m = \frac{\Lambda}{2\pi} > 0$)
- Inward radial velocity:
  $$u_{r, \text{sink}} = -\frac{m}{r} = -\frac{\Lambda}{2\pi r}, \quad u_{\theta, \text{sink}} = 0$$
- Velocity Potential:
  $$\phi_{\text{sink}}(r, \theta) = -m \ln r$$
- Streamfunction:
  $$\psi_{\text{sink}}(r, \theta) = -m \theta$$
- Complex Potential:
  $$W_{\text{sink}}(z) = -m \ln z$$

#### B. Free Vortex (Circulation $\Gamma = 2\pi k > 0$, Counter-Clockwise)
- Tangential swirl velocity:
  $$u_{r, \text{vortex}} = 0, \quad u_{\theta, \text{vortex}} = \frac{k}{r} = \frac{\Gamma}{2\pi r}$$
- Velocity Potential:
  $$\phi_{\text{vortex}}(r, \theta) = k \theta$$
- Streamfunction:
  $$\psi_{\text{vortex}}(r, \theta) = -k \ln r$$
- Complex Potential:
  $$W_{\text{vortex}}(z) = -i k \ln z$$

---

### 2.2 Combined Flow: Sink + Vortex (Bathtub Vortex)

#### Complex Potential $W(z)$:
Using linear superposition:
$$W(z) = W_{\text{sink}}(z) + W_{\text{vortex}}(z) = -(m + ik) \ln z$$

With $z = r e^{i\theta}$, we have $\ln z = \ln r + i\theta$:
$$W(z) = -(m + ik)(\ln r + i\theta) = -(m \ln r - k \theta) - i(k \ln r + m \theta)$$

Separating real and imaginary parts ($W(z) = \phi + i\psi$):
- **Velocity Potential:**
  $$\phi(r, \theta) = -m \ln r + k \theta$$
- **Streamfunction:**
  $$\psi(r, \theta) = -k \ln r - m \theta$$

---

### 2.3 Velocity Field

In polar coordinates:
$$u_r = \frac{\partial \phi}{\partial r} = \frac{1}{r}\frac{\partial \psi}{\partial \theta} = -\frac{m}{r}$$
$$u_\theta = \frac{1}{r}\frac{\partial \phi}{\partial \theta} = -\frac{\partial \psi}{\partial r} = \frac{k}{r}$$

Transforming to Cartesian coordinates $(x, y)$ where $r = \sqrt{x^2 + y^2}$, $\cos\theta = \frac{x}{r}$, and $\sin\theta = \frac{y}{r}$:
$$u = u_r \cos\theta - u_\theta \sin\theta = \left(-\frac{m}{r}\right)\left(\frac{x}{r}\right) - \left(\frac{k}{r}\right)\left(\frac{y}{r}\right) = -\frac{mx + ky}{x^2 + y^2}$$
$$v = u_r \sin\theta + u_\theta \cos\theta = \left(-\frac{m}{r}\right)\left(\frac{y}{r}\right) + \left(\frac{k}{r}\right)\left(\frac{x}{r}\right) = \frac{kx - my}{x^2 + y^2}$$

Velocity magnitude (speed):
$$|\vec{V}| = \sqrt{u_r^2 + u_\theta^2} = \sqrt{u^2 + v^2} = \frac{\sqrt{m^2 + k^2}}{r}$$

---

### 2.4 Streamline Equation: Logarithmic Spiral

A streamline satisfies $\psi(r, \theta) = \text{constant} = C$:
$$-k \ln r - m \theta = C$$
$$\ln r = -\frac{m}{k}\theta - \frac{C}{k}$$
$$r(\theta) = r_0 \exp\left(-\frac{m}{k}\theta\right)$$

This is the equation of an **equiangular (logarithmic) spiral**:
- As the fluid swirls counter-clockwise ($\theta$ increases), its radius $r$ decays exponentially toward the central sink at $r = 0$.
- The angle $\alpha$ between the streamline tangent and the radial direction is strictly constant everywhere:
  $$\tan\alpha = \left|\frac{r \, d\theta}{dr}\right| = \left|\frac{u_\theta}{u_r}\right| = \frac{k}{m} = \frac{\Gamma}{\Lambda} = \text{constant}$$
  - If $k \gg m$: swirl dominates, tight spiral with many turns before draining.
  - If $m \gg k$: suction dominates, nearly radial paths.
  - If $k \approx m$: balanced spiral drain.

---

### 2.5 Orthogonality of the Flow Net

Along an equipotential line ($\phi = \text{constant}$):
$$-m \ln r + k \theta = C_\phi \implies \ln r = \frac{k}{m}\theta - \frac{C_\phi}{m} \implies r(\theta) = r_0 \exp\left(\frac{k}{m}\theta\right)$$

Calculating the slopes in polar coordinates:
$$\left.\frac{dr}{r\,d\theta}\right|_{\psi} = -\frac{m}{k}, \qquad \left.\frac{dr}{r\,d\theta}\right|_{\phi} = +\frac{k}{m}$$

Their product is:
$$\left(-\frac{m}{k}\right) \times \left(\frac{k}{m}\right) = -1$$

Thus, **streamlines and equipotential lines are mutually orthogonal at all points**, forming an ideal conformal flow net.

---

### 2.6 Free Surface Depression (The Bathtub Funnel)

For a liquid with a free surface open to the atmosphere ($p = p_{\text{atm}}$), Bernoulli's equation along the free surface yields:
$$g z(r) + \frac{1}{2} |\vec{V}(r)|^2 = g z_\infty$$
$$z(r) - z_\infty = -\frac{|\vec{V}|^2}{2g} = -\frac{m^2 + k^2}{2g r^2}$$

This quadratic drop with $1/r^2$ forms the characteristic steep hollow vortex core (the "whirlpool funnel") seen above the drain in a bathtub.

---

## 3. Visualizations Generated

The script generates three high-resolution figures:

1. **`bathtub_vortex_streamlines.png`**:
   - Complete streamline pattern plotted using `matplotlib.pyplot.streamplot`.
   - Colormapped by local flow speed $|\vec{V}|$ from deep violet (slow far-field) to bright yellow/orange (fast near-field).
   - Overlaid analytical logarithmic spirals: $r(\theta) = r_0 e^{-(m/k)\theta}$.
   - Indication of the central singularity $(0, 0)$ representing the drain.

2. **`bathtub_vortex_flow_net.png`**:
   - Orthogonal flow net covering all four quadrants.
   - Streamlines ($\psi = \text{const}$, blue solid curves with inward direction arrows).
   - Equipotential lines ($\phi = \text{const}$, red dashed curves).
   - Demonstrates that $\nabla\phi \cdot \nabla\psi = 0$ everywhere in the flow field.

3. **`bathtub_vortex_analysis.png`**:
   - **(a) Streamlines & Speed Distribution:** Streamlines overlaid on continuous velocity magnitude contours.
   - **(b) Velocity Quiver Plot:** Vector grid illustrating local flow direction and magnitude scaling.
   - **(c) Velocity Decay Profile:** Explicit curves comparing $|u_r| = m/r$, $u_\theta = k/r$, and $|\vec{V}| = \sqrt{m^2+k^2}/r$ vs radial distance $r$.
   - **(d) 3D Bathtub Vortex Funnel:** 3D surface plot showing the free surface depression $z(r) = -\frac{m^2+k^2}{2gr^2}$ formed by the vortex core.

---

## 4. File Structure

```text
assignment_1/
├── streamline_plot.py               # Main Python simulation and plotting script
├── requirements.txt                 # Dependencies (numpy, matplotlib, scipy)
├── bathtub_vortex_streamlines.png   # Streamline plot with velocity colormap
├── bathtub_vortex_flow_net.png      # Orthogonal flow net (streamlines + equipotential)
├── bathtub_vortex_analysis.png      # 4-panel comprehensive analysis figure
└── README.md                        # Documentation and mathematical formulation
```

---

## 5. Requirements & Installation

The code requires Python 3.8+ and standard scientific packages:
- `numpy`
- `matplotlib`
- `scipy`

Install dependencies using:
```bash
pip install -r requirements.txt
```

---

## 6. How to Run the Code

To run with default parameters ($m = 1.0$, $k = 1.5$):
```bash
python3 streamline_plot.py
```

### Command-Line Arguments
The script supports flexible command-line customization:
```bash
python3 streamline_plot.py -h
```

| Argument | Flag | Default | Description |
| :--- | :--- | :--- | :--- |
| `--sink-strength` | `-m` | `1.0` | Sink strength parameter $m = \Lambda / (2\pi)$ |
| `--vortex-strength`| `-k` | `1.5` | Vortex circulation parameter $k = \Gamma / (2\pi)$ |
| `--domain-limit` | `-L` | `3.5` | Half-width of square domain $[-L, L]$ |
| `--grid-size` | `-N` | `320` | Grid resolution per dimension |
| `--output-dir` | `-o` | `.` | Directory to save generated plot files |

#### Example Custom Runs:
- **Stronger Swirl:**
  ```bash
  python3 streamline_plot.py -m 0.5 -k 2.5
  ```
- **Stronger Suction:**
  ```bash
  python3 streamline_plot.py -m 2.0 -k 0.8
  ```
- **High Resolution Grid:**
  ```bash
  python3 streamline_plot.py -N 500 -L 4.0
  ```

---

## 7. Key Results & Conclusions

1. **Logarithmic Spiral Trajectories:** The combination of radial inward suction ($u_r \propto 1/r$) and tangential irrotational circulation ($u_\theta \propto 1/r$) produces logarithmic spiral streamlines whose pitch angle $\alpha = \arctan(k/m)$ remains strictly invariant with radial distance.
2. **Velocity Singularity:** Both velocity components diverge as $r \to 0$ ($|\vec{V}| \propto 1/r$), corresponding to the localized core singularity of potential flow theory. In viscous reality, a finite viscous core (such as a Rankine or Lamb-Oseen vortex) regularizes the center.
3. **Flow Net Orthogonality:** The analytical flow net rigorously demonstrates the mutual orthogonality of streamfunction contours $\psi = \text{const}$ and velocity potential contours $\phi = \text{const}$, confirming the conformal mapping properties of analytic complex potentials $W(z)$.
4. **Bathtub Whirlpool Funnel:** The $1/r^2$ Bernoulli depression captures the physics of free-surface drawdown into an open drain.
