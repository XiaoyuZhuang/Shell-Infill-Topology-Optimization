# Shell-Infill Topology Optimization (DSP & Erosion Approach)

[![MATLAB](https://img.shields.io/badge/MATLAB-Research-blue)](https://www.mathworks.com/)

> **Note to Researchers:**
> Reproduction of topology optimization papers can be an arduous and painful process. I have synthesized methods from multiple papers to build this Shell-Infill optimization framework. While this code is experimental and may not be perfectly stable in all scenarios, I am sharing it in the spirit of open source to save others time.
>
> **致研究者：**
> 复现拓扑优化论文往往是一个痛苦的过程。本项目参考了多篇文献，构建了这套基于腐蚀法的壳-填充优化框架。虽然代码在某些参数下可能不够稳定，但我希望通过开源分享，能为正在攻克类似问题的你提供一些参考和帮助。

## 📖 Introduction

This repository contains a MATLAB implementation for generating **Shell-Infill structures** (structures with a stiff outer coating and a compliant porous core).

The method relies on a **Double-Filter / Double-Projection** strategy:
1.  **Macro Filter (DSP):** Defines the base topology and controls the minimum length scale.
2.  **Erosion Filter:** "Erodes" the base topology to identify the inner core.
3.  **Boolean Subtraction:** The Shell is mathematically extracted by subtracting the Core from the Base.

## 🧠 Methodology & Principles

### 1. Geometric Representation (How the Shell is Extracted)

The core idea is to generate two geometric fields from a single design variable field $\mathbf{x}$.

* **Base Structure ($\mu$):**
    Generated using standard **Design Space Projection (DSP)** (or Heaviside Projection). This defines the overall shape of the structure.
    $$\mu = H(\tilde{x}, \beta, \eta_{dsp})$$

* **Eroded Core ($\phi$):**
    We take the Base Structure $\mu$ and apply a second filtering and projection step with a **higher threshold** (Erosion). This shrinks the boundaries of $\mu$.
    $$\phi = H(\tilde{\mu}, \beta, \eta_{erosion})$$
    * If $\eta_{erosion} > 0.5$, the structure shrinks (Erosion).
    * The difference between the original boundary and the eroded boundary determines the **Shell Thickness**.

* **The Shell:**
    The shell is simply the region existing in the Base but not in the Core:
    $$\text{Region}_{shell} = \mu - \phi$$

### 2. Material Interpolation Scheme

To optimize for a stiff shell and a compliant infill, we use a multi-phase SIMP (Solid Isotropic Material with Penalization) interpolation.

The Young's Modulus $E$ at any element is interpolated as:

$$E(\mathbf{x}) = E_{min} + (\underbrace{\mu - \phi}_{\text{Shell}})^p E_{shell} + (\underbrace{\phi}_{\text{Infill}})^p E_{infill}$$

* $E_{shell}$: Modulus of the stiff coating (e.g., 1.0).
* $E_{infill}$: Modulus of the internal lattice/porous material (e.g., 0.2).
* $p$: Penalization power (typically 3).

### 3. Workflow Visualized

1.  **Design Variable ($x$)** $\xrightarrow{\text{Filter } R_1}$ **Density** $\xrightarrow{\text{Proj}}$ **Base Topology ($\mu$)**
2.  **Base Topology ($\mu$)** $\xrightarrow{\text{Filter } R_2}$ **Eroded Field** $\xrightarrow{\text{Proj (High }\eta)}$ **Core ($\phi$)**
3.  **Physical Model:** Assign $E_{shell}$ to $(\mu-\phi)$ and $E_{infill}$ to $\phi$.
4.  **FEA & Sensitivity:** Solve $Ku=F$ and update $x$ using GDA (Gradient Descent with Adaptive step).

## ⚙️ Key Parameters

Because this method relies on geometric convolution, the **Filter Radii** are crucial:

* `R_dsp` (Filter 1): Controls the minimum feature size of the *entire* structure.
* `R_ero` (Filter 2): Controls the **Shell Thickness**. A larger `R_ero` results in a thicker shell (as the core erodes more).
* `eta_ero` (Erosion Threshold): Controls how aggressive the erosion is. Typically set around `0.7` to `0.8`.

## ⚠️ Disclaimer & Tips

* **Stability:** This is a non-linear problem with high non-convexity. If the optimization oscillates, try reducing the `step_size` or increasing the `beta` update interval.
* **Boundary Conditions:** The code includes a specific treatment for the left boundary (Symmetry) to prevent the shell from wrapping around the cut plane. You may need to adjust `pad_size` logic for different BCs.
* **Beta Continuation:** The projection sharpness (`beta`) increases gradually. If the shell is not distinct enough, check if `beta` has reached a sufficiently high value (e.g., 32 or 64).

## 🤝 Contribution

Issues and Pull Requests are welcome! If you find a way to stabilize the convergence further or implement a better optimizer (like MMA), please feel free to contribute.

## 📧 Contact

**Xiaoyu Zhuang**
* Dalian University of Technology
* Research Interest: Topology Optimization, Multi-scale Structures
