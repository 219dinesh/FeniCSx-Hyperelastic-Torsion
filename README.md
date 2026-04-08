# 3D Hyperelastic Torsion Simulation in FEniCSx

![FEniCSx](https://img.shields.io/badge/FEniCSx-DOLFINx-blue.svg)
![Python](https://img.shields.io/badge/Python-3.10+-yellow.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)

This repository contains a high-performance 3D computational mechanics simulation of a column undergoing large torsional deformation. The model is implemented using **FEniCSx (DOLFINx)** and utilizes a robust nonlinear solver with the MUMPS direct solver for optimal computational speed.

## Visualizing the Torsion
*![Torsion Animation](media/torsion_animation.gif)*

## Mathematical Formulation

The simulation models large strain hyperelasticity using a **Neo-Hookean** material model. 

### Kinematics
* Deformation Gradient: $F = I + \nabla u$
* Right Cauchy-Green Tensor: $C = F^T F$
* Invariants: $I_c = \text{tr}(C)$ and $J = \det(F)$

### Strain Energy Density Function
The hyperelastic behavior is governed by the following strain energy density function ($\psi$):

$$\psi = \frac{\mu}{2} (I_c - 3) - \mu \ln(J) + \frac{\lambda}{2} (\ln(J))^2$$

Where $\mu$ and $\lambda$ are the Lamé parameters derived from Young's Modulus ($E = 100$) and Poisson's ratio ($\nu = 0.3$).

### Weak Form of the Equilibrium Equation

The simulation solves for the displacement field $u \in V$ such that for all test functions $v \in \hat{V}$, the virtual work equation is satisfied:

$$\int_{\Omega} P : \nabla v \ dx - \int_{\Omega} B \cdot v \ dx - \int_{\partial \Omega_T} T \cdot v \ ds = 0$$

**Where:**
* **$P$**: The First Piola-Kirchhoff stress tensor.
* **$\nabla v$**: The gradient of the test function with respect to the reference configuration.
* **$B$**: Body forces per unit reference volume (e.g., gravity).
* **$T$**: Traction forces applied to the boundary surface $\partial \Omega_T$.

The first Piola-Kirchhoff stress tensor is defined as $P = \frac{\partial \psi}{\partial F}$. The weak form of the equilibrium equation (ignoring body forces and tractions for this specific problem) is formulated as:

$$\int_{\Omega} P : \nabla v \ dx = 0$$

## Installation and Setup

This code is designed to run in a FEniCSx environment, typically managed via Conda or Docker. It is highly recommended to run this within a **WSL (Windows Subsystem for Linux)** environment if you are on a Windows machine.
