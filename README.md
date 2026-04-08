# 3D Hyperelastic Torsion Simulation in FEniCSx

![FEniCSx](https://img.shields.io/badge/FEniCSx-DOLFINx-blue.svg)
![Python](https://img.shields.io/badge/Python-3.10+-yellow.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)

This repository contains a high-performance 3D computational mechanics simulation of a column undergoing large torsional deformation. The model is implemented using **FEniCSx (DOLFINx)** and utilizes a robust nonlinear solver with the MUMPS direct solver for optimal computational speed.

## Visualizing the Torsion
*(Add a GIF of your ParaView animation here! Save it in the media folder and link it like this: `![Torsion Animation](media/torsion_animation.gif)`)*

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

### Weak Form
The first Piola-Kirchhoff stress tensor is defined as $P = \frac{\partial \psi}{\partial F}$. The weak form of the equilibrium equation (ignoring body forces and tractions for this specific problem) is formulated as:

$$\int_{\Omega} P : \nabla v \, dx = 0$$

## Installation and Setup

This code is designed to run in a FEniCSx environment, typically managed via Conda or Docker. It is highly recommended to run this within a **WSL (Windows Subsystem for Linux)** environment if you are on a Windows machine.
