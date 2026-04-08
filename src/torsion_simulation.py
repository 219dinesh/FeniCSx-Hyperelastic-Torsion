from dolfinx import default_scalar_type
from dolfinx.fem.petsc import NonlinearProblem
import numpy as np
import ufl
from mpi4py import MPI
from dolfinx import fem, mesh, io


#--------------------------------------------
# GEOMETRY AND MESH SETUP
#--------------------------------------------

Length = 1.0
Width = 1.0
Height = 5.0

# Create a 3D Column
domain = mesh.create_box(MPI.COMM_WORLD, [[0.0, 0.0, 0.0], [Length, Width, Height]], [6, 6, 40], mesh.CellType.hexahedron)

V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))

#--------------------------------------------
# BOUNDARY CONDITIONS
#--------------------------------------------

# Bottom Boundary (z=0): Fixed
def bottom(x):
     return np.isclose(x[2], 0.0)
 
fdim = domain.topology.dim - 1
bottom_dofs = fem.locate_dofs_geometrical(V, bottom)

u_bc = np.array([0.0, 0.0, 0.0], dtype=default_scalar_type)

bc_bottom = fem.dirichletbc(u_bc, bottom_dofs, V)

# Top Boundary (z = Height) : Twisted

theta = fem.Constant(domain, default_scalar_type(0.0))

x = ufl.SpatialCoordinate(domain)
x_new =  ufl.as_vector((
    x[0] * (ufl.cos(theta)-1) - x[1] * ufl.sin(theta),
    x[0] * ufl.sin(theta) + x[1]*(ufl.cos(theta) - 1),
    0.0
))

expr_top = fem.Expression(x_new, V.element.interpolation_points)

def top(x):
    return np.isclose(x[2], Height)

u_top = fem.Function(V)

top_dofs = fem.locate_dofs_geometrical(V, top)

bc_top = fem.dirichletbc(u_top, top_dofs)

bcs = [bc_bottom, bc_top]

#--------------------------------------------
#  BODY FORCES AND TRACTION
#--------------------------------------------

B = fem.Constant(domain, default_scalar_type((0, 0, 0)))
T = fem.Constant(domain, default_scalar_type((0, 0, 0)))

#--------------------------------------------
# TEST AND SOLUTION FUNCTIONS
#--------------------------------------------

v = ufl.TestFunction(V)
u = fem.Function(V)

#--------------------------------------------
# KINEMATICS AND MATERIAL LAW
#--------------------------------------------

# Spatial Dimension
d =  len(u)

# Identity Tensor
I = ufl.variable(ufl.Identity(d))

# Deformation gradient
F = ufl.variable(I + ufl.grad(u))

# Right Cauchy Green tensor
C = ufl.variable(F.T * F)

# Invariants of deformation tensors
Ic = ufl.variable(ufl.tr(C))
J = ufl.variable(ufl.det(F))

# Material Parameters

E = default_scalar_type(1.0e2)
nu = default_scalar_type(0.3)
mu = fem.Constant(domain, E/(2 * (1 + nu)))
lmbda = fem.Constant(domain, E * nu / ((1 + nu) * (1 - 2 * nu)))

# Strain Energy Density function for Neo Hookean model

psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2

# First Piola- Kirchoff Stress

P = ufl.diff(psi, F)

#--------------------------------------------
# WEAK FORM OF EQUILIBRIUM EQUATION
#--------------------------------------------

residual = (
    ufl.inner(ufl.grad(v), P) * ufl.dx - ufl.inner(v, B) * ufl.dx - ufl.inner(v, T) * ufl.ds
)

#--------------------------------------------
# NON-LINEAR SOLVER SETUP
#--------------------------------------------

petsc_options = {
    "snes_type": "newtonls",
    "snes_linesearch_type": "none",
    "snes_monitor": None,
    "snes_atol": 1e-8,
    "snes_rtol": 1e-8,
    "snes_stol": 1e-8,
    "ksp_type": "preonly",
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",    
}

problem = NonlinearProblem(
    residual,
    u,
    bcs=bcs,
    petsc_options=petsc_options,
    petsc_options_prefix="Torsion",
)

with io.XDMFFile(domain.comm, "deformation_torsion.xdmf", "w") as xdmf :
    xdmf.write_mesh(domain)
    
    V_out = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))
    u_out = fem.Function(V_out)
    
    Target_angle = np.pi * 6
    steps = 100
    
    for step in range (1, steps + 1):
        # Calculate Current Angle
        current_angle = step * (Target_angle/steps)
        theta.value = current_angle
        u_top.interpolate(expr_top)
        problem.solve()
        
        u_out.interpolate(u)
        
        xdmf.write_function(u_out, current_angle)

print("Simulation Completed")
