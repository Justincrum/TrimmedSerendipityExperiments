#Firedrake script for solving pde
#on reentrant corner mesh.

from firedrake import *
import numpy as np

mesh = Mesh('corner_mesh.msh')

#mesh = UnitSquareMesh(8, 8, quadrilateral=True)

V = FunctionSpace(mesh, "S", 2)

u = TrialFunction(V)
v = TestFunction(V)

n = FacetNormal(mesh)

f = Constant(0.0)
a = inner(grad(u), grad(v))*dx - inner(dot(grad(u), n), v)*ds
L = inner(v,f)*dx

#set up the solution
x, y = SpatialCoordinate(mesh)
r = (x**2 + y**2)**(0.5)
alpha =  2/3
theta = atan(y / x) + pi

wex = r**alpha * sin(alpha * theta)


BC0 = DirichletBC(V, 0.0, [8])
BC1 = DirichletBC(V, wex, [9])


w = Function(V)

params = {"snes_type": "newtonls",
                  "snes_linesearch_type": "basic",
                  "snes_monitor": None,
                  "snes_converged_reason": None,
                  "mat_type": "aij",
                  "snes_max_it": 10,
                  "snes_lag_jacobian": -2,
                  "snes_lag_preconditioner": -2,
                  "ksp_type": "preonly",
                  "ksp_converged_reason": None,
                  "ksp_monitor_true_residual": None,
                  "pc_type": "lu",
                  "snes_rtol": 1e-12,
                  "snes_atol": 1e-20,
                  "pc_factor_mat_solver_type": "mumps",
                  "mat_mumps_icntl_14": "1000"}


solve(a == L, w, bcs=[BC0, BC1], solver_parameters=params)

err = norms.errornorm(wex, w)
print(err)
outfile = File("reentrant_corner.pvd")
outfile.write(w)

