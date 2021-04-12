from firedrake import *

# load the mesh generated with Gmsh
mesh = Mesh('immersed_domain.msh')

# define the space of linear Lagrangian finite elements
V = FunctionSpace(mesh, "S", 1)

# define the trial function u and the test function v
u = TrialFunction(V)
v = TestFunction(V)

#x, y = SpatialCoordinate(mesh)

# define the bilinear form of the problem under consideration
# to specify the domain of integration, the surface tag is specified in brackets after dx
# in this example, 3 is the tag of the rectangle without the disc, and 4 is the disc tag
a = 2*dot(grad(v), grad(u))*dx(4) + dot(grad(v), grad(u))*dx(3) + v*u*dx

# define the linear form of the problem under consideration
# to specify the boundary of the boundary integral, the boundary tag is specified after dS
# note the use of dS due to 13 not being an external boundary
# Since the dS integral is an interior one, we must restrict the
# test function: since the space is continuous, we arbitrarily pick
# the '+' side.
L = Constant(5.) * v * dx + Constant(3.)*v('+')*dS(13)

# set homogeneous Dirichlet boundary conditions on the rectangle boundaries
# the tag  11 referes to the horizontal edges, the tag 12 refers to the vertical edges
DirBC = DirichletBC(V, 0, [11, 12])

# define u to contain the solution to the problem under consideration
u = Function(V)

# solve the variational problem
solve(a == L, u, bcs=DirBC, solver_parameters={'ksp_type': 'cg'})

outfile = File("immersed_output.pvd")
outfile.write(u)

import matplotlib.pyplot as plt
import numpy as np

fig, axes = plt.subplots()
#mesh = UnitSquareMesh(10,10)
#V = FunctionSpace(mesh, "CG", 1)
#x = SpatialCoordinate(mesh)
#u = Function(V)
#u.interpolate(x[0]+x[1])

levels = np.linspace(0, 6, 51)
contours = tricontourf(u, levels=levels, axes=axes, cmap="inferno")
axes.set_aspect("equal")
fig.colorbar(contours)
fig.show()
plt.show()
