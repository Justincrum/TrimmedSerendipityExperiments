from firedrake import *
import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(
        description="Allows for input of order and mesh refinement.")
parser.add_argument("-O", "--Order", 
        type=int, help="Input the order of the polynomials.")
args = parser.parse_args()

for n in range(args.Order, args.Order + 1):

    ###Setting up the mesh.
    polyDegree = n
    mesh = Mesh('trapezoids.msh')
    pickASpace = FunctionSpace(mesh, "SminusCurl", polyDegree)
    dofs = pickASpace.dim()
        
    ###Solve the problem.
    x, y = SpatialCoordinate(mesh)
    uex = sin(pi*x)*sin(pi*y)
    sigmaex = grad(uex)
    params={"mat_type": "aij", 
                "ksp_type": "cg", 
                "pc_type": "bjacobi", 
                "sub_pc_type": "ilu", 
                "ksp_rtol": 1e-10}
    #For scalar elements.
    #err = errornorm(uex, project(uex, pickASpace, 
    #                solver_parameters=params))
    err = errornorm(sigmaex, project(sigmaex, pickASpace, 
                       solver_parameters=params)) #For vector elements.
    currentData = [polyDegree, dofs, err]
    u = project(sigmaex, pickASpace, solver_parameters=params)
print(currentData)


#import matplotlib.pyplot as plt
#import numpy as np

#fig, axes = plt.subplots()
#mesh = UnitSquareMesh(10,10)
#V = FunctionSpace(mesh, "CG", 1)
#x = SpatialCoordinate(mesh)
#u = Function(V)
#u.interpolate(x[0]+x[1])

#levels = np.linspace(0, 6, 51)
#contours = tricontourf(u, levels=levels, axes=axes, cmap="inferno")
#axes.set_aspect("equal")
#fig.colorbar(contours)
#fig.show()
#plt.show()
