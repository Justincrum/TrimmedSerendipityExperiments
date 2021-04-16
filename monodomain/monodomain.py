#Solving monodomain equation in Firedrake
#Justin Crum
#4/16/21

#To do:
#1.  Implement iIon as in Ten Tusscher
#2.  Implement the time stepping (currently left blank, see the Firedrake example Josh sent
#    as well as the Niederer paper for proper time stepping values to use)
#3.  Check that sigma is set up properly (currently using Vincent paper values)
#4.  Once I set up time stepping, I'll also need to set up initial conditions
#    including what the activiation in a corner will look like.

from firedrake import *
import argparse

parser = argparse.ArgumentParser(
        description="Allows for input of order and mesh refinement.")
parser.add_argument("-O", "--Order", 
        type=int, help="Input the order of the polynomials.")
args = parser.parse_args()

#Create an extruded mesh with geometry 20x7x3, with cube elements of size 0.1x0.1x0.1
msh = RectangleMesh(200, 70, 20, 7, quadrilateral=True)
mesh = ExtrudedMesh(msh, layers=30, layer_height=0.1)
polyOrder = args.Order

#Set up the function space and test/trial functions.
h1Space = FunctionSpace(mesh, "S", polyOrder)
u = TrialFunction(h1Space)
v = TestFunction(h1Space)

#Model parameters, taken from the literature.
chi = Constant(140)    
iIon = Constant(1)   #This won't actually be a constant later.
capacitance = Constant(0.01)
sigma =  as_vector([0.133, 0.0176, 0.0176])

#Time step for dU/dT -> dUt
dUt = 

#Problem setup.
a = (chi * v * capacitance * dUt + chi * v * iIon)*dx
L = -(grad(v) * sigma * grad(u))*dx

#No need to set up explicit boundary conditions.  We have already
#applied the Neumann BC of \grad u \cdot n = 0 by getting rid of 
#the boundary integral including a \grad u \cdot n factor.

#Set up solver parameters and solve the problem.
w = Function(h1Space)

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
          "mat_mumps_icntl_14": "200",
          "mat_mumps_icntl_11": "2"}

solve(a == L, w, solver_params=params)



