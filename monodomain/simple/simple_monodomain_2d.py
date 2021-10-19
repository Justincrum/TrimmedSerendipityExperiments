#Solving monodomain equation in Firedrake
#Justin Crum
#4/16/21

#To do:
#1.  Implement iIon as in Ten Tusscher
#3.  Need to set up initial conditions including what the activiation in
#    a corner will look like.

from firedrake import *
from irksome import GaussLegendre, Dt, TimeStepper, RadauIIA
import argparse

# parser = argparse.ArgumentParser(
#         description="Allows for input of order and mesh refinement.")
# parser.add_argument("-O", "--Order", 
#         type=int, help="Input the order of the polynomials.")
# args = parser.parse_args()


#Create an extruded mesh with geometry 20x7x3, with cube elements of size 0.1x0.1x0.1
#mesh = Mesh('simple_rectangle_mesh.msh')
mesh = RectangleMesh(20, 7, 20, 7, quadrilateral=True)
#mesh = ExtrudedMesh(msh, layers=3)#, layer_height=1)
#polyOrder = args.Order
polyOrder = 2

#Set up the function space and test/trial functions.
h1Space = FunctionSpace(mesh, "Q", polyOrder)
uCurr = Function(h1Space, name="u")

#u = TrialFunction(h1Space)
v = TestFunction(h1Space)

uCurr.assign(-8.23)

#Set the initial conditions.
x, y = SpatialCoordinate(mesh)

#Set up an output file.  I think I need to set up initial conditions before this
#if I want to use to IC's other than 0.
#outfile = File("monodomain.pvd")
#outfile.write(uCurr)

#Model parameters, taken from the literature.
A = Constant(0.04)
vRest = Constant(-85.23)
vPeak = Constant(40.)
vTh = Constant(-65.)
chi = Constant(140.)    

capacitance = Constant(0.01)
sigma =  as_matrix([[0.133, 0.0], [0.0, 0.0176]])

#Set up Irksome to be used.
butcher_tableau = RadauIIA(1)
ns = butcher_tableau.num_stages

#Time step for dU/dT -> dUt
tFinal = 0.2
dt = Constant(0.0025)
t = Constant(0.0)
step = 0

#print(uCurr((x,y)))


#appStim_On = 50,000 * Function(h1Space).assign(conditional(max(x, y) < 1.5, 1, 0))
#appStimulus = Function(h1Space)

#appStimulus = Function(h1Space).assign(conditional(t <= 0.002, appStim_On, 0))

App = conditional(And(And(x < 1.5, y < 1.5), t < 2), Constant(75), Constant(0.0))
iIon = A**2 * (uCurr - vRest) * (uCurr - vTh) * (uCurr - vPeak) + App

#iIonWApp = A**2 * (uCurr - vRest) * (uCurr - vTh) * (uCurr - vPeak) + 112.5
#iIon = A**2 * (uCurr - vRest) * (uCurr - vTh) * (uCurr - vPeak)
"""
params = {"snes_type": "newtonls",
        "snes_linesearch_type": "basic",
        "snes_monitor": None,
        "snes_converged_reason": None,
        "mat_type": "aij",
        "snes_max_it": 20,
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
"""

params = {"mat_type": "aij",
           "ksp_type": "preonly", #gmres or cg
           "pc_type": "lu", #gamg
           #"snes_lag_jacobian": -2,
            #"snes_lag_preconditioner": -2,
           "snes_atol" : 1e-10,
           "snes_monitor" : None
           }


#F = inner(Dt(uCurr), v)*dx + 1 / (chi * capacitance) * inner(grad(v), sigma * grad(uCurr))*dx -1 / capacitance * inner(iIon, v)*dx(3) - 1 / capacitance * inner(iIonWApp, v)*dx(4)
#F2 = inner(Dt(uCurr), v)*dx + 1 / (chi * capacitance) * inner(grad(v), sigma * grad(uCurr))*dx -1 / capacitance * inner(iIon, v)*dx(3) - 1 / capacitance * inner(iIon, v)*dx(4)

F = inner(Dt(uCurr), v)*dx + 1 / (chi * capacitance) * inner(grad(v), sigma * grad(uCurr))*dx + 1 / capacitance * inner(iIon, v)*dx

stepper = TimeStepper(F, butcher_tableau, t, dt, uCurr,
                      solver_parameters=params)
#stepper2 = TimeStepper(F2, butcher_tableau, t, dt, uCurr,
#                      solver_parameters=params)

#while (float(t)< tFinal):
for j in range(100):
    print("Time step", j)
    #if (float(t) + float(dt) > tFinal):
        #dt.assign(tFinal - float(t))
    stepper.advance()
    t.assign(float(t) + float(dt)) 

outfile = File("simple_monodomain.pvd")
outfile.write(uCurr)




        #Problem setup.
        #a = (v * u)*dx + (dt / (chi * capacitance) * dot(grad(v), sigma * grad(u)))*dx
        #L = (v * uCurr)*dx - (dt / capacitance * iIon * v)*dx
        #F = inner(Dt(u), v)*dx + 1 / (chi * capacitance) * inner(grad(v), sigma * grad(u))*dx -1 / capacitance * inner(iIon, v)*dx 
        #iIon = A**2 * (uCurr - vRest) * (uCurr - vTh) * (uCurr - vPeak)
        #No need to set up explicit boundary conditions.  We have already
        #applied the Neumann BC of \grad u \cdot n = 0 by getting rid of 
        #the boundary integral including a \grad u \cdot n factor.
"""
        #Set up solver parameters and solve.
#        params = {"snes_type": "newtonls",
##                "snes_linesearch_type": "basic",
 #               "snes_monitor": None, 
 #               "snes_converged_reason": None,
 #               "mat_type": "aij",
 #               "snes_max_it": 10,
 #               "snes_lag_jacobian": -2,
 #               "snes_lag_preconditioner": -2,
 #               "ksp_type": "preonly",
 #               "ksp_converged_reason": None,
 #               "ksp_monitor_true_residual": None,
 #               "pc_type": "lu",
 #               "snes_rtol": 1e-12,
 #               "snes_atol": 1e-20,
 #               "pc_factor_mat_solver_type": "mumps",
 #               "mat_mumps_icntl_14": "200",
 #               "mat_mumps_icntl_11": "2"}

        #stepper = TimeStepper(F, butcher_tableau, t, dt, u,
        #              solver_parameters=params)

        #solve(a == L, uCurr, solver_parameters=params)
        #t += dt
#        if step % 5 == 0:
#                outfile.write(uCurr, time=t)

"""