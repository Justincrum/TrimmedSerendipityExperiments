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
import numpy as np

# parser = argparse.ArgumentParser(
#         description="Allows for input of order and mesh refinement.")
# parser.add_argument("-O", "--Order", 
#         type=int, help="Input the order of the polynomials.")
# args = parser.parse_args()

def alphaM(potential):
    val = 0.1 * (25 - potential) / (exp((25 - potential) / 10) - 1)
    return val

def alphaH(potential):
    val = 0.07 * exp(-potential / 20)
    return val

def alphaN(potential):
    val = 0.01 * (10 - potential) / (exp((10 - potential) / 10) - 1)
    return val

def betaM(potential):
    val = 4 * exp(-potential / 18)
    return val

def betaH(potential):
    val = 1 / (exp((30 - potential) / 10) + 1)
    return val

def betaN(potential):
    val = 0.125 * exp(-potential / 80)
    return val


#Create an extruded mesh with geometry 20x7x3, with cube elements of size 0.1x0.1x0.1
#mesh = Mesh('simple_rectangle_mesh.msh')
mesh = RectangleMesh(20, 7, 20, 7, quadrilateral=True)
#mesh = ExtrudedMesh(msh, layers=3)#, layer_height=1)
#polyOrder = args.Order
polyOrder = 1

#Set up the function space and test/trial functions.
V = FunctionSpace(mesh, "Q", polyOrder)
Z = V * V * V *V

x, y = SpatialCoordinate(mesh)
t = Constant(0.0)

InitialCond = conditional(And(And(x < 10.0, y < 10.0), t < 0.002), 2.0, 0.0)
temp = Function(V).project(InitialCond)
uu = Function(Z)
uu.sub(0).project(InitialCond)
uu.sub(1).interpolate(Constant(0.0))
uu.sub(2).project(Constant(0.0))
uu.sub(3).project(Constant(0.0))


uCurr, m, n, h = split(uu)
vu, vm, vn, vh = TestFunctions(Z)

uCurr, m, n, h = uu.split()

outfile = File("HHmonodomain_2d.pvd")
outfile.write(uCurr)


"""
#uCurr.assign(-8.23)
x, y = SpatialCoordinate(mesh)

#Model parameters, taken from the literature.
vRest = Constant(-85.23)
vPeak = Constant(40.)
vTh = Constant(-65.)

chi = Constant(140.)    
GNa = 120
GK = 36
GL = 0.3
VNa = 115
VK = -12
VL = 10.6
capacitance = Constant(1.0)
sigma =  as_matrix([[0.133, 0.0], [0.0, 0.0176]])

#Set up Irksome to be used.
butcher_tableau = RadauIIA(1)
ns = butcher_tableau.num_stages

tFinal = 0.2
dt = Constant(0.0025)
t = Constant(0.0)
step = 0

App = conditional(And(And(x < 1.5, y < 1.5), t < 2), Constant(75), Constant(0.0))

iNa = GNa * m**3 * h * (uCurr - VNa)
iK = GK * n**4 * (uCurr - VK)
iL = GL * (uCurr - VL)

iIon = iNa + iK + iL + App

params = {"mat_type": "aij",
           "ksp_type": "preonly", #gmres or cg
           "pc_type": "lu", #gamg
           #"snes_lag_jacobian": -2,
            #"snes_lag_preconditioner": -2,
           "snes_atol" : 1e-10,
           "snes_monitor" : None
           }


Fu = inner(Dt(uCurr), vu)*dx + 1 / (chi * capacitance) * inner(grad(vu), sigma * grad(uCurr))*dx + 1 / capacitance * inner(iIon, vu)*dx
Fm = inner(Dt(m), vm)*dx - inner(alphaM(uCurr) * (1 - m), vm)*dx - inner(betaM(uCurr) * m, vm)*dx
Fn = inner(Dt(n), vn)*dx - inner(alphaN(uCurr) * (1 - n), vn)*dx - inner(betaN(uCurr) * n, vn)*dx
Fh = inner(Dt(h), vh)*dx - inner(alphaH(uCurr) * (1 - h), vh)*dx - inner(betaH(uCurr) * h, vh)*dx

F = Fu + Fm + Fn + Fh

stepper = TimeStepper(F, butcher_tableau, t, dt, uu,
                      solver_parameters=params)

#while (float(t)< tFinal):
for j in range(150):
    print("Time step", j)
    #if (float(t) + float(dt) > tFinal):
        #dt.assign(tFinal - float(t))
    stepper.advance()
    t.assign(float(t) + float(dt)) 

uFinal, mFinal, nFinal, hFinal = uu.split()
outfile = File("HHmonodomain_2d.pvd")
outfile.write(uFinal)

"""


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