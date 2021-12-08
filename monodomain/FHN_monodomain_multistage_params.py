# Solving monodomain equation in Firedrake
# Justin Crum
# 4/16/21

# To do:
# 1.  Implement iIon as in Ten Tusscher
# 3.  Need to set up initial conditions including what the activiation in
#    a corner will look like.

from firedrake import (And, Constant, File, Function, FunctionSpace,
                       RectangleMesh, SpatialCoordinate, TestFunctions,
                       as_matrix, conditional, dx, grad, inner, interpolate,
                       split)
from firedrake.petsc import PETSc
from irksome import Dt, RadauIIA, TimeStepper
from mpi4py import MPI


PETSc.Log().begin()
# def get_time(event, comm=COMM_WORLD):
#     return comm.allreduce(PETSc.Log.Event(event).getPerfInfo()["time"], op=MPI.SUM) / comm.size


# mesh = Mesh('simple_rectangle_mesh.msh')
mesh = RectangleMesh(100, 100, 70, 70, quadrilateral=True)
# mesh = ExtrudedMesh(msh, layers=3)#, layer_height=1)
polyOrder = 3

# Set up the function space and test/trial functions.
V = FunctionSpace(mesh, "Q", polyOrder)
Z = V * V

x, y = SpatialCoordinate(mesh)
dt = Constant(0.1)
t = Constant(0.0)

InitialPotential = conditional(x < 3.5, Constant(2.0), Constant(-1.28791))
InitialCell = conditional(And(And(31 <= x, x < 39), And(0 <= y, y < 35)),
                          Constant(2.0), Constant(-0.5758))

eps = Constant(0.1)
beta = Constant(1.0)
gamma = Constant(0.5)

chi = Constant(1.0)
capacitance = Constant(1.0)

sigma1 = 1.0
sigma2 = 1.0

sigma = as_matrix([[sigma1, 0.0], [0.0, sigma2]])


# Set up Irksome to be used.
butcher_tableau = RadauIIA(2)
ns = butcher_tableau.num_stages

uu = Function(Z)
vu, vc = TestFunctions(Z)
uu.sub(0).interpolate(InitialPotential)
uu.sub(1).interpolate(InitialCell)

(uCurr, cCurr) = split(uu)

Fu = inner(chi * capacitance * Dt(uCurr), vu)*dx \
    + inner(grad(uCurr), sigma * grad(vu))*dx \
    + inner((chi/eps) * (-uCurr + (uCurr**3 / 3) + cCurr), vu)*dx
Fc = inner(Dt(cCurr), vc)*dx - inner(eps * uCurr, vc)*dx \
    - inner(beta * eps, vc)*dx + inner(gamma * eps * cCurr, vc)*dx
F = Fu + Fc


params = {"snes_atol": 1e-10,
          "snes_monitor": None,
          "mat_type": "matfree",
          "pmat_type": "matfree",
          "ksp_type": "fgmres",  # gmres or cg preonly
          "ksp_monitor": None,
          "pc_type": "python",
          "pc_python_type": "irksome.RanaLD",
          "aux": {
              "mat_type": "matfree",
              "pc_type": "fieldsplit",
              "pc_fieldsplit_type": "multiplicative"}}

per_stage = {
    "ksp_type": "preonly",
    "pc_type": "fieldsplit",
    "pc_fieldsplit_type": "additive",
    "fieldsplit_0": {
        "ksp_type": "preonly",
        "pc_type": "python",
        "pc_python_type": "firedrake.AssembledPC",
        "assembled_pc_type": "gamg"
    },
    "fieldsplit_1": {
        "ksp_type": "preonly",
        "pc_type": "python",
        "pc_python_type": "firedrake.AssembledPC",
        "assembled_pc_type": "gamg"
    }}

for s in range(ns):
    params["aux"][f"pc_fieldsplit_{s}_fields"] = f"{2*s},{2*s+1}"
    params["aux"][f"fieldsplit_{s}"] = per_stage


# params = {"snes_atol": 1e-10,
#           "snes_monitor": None,
#           "mat_type": "aij",
#           "ksp_type": "fgmres",  # gmres or cg preonly
#           "ksp_monitor": None,
#           "pc_type": "python",
#           "pc_python_type": "irksome.RanaLD",
#           "aux": {
#               "pc_type": "fieldsplit",
#               "pc_fieldsplit_type": "multiplicative"}}

# per_stage = {
#     "ksp_type": "preonly",
#     "pc_type": "fieldsplit",
#     "pc_fieldsplit_type": "additive",
#     "fieldsplit_0": {
#         "ksp_type": "preonly",
#         "pc_type": "gamg",
#     },
#     "fieldsplit_1": {
#         "ksp_type": "preonly",
#         "pc_type": "icc",
#     }}

# for s in range(ns):
#     params["aux"][f"pc_fieldsplit_{s}_fields"] = f"{2*s},{2*s+1}"
#     params["aux"][f"fieldsplit_{s}"] = per_stage


# stepper = TimeStepper(F, butcher_tableau, t, dt, uu,
#                       solver_parameters=params)

uFinal, cFinal = uu.split()
outfile1 = File("FHN_results/FHN_2d_u.pvd")
outfile2 = File("FHN_results/FHN_2d_c.pvd")
outfile1.write(uFinal, time=0)
outfile2.write(cFinal, time=0)

with PETSc.Log.Event("MonodomainSolveTime"):
    # def get_time(event, comm=COMM_WORLD):
    #     return comm.allreduce(PETSc.Log.Event(event).getPerfInfo()["time"], op=MPI.SUM) / comm.size

    stepper = TimeStepper(F, butcher_tableau, t, dt, uu,
                      solver_parameters=params)
    
    snes = stepper.solver.snes
    newt = 0
    ksp = 0

    with PETSc.Log.Stage("SolvingStage"):

        for j in range(301):
            stepper.advance()
            t.assign(float(t) + float(dt))
            newt += snes.getIterationNumber()
            ksp += snes.getLinearSolveIterations()

            #snes = get_time("SNESSolve")
            # ksp = get_time("KSPSolve")
            # pcsetup = get_time("PCSetUp")
            # pcapply = get_time("PCApply")
            # jac = get_time("SNESJacobianEval")
            # residual = get_time("SNESFunctionEval")

            # uCurr, cCurr = split(uu)
            if (j % 5 == 0):
                print("Time step", j)
                outfile1.write(uFinal, time=j * float(dt))
                outfile2.write(cFinal, time=j * float(dt))

        snestime = PETSc.Log.Event("SNESSolve").getPerfInfo()["time"]
        ksptime = PETSc.Log.Event("KSPSolve").getPerfInfo()["time"]

MonodomainRunTime = PETSc.Log.Event("MonodomainSolveTime").getPerfInfo()["time"]
print("Average number of newton iterations:", newt / 301)
print("Average number of gmres iterations:", ksp / newt)
print("Total newton iteration time was:", snestime)
print("Total linear solving time was:", ksptime)
print("Total monodomain equation solve time was:", MonodomainRunTime)




