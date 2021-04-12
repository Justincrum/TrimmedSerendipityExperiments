#
#file:  Mixed_Poisson_2d.py
#author:  Justin Crum
#date: 3/19/21
#
"""
Copyright 2021 Justin Crum

Permission is hereby granted, free of charge, to any person obtaining 
a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation 
the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the Software 
is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR 
THE USE OR OTHER DEALINGS IN THE SOFTWARE.


"""
from firedrake import *
import argparse
from firedrake.petsc import PETSc

parser = argparse.ArgumentParser(
        description="Allows for input of order and mesh refinement.")
parser.add_argument("-O", "--Order", 
        type=int, help="Input the order of the polynomials.")
args = parser.parse_args()

for n in range(args.Order, args.Order + 1):

        ###Mesh set up.
        polyDegree = n
        mesh = Mesh('trapezoids.msh')

        ###Function space set up.
        #Could use SminusDiv (hdiv) + DPC (L^2) or
        #RTCF (hdiv) + DQ (L^2).
        hDivSpace = FunctionSpace(mesh, "RTCF", polyDegree)
        l2Space = FunctionSpace(mesh, "DQ", polyDegree - 1)
        mixedSpace = hDivSpace * l2Space
        dofs = mixedSpace.dim()

        ###Problem set up.
        sigma, u = TrialFunctions(mixedSpace)
        tau, v = TestFunctions(mixedSpace)

        x, y = SpatialCoordinate(mesh)
        uex = sin(pi*x)*sin(pi*y)
        sigmaex = grad(uex)
        f = -div(grad(uex))

        a = (dot(sigma, tau) + div(tau)*u + div(sigma)*v)*dx
        l = -f*v*dx
        w = Function(mixedSpace)

        ###Solver parameters and solving the problem.
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
        PETSc.Log.begin()
        with PETSc.Log.Event("Solve"):
            solve(a == l, w, solver_parameters=params)
        time = PETSc.Log.Event("Solve").getPerfInfo()["time"] 
        sigma, u = w.split()

        ###Data collection and printing.
        errVal = norms.errornorm(uex, u)
        sigErrVal = norms.errornorm(sigmaex, sigma)
        info = [polyDegree, dofs, errVal, sigErrVal, time]
        PETSc.Sys.Print(info)
