#Script to process vtu files and produce an activation time field.
import vtk
import argparse
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk
import numpy as np

timesteps = 60
dt = 0.5
number_of_points = 10201
reader = vtk.vtkXMLUnstructuredGridReader()
Potential_Array = np.zeros([number_of_points,0])

#Create array of 10201 by timesteps.
for j in range(timesteps):
    file_name = "FHN_results/FHN_2d_u_" + str(j + 1) + ".vtu"
    reader.SetFileName(file_name)
    reader.Update()
    output = reader.GetOutput()
    potential = output.GetPointData().GetArray("function_3[0]")
    N = potential.GetNumberOfTuples()
    potentialNP = vtk_to_numpy(potential)
    Potential_Array = np.column_stack((Potential_Array, potentialNP))

#Search array and determine activation time.
activation_time = np.zeros([number_of_points,])
curr_point = 0
curr_step = 0
while curr_point < number_of_points:
    if Potential_Array[curr_point,curr_step] >= 0.0:
        activation_time[curr_point] = curr_step * dt
        curr_point += 1
        curr_step = 0
    else:
        curr_step += 1
        if curr_step == timesteps:
            activation_time[curr_point] = 1000
            curr_point += 1
            curr_step = 0

vtu_array = numpy_to_vtk(activation_time)


## Code to write the array out on the mesh

vtu_array.SetName("Activations")

reader.SetFileName("FHN_results/FHN_2d_u_1.vtu")
reader.Update()
baseMesh = reader.GetOutput()


outputMesh = vtk.vtkUnstructuredGrid()
outputMesh.SetPoints(baseMesh.GetPoints())
outputMesh.SetCells(vtk.VTK_QUAD, baseMesh.GetCells())
outputMesh.GetPointData().AddArray(vtu_array)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName("activations.vtu")
writer.SetInputData(outputMesh)
writer.Write()




