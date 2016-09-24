'''
Copyright (c) 2016, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory
Written by Harsh Bhatia (bhatia4@llnl.gov).
CODE-701045.
All rights reserved.

This file is part of JacobiSetComputation v1.0.

For details, see https://github.com/bhatiaharsh/JacobiSetComputation.
For more details on the Licence, please read LICENCE file.
'''

import vtk
from sys import argv, exit

if len(argv) != 2:
    print "Usage: %s <input_file>" % argv[0]
    exit(0)

## details of inp
infilename = argv[1]
outfilename = infilename[:-3]+"vtk"

# ---------------------------------
print " Reading " + infilename + "..."
f = open(infilename, 'r')

f.readline() # comment line
ne = int(f.readline())

print " Expecting", ne, "edges"

vpoints = vtk.vtkPoints()
vplines = vtk.vtkCellArray()

for i in range(0, ne):

    #currNoPoints = vpoints.GetNumberOfPoints()
    t = f.readline().split()

    # add two points
    vpoints.InsertNextPoint( float(t[1]), float(t[2]), float(t[3]) )
    vpoints.InsertNextPoint( float(t[5]), float(t[6]), float(t[7]) )

    # add an edge (polyline with 2 points)
    vpline = vtk.vtkPolyLine()
    vpline.GetPointIds().SetNumberOfIds(2)
    vpline.GetPointIds().SetId(0, 2*i)
    vpline.GetPointIds().SetId(1, 2*i+1)

    vplines.InsertNextCell(vpline)

# ---------------------------------
print " Writing " + outfilename + "..."

vpolyData = vtk.vtkPolyData()
vpolyData.SetPoints(vpoints)
vpolyData.SetLines(vplines)

vtkWriter = vtk.vtkPolyDataWriter()
vtkWriter.SetFileName(outfilename)
vtkWriter.SetFileTypeToBinary()
vtkWriter.SetInputData(vpolyData)
vtkWriter.Write()
