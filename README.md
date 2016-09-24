# Jacobi Set Computation v1.0

Copyright (c) 2016, Lawrence Livermore National Security, LLC.

Produced at the Lawrence Livermore National Laboratory.

Written by Harsh Bhatia (bhatia4@llnl.gov).

CODE-701045.

All rights reserved.

### Purpose & Scope

This software is a comamnd line utility to compute the Jacobi set of two piecewise-linear scalar functions defined on a 2D manifold (triangulated mesh in 3D). The algorithm was proposed by Edelsbrunner [1], and was described in more detail by Natrajan [2]. This software implements the algorithm given by [1, 2] using Simulation of Simplicity [3], which uses symbolic perturbations to resolve degenerate geometric cases.

[1] H. Edelsbrunner, J. Harer. Jacobi sets of multiple Morse functions. F. Cucker, R. DeVore, P. Olver, E. Süli (Eds.), Foundations of Computational Mathematics, Minneapolis 2002, Cambridge University Press (2002), pp. 37–57.

[2] Vijay Natarajan. Topological analysis of scalar functions for scientific data visualization. Ph.D. Thesis, Department of Computer Science, Duke University, 2004.

[3] H. Edelsbrunner, E.P. Mücke, Simulation of simplicity: a technique to cope with degenerate
cases in geometric algorithms. ACM Trans. Graph. 9, (1990), pp. 66–104.

### Dependencies

The TrajectoryExplorer has the following software dependencies:

1. `trimesh2` -- An open-source library to rerpesent and manipulate triangle meshes. <http://gfx.cs.princeton.edu/proj/trimesh2>

2. `Simulation of Simplicity` -- The code for SoS was provided by Ernst Mücke (<http://www.geom.uiuc.edu/software/cglist/GeomDir/>) as a sub-module for a different software. The original SoS code was written in C in the 1990's, and requires some modifications for compilation with modern C++ compilers. A **patch** has been provided to fix the code and build it (`patch_SOS.txt`)

These dependencies can be autmatically downloaded and installed using the `install_deps.sh` script.

### Installation & Execution

The JacobiSetComputation tool can be installed by first building the dependencies (`trimesh2` and `SimulationOfSimplicity`) as static libraries, and then using the cmake system.

```
$ pwd
JacobiSetComputation
$ sh install_deps.sh
$ mkdir build && cd build
$ cmake ../
$ make
```

The program requires a triangle mesh, and computes the Jacobi set of x-coordinate and z-coordinate of the mesh. Modify `src/main.cpp` to add the funcionality to read pre-defined functions (f and g) in a desired format. The Jacobi set is written as a collection of edges in a file named `<inputMesh_jacobi.txt>`.

```
$ pwd
JacobiSetComputation
$ ./build/JacobiSetComputation al.obj 
Reading al.obj... Done.
Computing edges... Done. in 10.597000 msec.
Finding vertex neighbors... Done.
Loaded mesh: 3618 verts, 10728 edges, 7152 faces
 -------------- Creating SOS Matrix ........... (#fix=15.14) scale = 0.000000000000010
SoS: matrix[3618,2] @ 5 Lia digits; lia_length (10); 0.331 Mb.
Computing Jacobi set... Done! Found 1246 edges in the Jacobi set!
Writing Jacobi set to file al_jacobi.txt... Done!
```

A python utility is provided to convert the Jacobi set into vtk polydata, which can be loaded into Paraview or VisIt.

```
$ pwd
JacobiSetComputation
$ python utilities/jstovtk.py al_jacobi.txt
 Reading al_jacobi.txt...
 Expecting 1246 edges
 Writing al_jacobi.vtk...
```