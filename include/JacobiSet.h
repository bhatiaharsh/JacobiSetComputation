/**
Copyright (c) 2016, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory
Written by Harsh Bhatia (bhatia4@llnl.gov).
CODE-701045.
All rights reserved.

This file is part of JacobiSetComputation v1.0.

For details, see https://github.com/bhatiaharsh/JacobiSetComputation.
For more details on the Licence, please read LICENCE file.
*/

#ifndef _JACOBISET_H_
#define _JACOBISET_H_

#include <vector>
#include "TriMeshJ.h"

//! Store properties of an edge in Jacobi set
struct JacobiEdge {

    int fv1, fv2;       // the function vertices
    bool min_f, min_g;  // a restricted min of f or g
    bool pos_align;     // alignment

    JacobiEdge(){}
    JacobiEdge(int v1, int v2, bool minf, bool ming, bool pa){

        fv1 = v1;   fv2 = v2;
        pos_align = pa;
        min_f = minf;      min_g = ming;
    }
};


//! Compute and store the Jacobi set of 2 functions
class JacobiSet {

private:
    const trimesh::TriMeshJ *mesh;
    const std::vector<double> *f;
    const std::vector<double> *g;

    std::vector<JacobiEdge> jacobiedges;

public:
    JacobiSet(const trimesh::TriMeshJ *mesh_, const std::vector<double> *f_, const std::vector<double> *g_) :
        mesh(mesh_), f(f_), g(g_) {

#ifdef USE_SOS
        createSoS();
#endif
    }

    bool createSoS();

    //! sort three values
    static int sort(int &a, int &b, int &c);

    //! compare two floating point values
    static bool are_equal(double a, double b, int precision = 12);

    //! POS1 and POS2 functions from Vijay's thesis Pg 111
    bool POS1(int a, int b, int v) const;
    bool POS2(int a, int b, const std::vector<double> &F) const ;

    bool is_lowerLink(int i, int j, int k, const std::vector<double> *F) const ;

    //! compute alignment of a Jacobi edge
    bool alignment(int v1, int v2) const;

    //! compute Jacobi set
    void compute();

    //! write the Jacobi set
    void write(const std::string &filename) const;
};


#endif
