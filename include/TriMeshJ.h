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


#ifndef _TRIMESHJ_H_
#define _TRIMESHJ_H_

#include <set>
#include "TriMesh.h"

//! extend class TriMesh to add additional functionalities useful for JS computation
//! in particular, we need a datastructure to store all edges

namespace trimesh {

typedef std::pair<int,int> Edge;                    // (vert,vert) pair

class TriMeshJ : public TriMesh {

public:
    TriMeshJ(const std::string &filename)
        : TriMesh ( *read(filename) )
    {}

    // store a set of edges
    std::set<Edge> edges;

    Edge create_edge(const int &va, const int &vb) const{
        int v1 = (va < vb) ? va : vb;
        int v2 = (va > vb) ? va : vb;

        return Edge(v1,v2);
    }

    void need_edges();
    std::pair<int, int> get_e_link(const Edge &edge) const;
};

}

#endif
