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


#include "TriMeshJ.h"
#include <timestamp.h>

using namespace trimesh;

// ================================================================
// link of an edge is a pair of two vertices.
    // requires neighborhood information
std::pair<int, int> TriMeshJ::get_e_link(const Edge &edge) const {

    int v1 = edge.first;
    int v2 = edge.second;

    const std::vector <int> &n1 = neighbors[v1];
    const std::vector <int> &n2 = neighbors[v2];

    int va = -1;
    int vb = -1;

    for(int i = 0; i < n1.size(); i++){

        std::vector<int>::const_iterator iter = find(n2.begin(), n2.end(), n1[i]);

        if( iter != n2.end() ){

            if( va == -1 ){
                va = n1[i];
            }
            else if( vb == -1 ){
                vb = n1[i];
                break;
            }
        }
    }
    return std::pair<int, int>(va, vb);
}

// ================================================================
// compute edges of the mesh
void TriMeshJ::need_edges(){

    if (!edges.empty())
        return;

    need_faces();
    int nf = faces.size();

    if( nf == 0 ){
        dprintf("Edges can't be computed without faces.\n");
        return;
    }

    timestamp st = now();
    dprintf("Computing edges... ");
    fflush(stdout);

    edges.clear();

    for(int i = 0; i < nf; i++) {
    for(int j = 0; j < 3; j++ ){
        edges.insert( create_edge(faces[i][(j+1)%3], faces[i][(j+2)%3]) );
    }
    }
    double et = now()-st;
    if(et < 1.0)    dprintf("Done. in %f msec.\n",1000.00*et);
    else            dprintf("Done. in %f sec.\n",et);
}
