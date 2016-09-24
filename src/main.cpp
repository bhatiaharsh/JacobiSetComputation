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


#include <iostream>
#include "TriMeshJ.h"
#include "JacobiSet.h"

using namespace std;
using namespace trimesh;

// ---------------------------------------------------------------

template <typename T>
T* read_binary(const std::string &filename, size_t &sz) {

    // open the file
    FILE *datafile = fopen(filename.c_str(), "rb");
    if(!datafile){
        printf(" Unable to open data file %s\n", filename.c_str());
        return 0;
    }

    printf(" Reading binary file %s...", filename.c_str());
    fflush(stdout);

    // find the size of the file
    fseek(datafile, 0, SEEK_END);
    size_t fsz = ftell(datafile);
    rewind(datafile);

    sz = fsz / sizeof(T);

    if( fsz % sizeof(T) != 0 ){
        printf("\n\t - Invalid number of values in file %s. Size of file = %'ld, Size of each value = %'lu.\n",
               filename.c_str(), fsz, sizeof(T));
        fclose(datafile);
        return 0;
    }

    // now read the data
    T* data = new T[sz];
    size_t rd_sz = fread(data, sizeof(T), sz, datafile);
    if(rd_sz != sz){
        printf("\n\t - Expected %'ld, but read %'ld values!\n", sz, rd_sz);
    }

    // return
    fclose(datafile);
    printf(" Done! Read %'lu values!\n", sz);
    return data;
}

int main(int argc, char *argv[]){

    if(argc != 2) {
        printf(" Usage: %s <meshfile>\n", argv[0]);
        return 1;
    }

    std::string filename (argv[1]);

    trimesh::TriMeshJ *mesh = new trimesh::TriMeshJ(filename);

    mesh->need_edges();
    mesh->need_neighbors();
    printf("Loaded mesh: %d verts, %d edges, %d faces\n", mesh->vertices.size(), mesh->edges.size(), mesh->faces.size());

    int nv = mesh->vertices.size();

    // create test functions
    std::vector<double> f, g;
    size_t fsz, gsz;
#if 0
    double *af = read_binary<double>("f.raw", fsz);
    double *ag = read_binary<double>("g.raw", gsz);

    if(fsz != gsz || fsz != nv){
        printf(" Invalid data. sizes: mesh %d, f %d, g %d\n", nv, fsz, gsz);
        return 1;
    }

    f.insert(f.end(), af, af+fsz);
    g.insert(g.end(), ag, ag+gsz);

#else

    f.resize(nv);
    g.resize(nv);
    for(int i = 0; i < nv; i++) {

        f[i] = mesh->vertices[i][1];
        g[i] = mesh->vertices[i][0];
    }

    FILE *datafile = fopen("f.raw", "wb");
    fwrite(f.data(), sizeof(double), f.size(), datafile);
    fclose(datafile);

    datafile = fopen("g.raw", "wb");
    fwrite(g.data(), sizeof(double), g.size(), datafile);
    fclose(datafile);

    //exit(1);
#endif
    // now compute Jacobi set
    JacobiSet js(mesh, &f, &g);
    js.compute();

    std::string outfile = filename.substr(0, filename.find_last_of("."));
    outfile = outfile + "_jacobi.txt";
    js.write(outfile);
}
