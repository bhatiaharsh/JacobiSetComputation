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
#include "JacobiSet.h"
#include <iostream>
#include <fstream>

#ifdef USE_SOS
    #include "sos_utils.h"
#endif

// ==========================================================================================
bool JacobiSet::createSoS() {

    size_t nv = mesh->vertices.size();

#ifndef USE_SOS
    return true;
#else
    printf(" -------------- Creating SOS Matrix ........... ");

    static int lia_count = 1;

    if( !sos_is_down() ){
        sos_shutdown();
        lia_count++;
    }

    // initialize SoS just once
    SoS_metadata sm;

    sm.name = NULL;
    sm.title = NULL;
    sm.lines = 0;

    // p_x = (f_x, g_x)
    sm.data_size = nv;          // no of values
    sm.data_dim  = 2;           // dim of points
    //sm.simp_size = H_U.size();     // no of simplices
    //sm.simp_dim  = 1;        // dim of simplices

    sm.fix_w = 15;
    sm.fix_a = 14;

    sm.scale = 1.0;     //ie, no scaling
    sm.decimals = 10;   //ie, int coordinates
    sm.has_weights = 0; //ie, unknown

    // ---
    // test the metadata

    // 0 <= a < w < MAX
    if (!( (0 <= sm.fix_a) && (sm.fix_a < sm.fix_w) && (sm.fix_w <= MAX_DECIMALS) ) )
        printf("\n SoS Metadata Error: invalid #fix=%d.%d\n", sm.fix_w, sm.fix_a );

    if( sm.data_size == 0 || sm.data_dim == 0 )
        printf("\n SoS Metadata Error: invalid data. size = %d, dim = %d\n", sm.data_size, sm.data_dim );

    if (sm.fix_w != 0){

        print ("(#fix=%d.%d) ", sm.fix_w, sm.fix_a);
        sm.decimals = sm.fix_w;
        sm.scale = sm.scale / powf(10, (double) sm.fix_a);
        //sm.scale = sm.scale / exp10 ((double) sm.fix_a);
        printf("scale = %.15f\n", sm.scale);
    }
    // ---

    //printf("initialize now!\n");

    sos_matrix (sm.data_size, sm.data_dim, sm.scale, Lia_DIGITS (5 * sm.decimals + 3),  Lia_DIGITS (2 * sm.decimals + 1));

    //printf(" lia_limit now!!\n");
    lia_stack_limit( lia_count*((sm.data_size + 1) * sm.data_dim));

    // --------------------
    for( int v = 0; v < nv; v++ ){

        //Printf(" loading element %d\n", v);
        SoSUtils::ffp_param_push2 (v+1, 1, SoSUtils::float_to_fixed(f->at(v), sm.fix_a), sm.fix_w, sm.fix_a);
        SoSUtils::ffp_param_push2 (v+1, 2, SoSUtils::float_to_fixed(g->at(v), sm.fix_a), sm.fix_w, sm.fix_a);
    }
    return true;
#endif
}


// ==========================================================================================
// simple utility functions


// compare two floating point valies
bool JacobiSet::are_equal(double a, double b, int precision) {

    double error = pow(10.0, double(-1.0*precision));
    return (fabs(a-b) < error);
}

// Sort three values, and return the number of swaps needed for sorting
int JacobiSet::sort(int &a, int &b, int &c) {

#ifdef USE_SOS
    return basic_isort3 (&a, &b, &c);
#else
    int swaps = 0;
    if (a > b){
        std::swap (a, b);
        swaps++;
    }
    if (b > c){
        std::swap (b, c);
        swaps++;
        if (a > b){
            std::swap (a, b);
            swaps++;
        }
    }
    return (swaps);
#endif
}

// ------------------------------------------------------------------------------------------
// Compute the alignment of an edge in Jacobi set
bool JacobiSet::alignment(int v1, int v2) const {

#ifdef USE_SOS
    return ( sos_smaller(v1+1, 1, v2+1, 1) == sos_smaller(v1+1, 2, v2+1, 2) );
#else
    return ( (f->at(v2) > f->at(v1)) == (g->at(v2) > g->at(v1) ) );
#endif
}

// ------------------------------------------------------------------------------------------
// POS1 function -- explained on Pg 111 in Vijay's thesis.
bool JacobiSet::POS1(int a, int b, int v) const{

#ifdef USE_SOS
   return sos_lambda3(a+1,b+1,v+1)->signum < 0;
#else
    // Algo in Vijay's Thesis. Pg 111
    double X = (f->at(a) * (g->at(v)-g->at(b))) +
               (f->at(b) * (g->at(a)-g->at(v))) +
               (f->at(v) * (g->at(b)-g->at(a)));

    if( !are_equal(X, 0.0f) )              return (X > 0);
    if( !are_equal(g->at(v), g->at(b)) )   return (g->at(v) > g->at(b));
    if( !are_equal(g->at(a), g->at(v)) )   return (g->at(a) > g->at(v));
    if( !are_equal(f->at(b), f->at(v)) )   return (f->at(b) > f->at(v));
    if( !are_equal(f->at(v), f->at(a)) )   return (f->at(v) > f->at(a));

    printf(" --- JacobiSet::POS1 failed!\n");
    return true;
#endif
}

// ------------------------------------------------------------------------------------------
// POS2 function -- explained on Pg 111 in Vijay's thesis.
bool JacobiSet::POS2(int a, int b, const std::vector<double> &F) const {

#ifdef USE_SOS
    int sos_index = (&F == f) ? 1 : (&F == g) ? 2 : 0;
    if(sos_index == 0)
        throw std::invalid_argument("JacobiSet::is_smaller -- invalid scalar field!");

    return !sos_smaller(a+1, sos_index, b+1, sos_index);
#else
    return are_equal(F[a], F[b]) ? (a < b) : (F[a] < F[b]);
#endif
}

// ==========================================================================================
// Functions for computing Jacobi Sets.

/* ------------------------------------------------------------------------------------------
    Test if vertex 'k' lies in the lower link of edge (i,j)
    Implementation of Pg 111 of Vijay's thesis.
    POS1 test is same for all inputs of the triangle (ijk), but POS2 is not
*/
bool JacobiSet::is_lowerLink(int i, int j, int k, const std::vector<double> *F) const{

    int a = i;      int b = j;      int v = k;
    sort(a, b, v);

    int index_i = (i==a) ? 0 : ((i==b) ? 1 : 2);
    int index_j = (j==a) ? 0 : ((j==b) ? 1 : 2);

    bool is_POS1 = POS1(a,b,v);
    bool is_POS2 = false;

    if(index_j == (index_i+1)%3)	is_POS2 = POS2(i, j, *F);
    else							is_POS2 = POS2(j, i, *F);

    return (is_POS1 == is_POS2);
}

// ------------------------------------------------------------------------------------------
// Compute the Jacobi set for checking all edges of the mesh for criticality
void JacobiSet::compute() {

    if( mesh->vertices.empty() || mesh->edges.empty() ){
        printf(" Invalid mesh. %d vertices and %d edges!\n", mesh->vertices.size(), mesh->edges.size());
        return;
    }

    if( f->size() != g->size() || f->size() != mesh->vertices.size() ) {
        printf(" Size mismatch. %d vertices, but functions contain %d and %d values!\n", mesh->vertices.size(), f->size(), g->size());
        return;
    }

    printf("Computing Jacobi set...");
    fflush(stdout);

    // check all edges in the mesh
    for(std::set<trimesh::Edge>::const_iterator iter = mesh->edges.begin(); iter != mesh->edges.end(); iter++){

        // the vertices of this edge
        int e1 = iter->first;       int e2 = iter->second;

        // the vertices in the link
        std::pair<int,int> lnk = mesh->get_e_link(*iter);
        int v1 = lnk.first;          int v2 = lnk.second;

        // if this is a boundary edge, skip
        if(v1 == -1 || v2 == -1)
            continue;

        // for all interior edges, check the criticality wrt f and g
        bool f_lower_v1 = is_lowerLink(e1, e2, v1, g);
        bool f_lower_v2 = is_lowerLink(e1, e2, v2, g);
        bool g_lower_v1 = is_lowerLink(e1, e2, v1, f);
        bool g_lower_v2 = is_lowerLink(e1, e2, v2, f);

        bool f_crit = (f_lower_v1 == f_lower_v2);
        bool g_crit = (g_lower_v1 == g_lower_v2);

        if(f_crit != g_crit){
            std::cout<<" Warning: JacobiSet::compute() -- mismatch: "<<e1<<", "<<e2<<" :: "<<f_lower_v1<<" - "<<f_lower_v2<<" :: "<<g_lower_v1<<", "<<g_lower_v2<< std::endl;
        }

        if(f_crit){
            jacobiedges.push_back(JacobiEdge(e1, e2, !f_lower_v1, g_lower_v1, alignment(e1, e2)));
        }
    }
    printf(" Done! Found %d edges in the Jacobi set!\n", jacobiedges.size());
}

// ------------------------------------------------------------------------------------------
// write the Jacobi set
void JacobiSet::write(const std::string &filename) const {

    std::ofstream ofile (filename.c_str());
    if(!ofile.is_open()) {
        std::cout << " Could not open file " << filename << " to writing!\n";
        return;
    }

    printf("Writing Jacobi set to file %s...", filename.c_str());
    fflush(stdout);

    ofile << "JacobiSet\n";
    ofile << jacobiedges.size() << std::endl;

    for(int i = 0; i < jacobiedges.size(); i++) {
        const JacobiEdge &je = jacobiedges.at(i);
        const trimesh::point &p1 = mesh->vertices[je.fv1];
        const trimesh::point &p2 = mesh->vertices[je.fv2];
        ofile << je.fv1 << " " << p1[0] << " " << p1[1] << " " << p1[2] << " " <<
                 je.fv2 << " " << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
    }
    ofile.close();

    printf(" Done!\n");
}
