/*
============================================================================================
   This is the main "gateway" routine that interfaces MATLAB with a custom
   FEM Degree-of-Freedom (DoF) allocator for a particular finite element.

   NOTE: this code is part of the FELICITY package:
   Finite ELement Implementation and Computational Interface Tool for You

   OUTPUTS
   -------
   Matrix where each row is the local-to-global DoF map for an element in the mesh.
   (see code below for more info)

   NOTE: portions of this code are automatically generated!

   WARNING!: Make sure all inputs to this mex function are NOT sparse matrices!!!
             Only FULL matrices are allowed as inputs!!!

   Copyright (c) 10-16-2016,  Shawn W. Walker
============================================================================================
*/

// include any libraries you need here
#include <algorithm>
#include <mex.h> // <-- This one is required

/*------------ BEGIN: Auto Generate ------------*/
// define  input indices
#define PRHS_Tri_List                                                         0

// define output indices
#define PLHS_nedelec_1stkind_deg1_dim2_DoF_Map                                0
/*------------   END: Auto Generate ------------*/

/* include classes and other sub-routines */
#include "Misc_Files.h"
#include "Triangle_Edge_Search.cc"
#include "nedelec_1stkind_deg1_dim2_DoF_Allocator.cc"


// note: 'prhs' represents the Right-Hand-Side arguments from MATLAB (inputs)
//       'plhs' represents the  Left-Hand-Side arguments from MATLAB (outputs)


/***************************************************************************************/
// define the "gateway" function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*------------ BEGIN: Auto Generate ------------*/
    /* BEGIN: Error Checking */
    if (nrhs!=1)
        {
        printf("ERROR: 1 inputs required!\n");
        printf("\n");
        printf("      INPUTS                                                         ORDER \n");
        printf("      -------------------                                            ----- \n");
        printf("      Tri_List                                                         0 \n");
        printf("\n");
        printf("      OUTPUTS                                                        ORDER \n");
        printf("      -------------------                                            ----- \n");
        printf("      nedelec_1stkind_deg1_dim2_DoF_Map                                0 \n");
        printf("\n");
        mexErrMsgTxt("Check the number of input arguments!");
        }
    if ((nlhs < 1)||(nlhs > 1)) mexErrMsgTxt("1~1 outputs are needed!");
    /* END: Error Checking */
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // declare triangle edge searcher object
    TRIANGLE_EDGE_SEARCH    Tri_Edge_Search;
    Tri_Edge_Search.Setup_Mesh(prhs[PRHS_Tri_List]);

    // declare DoF allocator
    nedelec_1stkind_deg1_dim2_DoF_Allocator      nedelec_1stkind_deg1_dim2_DoF_Obj;
    // do some initial allocation
    plhs[PLHS_nedelec_1stkind_deg1_dim2_DoF_Map] = nedelec_1stkind_deg1_dim2_DoF_Obj.Init_DoF_Map(Tri_Edge_Search.Num_Tri);
    // create the DoFmap
    nedelec_1stkind_deg1_dim2_DoF_Obj.Fill_DoF_Map(&Tri_Edge_Search);

    /*------------   END: Auto Generate ------------*/

    // output the Euler characteristic of the domain
    const int CHI = Tri_Edge_Search.Num_Unique_Vertices - Tri_Edge_Search.Num_Unique_Edges + Tri_Edge_Search.Num_Tri;
    const bool INVALID = (Tri_Edge_Search.Num_Unique_Vertices==0) || (Tri_Edge_Search.Num_Unique_Edges==0) || (Tri_Edge_Search.Num_Tri==0);
    if (!INVALID)
        {
        mexPrintf("Euler Characteristic of the Domain: CHI = V - E + F\n");
        mexPrintf("                                     %d = %d - %d + %d \n",CHI,Tri_Edge_Search.Num_Unique_Vertices,Tri_Edge_Search.Num_Unique_Edges,Tri_Edge_Search.Num_Tri);
        }
}

/***/
