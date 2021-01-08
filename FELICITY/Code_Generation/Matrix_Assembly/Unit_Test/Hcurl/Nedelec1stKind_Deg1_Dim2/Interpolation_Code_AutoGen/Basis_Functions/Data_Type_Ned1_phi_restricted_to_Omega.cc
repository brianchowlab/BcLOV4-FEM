/*
============================================================================================
   This class contains data about a given FE function space, and methods for computing
   transformations of the local basis functions.

   This code references the header files:

   matrix_vector_defn.h
   matrix_vector_ops.h
   geometric_computations.h
   basis_function_computations.h


   NOTE: portions of this code are automatically generated!

   Copyright (c) 01-15-2018,  Shawn W. Walker
============================================================================================
*/

/*------------ BEGIN: Auto Generate ------------*/
// define the name of the FE basis function (should be the same as the filename of this file)
#define SpecificFUNC        Data_Type_Ned1_phi_restricted_to_Omega
#define SpecificFUNC_str   "Data_Type_Ned1_phi_restricted_to_Omega"

// set the type of function space
#define SPACE_type  "CG - nedelec_1stkind_deg1_dim2"
// set the name of function space
#define SPACE_name  "Ned1"

// set the Subdomain topological dimension
#define SUB_TD  2
// set the Domain of Integration (DoI) topological dimension
#define DOI_TD  2
// set the geometric dimension
#define GD  2
// set the number of cartesian tuple components (m*n) = 1 * 1
#define NC  1
// NOTE: the (i,j) tuple component is accessed by the linear index k = i + (j-1)*m
// set the number of quad points
#define NQ  1
// set the number of basis functions
#define NB  3
/*------------   END: Auto Generate ------------*/

/* C++ (Specific) FE Function class definition */
class SpecificFUNC: public ABSTRACT_FEM_Function_Class // derive from base class
{
public:
    int*     Elem_DoF[NB];    // element DoF list

    // data structure containing information on the function evaluations.
    // Note: this data is evaluated at several quadrature points!
    // vector valued H(curl) basis functions
    VEC_2x1 Func_vv_Value[NB][NQ];

    // constructor
    SpecificFUNC ();
    ~SpecificFUNC (); // destructor
    void Setup_Function_Space(const mxArray*);
    void Get_Local_to_Global_DoFmap(const int&, int*) const;
                   // need the "const" to ENSURE that nothing in this object will change!
    void Transform_Basis_Functions();
    const CLASS_geom_Omega_embedded_in_Omega_restricted_to_Omega*  Mesh;

private:
    void Get_Basis_Sign_Change(const int Mesh_Vertex[SUB_TD+1], double Basis_Sign[NB]) const;
    void Map_Basis_p1();
};

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* constructor */
SpecificFUNC::SpecificFUNC () :
ABSTRACT_FEM_Function_Class () // call the base class constructor
{
    Name       = (char*) SpecificFUNC_str;
    Type       = (char*) SPACE_type;
    Space_Name = (char*) SPACE_name;
    Sub_TopDim = SUB_TD;
    DoI_TopDim = DOI_TD;
    GeoDim     = GD;
    Num_Basis  = NB;
    Num_Comp   = NC;
    Num_QP     = NQ;
    Mesh       = NULL;

    // init DoF information to NULL
    for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
        Elem_DoF[basis_i] = NULL;

    // init everything to zero
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
    for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
        Func_vv_Value[basis_i][qp_i].Set_To_Zero();
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/


/***************************************************************************************/
/* DE-structor */
SpecificFUNC::~SpecificFUNC ()
{
}
/***************************************************************************************/


/***************************************************************************************/
/* put incoming function data from MATLAB into a nice struct  */
void SpecificFUNC::Setup_Function_Space(const mxArray* Elem)          // inputs
{
    Init_Function_Space(Elem);

    // split up the columns of the element data
    Elem_DoF[0] = (int *) mxGetPr(Elem);
    for (int basis_i = 1; (basis_i < Num_Basis); basis_i++)
        Elem_DoF[basis_i] = Elem_DoF[basis_i-1] + Num_Elem;
}
/***************************************************************************************/


/***************************************************************************************/
/* get the local DoFs on the given element.
   Note: elem_index is in the   C-style (i.e. 0 <= elem_index <= Num_Elem - 1),
         Indices is in the MATLAB-style (i.e. 1 <= Indices[:] <= max(Elem_DoF)). */
void SpecificFUNC::Get_Local_to_Global_DoFmap(const int& elem_index, int* Indices) const  // inputs
{
    /* error check: */
    if (elem_index < 0)
        {
        mexPrintf("ERROR: Given element index #%d is not positive. It must be > 0!\n",elem_index+1);
        mexPrintf("ERROR: There is an issue with the Finite Element Space = %s!\n",Space_Name);
        mexErrMsgTxt("ERROR: Make sure your inputs are valid!");
        }
    else if (elem_index >= Num_Elem)
        {
        mexPrintf("ERROR: Given element index #%d exceeds the number of elements in the finite element (FE) space.\n",elem_index+1);
		mexPrintf("It must be <= %d!  OR  Your FE space DoFmap is not defined correctly!\n",Num_Elem);
		mexPrintf("   For example, the number of rows in DoFmap should *equal*\n");
		mexPrintf("       the number of mesh elements in the (sub)-domain.\n");
        mexPrintf("ERROR: There is an issue with the Finite Element Space = %s!\n",Space_Name);
        mexErrMsgTxt("ERROR: Make sure your inputs are valid!");
        }

    // get local to global index map for the current element
    for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
        {
        int DoF_index = Elem_DoF[basis_i][elem_index] - 1; // shifted for C - style indexing
        Indices[basis_i] = DoF_index;
        }
}
/***************************************************************************************/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* get sign changes to apply to basis functions to account for changes in orientation
  (with respect to the "standard" reference orientation) on the current mesh element.  */
void SpecificFUNC::Get_Basis_Sign_Change(const int Mesh_Vertex[SUB_TD+1], double Basis_Sign[NB]) const
{
    double Mesh_Orient[3] = {1.0, 1.0, 1.0}; // into to all positive

    /* edge orientation is chosen in the following way:
       Each edge points from the vertex of lower (global) index to the vertex of higher (global) index;
       On a triangle, a positive orientation means:
       edge #0: Mesh_Vertex[1] < Mesh_Vertex[2]
       edge #1: Mesh_Vertex[0] < Mesh_Vertex[2]
       edge #2: Mesh_Vertex[0] < Mesh_Vertex[1] */

    // flip orientation if global vertex indices do not satisfy the above
    if (Mesh_Vertex[1] > Mesh_Vertex[2]) Mesh_Orient[0] = -1.0;
    if (Mesh_Vertex[0] > Mesh_Vertex[2]) Mesh_Orient[1] = -1.0;
    if (Mesh_Vertex[0] > Mesh_Vertex[1]) Mesh_Orient[2] = -1.0;

    // get edge orientation "signature" (takes values from 0 to 7)
    const int Edge_Orientation_Signature = (int) ( (Mesh_Orient[0] < 0.0) * 1 +
                                                   (Mesh_Orient[1] < 0.0) * 2 +
                                                   (Mesh_Orient[2] < 0.0) * 4 );

    // BEGIN: determine which basis functions must change their signs
    /* 2-D triangle element has 3 edges */
    switch (Edge_Orientation_Signature)
    {
    case 1: // (-1, 1, 1) edge orientation
        Basis_Sign[0] = -1.0;
        break;
    case 2: // (1, -1, 1) edge orientation
        Basis_Sign[1] = -1.0;
        break;
    case 3: // (-1, -1, 1) edge orientation
        Basis_Sign[0] = -1.0;
        Basis_Sign[1] = -1.0;
        break;
    case 4: // (1, 1, -1) edge orientation
        Basis_Sign[2] = -1.0;
        break;
    case 5: // (-1, 1, -1) edge orientation
        Basis_Sign[0] = -1.0;
        Basis_Sign[2] = -1.0;
        break;
    case 6: // (1, -1, -1) edge orientation
        Basis_Sign[1] = -1.0;
        Basis_Sign[2] = -1.0;
        break;
    case 7: // (-1, -1, -1) edge orientation
        Basis_Sign[0] = -1.0;
        Basis_Sign[1] = -1.0;
        Basis_Sign[2] = -1.0;
        break;
    default: ; // Edge_Orientation_Signature==0 (1, 1, 1) edge orientation
        // do nothing; everything is positive already
    }
    // END: determine which basis functions must change their signs
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* compute the correct local transformation */
void SpecificFUNC::Transform_Basis_Functions()
{
    Map_Basis_p1();
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/


/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* map basis functions from the standard reference element
       to an actual element in the Domain.     */
void SpecificFUNC::Map_Basis_p1()
{

/*------------ BEGIN: Auto Generate ------------*/
    /* call sub-routine to determine sign flips for 2-D H(curl) basis functions */
    double Basis_Sign[NB] = {1.0, 1.0, 1.0}; // init to all positive (no sign change)
    // retrieve global vertex indices on the current cell
    int Vtx_Ind[SUB_TD+1];
    Mesh->Get_Current_Cell_Vertex_Indices(Vtx_Ind);
    Get_Basis_Sign_Change(Vtx_Ind, Basis_Sign);
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, nedelec_1stkind_deg1_dim2
// the Subdomain             has topological dimension = 2
// the Domain of Integration has topological dimension = 2
// geometric dimension = 2

    // get "Val" of basis functions
    VEC_2x1 phi_Val[NB];

    // convenience variables
    const double& x_hat = local_coord[0];
    const double& y_hat = local_coord[1];
    const double& z_hat = local_coord[2];

    phi_Val[0].Set_Equal_To(-y_hat, x_hat);
    phi_Val[1].Set_Equal_To(y_hat, -x_hat+1.0);
    phi_Val[2].Set_Equal_To(-y_hat+1.0, x_hat);

/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
    /*** compute basis function quantities ***/
    // map basis vectors over (indexing is in the C style)
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        // evaluate for each basis function
        for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
            {
            VEC_2x1 vv_temp;
            // multiply by sign change
            Scalar_Mult_Vector(phi_Val[basis_i], Basis_Sign[basis_i], vv_temp);
            // multiply by inverse Jacobian matrix transpose
            Mat_Transpose_Vec(Mesh->Map_PHI_Inv_Grad[0], vv_temp, Func_vv_Value[basis_i][qp_i]);
            }
        }
/*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

// remove those macros!
#undef SpecificFUNC
#undef SpecificFUNC_str

#undef SPACE_type
#undef SPACE_name
#undef SUB_TD
#undef DOI_TD
#undef GD
#undef NC
#undef NB
#undef NQ

/***/
