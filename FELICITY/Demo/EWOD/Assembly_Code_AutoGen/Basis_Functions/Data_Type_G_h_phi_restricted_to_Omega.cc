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
#define SpecificFUNC        Data_Type_G_h_phi_restricted_to_Omega
#define SpecificFUNC_str   "Data_Type_G_h_phi_restricted_to_Omega"

// set the type of function space
#define SPACE_type  "CG - lagrange_deg1_dim2"
// set the name of function space
#define SPACE_name  "G_h"

// set the Subdomain topological dimension
#define SUB_TD  2
// set the Domain of Integration (DoI) topological dimension
#define DOI_TD  2
// set the geometric dimension
#define GD  2
// set the number of cartesian tuple components (m*n) = 2 * 1
#define NC  2
// NOTE: the (i,j) tuple component is accessed by the linear index k = i + (j-1)*m
// set the number of quad points
#define NQ  6
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
    // local function evaluated at a quadrature point in reference element
    // (this is a pointer because it will change depending on the local mesh entity)
    SCALAR  (*Func_f_Value)[NB][NQ];
    // intrinsic gradient of basis function
    VEC_2x1 Func_f_Grad[NB][NQ];

    // constructor
    SpecificFUNC ();
    ~SpecificFUNC (); // destructor
    void Setup_Function_Space(const mxArray*);
    void Get_Local_to_Global_DoFmap(const int&, int*) const;
                   // need the "const" to ENSURE that nothing in this object will change!
    void Transform_Basis_Functions();
    const CLASS_geom_Omega_embedded_in_Omega_restricted_to_Omega*  Mesh;

private:
    void Map_Basis_p1();
    SCALAR Value_p1[NB][NQ];
    void Basis_Value_p1(SCALAR V[NB][NQ]);
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
        Func_f_Grad[basis_i][qp_i].Set_To_Zero();
    // init basis function values on local mesh entities
    Basis_Value_p1(Value_p1);
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
// Local Element defined on Subdomain: CG, lagrange_deg1_dim2
// the Subdomain             has topological dimension = 2
// the Domain of Integration has topological dimension = 2
// geometric dimension = 2
// Number of Quadrature Points = 6

    // get "Grad" of basis functions
    VEC_2x1 phi_Grad[NQ][NB];

    phi_Grad[0][0].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00);
    phi_Grad[0][1].Set_Equal_To(1.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Grad[0][2].Set_Equal_To(0.00000000000000000E+00, 1.00000000000000000E+00);
    phi_Grad[1][0].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00);
    phi_Grad[1][1].Set_Equal_To(1.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Grad[1][2].Set_Equal_To(0.00000000000000000E+00, 1.00000000000000000E+00);
    phi_Grad[2][0].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00);
    phi_Grad[2][1].Set_Equal_To(1.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Grad[2][2].Set_Equal_To(0.00000000000000000E+00, 1.00000000000000000E+00);
    phi_Grad[3][0].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00);
    phi_Grad[3][1].Set_Equal_To(1.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Grad[3][2].Set_Equal_To(0.00000000000000000E+00, 1.00000000000000000E+00);
    phi_Grad[4][0].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00);
    phi_Grad[4][1].Set_Equal_To(1.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Grad[4][2].Set_Equal_To(0.00000000000000000E+00, 1.00000000000000000E+00);
    phi_Grad[5][0].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00);
    phi_Grad[5][1].Set_Equal_To(1.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Grad[5][2].Set_Equal_To(0.00000000000000000E+00, 1.00000000000000000E+00);


// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.16847572980458958E-01, 9.15762135097710067E-02}, \
//     {9.15762135097710067E-02, 8.16847572980458958E-01}, \
//     {9.15762135097710067E-02, 9.15762135097710067E-02}, \
//     {1.08103018168069998E-01, 4.45948490915965001E-01}, \
//     {4.45948490915965001E-01, 1.08103018168069998E-01}, \
//     {4.45948490915965001E-01, 4.45948490915965001E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     5.49758718276609978E-02, \
//     5.49758718276609978E-02, \
//     5.49758718276609978E-02, \
//     1.11690794839005500E-01, \
//     1.11690794839005500E-01, \
//     1.11690794839005500E-01  \
//     };
/*------------   END: Auto Generate ------------*/
    // set basis function values to the correct mesh entity
    Func_f_Value = &Value_p1; // point to the correct mesh entity

/*------------ BEGIN: Auto Generate ------------*/
    /*** compute basis function quantities ***/
    // map gradient from local to global coordinates (indexing is in the C style)
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        // multiply local gradient by transpose of Inv_Grad(PHI) matrix
        for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
            {
            Mat_Transpose_Vec(Mesh->Map_PHI_Inv_Grad[0], phi_Grad[qp_i][basis_i], Func_f_Grad[basis_i][qp_i]);
            }
        }
/*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* evaluate basis functions (no derivatives!) on the local reference element. */
void SpecificFUNC::Basis_Value_p1(SCALAR BF_V[NB][NQ])
{

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, lagrange_deg1_dim2
// the Subdomain             has topological dimension = 2
// the Domain of Integration has topological dimension = 2
// geometric dimension = 2
// Number of Quadrature Points = 6

    // get "Val" of basis functions
    SCALAR phi_Val[NQ][NB];

    phi_Val[0][0].Set_Equal_To(9.15762135097700353E-02);
    phi_Val[0][1].Set_Equal_To(8.16847572980458958E-01);
    phi_Val[0][2].Set_Equal_To(9.15762135097710067E-02);
    phi_Val[1][0].Set_Equal_To(9.15762135097700353E-02);
    phi_Val[1][1].Set_Equal_To(9.15762135097710067E-02);
    phi_Val[1][2].Set_Equal_To(8.16847572980458958E-01);
    phi_Val[2][0].Set_Equal_To(8.16847572980457959E-01);
    phi_Val[2][1].Set_Equal_To(9.15762135097710067E-02);
    phi_Val[2][2].Set_Equal_To(9.15762135097710067E-02);
    phi_Val[3][0].Set_Equal_To(4.45948490915965001E-01);
    phi_Val[3][1].Set_Equal_To(1.08103018168069998E-01);
    phi_Val[3][2].Set_Equal_To(4.45948490915965001E-01);
    phi_Val[4][0].Set_Equal_To(4.45948490915965001E-01);
    phi_Val[4][1].Set_Equal_To(4.45948490915965001E-01);
    phi_Val[4][2].Set_Equal_To(1.08103018168069998E-01);
    phi_Val[5][0].Set_Equal_To(1.08103018168069998E-01);
    phi_Val[5][1].Set_Equal_To(4.45948490915965001E-01);
    phi_Val[5][2].Set_Equal_To(4.45948490915965001E-01);


// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.16847572980458958E-01, 9.15762135097710067E-02}, \
//     {9.15762135097710067E-02, 8.16847572980458958E-01}, \
//     {9.15762135097710067E-02, 9.15762135097710067E-02}, \
//     {1.08103018168069998E-01, 4.45948490915965001E-01}, \
//     {4.45948490915965001E-01, 1.08103018168069998E-01}, \
//     {4.45948490915965001E-01, 4.45948490915965001E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     5.49758718276609978E-02, \
//     5.49758718276609978E-02, \
//     5.49758718276609978E-02, \
//     1.11690794839005500E-01, \
//     1.11690794839005500E-01, \
//     1.11690794839005500E-01  \
//     };
/*------------   END: Auto Generate ------------*/
    // copy function evaluations over (indexing is in the C style)
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        // evaluate for each basis function
        for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
            {
            BF_V[basis_i][qp_i].a = phi_Val[qp_i][basis_i].a;
            }
        }
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
