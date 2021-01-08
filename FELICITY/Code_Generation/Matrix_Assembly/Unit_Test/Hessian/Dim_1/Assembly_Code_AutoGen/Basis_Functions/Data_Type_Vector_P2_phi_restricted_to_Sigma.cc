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
#define SpecificFUNC        Data_Type_Vector_P2_phi_restricted_to_Sigma
#define SpecificFUNC_str   "Data_Type_Vector_P2_phi_restricted_to_Sigma"

// set the type of function space
#define SPACE_type  "CG - lagrange_deg2_dim1"
// set the name of function space
#define SPACE_name  "Vector_P2"

// set the Subdomain topological dimension
#define SUB_TD  1
// set the Domain of Integration (DoI) topological dimension
#define DOI_TD  1
// set the geometric dimension
#define GD  1
// set the number of cartesian tuple components (m*n) = 2 * 1
#define NC  2
// NOTE: the (i,j) tuple component is accessed by the linear index k = i + (j-1)*m
// set the number of quad points
#define NQ  4
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
    // arc-length derivative of the function
    SCALAR Func_f_d_ds[NB][NQ];
    // 2nd arc-length derivative of the function
    SCALAR Func_f_d2_ds2[NB][NQ];
    // hessian of basis function
    MAT_1x1 Func_f_Hess[NB][NQ];

    // constructor
    SpecificFUNC ();
    ~SpecificFUNC (); // destructor
    void Setup_Function_Space(const mxArray*);
    void Get_Local_to_Global_DoFmap(const int&, int*) const;
                   // need the "const" to ENSURE that nothing in this object will change!
    void Transform_Basis_Functions();
    const CLASS_geom_Sigma_embedded_in_Sigma_restricted_to_Sigma*  Mesh;

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
        Func_f_d_ds[basis_i][qp_i].Set_To_Zero();
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
    for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
        Func_f_d2_ds2[basis_i][qp_i].Set_To_Zero();
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
    for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
        Func_f_Hess[basis_i][qp_i].Set_To_Zero();
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
// Local Element defined on Subdomain: CG, lagrange_deg2_dim1
// the Subdomain             has topological dimension = 1
// the Domain of Integration has topological dimension = 1
// geometric dimension = 1
// Number of Quadrature Points = 4

    // get "Grad" of basis functions
    VEC_1x1 phi_Grad[NQ][NB];
    // get "Hess" of basis functions
    MAT_1x1 phi_Hess[NQ][NB];
    // get "d_ds" of basis functions
    SCALAR phi_d_ds[NQ][NB];
    // get "d2_ds2" of basis functions
    SCALAR phi_d2_ds2[NQ][NB];

    phi_Grad[0][0].Set_Equal_To(-2.72227262318810492E+00);
    phi_Grad[0][1].Set_Equal_To(-7.22272623188105145E-01);
    phi_Grad[0][2].Set_Equal_To(3.44454524637621029E+00);
    phi_Grad[1][0].Set_Equal_To(-1.67996208716971274E+00);
    phi_Grad[1][1].Set_Equal_To(3.20037912830287263E-01);
    phi_Grad[1][2].Set_Equal_To(1.35992417433942547E+00);
    phi_Grad[2][0].Set_Equal_To(-3.20037912830287485E-01);
    phi_Grad[2][1].Set_Equal_To(1.67996208716971251E+00);
    phi_Grad[2][2].Set_Equal_To(-1.35992417433942503E+00);
    phi_Grad[3][0].Set_Equal_To(7.22272623188105367E-01);
    phi_Grad[3][1].Set_Equal_To(2.72227262318810537E+00);
    phi_Grad[3][2].Set_Equal_To(-3.44454524637621073E+00);

    phi_Hess[0][0].Set_Equal_To(4.00000000000000000E+00);
    phi_Hess[0][1].Set_Equal_To(4.00000000000000000E+00);
    phi_Hess[0][2].Set_Equal_To(-8.00000000000000000E+00);
    phi_Hess[1][0].Set_Equal_To(4.00000000000000000E+00);
    phi_Hess[1][1].Set_Equal_To(4.00000000000000000E+00);
    phi_Hess[1][2].Set_Equal_To(-8.00000000000000000E+00);
    phi_Hess[2][0].Set_Equal_To(4.00000000000000000E+00);
    phi_Hess[2][1].Set_Equal_To(4.00000000000000000E+00);
    phi_Hess[2][2].Set_Equal_To(-8.00000000000000000E+00);
    phi_Hess[3][0].Set_Equal_To(4.00000000000000000E+00);
    phi_Hess[3][1].Set_Equal_To(4.00000000000000000E+00);
    phi_Hess[3][2].Set_Equal_To(-8.00000000000000000E+00);

    phi_d_ds[0][0].Set_Equal_To(-2.72227262318810492E+00);
    phi_d_ds[0][1].Set_Equal_To(-7.22272623188105145E-01);
    phi_d_ds[0][2].Set_Equal_To(3.44454524637621029E+00);
    phi_d_ds[1][0].Set_Equal_To(-1.67996208716971274E+00);
    phi_d_ds[1][1].Set_Equal_To(3.20037912830287263E-01);
    phi_d_ds[1][2].Set_Equal_To(1.35992417433942547E+00);
    phi_d_ds[2][0].Set_Equal_To(-3.20037912830287485E-01);
    phi_d_ds[2][1].Set_Equal_To(1.67996208716971251E+00);
    phi_d_ds[2][2].Set_Equal_To(-1.35992417433942503E+00);
    phi_d_ds[3][0].Set_Equal_To(7.22272623188105367E-01);
    phi_d_ds[3][1].Set_Equal_To(2.72227262318810537E+00);
    phi_d_ds[3][2].Set_Equal_To(-3.44454524637621073E+00);

    phi_d2_ds2[0][0].Set_Equal_To(4.00000000000000000E+00);
    phi_d2_ds2[0][1].Set_Equal_To(4.00000000000000000E+00);
    phi_d2_ds2[0][2].Set_Equal_To(-8.00000000000000000E+00);
    phi_d2_ds2[1][0].Set_Equal_To(4.00000000000000000E+00);
    phi_d2_ds2[1][1].Set_Equal_To(4.00000000000000000E+00);
    phi_d2_ds2[1][2].Set_Equal_To(-8.00000000000000000E+00);
    phi_d2_ds2[2][0].Set_Equal_To(4.00000000000000000E+00);
    phi_d2_ds2[2][1].Set_Equal_To(4.00000000000000000E+00);
    phi_d2_ds2[2][2].Set_Equal_To(-8.00000000000000000E+00);
    phi_d2_ds2[3][0].Set_Equal_To(4.00000000000000000E+00);
    phi_d2_ds2[3][1].Set_Equal_To(4.00000000000000000E+00);
    phi_d2_ds2[3][2].Set_Equal_To(-8.00000000000000000E+00);


// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     6.94318442029737137E-02, \
//     3.30009478207571816E-01, \
//     6.69990521792428129E-01, \
//     9.30568155797026342E-01  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     1.73927422568727008E-01, \
//     3.26072577431272992E-01, \
//     3.26072577431272992E-01, \
//     1.73927422568726953E-01  \
//     };
/*------------   END: Auto Generate ------------*/
    // set basis function values to the correct mesh entity
    Func_f_Value = &Value_p1; // point to the correct mesh entity

/*------------ BEGIN: Auto Generate ------------*/
    /*** compute basis function quantities ***/
    // compute d/ds of the basis function on the true element in the domain (indexing is in the C style)
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        // multiply by 1 / Det(Jacobian)
        for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
            {
            Func_f_d_ds[basis_i][qp_i].a = phi_d_ds[qp_i][basis_i].a * Mesh->Map_Inv_Det_Jac[qp_i].a;
            }
        }
    // compute d^2/ds^2 of the basis function on the true element in the domain (indexing is in the C style)
    // note: this requires the 1st arc-length derivative of the function
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        // multiply 1 / Det(Jacobian)^2 by the difference of the hessian terms
        for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
            {
            // compute square of inverse jacobian determinant
            const double Inv_Det_Jac_Square = Mesh->Map_Inv_Det_Jac[qp_i].a * Mesh->Map_Inv_Det_Jac[qp_i].a;
            const double Diff_Hess_Terms = phi_d2_ds2[qp_i][basis_i].a - ( Func_f_d_ds[basis_i][qp_i].a * Mesh->Map_PHI_Hess[qp_i].m[0][0][0] );
            Func_f_d2_ds2[basis_i][qp_i].a = Inv_Det_Jac_Square * Diff_Hess_Terms;
            }
        }
    // map hessian from local to global coordinates (indexing is in the C style)
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        // formula involves factors of the inverse of the jacobian (or metric) of PHI
        for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
            {
            const double Inv_Det_Jac_Square = Mesh->Map_Inv_Det_Jac[qp_i].a * Mesh->Map_Inv_Det_Jac[qp_i].a;
            Func_f_Hess[basis_i][qp_i].m[0][0] = phi_Hess[qp_i][basis_i].m[0][0] * Inv_Det_Jac_Square + phi_Grad[qp_i][basis_i].v[0] * Mesh->Map_PHI_Hess_Inv_Map[qp_i].m[0][0][0];
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
// Local Element defined on Subdomain: CG, lagrange_deg2_dim1
// the Subdomain             has topological dimension = 1
// the Domain of Integration has topological dimension = 1
// geometric dimension = 1
// Number of Quadrature Points = 4

    // get "Val" of basis functions
    SCALAR phi_Val[NQ][NB];

    phi_Val[0][0].Set_Equal_To(8.01346029369930934E-01);
    phi_Val[0][1].Set_Equal_To(-5.97902822241216866E-02);
    phi_Val[0][2].Set_Equal_To(2.58444252854190815E-01);
    phi_Val[1][0].Set_Equal_To(2.27784076790952195E-01);
    phi_Val[1][1].Set_Equal_To(-1.12196966793904188E-01);
    phi_Val[1][2].Set_Equal_To(8.84412890002951979E-01);
    phi_Val[2][0].Set_Equal_To(-1.12196966793904160E-01);
    phi_Val[2][1].Set_Equal_To(2.27784076790952084E-01);
    phi_Val[2][2].Set_Equal_To(8.84412890002952090E-01);
    phi_Val[3][0].Set_Equal_To(-5.97902822241216450E-02);
    phi_Val[3][1].Set_Equal_To(8.01346029369931045E-01);
    phi_Val[3][2].Set_Equal_To(2.58444252854190593E-01);


// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     6.94318442029737137E-02, \
//     3.30009478207571816E-01, \
//     6.69990521792428129E-01, \
//     9.30568155797026342E-01  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     1.73927422568727008E-01, \
//     3.26072577431272992E-01, \
//     3.26072577431272992E-01, \
//     1.73927422568726953E-01  \
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
