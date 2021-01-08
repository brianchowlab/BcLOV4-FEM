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
#define SpecificFUNC        Data_Type_RT1_phi_restricted_to_Omega
#define SpecificFUNC_str   "Data_Type_RT1_phi_restricted_to_Omega"

// set the type of function space
#define SPACE_type  "CG - raviart_thomas_deg1_dim3"
// set the name of function space
#define SPACE_name  "RT1"

// set the Subdomain topological dimension
#define SUB_TD  3
// set the Domain of Integration (DoI) topological dimension
#define DOI_TD  3
// set the geometric dimension
#define GD  3
// set the number of cartesian tuple components (m*n) = 1 * 1
#define NC  1
// NOTE: the (i,j) tuple component is accessed by the linear index k = i + (j-1)*m
// set the number of quad points
#define NQ  10
// set the number of basis functions
#define NB  15
/*------------   END: Auto Generate ------------*/

/* C++ (Specific) FE Function class definition */
class SpecificFUNC: public ABSTRACT_FEM_Function_Class // derive from base class
{
public:
    int*     Elem_DoF[NB];    // element DoF list

    // data structure containing information on the function evaluations.
    // Note: this data is evaluated at several quadrature points!
    // vector valued H(div) basis functions
    VEC_3x1 Func_vv_Value[NB][NQ];
    // divergence of vector valued H(div) basis functions
    SCALAR Func_vv_Div[NB][NQ];

    // constructor
    SpecificFUNC ();
    ~SpecificFUNC (); // destructor
    void Setup_Function_Space(const mxArray*);
    void Get_Local_to_Global_DoFmap(const int&, int*) const;
                   // need the "const" to ENSURE that nothing in this object will change!
    void Transform_Basis_Functions();
    const CLASS_geom_Omega_embedded_in_Omega_restricted_to_Omega*  Mesh;

private:
    void Get_Basis_Sign_Change(const double Mesh_Orient[SUB_TD+1], double Basis_Sign[NB]) const;
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
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
    for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
        Func_vv_Div[basis_i][qp_i].Set_To_Zero();
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
void SpecificFUNC::Get_Basis_Sign_Change(const double Mesh_Orient[SUB_TD+1], double Basis_Sign[NB]) const
{
    // get facet (face) orientation "signature" (takes values from 0 to 15)
    const int Facet_Orientation_Signature = (int) ( (Mesh_Orient[0] < 0.0) * 1 +
                                                    (Mesh_Orient[1] < 0.0) * 2 +
                                                    (Mesh_Orient[2] < 0.0) * 4 +
                                                    (Mesh_Orient[3] < 0.0) * 8 );

    // BEGIN: determine which basis functions must change their signs
    /* 3-D tetrahedron element has 4 facets */
    switch (Facet_Orientation_Signature)
    {
    case 1: // (-1, 1, 1, 1) facet orientation
        Basis_Sign[0] = -1.0;
        Basis_Sign[1] = -1.0;
        Basis_Sign[2] = -1.0;
        break;
    case 2: // (1, -1, 1, 1) facet orientation
        Basis_Sign[3] = -1.0;
        Basis_Sign[4] = -1.0;
        Basis_Sign[5] = -1.0;
        break;
    case 3: // (-1, -1, 1, 1) facet orientation
        Basis_Sign[0] = -1.0;
        Basis_Sign[1] = -1.0;
        Basis_Sign[2] = -1.0;
        Basis_Sign[3] = -1.0;
        Basis_Sign[4] = -1.0;
        Basis_Sign[5] = -1.0;
        break;
    case 4: // (1, 1, -1, 1) facet orientation
        Basis_Sign[6] = -1.0;
        Basis_Sign[7] = -1.0;
        Basis_Sign[8] = -1.0;
        break;
    case 5: // (-1, 1, -1, 1) facet orientation
        Basis_Sign[0] = -1.0;
        Basis_Sign[1] = -1.0;
        Basis_Sign[2] = -1.0;
        Basis_Sign[6] = -1.0;
        Basis_Sign[7] = -1.0;
        Basis_Sign[8] = -1.0;
        break;
    case 6: // (1, -1, -1, 1) facet orientation
        Basis_Sign[3] = -1.0;
        Basis_Sign[4] = -1.0;
        Basis_Sign[5] = -1.0;
        Basis_Sign[6] = -1.0;
        Basis_Sign[7] = -1.0;
        Basis_Sign[8] = -1.0;
        break;
    case 7: // (-1, -1, -1, 1) facet orientation
        Basis_Sign[0] = -1.0;
        Basis_Sign[1] = -1.0;
        Basis_Sign[2] = -1.0;
        Basis_Sign[3] = -1.0;
        Basis_Sign[4] = -1.0;
        Basis_Sign[5] = -1.0;
        Basis_Sign[6] = -1.0;
        Basis_Sign[7] = -1.0;
        Basis_Sign[8] = -1.0;
        break;
    case 8: // (1, 1, 1, -1) facet orientation
        Basis_Sign[9] = -1.0;
        Basis_Sign[10] = -1.0;
        Basis_Sign[11] = -1.0;
        break;
    case 9: // (-1, 1, 1, -1) facet orientation
        Basis_Sign[0] = -1.0;
        Basis_Sign[1] = -1.0;
        Basis_Sign[2] = -1.0;
        Basis_Sign[9] = -1.0;
        Basis_Sign[10] = -1.0;
        Basis_Sign[11] = -1.0;
        break;
    case 10: // (1, -1, 1, -1) facet orientation
        Basis_Sign[3] = -1.0;
        Basis_Sign[4] = -1.0;
        Basis_Sign[5] = -1.0;
        Basis_Sign[9] = -1.0;
        Basis_Sign[10] = -1.0;
        Basis_Sign[11] = -1.0;
        break;
    case 11: // (-1, -1, 1, -1) facet orientation
        Basis_Sign[0] = -1.0;
        Basis_Sign[1] = -1.0;
        Basis_Sign[2] = -1.0;
        Basis_Sign[3] = -1.0;
        Basis_Sign[4] = -1.0;
        Basis_Sign[5] = -1.0;
        Basis_Sign[9] = -1.0;
        Basis_Sign[10] = -1.0;
        Basis_Sign[11] = -1.0;
        break;
    case 12: // (1, 1, -1, -1) facet orientation
        Basis_Sign[6] = -1.0;
        Basis_Sign[7] = -1.0;
        Basis_Sign[8] = -1.0;
        Basis_Sign[9] = -1.0;
        Basis_Sign[10] = -1.0;
        Basis_Sign[11] = -1.0;
        break;
    case 13: // (-1, 1, -1, -1) facet orientation
        Basis_Sign[0] = -1.0;
        Basis_Sign[1] = -1.0;
        Basis_Sign[2] = -1.0;
        Basis_Sign[6] = -1.0;
        Basis_Sign[7] = -1.0;
        Basis_Sign[8] = -1.0;
        Basis_Sign[9] = -1.0;
        Basis_Sign[10] = -1.0;
        Basis_Sign[11] = -1.0;
        break;
    case 14: // (1, -1, -1, -1) facet orientation
        Basis_Sign[3] = -1.0;
        Basis_Sign[4] = -1.0;
        Basis_Sign[5] = -1.0;
        Basis_Sign[6] = -1.0;
        Basis_Sign[7] = -1.0;
        Basis_Sign[8] = -1.0;
        Basis_Sign[9] = -1.0;
        Basis_Sign[10] = -1.0;
        Basis_Sign[11] = -1.0;
        break;
    case 15: // (-1, -1, -1, -1) facet orientation
        Basis_Sign[0] = -1.0;
        Basis_Sign[1] = -1.0;
        Basis_Sign[2] = -1.0;
        Basis_Sign[3] = -1.0;
        Basis_Sign[4] = -1.0;
        Basis_Sign[5] = -1.0;
        Basis_Sign[6] = -1.0;
        Basis_Sign[7] = -1.0;
        Basis_Sign[8] = -1.0;
        Basis_Sign[9] = -1.0;
        Basis_Sign[10] = -1.0;
        Basis_Sign[11] = -1.0;
        break;
    default: ; // Facet_Orientation_Signature==0 (1, 1, 1, 1) facet orientation
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
    /* call sub-routine to determine sign flips for H(div) basis functions */
    double Basis_Sign[NB] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; // init to all positive (no sign change)
    Get_Basis_Sign_Change(Mesh->Orientation, Basis_Sign);
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, raviart_thomas_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 3
// geometric dimension = 3
// Number of Quadrature Points = 10

    // get "Val" of basis functions
    VEC_3x1 phi_Val[NQ][NB];
    // get "Div" of basis functions
    SCALAR phi_Div[NQ][NB];

    phi_Val[0][0].Set_Equal_To(2.87223286114884235E+00, 1.59003371985847486E+00, 1.59003371985847486E+00);
    phi_Val[0][1].Set_Equal_To(-9.57410953716280710E-01, -1.10543712769036562E+00, -2.42298296084054510E-01);
    phi_Val[0][2].Set_Equal_To(-9.57410953716280710E-01, -2.42298296084054510E-01, -1.10543712769036562E+00);
    phi_Val[0][3].Set_Equal_To(1.59003371985847486E+00, -2.42298296084054621E-01, -2.42298296084054621E-01);
    phi_Val[0][4].Set_Equal_To(1.59003371985847464E+00, -2.42298296084054510E-01, -1.10543712769036562E+00);
    phi_Val[0][5].Set_Equal_To(1.59003371985847464E+00, -1.10543712769036562E+00, -2.42298296084054510E-01);
    phi_Val[0][6].Set_Equal_To(-9.57410953716281155E-01, 2.30514637749070150E+00, -2.42298296084054621E-01);
    phi_Val[0][7].Set_Equal_To(2.87223286114884235E+00, -6.05230030086579163E+00, 1.59003371985847486E+00);
    phi_Val[0][8].Set_Equal_To(-9.57410953716280710E-01, 2.30514637749070062E+00, -1.10543712769036562E+00);
    phi_Val[0][9].Set_Equal_To(-9.57410953716281155E-01, -2.42298296084054621E-01, 2.30514637749070150E+00);
    phi_Val[0][10].Set_Equal_To(-9.57410953716280710E-01, -1.10543712769036562E+00, 2.30514637749070062E+00);
    phi_Val[0][11].Set_Equal_To(2.87223286114884235E+00, 1.59003371985847486E+00, -6.05230030086579163E+00);
    phi_Val[0][12].Set_Equal_To(9.81269020585914298E+00, -1.83233201594252937E+00, -1.83233201594252937E+00);
    phi_Val[0][13].Set_Equal_To(-4.73313541857348988E-16, 4.31569415803155554E+00, -1.19784575537172798E-16);
    phi_Val[0][14].Set_Equal_To(-4.73313541857348988E-16, -1.19784575537172798E-16, 4.31569415803155554E+00);
    phi_Val[1][0].Set_Equal_To(-1.10543712769036562E+00, -2.42298296084054510E-01, -2.42298296084054510E-01);
    phi_Val[1][1].Set_Equal_To(-2.42298296084054510E-01, -1.10543712769036562E+00, -2.42298296084054510E-01);
    phi_Val[1][2].Set_Equal_To(-2.42298296084054510E-01, -2.42298296084054510E-01, -1.10543712769036562E+00);
    phi_Val[1][3].Set_Equal_To(-6.05230030086579163E+00, 1.59003371985847464E+00, 1.59003371985847464E+00);
    phi_Val[1][4].Set_Equal_To(2.30514637749070062E+00, -2.42298296084054510E-01, -1.10543712769036562E+00);
    phi_Val[1][5].Set_Equal_To(2.30514637749070062E+00, -1.10543712769036562E+00, -2.42298296084054510E-01);
    phi_Val[1][6].Set_Equal_To(1.59003371985847464E+00, -6.05230030086579163E+00, 1.59003371985847464E+00);
    phi_Val[1][7].Set_Equal_To(-1.10543712769036562E+00, 2.30514637749070062E+00, -2.42298296084054510E-01);
    phi_Val[1][8].Set_Equal_To(-2.42298296084054510E-01, 2.30514637749070062E+00, -1.10543712769036562E+00);
    phi_Val[1][9].Set_Equal_To(1.59003371985847464E+00, 1.59003371985847464E+00, -6.05230030086579163E+00);
    phi_Val[1][10].Set_Equal_To(-2.42298296084054510E-01, -1.10543712769036562E+00, 2.30514637749070062E+00);
    phi_Val[1][11].Set_Equal_To(-1.10543712769036562E+00, -2.42298296084054510E-01, 2.30514637749070062E+00);
    phi_Val[1][12].Set_Equal_To(6.14802617397408468E+00, 1.83233201594252915E+00, 1.83233201594252915E+00);
    phi_Val[1][13].Set_Equal_To(1.83233201594252915E+00, 6.14802617397408468E+00, 1.83233201594252915E+00);
    phi_Val[1][14].Set_Equal_To(1.83233201594252915E+00, 1.83233201594252915E+00, 6.14802617397408468E+00);
    phi_Val[2][0].Set_Equal_To(-1.10543712769036562E+00, -2.42298296084054510E-01, -9.57410953716280710E-01);
    phi_Val[2][1].Set_Equal_To(-2.42298296084054510E-01, -1.10543712769036562E+00, -9.57410953716280710E-01);
    phi_Val[2][2].Set_Equal_To(1.59003371985847486E+00, 1.59003371985847486E+00, 2.87223286114884235E+00);
    phi_Val[2][3].Set_Equal_To(2.30514637749070150E+00, -2.42298296084054621E-01, -9.57410953716281155E-01);
    phi_Val[2][4].Set_Equal_To(-6.05230030086579163E+00, 1.59003371985847486E+00, 2.87223286114884235E+00);
    phi_Val[2][5].Set_Equal_To(2.30514637749070062E+00, -1.10543712769036562E+00, -9.57410953716280710E-01);
    phi_Val[2][6].Set_Equal_To(-2.42298296084054621E-01, 2.30514637749070150E+00, -9.57410953716281155E-01);
    phi_Val[2][7].Set_Equal_To(-1.10543712769036562E+00, 2.30514637749070062E+00, -9.57410953716280710E-01);
    phi_Val[2][8].Set_Equal_To(1.59003371985847486E+00, -6.05230030086579163E+00, 2.87223286114884235E+00);
    phi_Val[2][9].Set_Equal_To(-2.42298296084054621E-01, -2.42298296084054621E-01, 1.59003371985847486E+00);
    phi_Val[2][10].Set_Equal_To(-2.42298296084054510E-01, -1.10543712769036562E+00, 1.59003371985847464E+00);
    phi_Val[2][11].Set_Equal_To(-1.10543712769036562E+00, -2.42298296084054510E-01, 1.59003371985847464E+00);
    phi_Val[2][12].Set_Equal_To(4.31569415803155554E+00, -1.19784575537172798E-16, -4.73313541857348988E-16);
    phi_Val[2][13].Set_Equal_To(-1.19784575537172798E-16, 4.31569415803155554E+00, -4.73313541857348988E-16);
    phi_Val[2][14].Set_Equal_To(-1.83233201594252937E+00, -1.83233201594252937E+00, 9.81269020585914298E+00);
    phi_Val[3][0].Set_Equal_To(-1.10543712769036562E+00, -9.57410953716280710E-01, -2.42298296084054510E-01);
    phi_Val[3][1].Set_Equal_To(1.59003371985847486E+00, 2.87223286114884235E+00, 1.59003371985847486E+00);
    phi_Val[3][2].Set_Equal_To(-2.42298296084054510E-01, -9.57410953716280710E-01, -1.10543712769036562E+00);
    phi_Val[3][3].Set_Equal_To(2.30514637749070150E+00, -9.57410953716281155E-01, -2.42298296084054621E-01);
    phi_Val[3][4].Set_Equal_To(2.30514637749070062E+00, -9.57410953716280710E-01, -1.10543712769036562E+00);
    phi_Val[3][5].Set_Equal_To(-6.05230030086579163E+00, 2.87223286114884235E+00, 1.59003371985847486E+00);
    phi_Val[3][6].Set_Equal_To(-2.42298296084054621E-01, 1.59003371985847486E+00, -2.42298296084054621E-01);
    phi_Val[3][7].Set_Equal_To(-1.10543712769036562E+00, 1.59003371985847464E+00, -2.42298296084054510E-01);
    phi_Val[3][8].Set_Equal_To(-2.42298296084054510E-01, 1.59003371985847464E+00, -1.10543712769036562E+00);
    phi_Val[3][9].Set_Equal_To(-2.42298296084054621E-01, -9.57410953716281155E-01, 2.30514637749070150E+00);
    phi_Val[3][10].Set_Equal_To(1.59003371985847486E+00, 2.87223286114884235E+00, -6.05230030086579163E+00);
    phi_Val[3][11].Set_Equal_To(-1.10543712769036562E+00, -9.57410953716280710E-01, 2.30514637749070062E+00);
    phi_Val[3][12].Set_Equal_To(4.31569415803155554E+00, -4.73313541857348988E-16, -1.19784575537172798E-16);
    phi_Val[3][13].Set_Equal_To(-1.83233201594252937E+00, 9.81269020585914298E+00, -1.83233201594252937E+00);
    phi_Val[3][14].Set_Equal_To(-1.19784575537172798E-16, -4.73313541857348988E-16, 4.31569415803155554E+00);
    phi_Val[4][0].Set_Equal_To(0.00000000000000000E+00, -3.00000000000000000E+00, -3.00000000000000000E+00);
    phi_Val[4][1].Set_Equal_To(0.00000000000000000E+00, 1.50000000000000000E+00, 4.50000000000000000E+00);
    phi_Val[4][2].Set_Equal_To(0.00000000000000000E+00, 4.50000000000000000E+00, 1.50000000000000000E+00);
    phi_Val[4][3].Set_Equal_To(6.00000000000000000E+00, -3.00000000000000000E+00, -3.00000000000000000E+00);
    phi_Val[4][4].Set_Equal_To(-6.00000000000000000E+00, 4.50000000000000000E+00, 1.50000000000000000E+00);
    phi_Val[4][5].Set_Equal_To(-6.00000000000000000E+00, 1.50000000000000000E+00, 4.50000000000000000E+00);
    phi_Val[4][6].Set_Equal_To(0.00000000000000000E+00, 3.00000000000000000E+00, -3.00000000000000000E+00);
    phi_Val[4][7].Set_Equal_To(0.00000000000000000E+00, 3.00000000000000000E+00, -3.00000000000000000E+00);
    phi_Val[4][8].Set_Equal_To(0.00000000000000000E+00, -1.50000000000000000E+00, 1.50000000000000000E+00);
    phi_Val[4][9].Set_Equal_To(0.00000000000000000E+00, -3.00000000000000000E+00, 3.00000000000000000E+00);
    phi_Val[4][10].Set_Equal_To(0.00000000000000000E+00, 1.50000000000000000E+00, -1.50000000000000000E+00);
    phi_Val[4][11].Set_Equal_To(0.00000000000000000E+00, -3.00000000000000000E+00, 3.00000000000000000E+00);
    phi_Val[4][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[4][13].Set_Equal_To(0.00000000000000000E+00, 7.50000000000000000E+00, -7.50000000000000000E+00);
    phi_Val[4][14].Set_Equal_To(0.00000000000000000E+00, -7.50000000000000000E+00, 7.50000000000000000E+00);
    phi_Val[5][0].Set_Equal_To(1.50000000000000000E+00, 0.00000000000000000E+00, 4.50000000000000000E+00);
    phi_Val[5][1].Set_Equal_To(-3.00000000000000000E+00, 0.00000000000000000E+00, -3.00000000000000000E+00);
    phi_Val[5][2].Set_Equal_To(4.50000000000000000E+00, 0.00000000000000000E+00, 1.50000000000000000E+00);
    phi_Val[5][3].Set_Equal_To(3.00000000000000000E+00, 0.00000000000000000E+00, -3.00000000000000000E+00);
    phi_Val[5][4].Set_Equal_To(-1.50000000000000000E+00, 0.00000000000000000E+00, 1.50000000000000000E+00);
    phi_Val[5][5].Set_Equal_To(3.00000000000000000E+00, 0.00000000000000000E+00, -3.00000000000000000E+00);
    phi_Val[5][6].Set_Equal_To(-3.00000000000000000E+00, 6.00000000000000000E+00, -3.00000000000000000E+00);
    phi_Val[5][7].Set_Equal_To(1.50000000000000000E+00, -6.00000000000000000E+00, 4.50000000000000000E+00);
    phi_Val[5][8].Set_Equal_To(4.50000000000000000E+00, -6.00000000000000000E+00, 1.50000000000000000E+00);
    phi_Val[5][9].Set_Equal_To(-3.00000000000000000E+00, 0.00000000000000000E+00, 3.00000000000000000E+00);
    phi_Val[5][10].Set_Equal_To(-3.00000000000000000E+00, 0.00000000000000000E+00, 3.00000000000000000E+00);
    phi_Val[5][11].Set_Equal_To(1.50000000000000000E+00, 0.00000000000000000E+00, -1.50000000000000000E+00);
    phi_Val[5][12].Set_Equal_To(7.50000000000000000E+00, 0.00000000000000000E+00, -7.50000000000000000E+00);
    phi_Val[5][13].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[5][14].Set_Equal_To(-7.50000000000000000E+00, 0.00000000000000000E+00, 7.50000000000000000E+00);
    phi_Val[6][0].Set_Equal_To(1.50000000000000000E+00, 4.50000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[6][1].Set_Equal_To(4.50000000000000000E+00, 1.50000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[6][2].Set_Equal_To(-3.00000000000000000E+00, -3.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[6][3].Set_Equal_To(3.00000000000000000E+00, -3.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[6][4].Set_Equal_To(3.00000000000000000E+00, -3.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[6][5].Set_Equal_To(-1.50000000000000000E+00, 1.50000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[6][6].Set_Equal_To(-3.00000000000000000E+00, 3.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[6][7].Set_Equal_To(1.50000000000000000E+00, -1.50000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[6][8].Set_Equal_To(-3.00000000000000000E+00, 3.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[6][9].Set_Equal_To(-3.00000000000000000E+00, -3.00000000000000000E+00, 6.00000000000000000E+00);
    phi_Val[6][10].Set_Equal_To(4.50000000000000000E+00, 1.50000000000000000E+00, -6.00000000000000000E+00);
    phi_Val[6][11].Set_Equal_To(1.50000000000000000E+00, 4.50000000000000000E+00, -6.00000000000000000E+00);
    phi_Val[6][12].Set_Equal_To(7.50000000000000000E+00, -7.50000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[6][13].Set_Equal_To(-7.50000000000000000E+00, 7.50000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[6][14].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[7][0].Set_Equal_To(1.50000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[7][1].Set_Equal_To(-3.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[7][2].Set_Equal_To(-3.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[7][3].Set_Equal_To(-1.50000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[7][4].Set_Equal_To(3.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[7][5].Set_Equal_To(3.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[7][6].Set_Equal_To(4.50000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[7][7].Set_Equal_To(1.50000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[7][8].Set_Equal_To(-3.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[7][9].Set_Equal_To(4.50000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
    phi_Val[7][10].Set_Equal_To(-3.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
    phi_Val[7][11].Set_Equal_To(1.50000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
    phi_Val[7][12].Set_Equal_To(1.50000000000000000E+01, 0.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[7][13].Set_Equal_To(7.50000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[7][14].Set_Equal_To(7.50000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[8][0].Set_Equal_To(0.00000000000000000E+00, -3.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[8][1].Set_Equal_To(0.00000000000000000E+00, 1.50000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[8][2].Set_Equal_To(0.00000000000000000E+00, -3.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[8][3].Set_Equal_To(-6.00000000000000000E+00, 4.50000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[8][4].Set_Equal_To(6.00000000000000000E+00, -3.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[8][5].Set_Equal_To(-6.00000000000000000E+00, 1.50000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[8][6].Set_Equal_To(0.00000000000000000E+00, -1.50000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[8][7].Set_Equal_To(0.00000000000000000E+00, 3.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[8][8].Set_Equal_To(0.00000000000000000E+00, 3.00000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[8][9].Set_Equal_To(0.00000000000000000E+00, 4.50000000000000000E+00, -6.00000000000000000E+00);
    phi_Val[8][10].Set_Equal_To(0.00000000000000000E+00, 1.50000000000000000E+00, -6.00000000000000000E+00);
    phi_Val[8][11].Set_Equal_To(0.00000000000000000E+00, -3.00000000000000000E+00, 6.00000000000000000E+00);
    phi_Val[8][12].Set_Equal_To(0.00000000000000000E+00, 7.50000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[8][13].Set_Equal_To(0.00000000000000000E+00, 1.50000000000000000E+01, 0.00000000000000000E+00);
    phi_Val[8][14].Set_Equal_To(0.00000000000000000E+00, 7.50000000000000000E+00, 0.00000000000000000E+00);
    phi_Val[9][0].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -3.00000000000000000E+00);
    phi_Val[9][1].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -3.00000000000000000E+00);
    phi_Val[9][2].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 1.50000000000000000E+00);
    phi_Val[9][3].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 4.50000000000000000E+00);
    phi_Val[9][4].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 1.50000000000000000E+00);
    phi_Val[9][5].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -3.00000000000000000E+00);
    phi_Val[9][6].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 4.50000000000000000E+00);
    phi_Val[9][7].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, -3.00000000000000000E+00);
    phi_Val[9][8].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 1.50000000000000000E+00);
    phi_Val[9][9].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.50000000000000000E+00);
    phi_Val[9][10].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 3.00000000000000000E+00);
    phi_Val[9][11].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 3.00000000000000000E+00);
    phi_Val[9][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 7.50000000000000000E+00);
    phi_Val[9][13].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 7.50000000000000000E+00);
    phi_Val[9][14].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 1.50000000000000000E+01);

    phi_Div[0][0].Set_Equal_To(4.42116701036213300E+01);
    phi_Div[0][1].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[0][2].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[0][3].Set_Equal_To(-6.73722336787377962E+00);
    phi_Div[0][4].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[0][5].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[0][6].Set_Equal_To(-6.73722336787377962E+00);
    phi_Div[0][7].Set_Equal_To(4.42116701036213300E+01);
    phi_Div[0][8].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[0][9].Set_Equal_To(-6.73722336787377962E+00);
    phi_Div[0][10].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[0][11].Set_Equal_To(4.42116701036213300E+01);
    phi_Div[0][12].Set_Equal_To(-5.09488934714951114E+01);
    phi_Div[0][13].Set_Equal_To(-3.33066907387546962E-15);
    phi_Div[0][14].Set_Equal_To(-3.33066907387546962E-15);
    phi_Div[1][0].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[1][1].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[1][2].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[1][3].Set_Equal_To(4.42116701036213300E+01);
    phi_Div[1][4].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[1][5].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[1][6].Set_Equal_To(4.42116701036213300E+01);
    phi_Div[1][7].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[1][8].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[1][9].Set_Equal_To(4.42116701036213300E+01);
    phi_Div[1][10].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[1][11].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[1][12].Set_Equal_To(5.09488934714951043E+01);
    phi_Div[1][13].Set_Equal_To(5.09488934714951043E+01);
    phi_Div[1][14].Set_Equal_To(5.09488934714951043E+01);
    phi_Div[2][0].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[2][1].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[2][2].Set_Equal_To(4.42116701036213300E+01);
    phi_Div[2][3].Set_Equal_To(-6.73722336787377962E+00);
    phi_Div[2][4].Set_Equal_To(4.42116701036213300E+01);
    phi_Div[2][5].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[2][6].Set_Equal_To(-6.73722336787377962E+00);
    phi_Div[2][7].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[2][8].Set_Equal_To(4.42116701036213300E+01);
    phi_Div[2][9].Set_Equal_To(-6.73722336787377962E+00);
    phi_Div[2][10].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[2][11].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[2][12].Set_Equal_To(-3.33066907387546962E-15);
    phi_Div[2][13].Set_Equal_To(-3.33066907387546962E-15);
    phi_Div[2][14].Set_Equal_To(-5.09488934714951114E+01);
    phi_Div[3][0].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[3][1].Set_Equal_To(4.42116701036213300E+01);
    phi_Div[3][2].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[3][3].Set_Equal_To(-6.73722336787377962E+00);
    phi_Div[3][4].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[3][5].Set_Equal_To(4.42116701036213300E+01);
    phi_Div[3][6].Set_Equal_To(-6.73722336787377962E+00);
    phi_Div[3][7].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[3][8].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[3][9].Set_Equal_To(-6.73722336787377962E+00);
    phi_Div[3][10].Set_Equal_To(4.42116701036213300E+01);
    phi_Div[3][11].Set_Equal_To(-6.73722336787377696E+00);
    phi_Div[3][12].Set_Equal_To(-3.33066907387546962E-15);
    phi_Div[3][13].Set_Equal_To(-5.09488934714951114E+01);
    phi_Div[3][14].Set_Equal_To(-3.33066907387546962E-15);
    phi_Div[4][0].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[4][1].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[4][2].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[4][3].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[4][4].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[4][5].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[4][6].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[4][7].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[4][8].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[4][9].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[4][10].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[4][11].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[4][12].Set_Equal_To(0.00000000000000000E+00);
    phi_Div[4][13].Set_Equal_To(-6.00000000000000000E+01);
    phi_Div[4][14].Set_Equal_To(-6.00000000000000000E+01);
    phi_Div[5][0].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[5][1].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[5][2].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[5][3].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[5][4].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[5][5].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[5][6].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[5][7].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[5][8].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[5][9].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[5][10].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[5][11].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[5][12].Set_Equal_To(-6.00000000000000000E+01);
    phi_Div[5][13].Set_Equal_To(0.00000000000000000E+00);
    phi_Div[5][14].Set_Equal_To(-6.00000000000000000E+01);
    phi_Div[6][0].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[6][1].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[6][2].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[6][3].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[6][4].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[6][5].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[6][6].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[6][7].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[6][8].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[6][9].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[6][10].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[6][11].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[6][12].Set_Equal_To(-6.00000000000000000E+01);
    phi_Div[6][13].Set_Equal_To(-6.00000000000000000E+01);
    phi_Div[6][14].Set_Equal_To(0.00000000000000000E+00);
    phi_Div[7][0].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[7][1].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[7][2].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[7][3].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[7][4].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[7][5].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[7][6].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[7][7].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[7][8].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[7][9].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[7][10].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[7][11].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[7][12].Set_Equal_To(0.00000000000000000E+00);
    phi_Div[7][13].Set_Equal_To(6.00000000000000000E+01);
    phi_Div[7][14].Set_Equal_To(6.00000000000000000E+01);
    phi_Div[8][0].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[8][1].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[8][2].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[8][3].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[8][4].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[8][5].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[8][6].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[8][7].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[8][8].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[8][9].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[8][10].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[8][11].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[8][12].Set_Equal_To(6.00000000000000000E+01);
    phi_Div[8][13].Set_Equal_To(0.00000000000000000E+00);
    phi_Div[8][14].Set_Equal_To(6.00000000000000000E+01);
    phi_Div[9][0].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[9][1].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[9][2].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[9][3].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[9][4].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[9][5].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[9][6].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[9][7].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[9][8].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[9][9].Set_Equal_To(3.60000000000000000E+01);
    phi_Div[9][10].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[9][11].Set_Equal_To(-2.40000000000000000E+01);
    phi_Div[9][12].Set_Equal_To(6.00000000000000000E+01);
    phi_Div[9][13].Set_Equal_To(6.00000000000000000E+01);
    phi_Div[9][14].Set_Equal_To(0.00000000000000000E+00);


// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {5.68430584196844446E-01, 1.43856471934385194E-01, 1.43856471934385194E-01}, \
//     {1.43856471934385194E-01, 1.43856471934385194E-01, 1.43856471934385194E-01}, \
//     {1.43856471934385194E-01, 1.43856471934385194E-01, 5.68430584196844446E-01}, \
//     {1.43856471934385194E-01, 5.68430584196844446E-01, 1.43856471934385194E-01}, \
//     {0.00000000000000000E+00, 5.00000000000000000E-01, 5.00000000000000000E-01}, \
//     {5.00000000000000000E-01, 0.00000000000000000E+00, 5.00000000000000000E-01}, \
//     {5.00000000000000000E-01, 5.00000000000000000E-01, 0.00000000000000000E+00}, \
//     {5.00000000000000000E-01, 0.00000000000000000E+00, 0.00000000000000000E+00}, \
//     {0.00000000000000000E+00, 5.00000000000000000E-01, 0.00000000000000000E+00}, \
//     {0.00000000000000000E+00, 0.00000000000000000E+00, 5.00000000000000000E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     3.62941783134008988E-02, \
//     3.62941783134008988E-02, \
//     3.62941783134008988E-02, \
//     3.62941783134008988E-02, \
//     3.58165890217718337E-03, \
//     3.58165890217718337E-03, \
//     3.58165890217718337E-03, \
//     3.58165890217718337E-03, \
//     3.58165890217718337E-03, \
//     3.58165890217718337E-03  \
//     };
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
            VEC_3x1 vv_temp;
            // pre-multiply by sign change
            const double INV_DET_SIGN_FLIP = Basis_Sign[basis_i] * Mesh->Map_Inv_Det_Jac[0].a;
            // pre-multiply by 1/det(Jac)
            Scalar_Mult_Vector(phi_Val[qp_i][basis_i], INV_DET_SIGN_FLIP, vv_temp);
            Mat_Vec(Mesh->Map_PHI_Grad[0], vv_temp, Func_vv_Value[basis_i][qp_i]);
            }
        }
    // map divergence of basis vectors over (indexing is in the C style)
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        // evaluate for each basis function
        for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
            {
            // pre-multiply by 1/det(Jac) (and flip sign if necessary)
            Func_vv_Div[basis_i][qp_i].a = Basis_Sign[basis_i] * (phi_Div[qp_i][basis_i].a) * Mesh->Map_Inv_Det_Jac[0].a;
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
