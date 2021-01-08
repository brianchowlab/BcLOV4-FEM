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
#define SpecificFUNC        Data_Type_V_h_phi_restricted_to_Omega
#define SpecificFUNC_str   "Data_Type_V_h_phi_restricted_to_Omega"

// set the type of function space
#define SPACE_type  "CG - lagrange_deg2_dim2"
// set the name of function space
#define SPACE_name  "V_h"

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
#define NQ  37
// set the number of basis functions
#define NB  6
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
// Local Element defined on Subdomain: CG, lagrange_deg2_dim2
// the Subdomain             has topological dimension = 2
// the Domain of Integration has topological dimension = 2
// geometric dimension = 2
// Number of Quadrature Points = 37

    // get "Grad" of basis functions
    VEC_2x1 phi_Grad[NQ][NB];

    phi_Grad[0][0].Set_Equal_To(-3.33333333333333315E-01, -3.33333333333333315E-01);
    phi_Grad[0][1].Set_Equal_To(3.33333333333333315E-01, 0.00000000000000000E+00);
    phi_Grad[0][2].Set_Equal_To(0.00000000000000000E+00, 3.33333333333333315E-01);
    phi_Grad[0][3].Set_Equal_To(1.33333333333333326E+00, 1.33333333333333326E+00);
    phi_Grad[0][4].Set_Equal_To(-1.33333333333333326E+00, 0.00000000000000000E+00);
    phi_Grad[0][5].Set_Equal_To(0.00000000000000000E+00, -1.33333333333333326E+00);
    phi_Grad[1][0].Set_Equal_To(9.00551325848211293E-01, 9.00551325848211293E-01);
    phi_Grad[1][1].Set_Equal_To(2.80110265169642236E+00, 0.00000000000000000E+00);
    phi_Grad[1][2].Set_Equal_To(0.00000000000000000E+00, -9.00551325848211071E-01);
    phi_Grad[1][3].Set_Equal_To(9.94486741517888739E-02, 3.80110265169642236E+00);
    phi_Grad[1][4].Set_Equal_To(-9.94486741517888739E-02, -1.11022302462515654E-16);
    phi_Grad[1][5].Set_Equal_To(-3.70165397754463354E+00, -3.80110265169642236E+00);
    phi_Grad[2][0].Set_Equal_To(9.00551325848211293E-01, 9.00551325848211293E-01);
    phi_Grad[2][1].Set_Equal_To(-9.00551325848211071E-01, 0.00000000000000000E+00);
    phi_Grad[2][2].Set_Equal_To(0.00000000000000000E+00, 2.80110265169642236E+00);
    phi_Grad[2][3].Set_Equal_To(3.80110265169642236E+00, 9.94486741517888739E-02);
    phi_Grad[2][4].Set_Equal_To(-3.80110265169642236E+00, -3.70165397754463354E+00);
    phi_Grad[2][5].Set_Equal_To(-1.11022302462515654E-16, -9.94486741517888739E-02);
    phi_Grad[3][0].Set_Equal_To(-2.80110265169642236E+00, -2.80110265169642236E+00);
    phi_Grad[3][1].Set_Equal_To(-9.00551325848211071E-01, 0.00000000000000000E+00);
    phi_Grad[3][2].Set_Equal_To(0.00000000000000000E+00, -9.00551325848211071E-01);
    phi_Grad[3][3].Set_Equal_To(9.94486741517888739E-02, 9.94486741517888739E-02);
    phi_Grad[3][4].Set_Equal_To(-9.94486741517888739E-02, 3.70165397754463354E+00);
    phi_Grad[3][5].Set_Equal_To(3.70165397754463354E+00, -9.94486741517888739E-02);
    phi_Grad[4][0].Set_Equal_To(-6.56770170152329191E-01, -6.56770170152329191E-01);
    phi_Grad[4][1].Set_Equal_To(-3.13540340304658605E-01, 0.00000000000000000E+00);
    phi_Grad[4][2].Set_Equal_To(0.00000000000000000E+00, 6.56770170152329413E-01);
    phi_Grad[4][3].Set_Equal_To(1.65677017015232941E+00, 6.86459659695341395E-01);
    phi_Grad[4][4].Set_Equal_To(-1.65677017015232941E+00, -2.22044604925031308E-16);
    phi_Grad[4][5].Set_Equal_To(9.70310510456987796E-01, -6.86459659695341395E-01);
    phi_Grad[5][0].Set_Equal_To(-6.56770170152329191E-01, -6.56770170152329191E-01);
    phi_Grad[5][1].Set_Equal_To(6.56770170152329413E-01, 0.00000000000000000E+00);
    phi_Grad[5][2].Set_Equal_To(0.00000000000000000E+00, -3.13540340304658605E-01);
    phi_Grad[5][3].Set_Equal_To(6.86459659695341395E-01, 1.65677017015232941E+00);
    phi_Grad[5][4].Set_Equal_To(-6.86459659695341395E-01, 9.70310510456987796E-01);
    phi_Grad[5][5].Set_Equal_To(-2.22044604925031308E-16, -1.65677017015232941E+00);
    phi_Grad[6][0].Set_Equal_To(3.13540340304658827E-01, 3.13540340304658827E-01);
    phi_Grad[6][1].Set_Equal_To(6.56770170152329413E-01, 0.00000000000000000E+00);
    phi_Grad[6][2].Set_Equal_To(0.00000000000000000E+00, 6.56770170152329413E-01);
    phi_Grad[6][3].Set_Equal_To(1.65677017015232941E+00, 1.65677017015232941E+00);
    phi_Grad[6][4].Set_Equal_To(-1.65677017015232941E+00, -9.70310510456988240E-01);
    phi_Grad[6][5].Set_Equal_To(-9.70310510456988240E-01, -1.65677017015232941E+00);
    phi_Grad[7][0].Set_Equal_To(7.88244873543808566E-02, 7.88244873543808566E-02);
    phi_Grad[7][1].Set_Equal_To(1.15764897470876171E+00, 0.00000000000000000E+00);
    phi_Grad[7][2].Set_Equal_To(0.00000000000000000E+00, -7.88244873543808566E-02);
    phi_Grad[7][3].Set_Equal_To(9.21175512645619143E-01, 2.15764897470876171E+00);
    phi_Grad[7][4].Set_Equal_To(-9.21175512645619143E-01, 0.00000000000000000E+00);
    phi_Grad[7][5].Set_Equal_To(-1.23647346206314257E+00, -2.15764897470876171E+00);
    phi_Grad[8][0].Set_Equal_To(7.88244873543808566E-02, 7.88244873543808566E-02);
    phi_Grad[8][1].Set_Equal_To(-7.88244873543808566E-02, 0.00000000000000000E+00);
    phi_Grad[8][2].Set_Equal_To(0.00000000000000000E+00, 1.15764897470876171E+00);
    phi_Grad[8][3].Set_Equal_To(2.15764897470876171E+00, 9.21175512645619143E-01);
    phi_Grad[8][4].Set_Equal_To(-2.15764897470876171E+00, -1.23647346206314257E+00);
    phi_Grad[8][5].Set_Equal_To(0.00000000000000000E+00, -9.21175512645619143E-01);
    phi_Grad[9][0].Set_Equal_To(-1.15764897470876171E+00, -1.15764897470876171E+00);
    phi_Grad[9][1].Set_Equal_To(-7.88244873543808566E-02, 0.00000000000000000E+00);
    phi_Grad[9][2].Set_Equal_To(0.00000000000000000E+00, -7.88244873543808566E-02);
    phi_Grad[9][3].Set_Equal_To(9.21175512645619143E-01, 9.21175512645619143E-01);
    phi_Grad[9][4].Set_Equal_To(-9.21175512645619143E-01, 1.23647346206314257E+00);
    phi_Grad[9][5].Set_Equal_To(1.23647346206314257E+00, -9.21175512645619143E-01);
    phi_Grad[10][0].Set_Equal_To(5.44320073353064982E-01, 5.44320073353064982E-01);
    phi_Grad[10][1].Set_Equal_To(2.08864014670613019E+00, 0.00000000000000000E+00);
    phi_Grad[10][2].Set_Equal_To(0.00000000000000000E+00, -5.44320073353065204E-01);
    phi_Grad[10][3].Set_Equal_To(4.55679926646934852E-01, 3.08864014670613019E+00);
    phi_Grad[10][4].Set_Equal_To(-4.55679926646934852E-01, 1.11022302462515654E-16);
    phi_Grad[10][5].Set_Equal_To(-2.63296022005919506E+00, -3.08864014670613019E+00);
    phi_Grad[11][0].Set_Equal_To(5.44320073353064982E-01, 5.44320073353064982E-01);
    phi_Grad[11][1].Set_Equal_To(-5.44320073353065204E-01, 0.00000000000000000E+00);
    phi_Grad[11][2].Set_Equal_To(0.00000000000000000E+00, 2.08864014670613019E+00);
    phi_Grad[11][3].Set_Equal_To(3.08864014670613019E+00, 4.55679926646934852E-01);
    phi_Grad[11][4].Set_Equal_To(-3.08864014670613019E+00, -2.63296022005919506E+00);
    phi_Grad[11][5].Set_Equal_To(1.11022302462515654E-16, -4.55679926646934852E-01);
    phi_Grad[12][0].Set_Equal_To(-2.08864014670613019E+00, -2.08864014670613019E+00);
    phi_Grad[12][1].Set_Equal_To(-5.44320073353065204E-01, 0.00000000000000000E+00);
    phi_Grad[12][2].Set_Equal_To(0.00000000000000000E+00, -5.44320073353065204E-01);
    phi_Grad[12][3].Set_Equal_To(4.55679926646934852E-01, 4.55679926646934852E-01);
    phi_Grad[12][4].Set_Equal_To(-4.55679926646934852E-01, 2.63296022005919550E+00);
    phi_Grad[12][5].Set_Equal_To(2.63296022005919550E+00, -4.55679926646934852E-01);
    phi_Grad[13][0].Set_Equal_To(-9.81829200100329369E-01, -9.81829200100329369E-01);
    phi_Grad[13][1].Set_Equal_To(-9.63658400200658627E-01, 0.00000000000000000E+00);
    phi_Grad[13][2].Set_Equal_To(0.00000000000000000E+00, 9.81829200100329258E-01);
    phi_Grad[13][3].Set_Equal_To(1.98182920010032926E+00, 3.63415997993414147E-02);
    phi_Grad[13][4].Set_Equal_To(-1.98182920010032926E+00, 6.93889390390722838E-17);
    phi_Grad[13][5].Set_Equal_To(1.94548760030098800E+00, -3.63415997993414147E-02);
    phi_Grad[14][0].Set_Equal_To(-9.81829200100329369E-01, -9.81829200100329369E-01);
    phi_Grad[14][1].Set_Equal_To(9.81829200100329258E-01, 0.00000000000000000E+00);
    phi_Grad[14][2].Set_Equal_To(0.00000000000000000E+00, -9.63658400200658627E-01);
    phi_Grad[14][3].Set_Equal_To(3.63415997993414147E-02, 1.98182920010032926E+00);
    phi_Grad[14][4].Set_Equal_To(-3.63415997993414147E-02, 1.94548760030098800E+00);
    phi_Grad[14][5].Set_Equal_To(6.93889390390722838E-17, -1.98182920010032926E+00);
    phi_Grad[15][0].Set_Equal_To(9.63658400200658516E-01, 9.63658400200658516E-01);
    phi_Grad[15][1].Set_Equal_To(9.81829200100329258E-01, 0.00000000000000000E+00);
    phi_Grad[15][2].Set_Equal_To(0.00000000000000000E+00, 9.81829200100329258E-01);
    phi_Grad[15][3].Set_Equal_To(1.98182920010032926E+00, 1.98182920010032926E+00);
    phi_Grad[15][4].Set_Equal_To(-1.98182920010032926E+00, -1.94548760030098777E+00);
    phi_Grad[15][5].Set_Equal_To(-1.94548760030098777E+00, -1.98182920010032926E+00);
    phi_Grad[16][0].Set_Equal_To(-8.75445419388226065E-01, -8.75445419388226065E-01);
    phi_Grad[16][1].Set_Equal_To(-7.50890838776452019E-01, 0.00000000000000000E+00);
    phi_Grad[16][2].Set_Equal_To(0.00000000000000000E+00, 8.75445419388225954E-01);
    phi_Grad[16][3].Set_Equal_To(1.87544541938822595E+00, 2.49109161223547981E-01);
    phi_Grad[16][4].Set_Equal_To(-1.87544541938822595E+00, 1.11022302462515654E-16);
    phi_Grad[16][5].Set_Equal_To(1.62633625816467808E+00, -2.49109161223547981E-01);
    phi_Grad[17][0].Set_Equal_To(-8.75445419388226065E-01, -8.75445419388226065E-01);
    phi_Grad[17][1].Set_Equal_To(8.75445419388225954E-01, 0.00000000000000000E+00);
    phi_Grad[17][2].Set_Equal_To(0.00000000000000000E+00, -7.50890838776452019E-01);
    phi_Grad[17][3].Set_Equal_To(2.49109161223547981E-01, 1.87544541938822595E+00);
    phi_Grad[17][4].Set_Equal_To(-2.49109161223547981E-01, 1.62633625816467808E+00);
    phi_Grad[17][5].Set_Equal_To(1.11022302462515654E-16, -1.87544541938822595E+00);
    phi_Grad[18][0].Set_Equal_To(7.50890838776451908E-01, 7.50890838776451908E-01);
    phi_Grad[18][1].Set_Equal_To(8.75445419388225954E-01, 0.00000000000000000E+00);
    phi_Grad[18][2].Set_Equal_To(0.00000000000000000E+00, 8.75445419388225954E-01);
    phi_Grad[18][3].Set_Equal_To(1.87544541938822595E+00, 1.87544541938822595E+00);
    phi_Grad[18][4].Set_Equal_To(-1.87544541938822595E+00, -1.62633625816467786E+00);
    phi_Grad[18][5].Set_Equal_To(-1.62633625816467786E+00, -1.87544541938822595E+00);
    phi_Grad[19][0].Set_Equal_To(4.93531175311891612E-01, 4.93531175311891612E-01);
    phi_Grad[19][1].Set_Equal_To(-9.11694841385502386E-01, 0.00000000000000000E+00);
    phi_Grad[19][2].Set_Equal_To(0.00000000000000000E+00, 2.40522601669739400E+00);
    phi_Grad[19][3].Set_Equal_To(3.40522601669739400E+00, 8.83051586144976142E-02);
    phi_Grad[19][4].Set_Equal_To(-3.40522601669739400E+00, -2.89875719200928561E+00);
    phi_Grad[19][5].Set_Equal_To(4.18163666073610774E-01, -8.83051586144976142E-02);
    phi_Grad[20][0].Set_Equal_To(-2.40522601669739400E+00, -2.40522601669739400E+00);
    phi_Grad[20][1].Set_Equal_To(-9.11694841385502386E-01, 0.00000000000000000E+00);
    phi_Grad[20][2].Set_Equal_To(0.00000000000000000E+00, -4.93531175311891612E-01);
    phi_Grad[20][3].Set_Equal_To(5.06468824688108388E-01, 8.83051586144976142E-02);
    phi_Grad[20][4].Set_Equal_To(-5.06468824688108388E-01, 2.89875719200928561E+00);
    phi_Grad[20][5].Set_Equal_To(3.31692085808289638E+00, -8.83051586144976142E-02);
    phi_Grad[21][0].Set_Equal_To(4.93531175311891612E-01, 4.93531175311891612E-01);
    phi_Grad[21][1].Set_Equal_To(2.40522601669739400E+00, 0.00000000000000000E+00);
    phi_Grad[21][2].Set_Equal_To(0.00000000000000000E+00, -9.11694841385502386E-01);
    phi_Grad[21][3].Set_Equal_To(8.83051586144976142E-02, 3.40522601669739400E+00);
    phi_Grad[21][4].Set_Equal_To(-8.83051586144976142E-02, 4.18163666073610774E-01);
    phi_Grad[21][5].Set_Equal_To(-2.89875719200928561E+00, -3.40522601669739400E+00);
    phi_Grad[22][0].Set_Equal_To(9.11694841385502386E-01, 9.11694841385502386E-01);
    phi_Grad[22][1].Set_Equal_To(2.40522601669739400E+00, 0.00000000000000000E+00);
    phi_Grad[22][2].Set_Equal_To(0.00000000000000000E+00, -4.93531175311891612E-01);
    phi_Grad[22][3].Set_Equal_To(5.06468824688108388E-01, 3.40522601669739400E+00);
    phi_Grad[22][4].Set_Equal_To(-5.06468824688108388E-01, -4.18163666073610774E-01);
    phi_Grad[22][5].Set_Equal_To(-3.31692085808289638E+00, -3.40522601669739400E+00);
    phi_Grad[23][0].Set_Equal_To(-2.40522601669739400E+00, -2.40522601669739400E+00);
    phi_Grad[23][1].Set_Equal_To(-4.93531175311891612E-01, 0.00000000000000000E+00);
    phi_Grad[23][2].Set_Equal_To(0.00000000000000000E+00, -9.11694841385502386E-01);
    phi_Grad[23][3].Set_Equal_To(8.83051586144976142E-02, 5.06468824688108388E-01);
    phi_Grad[23][4].Set_Equal_To(-8.83051586144976142E-02, 3.31692085808289638E+00);
    phi_Grad[23][5].Set_Equal_To(2.89875719200928561E+00, -5.06468824688108388E-01);
    phi_Grad[24][0].Set_Equal_To(9.11694841385502386E-01, 9.11694841385502386E-01);
    phi_Grad[24][1].Set_Equal_To(-4.93531175311891612E-01, 0.00000000000000000E+00);
    phi_Grad[24][2].Set_Equal_To(0.00000000000000000E+00, 2.40522601669739400E+00);
    phi_Grad[24][3].Set_Equal_To(3.40522601669739400E+00, 5.06468824688108388E-01);
    phi_Grad[24][4].Set_Equal_To(-3.40522601669739400E+00, -3.31692085808289638E+00);
    phi_Grad[24][5].Set_Equal_To(-4.18163666073610774E-01, -5.06468824688108388E-01);
    phi_Grad[25][0].Set_Equal_To(-1.67750025875551073E-01, -1.67750025875551073E-01);
    phi_Grad[25][1].Set_Equal_To(-9.25517908789916133E-01, 0.00000000000000000E+00);
    phi_Grad[25][2].Set_Equal_To(0.00000000000000000E+00, 1.75776788291436503E+00);
    phi_Grad[25][3].Set_Equal_To(2.75776788291436503E+00, 7.44820912100838811E-02);
    phi_Grad[25][4].Set_Equal_To(-2.75776788291436503E+00, -1.59001785703881393E+00);
    phi_Grad[25][5].Set_Equal_To(1.09326793466546723E+00, -7.44820912100838811E-02);
    phi_Grad[26][0].Set_Equal_To(-1.75776788291436503E+00, -1.75776788291436503E+00);
    phi_Grad[26][1].Set_Equal_To(-9.25517908789916133E-01, 0.00000000000000000E+00);
    phi_Grad[26][2].Set_Equal_To(0.00000000000000000E+00, 1.67750025875551101E-01);
    phi_Grad[26][3].Set_Equal_To(1.16775002587555110E+00, 7.44820912100838811E-02);
    phi_Grad[26][4].Set_Equal_To(-1.16775002587555110E+00, 1.59001785703881393E+00);
    phi_Grad[26][5].Set_Equal_To(2.68328579170428094E+00, -7.44820912100838811E-02);
    phi_Grad[27][0].Set_Equal_To(-1.67750025875551073E-01, -1.67750025875551073E-01);
    phi_Grad[27][1].Set_Equal_To(1.75776788291436503E+00, 0.00000000000000000E+00);
    phi_Grad[27][2].Set_Equal_To(0.00000000000000000E+00, -9.25517908789916133E-01);
    phi_Grad[27][3].Set_Equal_To(7.44820912100838811E-02, 2.75776788291436503E+00);
    phi_Grad[27][4].Set_Equal_To(-7.44820912100838811E-02, 1.09326793466546723E+00);
    phi_Grad[27][5].Set_Equal_To(-1.59001785703881393E+00, -2.75776788291436503E+00);
    phi_Grad[28][0].Set_Equal_To(9.25517908789916133E-01, 9.25517908789916133E-01);
    phi_Grad[28][1].Set_Equal_To(1.75776788291436503E+00, 0.00000000000000000E+00);
    phi_Grad[28][2].Set_Equal_To(0.00000000000000000E+00, 1.67750025875551101E-01);
    phi_Grad[28][3].Set_Equal_To(1.16775002587555110E+00, 2.75776788291436503E+00);
    phi_Grad[28][4].Set_Equal_To(-1.16775002587555110E+00, -1.09326793466546723E+00);
    phi_Grad[28][5].Set_Equal_To(-2.68328579170428139E+00, -2.75776788291436503E+00);
    phi_Grad[29][0].Set_Equal_To(-1.75776788291436503E+00, -1.75776788291436503E+00);
    phi_Grad[29][1].Set_Equal_To(1.67750025875551101E-01, 0.00000000000000000E+00);
    phi_Grad[29][2].Set_Equal_To(0.00000000000000000E+00, -9.25517908789916133E-01);
    phi_Grad[29][3].Set_Equal_To(7.44820912100838811E-02, 1.16775002587555110E+00);
    phi_Grad[29][4].Set_Equal_To(-7.44820912100838811E-02, 2.68328579170428094E+00);
    phi_Grad[29][5].Set_Equal_To(1.59001785703881393E+00, -1.16775002587555110E+00);
    phi_Grad[30][0].Set_Equal_To(9.25517908789916133E-01, 9.25517908789916133E-01);
    phi_Grad[30][1].Set_Equal_To(1.67750025875551101E-01, 0.00000000000000000E+00);
    phi_Grad[30][2].Set_Equal_To(0.00000000000000000E+00, 1.75776788291436503E+00);
    phi_Grad[30][3].Set_Equal_To(2.75776788291436503E+00, 1.16775002587555110E+00);
    phi_Grad[30][4].Set_Equal_To(-2.75776788291436503E+00, -2.68328579170428139E+00);
    phi_Grad[30][5].Set_Equal_To(-1.09326793466546723E+00, -1.16775002587555110E+00);
    phi_Grad[31][0].Set_Equal_To(-7.05026370958718451E-02, -7.05026370958718451E-02);
    phi_Grad[31][1].Set_Equal_To(-6.13974074831363015E-01, 0.00000000000000000E+00);
    phi_Grad[31][2].Set_Equal_To(0.00000000000000000E+00, 1.54347143773549123E+00);
    phi_Grad[31][3].Set_Equal_To(2.54347143773549123E+00, 3.86025925168636930E-01);
    phi_Grad[31][4].Set_Equal_To(-2.54347143773549123E+00, -1.47296880063961932E+00);
    phi_Grad[31][5].Set_Equal_To(6.84476711927234915E-01, -3.86025925168636930E-01);
    phi_Grad[32][0].Set_Equal_To(-1.54347143773549123E+00, -1.54347143773549123E+00);
    phi_Grad[32][1].Set_Equal_To(-6.13974074831363015E-01, 0.00000000000000000E+00);
    phi_Grad[32][2].Set_Equal_To(0.00000000000000000E+00, 7.05026370958719006E-02);
    phi_Grad[32][3].Set_Equal_To(1.07050263709587190E+00, 3.86025925168636930E-01);
    phi_Grad[32][4].Set_Equal_To(-1.07050263709587190E+00, 1.47296880063961932E+00);
    phi_Grad[32][5].Set_Equal_To(2.15744551256685435E+00, -3.86025925168636930E-01);
    phi_Grad[33][0].Set_Equal_To(-7.05026370958718451E-02, -7.05026370958718451E-02);
    phi_Grad[33][1].Set_Equal_To(1.54347143773549123E+00, 0.00000000000000000E+00);
    phi_Grad[33][2].Set_Equal_To(0.00000000000000000E+00, -6.13974074831363015E-01);
    phi_Grad[33][3].Set_Equal_To(3.86025925168636930E-01, 2.54347143773549123E+00);
    phi_Grad[33][4].Set_Equal_To(-3.86025925168636930E-01, 6.84476711927234915E-01);
    phi_Grad[33][5].Set_Equal_To(-1.47296880063961932E+00, -2.54347143773549123E+00);
    phi_Grad[34][0].Set_Equal_To(6.13974074831363126E-01, 6.13974074831363126E-01);
    phi_Grad[34][1].Set_Equal_To(1.54347143773549123E+00, 0.00000000000000000E+00);
    phi_Grad[34][2].Set_Equal_To(0.00000000000000000E+00, 7.05026370958719006E-02);
    phi_Grad[34][3].Set_Equal_To(1.07050263709587190E+00, 2.54347143773549123E+00);
    phi_Grad[34][4].Set_Equal_To(-1.07050263709587190E+00, -6.84476711927235026E-01);
    phi_Grad[34][5].Set_Equal_To(-2.15744551256685435E+00, -2.54347143773549123E+00);
    phi_Grad[35][0].Set_Equal_To(-1.54347143773549123E+00, -1.54347143773549123E+00);
    phi_Grad[35][1].Set_Equal_To(7.05026370958719006E-02, 0.00000000000000000E+00);
    phi_Grad[35][2].Set_Equal_To(0.00000000000000000E+00, -6.13974074831363015E-01);
    phi_Grad[35][3].Set_Equal_To(3.86025925168636930E-01, 1.07050263709587190E+00);
    phi_Grad[35][4].Set_Equal_To(-3.86025925168636930E-01, 2.15744551256685435E+00);
    phi_Grad[35][5].Set_Equal_To(1.47296880063961932E+00, -1.07050263709587190E+00);
    phi_Grad[36][0].Set_Equal_To(6.13974074831363126E-01, 6.13974074831363126E-01);
    phi_Grad[36][1].Set_Equal_To(7.05026370958719006E-02, 0.00000000000000000E+00);
    phi_Grad[36][2].Set_Equal_To(0.00000000000000000E+00, 1.54347143773549123E+00);
    phi_Grad[36][3].Set_Equal_To(2.54347143773549123E+00, 1.07050263709587190E+00);
    phi_Grad[36][4].Set_Equal_To(-2.54347143773549123E+00, -2.15744551256685435E+00);
    phi_Grad[36][5].Set_Equal_To(-6.84476711927235026E-01, -1.07050263709587190E+00);


// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {3.33333333333333315E-01, 3.33333333333333315E-01}, \
//     {9.50275662924105591E-01, 2.48621685379472185E-02}, \
//     {2.48621685379472185E-02, 9.50275662924105591E-01}, \
//     {2.48621685379472185E-02, 2.48621685379472185E-02}, \
//     {1.71614914923835349E-01, 4.14192542538082353E-01}, \
//     {4.14192542538082353E-01, 1.71614914923835349E-01}, \
//     {4.14192542538082353E-01, 4.14192542538082353E-01}, \
//     {5.39412243677190428E-01, 2.30293878161404786E-01}, \
//     {2.30293878161404786E-01, 5.39412243677190428E-01}, \
//     {2.30293878161404786E-01, 2.30293878161404786E-01}, \
//     {7.72160036676532546E-01, 1.13919981661733713E-01}, \
//     {1.13919981661733713E-01, 7.72160036676532546E-01}, \
//     {1.13919981661733713E-01, 1.13919981661733713E-01}, \
//     {9.08539994983535368E-03, 4.95457300025082314E-01}, \
//     {4.95457300025082314E-01, 9.08539994983535368E-03}, \
//     {4.95457300025082314E-01, 4.95457300025082314E-01}, \
//     {6.22772903058869953E-02, 4.68861354847056488E-01}, \
//     {4.68861354847056488E-01, 6.22772903058869953E-02}, \
//     {4.68861354847056488E-01, 4.68861354847056488E-01}, \
//     {2.20762896536244035E-02, 8.51306504174348500E-01}, \
//     {2.20762896536244035E-02, 1.26617206172027097E-01}, \
//     {8.51306504174348500E-01, 2.20762896536244035E-02}, \
//     {8.51306504174348500E-01, 1.26617206172027097E-01}, \
//     {1.26617206172027097E-01, 2.20762896536244035E-02}, \
//     {1.26617206172027097E-01, 8.51306504174348500E-01}, \
//     {1.86205228025209703E-02, 6.89441970728591258E-01}, \
//     {1.86205228025209703E-02, 2.91937506468887775E-01}, \
//     {6.89441970728591258E-01, 1.86205228025209703E-02}, \
//     {6.89441970728591258E-01, 2.91937506468887775E-01}, \
//     {2.91937506468887775E-01, 1.86205228025209703E-02}, \
//     {2.91937506468887775E-01, 6.89441970728591258E-01}, \
//     {9.65064812921592324E-02, 6.35867859433872806E-01}, \
//     {9.65064812921592324E-02, 2.67625659273967975E-01}, \
//     {6.35867859433872806E-01, 9.65064812921592324E-02}, \
//     {6.35867859433872806E-01, 2.67625659273967975E-01}, \
//     {2.67625659273967975E-01, 9.65064812921592324E-02}, \
//     {2.67625659273967975E-01, 6.35867859433872806E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.58698830328720659E-02, \
//     4.00389977778240098E-03, \
//     4.00389977778240098E-03, \
//     4.00389977778240098E-03, \
//     2.34344494909108220E-02, \
//     2.34344494909108220E-02, \
//     2.34344494909108220E-02, \
//     2.32954700919882456E-02, \
//     2.32954700919882456E-02, \
//     2.32954700919882456E-02, \
//     1.55084716568981915E-02, \
//     1.55084716568981915E-02, \
//     1.55084716568981915E-02, \
//     5.39580636831563706E-03, \
//     5.39580636831563706E-03, \
//     5.39580636831563706E-03, \
//     1.60977671212158106E-02, \
//     1.60977671212158106E-02, \
//     1.60977671212158106E-02, \
//     7.72291710535079161E-03, \
//     7.72291710535079161E-03, \
//     7.72291710535079161E-03, \
//     7.72291710535079161E-03, \
//     7.72291710535079161E-03, \
//     7.72291710535079161E-03, \
//     8.91149496158933144E-03, \
//     8.91149496158933144E-03, \
//     8.91149496158933144E-03, \
//     8.91149496158933144E-03, \
//     8.91149496158933144E-03, \
//     8.91149496158933144E-03, \
//     1.85193418406923126E-02, \
//     1.85193418406923126E-02, \
//     1.85193418406923126E-02, \
//     1.85193418406923126E-02, \
//     1.85193418406923126E-02, \
//     1.85193418406923126E-02  \
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
// Local Element defined on Subdomain: CG, lagrange_deg2_dim2
// the Subdomain             has topological dimension = 2
// the Domain of Integration has topological dimension = 2
// geometric dimension = 2
// Number of Quadrature Points = 37

    // get "Val" of basis functions
    SCALAR phi_Val[NQ][NB];

    phi_Val[0][0].Set_Equal_To(-1.11111111111111105E-01);
    phi_Val[0][1].Set_Equal_To(-1.11111111111111105E-01);
    phi_Val[0][2].Set_Equal_To(-1.11111111111111105E-01);
    phi_Val[0][3].Set_Equal_To(4.44444444444444420E-01);
    phi_Val[0][4].Set_Equal_To(4.44444444444444420E-01);
    phi_Val[0][5].Set_Equal_To(4.44444444444444420E-01);
    phi_Val[1][0].Set_Equal_To(-2.36259136891286067E-02);
    phi_Val[1][1].Set_Equal_To(8.55772008167591136E-01);
    phi_Val[1][2].Set_Equal_To(-2.36259136891286345E-02);
    phi_Val[1][3].Set_Equal_To(9.45036547565145379E-02);
    phi_Val[1][4].Set_Equal_To(2.47250969763716756E-03);
    phi_Val[1][5].Set_Equal_To(9.45036547565144269E-02);
    phi_Val[2][0].Set_Equal_To(-2.36259136891286067E-02);
    phi_Val[2][1].Set_Equal_To(-2.36259136891286345E-02);
    phi_Val[2][2].Set_Equal_To(8.55772008167591136E-01);
    phi_Val[2][3].Set_Equal_To(9.45036547565145379E-02);
    phi_Val[2][4].Set_Equal_To(9.45036547565144269E-02);
    phi_Val[2][5].Set_Equal_To(2.47250969763716756E-03);
    phi_Val[3][0].Set_Equal_To(8.55772008167591025E-01);
    phi_Val[3][1].Set_Equal_To(-2.36259136891286345E-02);
    phi_Val[3][2].Set_Equal_To(-2.36259136891286345E-02);
    phi_Val[3][3].Set_Equal_To(2.47250969763717016E-03);
    phi_Val[3][4].Set_Equal_To(9.45036547565145379E-02);
    phi_Val[3][5].Set_Equal_To(9.45036547565145379E-02);
    phi_Val[4][0].Set_Equal_To(-7.10816179497600759E-02);
    phi_Val[4][1].Set_Equal_To(-1.12711556875204857E-01);
    phi_Val[4][2].Set_Equal_To(-7.10816179497600342E-02);
    phi_Val[4][3].Set_Equal_To(2.84326471799040248E-01);
    phi_Val[4][4].Set_Equal_To(6.86221849176644527E-01);
    phi_Val[4][5].Set_Equal_To(2.84326471799040192E-01);
    phi_Val[5][0].Set_Equal_To(-7.10816179497600759E-02);
    phi_Val[5][1].Set_Equal_To(-7.10816179497600342E-02);
    phi_Val[5][2].Set_Equal_To(-1.12711556875204857E-01);
    phi_Val[5][3].Set_Equal_To(2.84326471799040248E-01);
    phi_Val[5][4].Set_Equal_To(2.84326471799040192E-01);
    phi_Val[5][5].Set_Equal_To(6.86221849176644527E-01);
    phi_Val[6][0].Set_Equal_To(-1.12711556875204844E-01);
    phi_Val[6][1].Set_Equal_To(-7.10816179497600342E-02);
    phi_Val[6][2].Set_Equal_To(-7.10816179497600342E-02);
    phi_Val[6][3].Set_Equal_To(6.86221849176644638E-01);
    phi_Val[6][4].Set_Equal_To(2.84326471799040137E-01);
    phi_Val[6][5].Set_Equal_To(2.84326471799040137E-01);
    phi_Val[7][0].Set_Equal_To(-1.24223337524164881E-01);
    phi_Val[7][1].Set_Equal_To(4.25188935805309032E-02);
    phi_Val[7][2].Set_Equal_To(-1.24223337524164881E-01);
    phi_Val[7][3].Set_Equal_To(4.96893350096659525E-01);
    phi_Val[7][4].Set_Equal_To(2.12141081274479809E-01);
    phi_Val[7][5].Set_Equal_To(4.96893350096659525E-01);
    phi_Val[8][0].Set_Equal_To(-1.24223337524164881E-01);
    phi_Val[8][1].Set_Equal_To(-1.24223337524164881E-01);
    phi_Val[8][2].Set_Equal_To(4.25188935805309032E-02);
    phi_Val[8][3].Set_Equal_To(4.96893350096659525E-01);
    phi_Val[8][4].Set_Equal_To(4.96893350096659525E-01);
    phi_Val[8][5].Set_Equal_To(2.12141081274479809E-01);
    phi_Val[9][0].Set_Equal_To(4.25188935805309032E-02);
    phi_Val[9][1].Set_Equal_To(-1.24223337524164881E-01);
    phi_Val[9][2].Set_Equal_To(-1.24223337524164881E-01);
    phi_Val[9][3].Set_Equal_To(2.12141081274479809E-01);
    phi_Val[9][4].Set_Equal_To(4.96893350096659525E-01);
    phi_Val[9][5].Set_Equal_To(4.96893350096659525E-01);
    phi_Val[10][0].Set_Equal_To(-8.79644572181142309E-02);
    phi_Val[10][1].Set_Equal_To(4.20302207804075623E-01);
    phi_Val[10][2].Set_Equal_To(-8.79644572181142170E-02);
    phi_Val[10][3].Set_Equal_To(3.51857828872456868E-01);
    phi_Val[10][4].Set_Equal_To(5.19110488872389919E-02);
    phi_Val[10][5].Set_Equal_To(3.51857828872456979E-01);
    phi_Val[11][0].Set_Equal_To(-8.79644572181142309E-02);
    phi_Val[11][1].Set_Equal_To(-8.79644572181142170E-02);
    phi_Val[11][2].Set_Equal_To(4.20302207804075623E-01);
    phi_Val[11][3].Set_Equal_To(3.51857828872456868E-01);
    phi_Val[11][4].Set_Equal_To(3.51857828872456979E-01);
    phi_Val[11][5].Set_Equal_To(5.19110488872389919E-02);
    phi_Val[12][0].Set_Equal_To(4.20302207804075678E-01);
    phi_Val[12][1].Set_Equal_To(-8.79644572181142170E-02);
    phi_Val[12][2].Set_Equal_To(-8.79644572181142170E-02);
    phi_Val[12][3].Set_Equal_To(5.19110488872389850E-02);
    phi_Val[12][4].Set_Equal_To(3.51857828872456868E-01);
    phi_Val[12][5].Set_Equal_To(3.51857828872456868E-01);
    phi_Val[13][0].Set_Equal_To(-4.50142772879343414E-03);
    phi_Val[13][1].Set_Equal_To(-8.92031096533841758E-03);
    phi_Val[13][2].Set_Equal_To(-4.50142772879345149E-03);
    phi_Val[13][3].Set_Equal_To(1.80057109151737713E-02);
    phi_Val[13][4].Set_Equal_To(9.81911744592577795E-01);
    phi_Val[13][5].Set_Equal_To(1.80057109151737713E-02);
    phi_Val[14][0].Set_Equal_To(-4.50142772879343414E-03);
    phi_Val[14][1].Set_Equal_To(-4.50142772879345149E-03);
    phi_Val[14][2].Set_Equal_To(-8.92031096533841758E-03);
    phi_Val[14][3].Set_Equal_To(1.80057109151737713E-02);
    phi_Val[14][4].Set_Equal_To(1.80057109151737713E-02);
    phi_Val[14][5].Set_Equal_To(9.81911744592577795E-01);
    phi_Val[15][0].Set_Equal_To(-8.92031096533843319E-03);
    phi_Val[15][1].Set_Equal_To(-4.50142772879345149E-03);
    phi_Val[15][2].Set_Equal_To(-4.50142772879345149E-03);
    phi_Val[15][3].Set_Equal_To(9.81911744592577684E-01);
    phi_Val[15][4].Set_Equal_To(1.80057109151738060E-02);
    phi_Val[15][5].Set_Equal_To(1.80057109151738060E-02);
    phi_Val[16][0].Set_Equal_To(-2.91994147090216226E-02);
    phi_Val[16][1].Set_Equal_To(-5.45203685301995436E-02);
    phi_Val[16][2].Set_Equal_To(-2.91994147090216469E-02);
    phi_Val[16][3].Set_Equal_To(1.16797658836086532E-01);
    phi_Val[16][4].Set_Equal_To(8.79323880276069780E-01);
    phi_Val[16][5].Set_Equal_To(1.16797658836086546E-01);
    phi_Val[17][0].Set_Equal_To(-2.91994147090216226E-02);
    phi_Val[17][1].Set_Equal_To(-2.91994147090216469E-02);
    phi_Val[17][2].Set_Equal_To(-5.45203685301995436E-02);
    phi_Val[17][3].Set_Equal_To(1.16797658836086532E-01);
    phi_Val[17][4].Set_Equal_To(1.16797658836086546E-01);
    phi_Val[17][5].Set_Equal_To(8.79323880276069780E-01);
    phi_Val[18][0].Set_Equal_To(-5.45203685301995644E-02);
    phi_Val[18][1].Set_Equal_To(-2.91994147090216469E-02);
    phi_Val[18][2].Set_Equal_To(-2.91994147090216469E-02);
    phi_Val[18][3].Set_Equal_To(8.79323880276069669E-01);
    phi_Val[18][4].Set_Equal_To(1.16797658836086587E-01);
    phi_Val[18][5].Set_Equal_To(1.16797658836086587E-01);
    phi_Val[19][0].Set_Equal_To(-9.45533723744078625E-02);
    phi_Val[19][1].Set_Equal_To(-2.11015645238829541E-02);
    phi_Val[19][2].Set_Equal_To(5.98139023924751601E-01);
    phi_Val[19][3].Set_Equal_To(7.51747558806693250E-02);
    phi_Val[19][4].Set_Equal_To(4.31160204618524501E-01);
    phi_Val[19][5].Set_Equal_To(1.11809524743453990E-02);
    phi_Val[20][0].Set_Equal_To(5.98139023924751601E-01);
    phi_Val[20][1].Set_Equal_To(-2.11015645238829541E-02);
    phi_Val[20][2].Set_Equal_To(-9.45533723744078625E-02);
    phi_Val[20][3].Set_Equal_To(1.11809524743453990E-02);
    phi_Val[20][4].Set_Equal_To(4.31160204618524501E-01);
    phi_Val[20][5].Set_Equal_To(7.51747558806693250E-02);
    phi_Val[21][0].Set_Equal_To(-9.45533723744078625E-02);
    phi_Val[21][1].Set_Equal_To(5.98139023924751601E-01);
    phi_Val[21][2].Set_Equal_To(-2.11015645238829541E-02);
    phi_Val[21][3].Set_Equal_To(7.51747558806693250E-02);
    phi_Val[21][4].Set_Equal_To(1.11809524743453990E-02);
    phi_Val[21][5].Set_Equal_To(4.31160204618524501E-01);
    phi_Val[22][0].Set_Equal_To(-2.11015645238829541E-02);
    phi_Val[22][1].Set_Equal_To(5.98139023924751601E-01);
    phi_Val[22][2].Set_Equal_To(-9.45533723744078625E-02);
    phi_Val[22][3].Set_Equal_To(4.31160204618524501E-01);
    phi_Val[22][4].Set_Equal_To(1.11809524743453990E-02);
    phi_Val[22][5].Set_Equal_To(7.51747558806693250E-02);
    phi_Val[23][0].Set_Equal_To(5.98139023924751601E-01);
    phi_Val[23][1].Set_Equal_To(-9.45533723744078625E-02);
    phi_Val[23][2].Set_Equal_To(-2.11015645238829541E-02);
    phi_Val[23][3].Set_Equal_To(1.11809524743453990E-02);
    phi_Val[23][4].Set_Equal_To(7.51747558806693250E-02);
    phi_Val[23][5].Set_Equal_To(4.31160204618524501E-01);
    phi_Val[24][0].Set_Equal_To(-2.11015645238829541E-02);
    phi_Val[24][1].Set_Equal_To(-9.45533723744078625E-02);
    phi_Val[24][2].Set_Equal_To(5.98139023924751601E-01);
    phi_Val[24][3].Set_Equal_To(4.31160204618524501E-01);
    phi_Val[24][4].Set_Equal_To(7.51747558806693250E-02);
    phi_Val[24][5].Set_Equal_To(1.11809524743453990E-02);
    phi_Val[25][0].Set_Equal_To(-1.21482491102343992E-01);
    phi_Val[25][1].Set_Equal_To(-1.79270750636425641E-02);
    phi_Val[25][2].Set_Equal_To(2.61218491275656106E-01);
    phi_Val[25][3].Set_Equal_To(5.13510797478669151E-02);
    phi_Val[25][4].Set_Equal_To(8.05095879158003402E-01);
    phi_Val[25][5].Set_Equal_To(2.17441159844601536E-02);
    phi_Val[26][0].Set_Equal_To(2.61218491275656106E-01);
    phi_Val[26][1].Set_Equal_To(-1.79270750636425641E-02);
    phi_Val[26][2].Set_Equal_To(-1.21482491102343992E-01);
    phi_Val[26][3].Set_Equal_To(2.17441159844601536E-02);
    phi_Val[26][4].Set_Equal_To(8.05095879158003402E-01);
    phi_Val[26][5].Set_Equal_To(5.13510797478669151E-02);
    phi_Val[27][0].Set_Equal_To(-1.21482491102343992E-01);
    phi_Val[27][1].Set_Equal_To(2.61218491275656106E-01);
    phi_Val[27][2].Set_Equal_To(-1.79270750636425641E-02);
    phi_Val[27][3].Set_Equal_To(5.13510797478669151E-02);
    phi_Val[27][4].Set_Equal_To(2.17441159844601536E-02);
    phi_Val[27][5].Set_Equal_To(8.05095879158003402E-01);
    phi_Val[28][0].Set_Equal_To(-1.79270750636425606E-02);
    phi_Val[28][1].Set_Equal_To(2.61218491275656106E-01);
    phi_Val[28][2].Set_Equal_To(-1.21482491102343992E-01);
    phi_Val[28][3].Set_Equal_To(8.05095879158003402E-01);
    phi_Val[28][4].Set_Equal_To(2.17441159844601467E-02);
    phi_Val[28][5].Set_Equal_To(5.13510797478669082E-02);
    phi_Val[29][0].Set_Equal_To(2.61218491275656106E-01);
    phi_Val[29][1].Set_Equal_To(-1.21482491102343992E-01);
    phi_Val[29][2].Set_Equal_To(-1.79270750636425641E-02);
    phi_Val[29][3].Set_Equal_To(2.17441159844601536E-02);
    phi_Val[29][4].Set_Equal_To(5.13510797478669151E-02);
    phi_Val[29][5].Set_Equal_To(8.05095879158003402E-01);
    phi_Val[30][0].Set_Equal_To(-1.79270750636425606E-02);
    phi_Val[30][1].Set_Equal_To(-1.21482491102343992E-01);
    phi_Val[30][2].Set_Equal_To(2.61218491275656106E-01);
    phi_Val[30][3].Set_Equal_To(8.05095879158003402E-01);
    phi_Val[30][4].Set_Equal_To(5.13510797478669082E-02);
    phi_Val[30][5].Set_Equal_To(2.17441159844601467E-02);
    phi_Val[31][0].Set_Equal_To(-1.24378672270315974E-01);
    phi_Val[31][1].Set_Equal_To(-7.78794794293714704E-02);
    phi_Val[31][2].Set_Equal_To(1.72788009888158040E-01);
    phi_Val[31][3].Set_Equal_To(2.45461478722961535E-01);
    phi_Val[31][4].Set_Equal_To(6.80698220368467943E-01);
    phi_Val[31][5].Set_Equal_To(1.03310442720099885E-01);
    phi_Val[32][0].Set_Equal_To(1.72788009888158012E-01);
    phi_Val[32][1].Set_Equal_To(-7.78794794293714704E-02);
    phi_Val[32][2].Set_Equal_To(-1.24378672270315974E-01);
    phi_Val[32][3].Set_Equal_To(1.03310442720099885E-01);
    phi_Val[32][4].Set_Equal_To(6.80698220368468054E-01);
    phi_Val[32][5].Set_Equal_To(2.45461478722961535E-01);
    phi_Val[33][0].Set_Equal_To(-1.24378672270315974E-01);
    phi_Val[33][1].Set_Equal_To(1.72788009888158040E-01);
    phi_Val[33][2].Set_Equal_To(-7.78794794293714704E-02);
    phi_Val[33][3].Set_Equal_To(2.45461478722961535E-01);
    phi_Val[33][4].Set_Equal_To(1.03310442720099885E-01);
    phi_Val[33][5].Set_Equal_To(6.80698220368467943E-01);
    phi_Val[34][0].Set_Equal_To(-7.78794794293714704E-02);
    phi_Val[34][1].Set_Equal_To(1.72788009888158040E-01);
    phi_Val[34][2].Set_Equal_To(-1.24378672270315974E-01);
    phi_Val[34][3].Set_Equal_To(6.80698220368468054E-01);
    phi_Val[34][4].Set_Equal_To(1.03310442720099871E-01);
    phi_Val[34][5].Set_Equal_To(2.45461478722961507E-01);
    phi_Val[35][0].Set_Equal_To(1.72788009888158012E-01);
    phi_Val[35][1].Set_Equal_To(-1.24378672270315974E-01);
    phi_Val[35][2].Set_Equal_To(-7.78794794293714704E-02);
    phi_Val[35][3].Set_Equal_To(1.03310442720099885E-01);
    phi_Val[35][4].Set_Equal_To(2.45461478722961535E-01);
    phi_Val[35][5].Set_Equal_To(6.80698220368468054E-01);
    phi_Val[36][0].Set_Equal_To(-7.78794794293714704E-02);
    phi_Val[36][1].Set_Equal_To(-1.24378672270315974E-01);
    phi_Val[36][2].Set_Equal_To(1.72788009888158040E-01);
    phi_Val[36][3].Set_Equal_To(6.80698220368468054E-01);
    phi_Val[36][4].Set_Equal_To(2.45461478722961507E-01);
    phi_Val[36][5].Set_Equal_To(1.03310442720099871E-01);


// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {3.33333333333333315E-01, 3.33333333333333315E-01}, \
//     {9.50275662924105591E-01, 2.48621685379472185E-02}, \
//     {2.48621685379472185E-02, 9.50275662924105591E-01}, \
//     {2.48621685379472185E-02, 2.48621685379472185E-02}, \
//     {1.71614914923835349E-01, 4.14192542538082353E-01}, \
//     {4.14192542538082353E-01, 1.71614914923835349E-01}, \
//     {4.14192542538082353E-01, 4.14192542538082353E-01}, \
//     {5.39412243677190428E-01, 2.30293878161404786E-01}, \
//     {2.30293878161404786E-01, 5.39412243677190428E-01}, \
//     {2.30293878161404786E-01, 2.30293878161404786E-01}, \
//     {7.72160036676532546E-01, 1.13919981661733713E-01}, \
//     {1.13919981661733713E-01, 7.72160036676532546E-01}, \
//     {1.13919981661733713E-01, 1.13919981661733713E-01}, \
//     {9.08539994983535368E-03, 4.95457300025082314E-01}, \
//     {4.95457300025082314E-01, 9.08539994983535368E-03}, \
//     {4.95457300025082314E-01, 4.95457300025082314E-01}, \
//     {6.22772903058869953E-02, 4.68861354847056488E-01}, \
//     {4.68861354847056488E-01, 6.22772903058869953E-02}, \
//     {4.68861354847056488E-01, 4.68861354847056488E-01}, \
//     {2.20762896536244035E-02, 8.51306504174348500E-01}, \
//     {2.20762896536244035E-02, 1.26617206172027097E-01}, \
//     {8.51306504174348500E-01, 2.20762896536244035E-02}, \
//     {8.51306504174348500E-01, 1.26617206172027097E-01}, \
//     {1.26617206172027097E-01, 2.20762896536244035E-02}, \
//     {1.26617206172027097E-01, 8.51306504174348500E-01}, \
//     {1.86205228025209703E-02, 6.89441970728591258E-01}, \
//     {1.86205228025209703E-02, 2.91937506468887775E-01}, \
//     {6.89441970728591258E-01, 1.86205228025209703E-02}, \
//     {6.89441970728591258E-01, 2.91937506468887775E-01}, \
//     {2.91937506468887775E-01, 1.86205228025209703E-02}, \
//     {2.91937506468887775E-01, 6.89441970728591258E-01}, \
//     {9.65064812921592324E-02, 6.35867859433872806E-01}, \
//     {9.65064812921592324E-02, 2.67625659273967975E-01}, \
//     {6.35867859433872806E-01, 9.65064812921592324E-02}, \
//     {6.35867859433872806E-01, 2.67625659273967975E-01}, \
//     {2.67625659273967975E-01, 9.65064812921592324E-02}, \
//     {2.67625659273967975E-01, 6.35867859433872806E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.58698830328720659E-02, \
//     4.00389977778240098E-03, \
//     4.00389977778240098E-03, \
//     4.00389977778240098E-03, \
//     2.34344494909108220E-02, \
//     2.34344494909108220E-02, \
//     2.34344494909108220E-02, \
//     2.32954700919882456E-02, \
//     2.32954700919882456E-02, \
//     2.32954700919882456E-02, \
//     1.55084716568981915E-02, \
//     1.55084716568981915E-02, \
//     1.55084716568981915E-02, \
//     5.39580636831563706E-03, \
//     5.39580636831563706E-03, \
//     5.39580636831563706E-03, \
//     1.60977671212158106E-02, \
//     1.60977671212158106E-02, \
//     1.60977671212158106E-02, \
//     7.72291710535079161E-03, \
//     7.72291710535079161E-03, \
//     7.72291710535079161E-03, \
//     7.72291710535079161E-03, \
//     7.72291710535079161E-03, \
//     7.72291710535079161E-03, \
//     8.91149496158933144E-03, \
//     8.91149496158933144E-03, \
//     8.91149496158933144E-03, \
//     8.91149496158933144E-03, \
//     8.91149496158933144E-03, \
//     8.91149496158933144E-03, \
//     1.85193418406923126E-02, \
//     1.85193418406923126E-02, \
//     1.85193418406923126E-02, \
//     1.85193418406923126E-02, \
//     1.85193418406923126E-02, \
//     1.85193418406923126E-02  \
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
