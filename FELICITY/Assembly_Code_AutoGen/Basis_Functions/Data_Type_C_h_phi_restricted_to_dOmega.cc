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
#define SpecificFUNC        Data_Type_C_h_phi_restricted_to_dOmega
#define SpecificFUNC_str   "Data_Type_C_h_phi_restricted_to_dOmega"

// set the type of function space
#define SPACE_type  "CG - lagrange_deg1_dim3"
// set the name of function space
#define SPACE_name  "C_h"

// set the Subdomain topological dimension
#define SUB_TD  3
// set the Domain of Integration (DoI) topological dimension
#define DOI_TD  2
// set the geometric dimension
#define GD  3
// set the number of cartesian tuple components (m*n) = 1 * 1
#define NC  1
// NOTE: the (i,j) tuple component is accessed by the linear index k = i + (j-1)*m
// set the number of quad points
#define NQ  12
// set the number of basis functions
#define NB  4
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

    // constructor
    SpecificFUNC ();
    ~SpecificFUNC (); // destructor
    void Setup_Function_Space(const mxArray*);
    void Get_Local_to_Global_DoFmap(const int&, int*) const;
                   // need the "const" to ENSURE that nothing in this object will change!
    void Transform_Basis_Functions();
    const CLASS_geom_Omega_embedded_in_Omega_restricted_to_dOmega*  Mesh;

private:
    void Map_Basis_p1();
    void Map_Basis_n1();
    void Map_Basis_p2();
    void Map_Basis_n2();
    void Map_Basis_p3();
    void Map_Basis_n3();
    void Map_Basis_p4();
    void Map_Basis_n4();
    SCALAR Value_p1[NB][NQ];
    SCALAR Value_n1[NB][NQ];
    SCALAR Value_p2[NB][NQ];
    SCALAR Value_n2[NB][NQ];
    SCALAR Value_p3[NB][NQ];
    SCALAR Value_n3[NB][NQ];
    SCALAR Value_p4[NB][NQ];
    SCALAR Value_n4[NB][NQ];
    void Basis_Value_p1(SCALAR V[NB][NQ]);
    void Basis_Value_n1(SCALAR V[NB][NQ]);
    void Basis_Value_p2(SCALAR V[NB][NQ]);
    void Basis_Value_n2(SCALAR V[NB][NQ]);
    void Basis_Value_p3(SCALAR V[NB][NQ]);
    void Basis_Value_n3(SCALAR V[NB][NQ]);
    void Basis_Value_p4(SCALAR V[NB][NQ]);
    void Basis_Value_n4(SCALAR V[NB][NQ]);
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
    // init basis function values on local mesh entities
    Basis_Value_p1(Value_p1);
    Basis_Value_n1(Value_n1);
    Basis_Value_p2(Value_p2);
    Basis_Value_n2(Value_n2);
    Basis_Value_p3(Value_p3);
    Basis_Value_n3(Value_n3);
    Basis_Value_p4(Value_p4);
    Basis_Value_n4(Value_n4);
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
    // read in the embedding info
    const int Entity_Index = Mesh->Domain->DoI_Entity_Index;

    // if orientation is positive
    if (Entity_Index > 0)
        {
        // pick the correct one
        if (Entity_Index==1) Map_Basis_p1();
        else if (Entity_Index==2) Map_Basis_p2();
        else if (Entity_Index==3) Map_Basis_p3();
        else if (Entity_Index==4) Map_Basis_p4();
        else mexErrMsgTxt("ERROR: Entity_Index outside valid range or is zero!");
        }
    else // else it is negative
        {
        // pick the correct one
        if (Entity_Index==-1) Map_Basis_n1();
        else if (Entity_Index==-2) Map_Basis_n2();
        else if (Entity_Index==-3) Map_Basis_n3();
        else if (Entity_Index==-4) Map_Basis_n4();
        else mexErrMsgTxt("ERROR: Entity_Index outside valid range or is zero!");
        }
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
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12



// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
//     };
/*------------   END: Auto Generate ------------*/
    // set basis function values to the correct mesh entity
    Func_f_Value = &Value_p1; // point to the correct mesh entity

/*------------ BEGIN: Auto Generate ------------*/
    /*** compute basis function quantities ***/
/*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* map basis functions from the standard reference element
       to an actual element in the Domain.     */
void SpecificFUNC::Map_Basis_n1()
{

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12



// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
//     };
/*------------   END: Auto Generate ------------*/
    // set basis function values to the correct mesh entity
    Func_f_Value = &Value_n1; // point to the correct mesh entity

/*------------ BEGIN: Auto Generate ------------*/
    /*** compute basis function quantities ***/
/*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* map basis functions from the standard reference element
       to an actual element in the Domain.     */
void SpecificFUNC::Map_Basis_p2()
{

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12



// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
//     };
/*------------   END: Auto Generate ------------*/
    // set basis function values to the correct mesh entity
    Func_f_Value = &Value_p2; // point to the correct mesh entity

/*------------ BEGIN: Auto Generate ------------*/
    /*** compute basis function quantities ***/
/*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* map basis functions from the standard reference element
       to an actual element in the Domain.     */
void SpecificFUNC::Map_Basis_n2()
{

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12



// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
//     };
/*------------   END: Auto Generate ------------*/
    // set basis function values to the correct mesh entity
    Func_f_Value = &Value_n2; // point to the correct mesh entity

/*------------ BEGIN: Auto Generate ------------*/
    /*** compute basis function quantities ***/
/*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* map basis functions from the standard reference element
       to an actual element in the Domain.     */
void SpecificFUNC::Map_Basis_p3()
{

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12



// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
//     };
/*------------   END: Auto Generate ------------*/
    // set basis function values to the correct mesh entity
    Func_f_Value = &Value_p3; // point to the correct mesh entity

/*------------ BEGIN: Auto Generate ------------*/
    /*** compute basis function quantities ***/
/*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* map basis functions from the standard reference element
       to an actual element in the Domain.     */
void SpecificFUNC::Map_Basis_n3()
{

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12



// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
//     };
/*------------   END: Auto Generate ------------*/
    // set basis function values to the correct mesh entity
    Func_f_Value = &Value_n3; // point to the correct mesh entity

/*------------ BEGIN: Auto Generate ------------*/
    /*** compute basis function quantities ***/
/*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* map basis functions from the standard reference element
       to an actual element in the Domain.     */
void SpecificFUNC::Map_Basis_p4()
{

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12



// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
//     };
/*------------   END: Auto Generate ------------*/
    // set basis function values to the correct mesh entity
    Func_f_Value = &Value_p4; // point to the correct mesh entity

/*------------ BEGIN: Auto Generate ------------*/
    /*** compute basis function quantities ***/
/*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* map basis functions from the standard reference element
       to an actual element in the Domain.     */
void SpecificFUNC::Map_Basis_n4()
{

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12



// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
//     };
/*------------   END: Auto Generate ------------*/
    // set basis function values to the correct mesh entity
    Func_f_Value = &Value_n4; // point to the correct mesh entity

/*------------ BEGIN: Auto Generate ------------*/
    /*** compute basis function quantities ***/
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
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12

    // get "Val" of basis functions
    SCALAR phi_Val[NQ][NB];

    phi_Val[0][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[0][1].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[0][2].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[0][3].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[1][1].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][2].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][3].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[2][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[2][1].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[2][2].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[2][3].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[3][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[3][1].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[3][2].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[3][3].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[4][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[4][1].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[4][2].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[4][3].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[5][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[5][1].Set_Equal_To(5.01426509658179920E-01);
    phi_Val[5][2].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[5][3].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[6][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[6][1].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[6][2].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[6][3].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[7][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[7][1].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[7][2].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[7][3].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[8][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[8][1].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[8][2].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[8][3].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[9][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[9][1].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[9][2].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[9][3].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[10][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[10][1].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[10][2].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[10][3].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[11][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[11][1].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[11][2].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[11][3].Set_Equal_To(3.10352451033785004E-01);


// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
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

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* evaluate basis functions (no derivatives!) on the local reference element. */
void SpecificFUNC::Basis_Value_n1(SCALAR BF_V[NB][NQ])
{

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12

    // get "Val" of basis functions
    SCALAR phi_Val[NQ][NB];

    phi_Val[0][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[0][1].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[0][2].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[0][3].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[1][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[1][1].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][2].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[1][3].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[2][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[2][1].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[2][2].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[2][3].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[3][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[3][1].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[3][2].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[3][3].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[4][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[4][1].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[4][2].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[4][3].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[5][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[5][1].Set_Equal_To(5.01426509658179920E-01);
    phi_Val[5][2].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[5][3].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[6][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[6][1].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[6][2].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[6][3].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[7][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[7][1].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[7][2].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[7][3].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[8][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[8][1].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[8][2].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[8][3].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[9][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[9][1].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[9][2].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[9][3].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[10][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[10][1].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[10][2].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[10][3].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[11][0].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[11][1].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[11][2].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[11][3].Set_Equal_To(5.31450498448160016E-02);


// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
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

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* evaluate basis functions (no derivatives!) on the local reference element. */
void SpecificFUNC::Basis_Value_p2(SCALAR BF_V[NB][NQ])
{

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12

    // get "Val" of basis functions
    SCALAR phi_Val[NQ][NB];

    phi_Val[0][0].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[0][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[0][2].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[0][3].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[1][0].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[1][2].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[1][3].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[2][0].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[2][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[2][2].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[2][3].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[3][0].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[3][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[3][2].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[3][3].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[4][0].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[4][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[4][2].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[4][3].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[5][0].Set_Equal_To(5.01426509658179920E-01);
    phi_Val[5][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[5][2].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[5][3].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[6][0].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[6][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[6][2].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[6][3].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[7][0].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[7][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[7][2].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[7][3].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[8][0].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[8][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[8][2].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[8][3].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[9][0].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[9][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[9][2].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[9][3].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[10][0].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[10][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[10][2].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[10][3].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[11][0].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[11][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[11][2].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[11][3].Set_Equal_To(5.31450498448160016E-02);


// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
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

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* evaluate basis functions (no derivatives!) on the local reference element. */
void SpecificFUNC::Basis_Value_n2(SCALAR BF_V[NB][NQ])
{

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12

    // get "Val" of basis functions
    SCALAR phi_Val[NQ][NB];

    phi_Val[0][0].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[0][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[0][2].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[0][3].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][0].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[1][2].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][3].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[2][0].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[2][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[2][2].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[2][3].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[3][0].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[3][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[3][2].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[3][3].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[4][0].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[4][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[4][2].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[4][3].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[5][0].Set_Equal_To(5.01426509658179920E-01);
    phi_Val[5][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[5][2].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[5][3].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[6][0].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[6][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[6][2].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[6][3].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[7][0].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[7][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[7][2].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[7][3].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[8][0].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[8][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[8][2].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[8][3].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[9][0].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[9][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[9][2].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[9][3].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[10][0].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[10][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[10][2].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[10][3].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[11][0].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[11][1].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[11][2].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[11][3].Set_Equal_To(3.10352451033785004E-01);


// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
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

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* evaluate basis functions (no derivatives!) on the local reference element. */
void SpecificFUNC::Basis_Value_p3(SCALAR BF_V[NB][NQ])
{

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12

    // get "Val" of basis functions
    SCALAR phi_Val[NQ][NB];

    phi_Val[0][0].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[0][1].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[0][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[0][3].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][0].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][1].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[1][3].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[2][0].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[2][1].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[2][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[2][3].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[3][0].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[3][1].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[3][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[3][3].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[4][0].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[4][1].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[4][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[4][3].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[5][0].Set_Equal_To(5.01426509658179920E-01);
    phi_Val[5][1].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[5][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[5][3].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[6][0].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[6][1].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[6][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[6][3].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[7][0].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[7][1].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[7][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[7][3].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[8][0].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[8][1].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[8][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[8][3].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[9][0].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[9][1].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[9][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[9][3].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[10][0].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[10][1].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[10][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[10][3].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[11][0].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[11][1].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[11][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[11][3].Set_Equal_To(3.10352451033785004E-01);


// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
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

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* evaluate basis functions (no derivatives!) on the local reference element. */
void SpecificFUNC::Basis_Value_n3(SCALAR BF_V[NB][NQ])
{

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12

    // get "Val" of basis functions
    SCALAR phi_Val[NQ][NB];

    phi_Val[0][0].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[0][1].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[0][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[0][3].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[1][0].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][1].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[1][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[1][3].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[2][0].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[2][1].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[2][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[2][3].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[3][0].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[3][1].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[3][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[3][3].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[4][0].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[4][1].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[4][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[4][3].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[5][0].Set_Equal_To(5.01426509658179920E-01);
    phi_Val[5][1].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[5][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[5][3].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[6][0].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[6][1].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[6][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[6][3].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[7][0].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[7][1].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[7][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[7][3].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[8][0].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[8][1].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[8][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[8][3].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[9][0].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[9][1].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[9][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[9][3].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[10][0].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[10][1].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[10][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[10][3].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[11][0].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[11][1].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[11][2].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[11][3].Set_Equal_To(5.31450498448160016E-02);


// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
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

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* evaluate basis functions (no derivatives!) on the local reference element. */
void SpecificFUNC::Basis_Value_p4(SCALAR BF_V[NB][NQ])
{

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12

    // get "Val" of basis functions
    SCALAR phi_Val[NQ][NB];

    phi_Val[0][0].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[0][1].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[0][2].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[0][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[1][0].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][1].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[1][2].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[2][0].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[2][1].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[2][2].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[2][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[3][0].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[3][1].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[3][2].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[3][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[4][0].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[4][1].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[4][2].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[4][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[5][0].Set_Equal_To(5.01426509658179920E-01);
    phi_Val[5][1].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[5][2].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[5][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[6][0].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[6][1].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[6][2].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[6][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[7][0].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[7][1].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[7][2].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[7][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[8][0].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[8][1].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[8][2].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[8][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[9][0].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[9][1].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[9][2].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[9][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[10][0].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[10][1].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[10][2].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[10][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[11][0].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[11][1].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[11][2].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[11][3].Set_Equal_To(0.00000000000000000E+00);


// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
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

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* evaluate basis functions (no derivatives!) on the local reference element. */
void SpecificFUNC::Basis_Value_n4(SCALAR BF_V[NB][NQ])
{

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, lagrange_deg1_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 2
// geometric dimension = 3
// Number of Quadrature Points = 12

    // get "Val" of basis functions
    SCALAR phi_Val[NQ][NB];

    phi_Val[0][0].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[0][1].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[0][2].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[0][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[1][0].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][1].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[1][2].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[1][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[2][0].Set_Equal_To(8.73821971016995991E-01);
    phi_Val[2][1].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[2][2].Set_Equal_To(6.30890144915020046E-02);
    phi_Val[2][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[3][0].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[3][1].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[3][2].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[3][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[4][0].Set_Equal_To(2.49286745170910956E-01);
    phi_Val[4][1].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[4][2].Set_Equal_To(5.01426509658179032E-01);
    phi_Val[4][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[5][0].Set_Equal_To(5.01426509658179920E-01);
    phi_Val[5][1].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[5][2].Set_Equal_To(2.49286745170910012E-01);
    phi_Val[5][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[6][0].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[6][1].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[6][2].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[6][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[7][0].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[7][1].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[7][2].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[7][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[8][0].Set_Equal_To(5.31450498448159947E-02);
    phi_Val[8][1].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[8][2].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[8][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[9][0].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[9][1].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[9][2].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[9][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[10][0].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[10][1].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[10][2].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[10][3].Set_Equal_To(0.00000000000000000E+00);
    phi_Val[11][0].Set_Equal_To(6.36502499121399001E-01);
    phi_Val[11][1].Set_Equal_To(5.31450498448160016E-02);
    phi_Val[11][2].Set_Equal_To(3.10352451033785004E-01);
    phi_Val[11][3].Set_Equal_To(0.00000000000000000E+00);


// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {8.73821971016995991E-01, 6.30890144915020046E-02}, \
//     {6.30890144915020046E-02, 8.73821971016995991E-01}, \
//     {6.30890144915020046E-02, 6.30890144915020046E-02}, \
//     {5.01426509658179032E-01, 2.49286745170910012E-01}, \
//     {2.49286745170910012E-01, 5.01426509658179032E-01}, \
//     {2.49286745170910012E-01, 2.49286745170910012E-01}, \
//     {6.36502499121399001E-01, 3.10352451033785004E-01}, \
//     {6.36502499121399001E-01, 5.31450498448160016E-02}, \
//     {3.10352451033785004E-01, 6.36502499121399001E-01}, \
//     {3.10352451033785004E-01, 5.31450498448160016E-02}, \
//     {5.31450498448160016E-02, 6.36502499121399001E-01}, \
//     {5.31450498448160016E-02, 3.10352451033785004E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     2.54224531851034996E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     5.83931378631894968E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02, \
//     4.14255378091870005E-02  \
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
