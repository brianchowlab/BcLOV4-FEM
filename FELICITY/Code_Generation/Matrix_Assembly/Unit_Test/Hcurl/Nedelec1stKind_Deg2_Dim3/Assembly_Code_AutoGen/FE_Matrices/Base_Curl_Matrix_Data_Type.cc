/*
============================================================================================
   This file contains an implementation of a derived C++ Class from the abstract base class
   in 'Local_FE_Matrix.cc'.

   NOTE: portions of this code are automatically generated!

   Copyright (c) 06-14-2016,  Shawn W. Walker
============================================================================================
*/

/*------------ BEGIN: Auto Generate ------------*/
// FE matrix is:
//
// \int_{Omega}  p_v1_t1*(vv_v1_t1_curl + vv_v2_t1_curl + vv_v3_t1_curl)
//
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
// define the name of the FEM matrix (should be the same as the filename of this file)
#define SpecificFEM        Base_Curl_Matrix_Data_Type
#define SpecificFEM_str   "Curl_Matrix"

// the row function space is Type = CG, Name = "nedelec_1stkind_deg2_dim3"

// set the number of cartesian tuple components (m*n) = 1 * 1
#define ROW_NC  1
// NOTE: the (i,j) tuple component is accessed by the linear index k = i + (j-1)*m
// set the number of basis functions on each element
#define ROW_NB  20

// the col function space is Type = CG, Name = "lagrange_deg1_dim3"

// set the number of cartesian tuple components (m*n) = 1 * 1
#define COL_NC  1
// NOTE: the (i,j) tuple component is accessed by the linear index k = i + (j-1)*m
// set the number of basis functions on each element
#define COL_NB  4
/*------------   END: Auto Generate ------------*/


/***************************************************************************************/
/* C++ (Specific) Local FE matrix class definition */
class SpecificFEM: public Base_FE_MATRIX_Class // derive from base class
{
public:

    // data structure for sub-matrices of the global FE matrix:
    // row/col offsets for inserting into global FE matrix
    int     Row_Shift[ROW_NC];
    int     Col_Shift[COL_NC];

    /*------------ BEGIN: Auto Generate ------------*/
    // access local mesh geometry info
    const CLASS_geom_Omega_embedded_in_Omega_restricted_to_Omega*  geom_Omega_embedded_in_Omega_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // pointers for accessing FE basis functions
    const Data_Type_Ned2_phi_restricted_to_Omega*  Ned2_phi_restricted_to_Omega;
    const Data_Type_P1_phi_restricted_to_Omega*  P1_phi_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // pointers for accessing coefficient functions
    const Data_Type_u_old_restricted_to_Omega*  u_old_restricted_to_Omega;
    const Data_Type_u_soln_restricted_to_Omega*  u_soln_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // data tables to hold local FE tensors (matrices)
    double  FE_Tensor_0[ROW_NB*COL_NB];
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // constructor
    SpecificFEM (const ABSTRACT_FEM_Function_Class*, const ABSTRACT_FEM_Function_Class*);
    ~SpecificFEM (); // destructor
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    void Tabulate_Tensor(const CLASS_geom_Omega_embedded_in_Omega_restricted_to_Omega&);
    /*------------   END: Auto Generate ------------*/
private:
};


/***************************************************************************************/
/* constructor */
SpecificFEM::SpecificFEM (const ABSTRACT_FEM_Function_Class* row_basis_func,
                          const ABSTRACT_FEM_Function_Class* col_basis_func) :
Base_FE_MATRIX_Class () // call the base class constructor
{
    // set the 'Name' of the form for this local FE matrix
    Form_Name = (char*) SpecificFEM_str; // this should be similar to the Class identifier
    // record the type of form
    Form_Type = bilinear;
    // record the number of sub-matrices (in local FE matrix)
    Num_Sub_Matrices = 1;

    // BEGIN: simple error check
    if (row_basis_func->Num_Comp!=ROW_NC)
        {
        mexPrintf("ERROR: The number of components for the row FE space is NOT correct!\n");
        mexPrintf("The FE matrix '%s' expects %d row components.\n",SpecificFEM_str,ROW_NC);
        mexPrintf("Actual number of row components is %d.\n",row_basis_func->Num_Comp);
        mexErrMsgTxt("Please report this error!\n");
        }
    if (col_basis_func->Num_Comp!=COL_NC)
        {
        mexPrintf("ERROR: The number of components for the column FE space is NOT correct!\n");
        mexPrintf("The FE matrix '%s' expects %d col components.\n",SpecificFEM_str,ROW_NC);
        mexPrintf("Actual number of col components is %d.\n",col_basis_func->Num_Comp);
        mexErrMsgTxt("Please report this error!\n");
        }
    if (row_basis_func->Num_Basis!=ROW_NB)
        {
        mexPrintf("ERROR: The number of basis functions in the row FEM space is NOT correct!\n");
        mexPrintf("The FEM matrix '%s' expects %d row basis functions.\n",SpecificFEM_str,ROW_NB);
        mexPrintf("Actual number of row components is %d.\n",row_basis_func->Num_Basis);
        mexErrMsgTxt("Please report this error!\n");
        }
    if (col_basis_func->Num_Basis!=COL_NB)
        {
        mexPrintf("ERROR: The number of basis functions in the column FEM space is NOT correct!\n");
        mexPrintf("The FEM matrix '%s' expects %d col basis functions.\n",SpecificFEM_str,COL_NB);
        mexPrintf("Actual number of row components is %d.\n",col_basis_func->Num_Basis);
        mexErrMsgTxt("Please report this error!\n");
        }
    //   END: simple error check

    // record the size of the global matrix
    global_num_row = row_basis_func->Num_Comp*row_basis_func->Num_Nodes;
    global_num_col = col_basis_func->Num_Comp*col_basis_func->Num_Nodes;

    /*------------ BEGIN: Auto Generate ------------*/
    // input information for offsetting the sub-matrices
    Row_Shift[0] = 0*row_basis_func->Num_Nodes;
    Col_Shift[0] = 0*col_basis_func->Num_Nodes;
    /*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/


/***************************************************************************************/
/* destructor: should not usually need to be modified */
SpecificFEM::~SpecificFEM ()
{
}
/***************************************************************************************/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* Tabulate the tensor for the local element contribution */
void SpecificFEM::Tabulate_Tensor(const CLASS_geom_Omega_embedded_in_Omega_restricted_to_Omega& Mesh)
{
    const unsigned int NQ = 10;

    // Compute element tensor using quadrature

    // Loop quadrature points for integral
    for (unsigned int j = 0; j < COL_NB; j++)
        {
        for (unsigned int i = 0; i < ROW_NB; i++)
            {
            double  A0_value = 0.0; // initialize
            for (unsigned int qp = 0; qp < NQ; qp++)
                {
                const double  integrand_0 = (*P1_phi_restricted_to_Omega->Func_f_Value)[j][qp].a*(Ned2_phi_restricted_to_Omega->Func_vv_Curl[i][qp].v[0]+Ned2_phi_restricted_to_Omega->Func_vv_Curl[i][qp].v[1]+Ned2_phi_restricted_to_Omega->Func_vv_Curl[i][qp].v[2]);
                A0_value += integrand_0 * Mesh.Map_Det_Jac_w_Weight[qp].a;
                }
            FE_Tensor_0[j*ROW_NB + i] = A0_value;
            }
        }
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

// remove those macros!
#undef SpecificFEM
#undef SpecificFEM_str

#undef ROW_NC
#undef ROW_NB
#undef COL_NC
#undef COL_NB

/***/

