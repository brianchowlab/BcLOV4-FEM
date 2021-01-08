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
// \int_{Omega}  uu_v1_t1_hess11*vv_v1_t1_hess11 + uu_v1_t1_hess12*vv_v1_t1_hess12 + uu_v1_t1_hess13*vv_v1_t1_hess13 + uu_v1_t1_hess21*vv_v1_t1_hess21 + uu_v1_t1_hess22*vv_v1_t1_hess22 + uu_v1_t1_hess23*vv_v1_t1_hess23 + uu_v1_t1_hess31*vv_v1_t1_hess31 + uu_v1_t1_hess32*vv_v1_t1_hess32 + uu_v1_t1_hess33*vv_v1_t1_hess33 + uu_v1_t2_hess11*vv_v1_t2_hess11 + uu_v1_t2_hess12*vv_v1_t2_hess12 + uu_v1_t2_hess13*vv_v1_t2_hess13 + uu_v1_t2_hess21*vv_v1_t2_hess21 + uu_v1_t2_hess22*vv_v1_t2_hess22 + uu_v1_t2_hess23*vv_v1_t2_hess23 + uu_v1_t2_hess31*vv_v1_t2_hess31 + uu_v1_t2_hess32*vv_v1_t2_hess32 + uu_v1_t2_hess33*vv_v1_t2_hess33 + uu_v1_t3_hess11*vv_v1_t3_hess11 + uu_v1_t3_hess12*vv_v1_t3_hess12 + uu_v1_t3_hess13*vv_v1_t3_hess13 + uu_v1_t3_hess21*vv_v1_t3_hess21 + uu_v1_t3_hess22*vv_v1_t3_hess22 + uu_v1_t3_hess23*vv_v1_t3_hess23 + uu_v1_t3_hess31*vv_v1_t3_hess31 + uu_v1_t3_hess32*vv_v1_t3_hess32 + uu_v1_t3_hess33*vv_v1_t3_hess33
//
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
// define the name of the FEM matrix (should be the same as the filename of this file)
#define SpecificFEM        Base_Vector_Hess_Matrix_Data_Type
#define SpecificFEM_str   "Vector_Hess_Matrix"

// the row function space is Type = CG, Name = "lagrange_deg2_dim3"

// set the number of cartesian tuple components (m*n) = 3 * 1
#define ROW_NC  3
// NOTE: the (i,j) tuple component is accessed by the linear index k = i + (j-1)*m
// set the number of basis functions on each element
#define ROW_NB  10

// the col function space is Type = CG, Name = "lagrange_deg2_dim3"

// set the number of cartesian tuple components (m*n) = 3 * 1
#define COL_NC  3
// NOTE: the (i,j) tuple component is accessed by the linear index k = i + (j-1)*m
// set the number of basis functions on each element
#define COL_NB  10
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
    const Data_Type_Vector_P2_phi_restricted_to_Omega*  Vector_P2_phi_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // pointers for accessing coefficient functions
    const Data_Type_Vec_restricted_to_Omega*  Vec_restricted_to_Omega;
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
    Num_Sub_Matrices = 3;

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
    Row_Shift[1] = 1*row_basis_func->Num_Nodes;
    Row_Shift[2] = 2*row_basis_func->Num_Nodes;
    Col_Shift[0] = 0*col_basis_func->Num_Nodes;
    Col_Shift[1] = 1*col_basis_func->Num_Nodes;
    Col_Shift[2] = 2*col_basis_func->Num_Nodes;
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
        for (unsigned int i = j; i < ROW_NB; i++)
            {
            double  A0_value = 0.0; // initialize
            for (unsigned int qp = 0; qp < NQ; qp++)
                {
                const double  integrand_0 = Vector_P2_phi_restricted_to_Omega->Func_f_Hess[j][qp].m[0][0]*Vector_P2_phi_restricted_to_Omega->Func_f_Hess[i][qp].m[0][0]+Vector_P2_phi_restricted_to_Omega->Func_f_Hess[j][qp].m[0][1]*Vector_P2_phi_restricted_to_Omega->Func_f_Hess[i][qp].m[0][1]+Vector_P2_phi_restricted_to_Omega->Func_f_Hess[j][qp].m[0][2]*Vector_P2_phi_restricted_to_Omega->Func_f_Hess[i][qp].m[0][2]+Vector_P2_phi_restricted_to_Omega->Func_f_Hess[j][qp].m[1][0]*Vector_P2_phi_restricted_to_Omega->Func_f_Hess[i][qp].m[1][0]+Vector_P2_phi_restricted_to_Omega->Func_f_Hess[j][qp].m[1][1]*Vector_P2_phi_restricted_to_Omega->Func_f_Hess[i][qp].m[1][1]+Vector_P2_phi_restricted_to_Omega->Func_f_Hess[j][qp].m[1][2]*Vector_P2_phi_restricted_to_Omega->Func_f_Hess[i][qp].m[1][2]+Vector_P2_phi_restricted_to_Omega->Func_f_Hess[j][qp].m[2][0]*Vector_P2_phi_restricted_to_Omega->Func_f_Hess[i][qp].m[2][0]+Vector_P2_phi_restricted_to_Omega->Func_f_Hess[j][qp].m[2][1]*Vector_P2_phi_restricted_to_Omega->Func_f_Hess[i][qp].m[2][1]+Vector_P2_phi_restricted_to_Omega->Func_f_Hess[j][qp].m[2][2]*Vector_P2_phi_restricted_to_Omega->Func_f_Hess[i][qp].m[2][2];
                A0_value += integrand_0 * Mesh.Map_Det_Jac_w_Weight[qp].a;
                }
            FE_Tensor_0[j*ROW_NB + i] = A0_value;
            }
        }

    // Copy the lower triangular entries to the upper triangular part (by symmetry)
    for (unsigned int j = 0; j < COL_NB; j++)
        {
        for (unsigned int i = j+1; i < ROW_NB; i++)
            {
            FE_Tensor_0[i*ROW_NB + j] = FE_Tensor_0[j*ROW_NB + i];
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

