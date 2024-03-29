/*
============================================================================================
   Header file for a C++ Class that contains methods for generic finite element assembly.

   NOTE: portions of this code are automatically generated!

   Copyright (c) 06-20-2012,  Shawn W. Walker
============================================================================================
*/

/*------------ BEGIN: Auto Generate ------------*/
// set the number of FEM matrices to assemble
#define NUM_FEM_MAT    3

// In MATLAB, the output FEM matrix entry list should look like:
//            FEM.MAT
//               .Type
//
// Here, we define the strings that makes these variable names
#define OUT_MAT_str      "MAT"
#define OUT_FEM_NAME_str "Type"
/*------------   END: Auto Generate ------------*/

/***************************************************************************************/
/*** C++ class ***/
class Generic_FEM_Assembly
{
public:
    //Generic_FEM_Assembly (); // constructor
    Generic_FEM_Assembly (const mxArray *[], const mxArray *); // constructor
    ~Generic_FEM_Assembly (); // DE-structor

    /*------------ BEGIN: Auto Generate ------------*/
    // create access routines
    const Data_Type_Sca_restricted_to_Sigma* Get_Sca_restricted_to_Sigma_ptr() const { return &(Sca_restricted_to_Sigma); }
    const Data_Type_Vec_restricted_to_Sigma* Get_Vec_restricted_to_Sigma_ptr() const { return &(Vec_restricted_to_Sigma); }

    void Setup_Data (const mxArray*[]);
    void Assemble_Matrices ();
    void Output_Matrices (mxArray*[]);
    void Init_Output_Matrices (mxArray*[]);
    void Output_Matrix (mwIndex, mxArray*, mxArray*, mxArray*);
    void Access_Previous_FEM_Matrix (const mxArray*, const char*, const int&, PTR_TO_SPARSE&);
    void Read_Sparse_Ptr (const mxArray*, const int&, PTR_TO_SPARSE&);
    void Clear_Sparse_Ptr (PTR_TO_SPARSE&);

private:
    // these variables are defined from inputs coming from MATLAB

    // classes for (sub)domain(s) and topological entities
    CLASS_Domain_Sigma_embedded_in_Sigma_restricted_to_Sigma    Domain_Sigma_embedded_in_Sigma_restricted_to_Sigma;

    // mesh geometry classes (can be higher order)
    CLASS_geom_Sigma_embedded_in_Sigma_restricted_to_Sigma   geom_Sigma_embedded_in_Sigma_restricted_to_Sigma;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // pointers to each specific FEM matrix
    Base_Small_Matrix_Data_Type*    Base_Matrix_Small_Matrix;
    Base_Vector_Hess_Matrix_Data_Type*    Base_Matrix_Vector_Hess_Matrix;
    Base_dsds_Matrix_Data_Type*    Base_Matrix_dsds_Matrix;

    Block_Assemble_Small_Matrix_Data_Type*    Block_Assemble_Matrix_Small_Matrix;
    PTR_TO_SPARSE    Sparse_Data_Small_Matrix;
    Block_Assemble_Vector_Hess_Matrix_Data_Type*    Block_Assemble_Matrix_Vector_Hess_Matrix;
    PTR_TO_SPARSE    Sparse_Data_Vector_Hess_Matrix;
    Block_Assemble_dsds_Matrix_Data_Type*    Block_Assemble_Matrix_dsds_Matrix;
    PTR_TO_SPARSE    Sparse_Data_dsds_Matrix;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing and evaulating FEM (local) basis functions
    Data_Type_Scalar_P2_phi_restricted_to_Sigma      Scalar_P2_phi_restricted_to_Sigma;
    Data_Type_Vector_P2_phi_restricted_to_Sigma      Vector_P2_phi_restricted_to_Sigma;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing CONSTANT basis functions
    Data_Type_CONST_ONE_phi      Small_Matrix_Sigma_col_constant_phi;
    Data_Type_CONST_ONE_phi      Small_Matrix_Sigma_row_constant_phi;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing external functions
    Data_Type_Sca_restricted_to_Sigma      Sca_restricted_to_Sigma;
    Data_Type_Vec_restricted_to_Sigma      Vec_restricted_to_Sigma;
    /*------------   END: Auto Generate ------------*/

};

/***/
