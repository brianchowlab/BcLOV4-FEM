/*
============================================================================================
   Header file for a C++ Class that contains methods for generic finite element assembly.

   NOTE: portions of this code are automatically generated!

   Copyright (c) 06-20-2012,  Shawn W. Walker
============================================================================================
*/

/*------------ BEGIN: Auto Generate ------------*/
// set the number of FEM matrices to assemble
#define NUM_FEM_MAT    8

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
    const Data_Type_my_f_restricted_to_Sigma* Get_my_f_restricted_to_Sigma_ptr() const { return &(my_f_restricted_to_Sigma); }
    const Data_Type_old_soln_restricted_to_Sigma* Get_old_soln_restricted_to_Sigma_ptr() const { return &(old_soln_restricted_to_Sigma); }

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
    Base_Body_Force_Matrix_Data_Type*    Base_Matrix_Body_Force_Matrix;
    Base_Curv_Matrix_Data_Type*    Base_Matrix_Curv_Matrix;
    Base_Mass_Matrix_Data_Type*    Base_Matrix_Mass_Matrix;
    Base_Small_Matrix_Data_Type*    Base_Matrix_Small_Matrix;
    Base_Stiff_Matrix_Data_Type*    Base_Matrix_Stiff_Matrix;
    Base_Tangent_Matrix_Data_Type*    Base_Matrix_Tangent_Matrix;
    Base_Vector_Mass_Matrix_Data_Type*    Base_Matrix_Vector_Mass_Matrix;
    Base_Weighted_Mass_Matrix_Data_Type*    Base_Matrix_Weighted_Mass_Matrix;

    Block_Assemble_Body_Force_Matrix_Data_Type*    Block_Assemble_Matrix_Body_Force_Matrix;
    PTR_TO_SPARSE    Sparse_Data_Body_Force_Matrix;
    Block_Assemble_Curv_Matrix_Data_Type*    Block_Assemble_Matrix_Curv_Matrix;
    PTR_TO_SPARSE    Sparse_Data_Curv_Matrix;
    Block_Assemble_Mass_Matrix_Data_Type*    Block_Assemble_Matrix_Mass_Matrix;
    PTR_TO_SPARSE    Sparse_Data_Mass_Matrix;
    Block_Assemble_Small_Matrix_Data_Type*    Block_Assemble_Matrix_Small_Matrix;
    PTR_TO_SPARSE    Sparse_Data_Small_Matrix;
    Block_Assemble_Stiff_Matrix_Data_Type*    Block_Assemble_Matrix_Stiff_Matrix;
    PTR_TO_SPARSE    Sparse_Data_Stiff_Matrix;
    Block_Assemble_Tangent_Matrix_Data_Type*    Block_Assemble_Matrix_Tangent_Matrix;
    PTR_TO_SPARSE    Sparse_Data_Tangent_Matrix;
    Block_Assemble_Vector_Mass_Matrix_Data_Type*    Block_Assemble_Matrix_Vector_Mass_Matrix;
    PTR_TO_SPARSE    Sparse_Data_Vector_Mass_Matrix;
    Block_Assemble_Weighted_Mass_Matrix_Data_Type*    Block_Assemble_Matrix_Weighted_Mass_Matrix;
    PTR_TO_SPARSE    Sparse_Data_Weighted_Mass_Matrix;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing and evaulating FEM (local) basis functions
    Data_Type_Scalar_P2_phi_restricted_to_Sigma      Scalar_P2_phi_restricted_to_Sigma;
    Data_Type_Vector_P1_phi_restricted_to_Sigma      Vector_P1_phi_restricted_to_Sigma;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing CONSTANT basis functions
    Data_Type_CONST_ONE_phi      Body_Force_Matrix_Sigma_col_constant_phi;
    Data_Type_CONST_ONE_phi      Curv_Matrix_Sigma_col_constant_phi;
    Data_Type_CONST_ONE_phi      Small_Matrix_Sigma_col_constant_phi;
    Data_Type_CONST_ONE_phi      Small_Matrix_Sigma_row_constant_phi;
    Data_Type_CONST_ONE_phi      Tangent_Matrix_Sigma_col_constant_phi;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing external functions
    Data_Type_my_f_restricted_to_Sigma      my_f_restricted_to_Sigma;
    Data_Type_old_soln_restricted_to_Sigma      old_soln_restricted_to_Sigma;
    /*------------   END: Auto Generate ------------*/

};

/***/