/*
============================================================================================
   Header file for a C++ Class that contains methods for generic finite element assembly.

   NOTE: portions of this code are automatically generated!

   Copyright (c) 06-20-2012,  Shawn W. Walker
============================================================================================
*/

/*------------ BEGIN: Auto Generate ------------*/
// set the number of FEM matrices to assemble
#define NUM_FEM_MAT    9

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
    const Data_Type_u_old_restricted_to_Omega* Get_u_old_restricted_to_Omega_ptr() const { return &(u_old_restricted_to_Omega); }
    const Data_Type_u_soln_restricted_to_Omega* Get_u_soln_restricted_to_Omega_ptr() const { return &(u_soln_restricted_to_Omega); }

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
    CLASS_Domain_Omega_embedded_in_Omega_restricted_to_Omega    Domain_Omega_embedded_in_Omega_restricted_to_Omega;

    // mesh geometry classes (can be higher order)
    CLASS_geom_Omega_embedded_in_Omega_restricted_to_Omega   geom_Omega_embedded_in_Omega_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // pointers to each specific FEM matrix
    Base_CurlCurl_Matrix_Data_Type*    Base_Matrix_CurlCurl_Matrix;
    Base_Curl_Matrix_Data_Type*    Base_Matrix_Curl_Matrix;
    Base_L2_Error_Sq_Data_Type*    Base_Matrix_L2_Error_Sq;
    Base_Mass_Matrix_Data_Type*    Base_Matrix_Mass_Matrix;
    Base_Ned2_Proj_Mat_1_Data_Type*    Base_Matrix_Ned2_Proj_Mat_1;
    Base_Ned2_Proj_Mat_2_Data_Type*    Base_Matrix_Ned2_Proj_Mat_2;
    Base_Ned2_Proj_Mat_3_Data_Type*    Base_Matrix_Ned2_Proj_Mat_3;
    Base_RHS_CC_Data_Type*    Base_Matrix_RHS_CC;
    Base_Small_Matrix_Data_Type*    Base_Matrix_Small_Matrix;

    Block_Assemble_CurlCurl_Matrix_Data_Type*    Block_Assemble_Matrix_CurlCurl_Matrix;
    PTR_TO_SPARSE    Sparse_Data_CurlCurl_Matrix;
    Block_Assemble_Curl_Matrix_Data_Type*    Block_Assemble_Matrix_Curl_Matrix;
    PTR_TO_SPARSE    Sparse_Data_Curl_Matrix;
    Block_Assemble_L2_Error_Sq_Data_Type*    Block_Assemble_Matrix_L2_Error_Sq;
    PTR_TO_SPARSE    Sparse_Data_L2_Error_Sq;
    Block_Assemble_Mass_Matrix_Data_Type*    Block_Assemble_Matrix_Mass_Matrix;
    PTR_TO_SPARSE    Sparse_Data_Mass_Matrix;
    Block_Assemble_Ned2_Proj_Mat_1_Data_Type*    Block_Assemble_Matrix_Ned2_Proj_Mat_1;
    PTR_TO_SPARSE    Sparse_Data_Ned2_Proj_Mat_1;
    Block_Assemble_Ned2_Proj_Mat_2_Data_Type*    Block_Assemble_Matrix_Ned2_Proj_Mat_2;
    PTR_TO_SPARSE    Sparse_Data_Ned2_Proj_Mat_2;
    Block_Assemble_Ned2_Proj_Mat_3_Data_Type*    Block_Assemble_Matrix_Ned2_Proj_Mat_3;
    PTR_TO_SPARSE    Sparse_Data_Ned2_Proj_Mat_3;
    Block_Assemble_RHS_CC_Data_Type*    Block_Assemble_Matrix_RHS_CC;
    PTR_TO_SPARSE    Sparse_Data_RHS_CC;
    Block_Assemble_Small_Matrix_Data_Type*    Block_Assemble_Matrix_Small_Matrix;
    PTR_TO_SPARSE    Sparse_Data_Small_Matrix;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing and evaulating FEM (local) basis functions
    Data_Type_Ned2_phi_restricted_to_Omega      Ned2_phi_restricted_to_Omega;
    Data_Type_P1_phi_restricted_to_Omega      P1_phi_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing CONSTANT basis functions
    Data_Type_CONST_ONE_phi      L2_Error_Sq_Omega_col_constant_phi;
    Data_Type_CONST_ONE_phi      L2_Error_Sq_Omega_row_constant_phi;
    Data_Type_CONST_ONE_phi      Ned2_Proj_Mat_1_Omega_col_constant_phi;
    Data_Type_CONST_ONE_phi      Ned2_Proj_Mat_2_Omega_col_constant_phi;
    Data_Type_CONST_ONE_phi      Ned2_Proj_Mat_3_Omega_col_constant_phi;
    Data_Type_CONST_ONE_phi      RHS_CC_Omega_col_constant_phi;
    Data_Type_CONST_ONE_phi      Small_Matrix_Omega_col_constant_phi;
    Data_Type_CONST_ONE_phi      Small_Matrix_Omega_row_constant_phi;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing external functions
    Data_Type_u_old_restricted_to_Omega      u_old_restricted_to_Omega;
    Data_Type_u_soln_restricted_to_Omega      u_soln_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

};

/***/
