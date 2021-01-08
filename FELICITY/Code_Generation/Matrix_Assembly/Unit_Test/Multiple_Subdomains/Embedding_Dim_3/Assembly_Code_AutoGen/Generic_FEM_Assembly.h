/*
============================================================================================
   Header file for a C++ Class that contains methods for generic finite element assembly.

   NOTE: portions of this code are automatically generated!

   Copyright (c) 06-20-2012,  Shawn W. Walker
============================================================================================
*/

/*------------ BEGIN: Auto Generate ------------*/
// set the number of FEM matrices to assemble
#define NUM_FEM_MAT    5

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
    CLASS_Domain_Gamma_embedded_in_Global_restricted_to_Gamma    Domain_Gamma_embedded_in_Global_restricted_to_Gamma;
    CLASS_Domain_Gamma_embedded_in_Global_restricted_to_Sigma    Domain_Gamma_embedded_in_Global_restricted_to_Sigma;
    CLASS_Domain_Omega_1_embedded_in_Global_restricted_to_Omega_1    Domain_Omega_1_embedded_in_Global_restricted_to_Omega_1;
    CLASS_Domain_Omega_2_embedded_in_Global_restricted_to_Omega_2    Domain_Omega_2_embedded_in_Global_restricted_to_Omega_2;
    CLASS_Domain_Omega_embedded_in_Global_restricted_to_Gamma    Domain_Omega_embedded_in_Global_restricted_to_Gamma;
    CLASS_Domain_Omega_embedded_in_Global_restricted_to_Omega    Domain_Omega_embedded_in_Global_restricted_to_Omega;
    CLASS_Domain_Omega_embedded_in_Global_restricted_to_Omega_1    Domain_Omega_embedded_in_Global_restricted_to_Omega_1;
    CLASS_Domain_Omega_embedded_in_Global_restricted_to_Omega_2    Domain_Omega_embedded_in_Global_restricted_to_Omega_2;
    CLASS_Domain_Omega_embedded_in_Global_restricted_to_Sigma    Domain_Omega_embedded_in_Global_restricted_to_Sigma;
    CLASS_Domain_Sigma_embedded_in_Global_restricted_to_Sigma    Domain_Sigma_embedded_in_Global_restricted_to_Sigma;

    // mesh geometry classes (can be higher order)
    CLASS_geom_Gamma_embedded_in_Global_restricted_to_Gamma   geom_Gamma_embedded_in_Global_restricted_to_Gamma;
    CLASS_geom_Gamma_embedded_in_Global_restricted_to_Sigma   geom_Gamma_embedded_in_Global_restricted_to_Sigma;
    CLASS_geom_Omega_1_embedded_in_Global_restricted_to_Omega_1   geom_Omega_1_embedded_in_Global_restricted_to_Omega_1;
    CLASS_geom_Omega_2_embedded_in_Global_restricted_to_Omega_2   geom_Omega_2_embedded_in_Global_restricted_to_Omega_2;
    CLASS_geom_Omega_embedded_in_Global_restricted_to_Gamma   geom_Omega_embedded_in_Global_restricted_to_Gamma;
    CLASS_geom_Omega_embedded_in_Global_restricted_to_Omega   geom_Omega_embedded_in_Global_restricted_to_Omega;
    CLASS_geom_Omega_embedded_in_Global_restricted_to_Omega_1   geom_Omega_embedded_in_Global_restricted_to_Omega_1;
    CLASS_geom_Omega_embedded_in_Global_restricted_to_Omega_2   geom_Omega_embedded_in_Global_restricted_to_Omega_2;
    CLASS_geom_Omega_embedded_in_Global_restricted_to_Sigma   geom_Omega_embedded_in_Global_restricted_to_Sigma;
    CLASS_geom_Sigma_embedded_in_Global_restricted_to_Sigma   geom_Sigma_embedded_in_Global_restricted_to_Sigma;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // pointers to each specific FEM matrix
    Base_Bdy_Mixed_Data_Type*    Base_Matrix_Bdy_Mixed;
    Base_C1_Data_Type*    Base_Matrix_C1;
    Base_Edge_Mixed_Data_Type*    Base_Matrix_Edge_Mixed;
    Base_Stiff_Matrix_Data_Type*    Base_Matrix_Stiff_Matrix;
    Base_Weight_Mass_Data_Type*    Base_Matrix_Weight_Mass;

    Block_Assemble_Bdy_Mixed_Data_Type*    Block_Assemble_Matrix_Bdy_Mixed;
    PTR_TO_SPARSE    Sparse_Data_Bdy_Mixed;
    Block_Assemble_C1_Data_Type*    Block_Assemble_Matrix_C1;
    PTR_TO_SPARSE    Sparse_Data_C1;
    Block_Assemble_Edge_Mixed_Data_Type*    Block_Assemble_Matrix_Edge_Mixed;
    PTR_TO_SPARSE    Sparse_Data_Edge_Mixed;
    Block_Assemble_Stiff_Matrix_Data_Type*    Block_Assemble_Matrix_Stiff_Matrix;
    PTR_TO_SPARSE    Sparse_Data_Stiff_Matrix;
    Block_Assemble_Weight_Mass_Data_Type*    Block_Assemble_Matrix_Weight_Mass;
    PTR_TO_SPARSE    Sparse_Data_Weight_Mass;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing and evaulating FEM (local) basis functions
    Data_Type_M_Space_phi_restricted_to_Gamma      M_Space_phi_restricted_to_Gamma;
    Data_Type_U_Space_phi_restricted_to_Gamma      U_Space_phi_restricted_to_Gamma;
    Data_Type_U_Space_phi_restricted_to_Omega      U_Space_phi_restricted_to_Omega;
    Data_Type_U_Space_phi_restricted_to_Omega_1      U_Space_phi_restricted_to_Omega_1;
    Data_Type_U_Space_phi_restricted_to_Omega_2      U_Space_phi_restricted_to_Omega_2;
    Data_Type_E_Space_phi_restricted_to_Sigma      E_Space_phi_restricted_to_Sigma;
    Data_Type_M_Space_phi_restricted_to_Sigma      M_Space_phi_restricted_to_Sigma;
    Data_Type_U_Space_phi_restricted_to_Sigma      U_Space_phi_restricted_to_Sigma;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing CONSTANT basis functions
    Data_Type_CONST_ONE_phi      C1_Sigma_col_constant_phi;
    Data_Type_CONST_ONE_phi      C1_Sigma_row_constant_phi;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing external functions
    Data_Type_old_soln_restricted_to_Sigma      old_soln_restricted_to_Sigma;
    /*------------   END: Auto Generate ------------*/

};

/***/
