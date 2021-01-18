/*
============================================================================================
   Header file for a C++ Class that contains methods for generic finite element assembly.

   NOTE: portions of this code are automatically generated!

   Copyright (c) 06-20-2012,  Shawn W. Walker
============================================================================================
*/

/*------------ BEGIN: Auto Generate ------------*/
// set the number of FEM matrices to assemble
#define NUM_FEM_MAT    7

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
    CLASS_Domain_Gamma_embedded_in_Omega_restricted_to_Gamma    Domain_Gamma_embedded_in_Omega_restricted_to_Gamma;
    CLASS_Domain_Omega_embedded_in_Omega_restricted_to_Gamma    Domain_Omega_embedded_in_Omega_restricted_to_Gamma;
    CLASS_Domain_Omega_embedded_in_Omega_restricted_to_Omega    Domain_Omega_embedded_in_Omega_restricted_to_Omega;

    // mesh geometry classes (can be higher order)
    CLASS_geom_Gamma_embedded_in_Omega_restricted_to_Gamma   geom_Gamma_embedded_in_Omega_restricted_to_Gamma;
    CLASS_geom_Omega_embedded_in_Omega_restricted_to_Gamma   geom_Omega_embedded_in_Omega_restricted_to_Gamma;
    CLASS_geom_Omega_embedded_in_Omega_restricted_to_Omega   geom_Omega_embedded_in_Omega_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // pointers to each specific FEM matrix
    Base_A_Data_Type*    Base_Matrix_A;
    Base_B_Data_Type*    Base_Matrix_B;
    Base_C_Data_Type*    Base_Matrix_C;
    Base_D_Data_Type*    Base_Matrix_D;
    Base_K_Data_Type*    Base_Matrix_K;
    Base_M_Data_Type*    Base_Matrix_M;
    Base_chi_Data_Type*    Base_Matrix_chi;

    Block_Assemble_A_Data_Type*    Block_Assemble_Matrix_A;
    PTR_TO_SPARSE    Sparse_Data_A;
    Block_Assemble_B_Data_Type*    Block_Assemble_Matrix_B;
    PTR_TO_SPARSE    Sparse_Data_B;
    Block_Assemble_C_Data_Type*    Block_Assemble_Matrix_C;
    PTR_TO_SPARSE    Sparse_Data_C;
    Block_Assemble_D_Data_Type*    Block_Assemble_Matrix_D;
    PTR_TO_SPARSE    Sparse_Data_D;
    Block_Assemble_K_Data_Type*    Block_Assemble_Matrix_K;
    PTR_TO_SPARSE    Sparse_Data_K;
    Block_Assemble_M_Data_Type*    Block_Assemble_Matrix_M;
    PTR_TO_SPARSE    Sparse_Data_M;
    Block_Assemble_chi_Data_Type*    Block_Assemble_Matrix_chi;
    PTR_TO_SPARSE    Sparse_Data_chi;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing and evaulating FEM (local) basis functions
    Data_Type_M_h_phi_restricted_to_Gamma      M_h_phi_restricted_to_Gamma;
    Data_Type_V_h_phi_restricted_to_Gamma      V_h_phi_restricted_to_Gamma;
    Data_Type_Y_h_phi_restricted_to_Gamma      Y_h_phi_restricted_to_Gamma;
    Data_Type_G_h_phi_restricted_to_Omega      G_h_phi_restricted_to_Omega;
    Data_Type_Q_h_phi_restricted_to_Omega      Q_h_phi_restricted_to_Omega;
    Data_Type_V_h_phi_restricted_to_Omega      V_h_phi_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing CONSTANT basis functions
    Data_Type_CONST_ONE_phi      chi_Gamma_col_constant_phi;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing external functions
    /*------------   END: Auto Generate ------------*/

};

/***/