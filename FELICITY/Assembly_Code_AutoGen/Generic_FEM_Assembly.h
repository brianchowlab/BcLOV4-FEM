/*
============================================================================================
   Header file for a C++ Class that contains methods for generic finite element assembly.

   NOTE: portions of this code are automatically generated!

   Copyright (c) 06-20-2012,  Shawn W. Walker
============================================================================================
*/

/*------------ BEGIN: Auto Generate ------------*/
// set the number of FEM matrices to assemble
#define NUM_FEM_MAT    2

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
    const Data_Type_u_M_h_coef_restricted_to_dOmega* Get_u_M_h_coef_restricted_to_dOmega_ptr() const { return &(u_M_h_coef_restricted_to_dOmega); }
    const Data_Type_v_M_h_coef_restricted_to_dOmega* Get_v_M_h_coef_restricted_to_dOmega_ptr() const { return &(v_M_h_coef_restricted_to_dOmega); }

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
    CLASS_Domain_Omega_embedded_in_Omega_restricted_to_dOmega    Domain_Omega_embedded_in_Omega_restricted_to_dOmega;
    CLASS_Domain_dOmega_embedded_in_Omega_restricted_to_dOmega    Domain_dOmega_embedded_in_Omega_restricted_to_dOmega;

    // mesh geometry classes (can be higher order)
    CLASS_geom_Omega_embedded_in_Omega_restricted_to_dOmega   geom_Omega_embedded_in_Omega_restricted_to_dOmega;
    CLASS_geom_dOmega_embedded_in_Omega_restricted_to_dOmega   geom_dOmega_embedded_in_Omega_restricted_to_dOmega;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // pointers to each specific FEM matrix
    Base_B_phi_Data_Type*    Base_Matrix_B_phi;
    Base_B_psi_Data_Type*    Base_Matrix_B_psi;

    Block_Assemble_B_phi_Data_Type*    Block_Assemble_Matrix_B_phi;
    PTR_TO_SPARSE    Sparse_Data_B_phi;
    Block_Assemble_B_psi_Data_Type*    Block_Assemble_Matrix_B_psi;
    PTR_TO_SPARSE    Sparse_Data_B_psi;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing and evaulating FEM (local) basis functions
    Data_Type_C_h_phi_restricted_to_dOmega      C_h_phi_restricted_to_dOmega;
    Data_Type_M_h_phi_restricted_to_dOmega      M_h_phi_restricted_to_dOmega;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing CONSTANT basis functions
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing external functions
    Data_Type_u_M_h_coef_restricted_to_dOmega      u_M_h_coef_restricted_to_dOmega;
    Data_Type_v_M_h_coef_restricted_to_dOmega      v_M_h_coef_restricted_to_dOmega;
    /*------------   END: Auto Generate ------------*/

};

/***/