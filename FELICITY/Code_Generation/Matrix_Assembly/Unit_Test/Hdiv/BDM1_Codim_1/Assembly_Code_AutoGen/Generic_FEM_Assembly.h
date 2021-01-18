/*
============================================================================================
   Header file for a C++ Class that contains methods for generic finite element assembly.

   NOTE: portions of this code are automatically generated!

   Copyright (c) 06-20-2012,  Shawn W. Walker
============================================================================================
*/

/*------------ BEGIN: Auto Generate ------------*/
// set the number of FEM matrices to assemble
#define NUM_FEM_MAT    1

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
    CLASS_Domain_Gamma_embedded_in_Box_restricted_to_Gamma    Domain_Gamma_embedded_in_Box_restricted_to_Gamma;
    CLASS_Domain_Omega_embedded_in_Box_restricted_to_Gamma    Domain_Omega_embedded_in_Box_restricted_to_Gamma;

    // mesh geometry classes (can be higher order)
    CLASS_geom_Gamma_embedded_in_Box_restricted_to_Gamma   geom_Gamma_embedded_in_Box_restricted_to_Gamma;
    CLASS_geom_Omega_embedded_in_Box_restricted_to_Gamma   geom_Omega_embedded_in_Box_restricted_to_Gamma;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // pointers to each specific FEM matrix
    Base_NV_Data_Type*    Base_Matrix_NV;

    Block_Assemble_NV_Data_Type*    Block_Assemble_Matrix_NV;
    PTR_TO_SPARSE    Sparse_Data_NV;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing and evaulating FEM (local) basis functions
    Data_Type_BDM1_phi_restricted_to_Gamma      BDM1_phi_restricted_to_Gamma;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing CONSTANT basis functions
    Data_Type_CONST_ONE_phi      NV_Gamma_col_constant_phi;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing external functions
    /*------------   END: Auto Generate ------------*/

};

/***/