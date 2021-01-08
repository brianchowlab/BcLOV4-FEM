/*
============================================================================================
   Header file for a C++ Class that contains methods for generic finite element assembly.

   NOTE: portions of this code are automatically generated!

   Copyright (c) 06-20-2012,  Shawn W. Walker
============================================================================================
*/

/*------------ BEGIN: Auto Generate ------------*/
// set the number of FEM matrices to assemble
#define NUM_FEM_MAT    6

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
    const Data_Type_c0_restricted_to_Omega* Get_c0_restricted_to_Omega_ptr() const { return &(c0_restricted_to_Omega); }
    const Data_Type_u0_restricted_to_Omega* Get_u0_restricted_to_Omega_ptr() const { return &(u0_restricted_to_Omega); }

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
    Base_Elliptic_Form_Data_Type*    Base_Matrix_Elliptic_Form;
    Base_Lin_A_Data_Type*    Base_Matrix_Lin_A;
    Base_Lin_B_Data_Type*    Base_Matrix_Lin_B;
    Base_Mixed_Form_A_Data_Type*    Base_Matrix_Mixed_Form_A;
    Base_Mixed_Form_B_Data_Type*    Base_Matrix_Mixed_Form_B;
    Base_Simple_Form_Data_Type*    Base_Matrix_Simple_Form;

    Block_Assemble_Elliptic_Form_Data_Type*    Block_Assemble_Matrix_Elliptic_Form;
    PTR_TO_SPARSE    Sparse_Data_Elliptic_Form;
    Block_Assemble_Lin_A_Data_Type*    Block_Assemble_Matrix_Lin_A;
    PTR_TO_SPARSE    Sparse_Data_Lin_A;
    Block_Assemble_Lin_B_Data_Type*    Block_Assemble_Matrix_Lin_B;
    PTR_TO_SPARSE    Sparse_Data_Lin_B;
    Block_Assemble_Mixed_Form_A_Data_Type*    Block_Assemble_Matrix_Mixed_Form_A;
    PTR_TO_SPARSE    Sparse_Data_Mixed_Form_A;
    Block_Assemble_Mixed_Form_B_Data_Type*    Block_Assemble_Matrix_Mixed_Form_B;
    PTR_TO_SPARSE    Sparse_Data_Mixed_Form_B;
    Block_Assemble_Simple_Form_Data_Type*    Block_Assemble_Matrix_Simple_Form;
    PTR_TO_SPARSE    Sparse_Data_Simple_Form;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing and evaulating FEM (local) basis functions
    Data_Type_Const_Space_phi_restricted_to_Omega      Const_Space_phi_restricted_to_Omega;
    Data_Type_Const_Space_Omega_TupleSize_2_1_phi_restricted_to_Omega      Const_Space_Omega_TupleSize_2_1_phi_restricted_to_Omega;
    Data_Type_Vector_P1_phi_restricted_to_Omega      Vector_P1_phi_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing CONSTANT basis functions
    Data_Type_CONST_ONE_phi      Lin_A_Omega_col_constant_phi;
    Data_Type_CONST_ONE_phi      Lin_B_Omega_col_constant_phi;
    Data_Type_CONST_ONE_phi      Lin_B_Omega_row_constant_phi;
    Data_Type_CONST_ONE_phi      Mixed_Form_A_Omega_col_constant_phi;
    Data_Type_CONST_ONE_phi      Mixed_Form_B_Omega_row_constant_phi;
    Data_Type_CONST_ONE_phi      Simple_Form_Omega_col_constant_phi;
    Data_Type_CONST_ONE_phi      Simple_Form_Omega_row_constant_phi;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing external functions
    Data_Type_c0_restricted_to_Omega      c0_restricted_to_Omega;
    Data_Type_u0_restricted_to_Omega      u0_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

};

/***/
