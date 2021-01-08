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
    const Data_Type_BC_Out_restricted_to_Outlet* Get_BC_Out_restricted_to_Outlet_ptr() const { return &(BC_Out_restricted_to_Outlet); }

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
    CLASS_Domain_Omega_embedded_in_Omega_restricted_to_Outlet    Domain_Omega_embedded_in_Omega_restricted_to_Outlet;
    CLASS_Domain_Outlet_embedded_in_Omega_restricted_to_Outlet    Domain_Outlet_embedded_in_Omega_restricted_to_Outlet;

    // mesh geometry classes (can be higher order)
    CLASS_geom_Omega_embedded_in_Omega_restricted_to_Omega   geom_Omega_embedded_in_Omega_restricted_to_Omega;
    CLASS_geom_Omega_embedded_in_Omega_restricted_to_Outlet   geom_Omega_embedded_in_Omega_restricted_to_Outlet;
    CLASS_geom_Outlet_embedded_in_Omega_restricted_to_Outlet   geom_Outlet_embedded_in_Omega_restricted_to_Outlet;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // pointers to each specific FEM matrix
    Base_BC_Matrix_Data_Type*    Base_Matrix_BC_Matrix;
    Base_Div_Pressure_Data_Type*    Base_Matrix_Div_Pressure;
    Base_Stress_Matrix_Data_Type*    Base_Matrix_Stress_Matrix;

    Block_Assemble_BC_Matrix_Data_Type*    Block_Assemble_Matrix_BC_Matrix;
    PTR_TO_SPARSE    Sparse_Data_BC_Matrix;
    Block_Assemble_Div_Pressure_Data_Type*    Block_Assemble_Matrix_Div_Pressure;
    PTR_TO_SPARSE    Sparse_Data_Div_Pressure;
    Block_Assemble_Stress_Matrix_Data_Type*    Block_Assemble_Matrix_Stress_Matrix;
    PTR_TO_SPARSE    Sparse_Data_Stress_Matrix;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing and evaulating FEM (local) basis functions
    Data_Type_Scalar_P1_phi_restricted_to_Omega      Scalar_P1_phi_restricted_to_Omega;
    Data_Type_Vector_P2_phi_restricted_to_Omega      Vector_P2_phi_restricted_to_Omega;
    Data_Type_Vector_P2_phi_restricted_to_Outlet      Vector_P2_phi_restricted_to_Outlet;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing CONSTANT basis functions
    Data_Type_CONST_ONE_phi      BC_Matrix_Outlet_col_constant_phi;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing external functions
    Data_Type_BC_Out_restricted_to_Outlet      BC_Out_restricted_to_Outlet;
    /*------------   END: Auto Generate ------------*/

};

/***/
