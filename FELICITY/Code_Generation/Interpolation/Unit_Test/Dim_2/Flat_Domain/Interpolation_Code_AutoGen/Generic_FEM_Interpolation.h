/*
============================================================================================
   Header file for a C++ Class that contains methods for generic finite element
   interpolation.

   NOTE: portions of this code are automatically generated!

   Copyright (c) 01-29-2013,  Shawn W. Walker
============================================================================================
*/

/*------------ BEGIN: Auto Generate ------------*/
// set the number of FEM interpolations to evaluate
#define NUM_FEM_INTERP    2

// In MATLAB, the output (INTERP) interpolation data should look like:
//            INTERP.DATA
//                  .Name
//
// Here, we define the strings that makes these variable names
#define OUT_DATA_str    "DATA"
#define OUT_NAME_str    "Name"
/*------------   END: Auto Generate ------------*/

/***************************************************************************************/
/*** C++ class ***/
class Generic_FEM_Interpolation
{
public:
    //Generic_FEM_Interpolation (); // constructor
    Generic_FEM_Interpolation (const mxArray *[]); // constructor
    ~Generic_FEM_Interpolation (); // DE-structor

    /*------------ BEGIN: Auto Generate ------------*/
    // create access routines
    const Data_Type_p_restricted_to_Omega* Get_p_restricted_to_Omega_ptr() const { return &(p_restricted_to_Omega); }
    const Data_Type_v_restricted_to_Omega* Get_v_restricted_to_Omega_ptr() const { return &(v_restricted_to_Omega); }

    void Setup_Data (const mxArray*[]);
    void Evaluate_Interpolations ();
    void Output_Interpolations (mxArray*[]);
    void Init_Output_Data (mxArray*[]);
    void Output_Single_Interpolation (mwIndex, mxArray*, mxArray*, mxArray*);

private:
    // these variables are defined from inputs coming from MATLAB

    // classes for (sub)domain(s) and topological entities
    CLASS_Domain_Omega_embedded_in_Omega_restricted_to_Omega    Domain_Omega_embedded_in_Omega_restricted_to_Omega;

    // classes for accessing interpolation points on subdomains
    Unstructured_Interpolation_Class    Omega_Interp_Data;

    // mesh geometry classes (can be higher order)
    CLASS_geom_Omega_embedded_in_Omega_restricted_to_Omega   geom_Omega_embedded_in_Omega_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // pointers to each specific FEM interpolation
    I_p*    Iobj_I_p;
    I_v*    Iobj_I_v;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing and evaulating FEM (local) basis functions
    Data_Type_Scalar_P2_phi_restricted_to_Omega      Scalar_P2_phi_restricted_to_Omega;
    Data_Type_Vector_P1_phi_restricted_to_Omega      Vector_P1_phi_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing external functions
    Data_Type_p_restricted_to_Omega      p_restricted_to_Omega;
    Data_Type_v_restricted_to_Omega      v_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

};

/***/
