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
#define NUM_FEM_INTERP    3

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
    const Data_Type_NN_restricted_to_Sigma* Get_NN_restricted_to_Sigma_ptr() const { return &(NN_restricted_to_Sigma); }

    void Setup_Data (const mxArray*[]);
    void Evaluate_Interpolations ();
    void Output_Interpolations (mxArray*[]);
    void Init_Output_Data (mxArray*[]);
    void Output_Single_Interpolation (mwIndex, mxArray*, mxArray*, mxArray*);

private:
    // these variables are defined from inputs coming from MATLAB

    // classes for (sub)domain(s) and topological entities
    CLASS_Domain_Gamma_embedded_in_Gamma_restricted_to_Sigma    Domain_Gamma_embedded_in_Gamma_restricted_to_Sigma;
    CLASS_Domain_Sigma_embedded_in_Gamma_restricted_to_Sigma    Domain_Sigma_embedded_in_Gamma_restricted_to_Sigma;

    // classes for accessing interpolation points on subdomains
    Unstructured_Interpolation_Class    Sigma_Interp_Data;

    // mesh geometry classes (can be higher order)
    CLASS_geom_Gamma_embedded_in_Gamma_restricted_to_Sigma   geom_Gamma_embedded_in_Gamma_restricted_to_Sigma;
    CLASS_geom_Sigma_embedded_in_Gamma_restricted_to_Sigma   geom_Sigma_embedded_in_Gamma_restricted_to_Sigma;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // pointers to each specific FEM interpolation
    I_N*    Iobj_I_N;
    I_N_dot_NN*    Iobj_I_N_dot_NN;
    I_N_dot_T*    Iobj_I_N_dot_T;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing and evaulating FEM (local) basis functions
    Data_Type_V_Space_phi_restricted_to_Sigma      V_Space_phi_restricted_to_Sigma;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // create variables for accessing external functions
    Data_Type_NN_restricted_to_Sigma      NN_restricted_to_Sigma;
    /*------------   END: Auto Generate ------------*/

};

/***/
