/*
============================================================================================
   Methods for a C++ Class that does generic finite element interpolation.

   NOTE: portions of this code are automatically generated!

   Copyright (c) 01-29-2013,  Shawn W. Walker
============================================================================================
*/

#define GFI Generic_FEM_Interpolation

/***************************************************************************************/
/* constructor */
GFI::GFI (const mxArray *prhs[])
{
    /*------------ BEGIN: Auto Generate ------------*/
    // setup inputs
    const mxArray *PRHS_EMPTY = NULL;
    Domain_Omega_embedded_in_Omega_restricted_to_Omega.Setup_Data(prhs[PRHS_EMPTY_2], prhs[PRHS_Omega_Mesh_DoFmap], PRHS_EMPTY);

    Omega_Interp_Data.Setup("Omega", 2, prhs[PRHS_Omega_Interp_Data]);

    geom_Omega_embedded_in_Omega_restricted_to_Omega.Setup_Mesh_Geometry(prhs[PRHS_Omega_Mesh_Vertices], prhs[PRHS_Omega_Mesh_DoFmap], prhs[PRHS_EMPTY_1]);
    geom_Omega_embedded_in_Omega_restricted_to_Omega.Domain = &Domain_Omega_embedded_in_Omega_restricted_to_Omega;

    Setup_Data(prhs);
    /*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/


/***************************************************************************************/
/* destructor */
GFI::~GFI ()
{
}
/***************************************************************************************/


/***************************************************************************************/
/* setup data into a nice struct for internal use */
void GFI::Setup_Data(const mxArray *prhs[]) // input from MATLAB
{
    /*------------ BEGIN: Auto Generate ------------*/
    // setup access to FE basis functions
    Ned1_phi_restricted_to_Omega.Setup_Function_Space(prhs[PRHS_Ned1_DoFmap]);
    Ned1_phi_restricted_to_Omega.Mesh = &geom_Omega_embedded_in_Omega_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // setup access to external FE functions
    v_restricted_to_Omega.Setup_Function_Space(prhs[PRHS_v_Values], &Ned1_phi_restricted_to_Omega);
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // setup the desired interpolations to evaluate
    Iobj_I_v = new I_v(Omega_Interp_Data.Num_Pts);
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // pass pointers around
    Iobj_I_v->geom_Omega_embedded_in_Omega_restricted_to_Omega = &geom_Omega_embedded_in_Omega_restricted_to_Omega;
    Iobj_I_v->v_restricted_to_Omega = &v_restricted_to_Omega;

    /*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/


/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* here is the main calling routine for interpolating */
void GFI::Evaluate_Interpolations ()
{
    unsigned int Num_DoE_Interp_Pts = 0;

    // BEGIN: evaluate interpolations over the Expression Domain: Omega
    Num_DoE_Interp_Pts = Omega_Interp_Data.Num_Pts;
    // loop through each point
    for (unsigned int DoE_Pt = 0; DoE_Pt < Num_DoE_Interp_Pts; DoE_Pt++)
        {
        // read the DoE *element* index from the interpolation point data
        const unsigned int DoE_Elem_Index_MATLAB_style = Omega_Interp_Data.Cell_Index[DoE_Pt];
        if (DoE_Elem_Index_MATLAB_style==0) // cell index is INVALID, so ignore!
            {
            mexPrintf("Interpolation cell index is *invalid* for point index #%d.\n",DoE_Pt+1);
            mexPrintf("    No interpolation will be done at this point!\n");
            }
        else // only compute if the cell index is valid!
            {
            unsigned int DoE_Elem_Index = DoE_Elem_Index_MATLAB_style - 1; // need to offset for C-style indexing
            Domain_Omega_embedded_in_Omega_restricted_to_Omega.Read_Embed_Data(DoE_Elem_Index);

            // get the local simplex transformation
            // copy local interpolation coordinates
            Omega_Interp_Data.Copy_Local_X(DoE_Pt,geom_Omega_embedded_in_Omega_restricted_to_Omega.local_coord);
            geom_Omega_embedded_in_Omega_restricted_to_Omega.Compute_Local_Transformation();

            // perform pre-computations with FE basis functions
            // NOTE: this must come before the external FE coefficient functions
            // copy local interpolation coordinates
            Omega_Interp_Data.Copy_Local_X(DoE_Pt,Ned1_phi_restricted_to_Omega.local_coord);
            Ned1_phi_restricted_to_Omega.Transform_Basis_Functions();

            // perform pre-computations with external FE coefficient functions
            v_restricted_to_Omega.Compute_Func();

            // loop through the desired FEM interpolations
            Iobj_I_v->Eval_All_Interpolations(DoE_Pt);
            }
        }
    // END: evaluate interpolations over the Expression Domain: Omega

}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* this outputs all FEM interpolations as MATLAB vectors (in cell arrays) */
void GFI::Output_Interpolations (mxArray* plhs[])
{
    Output_Single_Interpolation(0, Iobj_I_v->mxInterp_Data, mxCreateString(Iobj_I_v->Name), plhs[0]);
    delete(Iobj_I_v);

}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/


/***************************************************************************************/
/* setup INTERP data to be output to MATLAB */
#define NUM_Interp_Fieldnames (sizeof(Interp_Fieldnames)/sizeof(*Interp_Fieldnames))
void GFI::Init_Output_Data (mxArray* plhs[])               // output
{
    // declare constant arrays (see 'Generic_FEM_Interpolation.h')
    const char *Interp_Fieldnames[] = {OUT_DATA_str, OUT_NAME_str};

    // declare parameters for outputing structures to MATLAB
    mwSize Interp_dims[2] = {1, 1}; // just initialize to a 1x1 struct

    // set the number of MATLAB structs to create
    Interp_dims[1] = NUM_FEM_INTERP;
    /*** setup LHS argument of MATLAB calling function (i.e. the output structs) ***/
    //                            2, 1xN struct,    X sub-fields,  field names
    plhs[0] = mxCreateStructArray(2, Interp_dims, NUM_Interp_Fieldnames, Interp_Fieldnames);
}
/***************************************************************************************/


/***************************************************************************************/
/* output single FEM interpolation data to MATLAB */
void GFI::Output_Single_Interpolation(mwIndex index, mxArray* Interp_ptr, mxArray* Interp_Name,   // input
                                      mxArray* mxOUT) // output
{
	// point the output MATLAB structure fields to the correct data blocks!
	mxSetFieldByNumber(mxOUT,index,mxGetFieldNumber(mxOUT,OUT_DATA_str), Interp_ptr);
	mxSetFieldByNumber(mxOUT,index,mxGetFieldNumber(mxOUT,OUT_NAME_str), Interp_Name);
}
/***************************************************************************************/

#undef GFI

/***/
