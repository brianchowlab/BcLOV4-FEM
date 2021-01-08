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
    Domain_Sigma_embedded_in_Sigma_restricted_to_Sigma.Setup_Data(prhs[PRHS_EMPTY_2], prhs[PRHS_Sigma_Mesh_DoFmap], PRHS_EMPTY);

    Sigma_Interp_Data.Setup("Sigma", 1, prhs[PRHS_Sigma_Interp_Data]);

    geom_Sigma_embedded_in_Sigma_restricted_to_Sigma.Setup_Mesh_Geometry(prhs[PRHS_Sigma_Mesh_Vertices], prhs[PRHS_Sigma_Mesh_DoFmap], prhs[PRHS_EMPTY_1]);
    geom_Sigma_embedded_in_Sigma_restricted_to_Sigma.Domain = &Domain_Sigma_embedded_in_Sigma_restricted_to_Sigma;

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
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // setup access to external FE functions
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // setup the desired interpolations to evaluate
    Iobj_I_kappa = new I_kappa(Sigma_Interp_Data.Num_Pts);
    Iobj_I_shape = new I_shape(Sigma_Interp_Data.Num_Pts);
    Iobj_I_tangent = new I_tangent(Sigma_Interp_Data.Num_Pts);
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // pass pointers around
    Iobj_I_kappa->geom_Sigma_embedded_in_Sigma_restricted_to_Sigma = &geom_Sigma_embedded_in_Sigma_restricted_to_Sigma;

    Iobj_I_shape->geom_Sigma_embedded_in_Sigma_restricted_to_Sigma = &geom_Sigma_embedded_in_Sigma_restricted_to_Sigma;

    Iobj_I_tangent->geom_Sigma_embedded_in_Sigma_restricted_to_Sigma = &geom_Sigma_embedded_in_Sigma_restricted_to_Sigma;

    /*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/


/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* here is the main calling routine for interpolating */
void GFI::Evaluate_Interpolations ()
{
    unsigned int Num_DoE_Interp_Pts = 0;

    // BEGIN: evaluate interpolations over the Expression Domain: Sigma
    Num_DoE_Interp_Pts = Sigma_Interp_Data.Num_Pts;
    // loop through each point
    for (unsigned int DoE_Pt = 0; DoE_Pt < Num_DoE_Interp_Pts; DoE_Pt++)
        {
        // read the DoE *element* index from the interpolation point data
        const unsigned int DoE_Elem_Index_MATLAB_style = Sigma_Interp_Data.Cell_Index[DoE_Pt];
        if (DoE_Elem_Index_MATLAB_style==0) // cell index is INVALID, so ignore!
            {
            mexPrintf("Interpolation cell index is *invalid* for point index #%d.\n",DoE_Pt+1);
            mexPrintf("    No interpolation will be done at this point!\n");
            }
        else // only compute if the cell index is valid!
            {
            unsigned int DoE_Elem_Index = DoE_Elem_Index_MATLAB_style - 1; // need to offset for C-style indexing
            Domain_Sigma_embedded_in_Sigma_restricted_to_Sigma.Read_Embed_Data(DoE_Elem_Index);

            // get the local simplex transformation
            // copy local interpolation coordinates
            Sigma_Interp_Data.Copy_Local_X(DoE_Pt,geom_Sigma_embedded_in_Sigma_restricted_to_Sigma.local_coord);
            geom_Sigma_embedded_in_Sigma_restricted_to_Sigma.Compute_Local_Transformation();


            // loop through the desired FEM interpolations
            Iobj_I_kappa->Eval_All_Interpolations(DoE_Pt);
            Iobj_I_shape->Eval_All_Interpolations(DoE_Pt);
            Iobj_I_tangent->Eval_All_Interpolations(DoE_Pt);
            }
        }
    // END: evaluate interpolations over the Expression Domain: Sigma

}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* this outputs all FEM interpolations as MATLAB vectors (in cell arrays) */
void GFI::Output_Interpolations (mxArray* plhs[])
{
    Output_Single_Interpolation(0, Iobj_I_kappa->mxInterp_Data, mxCreateString(Iobj_I_kappa->Name), plhs[0]);
    delete(Iobj_I_kappa);

    Output_Single_Interpolation(1, Iobj_I_shape->mxInterp_Data, mxCreateString(Iobj_I_shape->Name), plhs[0]);
    delete(Iobj_I_shape);

    Output_Single_Interpolation(2, Iobj_I_tangent->mxInterp_Data, mxCreateString(Iobj_I_tangent->Name), plhs[0]);
    delete(Iobj_I_tangent);

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
