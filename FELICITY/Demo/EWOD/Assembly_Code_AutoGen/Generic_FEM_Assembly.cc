/*
============================================================================================
   Methods for a C++ Class that does generic finite element assembly.

   NOTE: portions of this code are automatically generated!

   Copyright (c) 04-08-2010,  Shawn W. Walker
============================================================================================
*/

#define GFA Generic_FEM_Assembly

/***************************************************************************************/
/* constructor */
GFA::GFA (const mxArray *prhs[], const mxArray *Subset_Elem)
{
    /*------------ BEGIN: Auto Generate ------------*/
    // setup inputs
    Domain_Gamma_embedded_in_Omega_restricted_to_Gamma.Setup_Data(prhs[PRHS_Omega_Mesh_Subdomains], prhs[PRHS_Omega_Mesh_DoFmap], Subset_Elem);
    Domain_Omega_embedded_in_Omega_restricted_to_Gamma.Setup_Data(prhs[PRHS_Omega_Mesh_Subdomains], prhs[PRHS_Omega_Mesh_DoFmap], Subset_Elem);
    Domain_Omega_embedded_in_Omega_restricted_to_Omega.Setup_Data(prhs[PRHS_Omega_Mesh_Subdomains], prhs[PRHS_Omega_Mesh_DoFmap], Subset_Elem);

    geom_Gamma_embedded_in_Omega_restricted_to_Gamma.Setup_Mesh_Geometry(prhs[PRHS_Omega_Mesh_Vertices], prhs[PRHS_Omega_Mesh_DoFmap], prhs[PRHS_Omega_Mesh_Orient]);
    geom_Gamma_embedded_in_Omega_restricted_to_Gamma.Domain = &Domain_Gamma_embedded_in_Omega_restricted_to_Gamma;
    geom_Omega_embedded_in_Omega_restricted_to_Gamma.Setup_Mesh_Geometry(prhs[PRHS_Omega_Mesh_Vertices], prhs[PRHS_Omega_Mesh_DoFmap], prhs[PRHS_Omega_Mesh_Orient]);
    geom_Omega_embedded_in_Omega_restricted_to_Gamma.Domain = &Domain_Omega_embedded_in_Omega_restricted_to_Gamma;
    geom_Omega_embedded_in_Omega_restricted_to_Omega.Setup_Mesh_Geometry(prhs[PRHS_Omega_Mesh_Vertices], prhs[PRHS_Omega_Mesh_DoFmap], prhs[PRHS_Omega_Mesh_Orient]);
    geom_Omega_embedded_in_Omega_restricted_to_Omega.Domain = &Domain_Omega_embedded_in_Omega_restricted_to_Omega;

    Setup_Data(prhs);
    /*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/


/***************************************************************************************/
/* destructor */
GFA::~GFA ()
{
    // clear it
    mxFree(Sparse_Data_A.name);
    mxFree(Sparse_Data_B.name);
    mxFree(Sparse_Data_C.name);
    mxFree(Sparse_Data_D.name);
    mxFree(Sparse_Data_K.name);
    mxFree(Sparse_Data_M.name);
    mxFree(Sparse_Data_chi.name);
}
/***************************************************************************************/


/***************************************************************************************/
/* setup matrix data into a nice struct for internal use */
void GFA::Setup_Data(const mxArray *prhs[]) // input from MATLAB
{
    // access previously assembled matrices (if they exist)
    /*------------ BEGIN: Auto Generate ------------*/
    Access_Previous_FEM_Matrix(prhs[PRHS_OLD_FEM], "A", 0, Sparse_Data_A);
    Access_Previous_FEM_Matrix(prhs[PRHS_OLD_FEM], "B", 1, Sparse_Data_B);
    Access_Previous_FEM_Matrix(prhs[PRHS_OLD_FEM], "C", 2, Sparse_Data_C);
    Access_Previous_FEM_Matrix(prhs[PRHS_OLD_FEM], "D", 3, Sparse_Data_D);
    Access_Previous_FEM_Matrix(prhs[PRHS_OLD_FEM], "K", 4, Sparse_Data_K);
    Access_Previous_FEM_Matrix(prhs[PRHS_OLD_FEM], "M", 5, Sparse_Data_M);
    Access_Previous_FEM_Matrix(prhs[PRHS_OLD_FEM], "chi", 6, Sparse_Data_chi);
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // setup access to FE basis functions
    M_h_phi_restricted_to_Gamma.Setup_Function_Space(prhs[PRHS_M_h_DoFmap]);
    M_h_phi_restricted_to_Gamma.Mesh = &geom_Gamma_embedded_in_Omega_restricted_to_Gamma;
    V_h_phi_restricted_to_Gamma.Setup_Function_Space(prhs[PRHS_V_h_DoFmap]);
    V_h_phi_restricted_to_Gamma.Mesh = &geom_Omega_embedded_in_Omega_restricted_to_Gamma;
    Y_h_phi_restricted_to_Gamma.Setup_Function_Space(prhs[PRHS_Y_h_DoFmap]);
    Y_h_phi_restricted_to_Gamma.Mesh = &geom_Gamma_embedded_in_Omega_restricted_to_Gamma;
    G_h_phi_restricted_to_Omega.Setup_Function_Space(prhs[PRHS_G_h_DoFmap]);
    G_h_phi_restricted_to_Omega.Mesh = &geom_Omega_embedded_in_Omega_restricted_to_Omega;
    Q_h_phi_restricted_to_Omega.Setup_Function_Space(prhs[PRHS_Q_h_DoFmap]);
    Q_h_phi_restricted_to_Omega.Mesh = &geom_Omega_embedded_in_Omega_restricted_to_Omega;
    V_h_phi_restricted_to_Omega.Setup_Function_Space(prhs[PRHS_V_h_DoFmap]);
    V_h_phi_restricted_to_Omega.Mesh = &geom_Omega_embedded_in_Omega_restricted_to_Omega;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // setup correct number of components for CONSTANT basis functions
    chi_Gamma_col_constant_phi.Num_Comp = 1;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // setup access to external FE functions
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // setup the base matrices to compute
    Base_Matrix_A = new Base_A_Data_Type(&G_h_phi_restricted_to_Omega, &G_h_phi_restricted_to_Omega);
    Base_Matrix_B = new Base_B_Data_Type(&Q_h_phi_restricted_to_Omega, &V_h_phi_restricted_to_Omega);
    Base_Matrix_C = new Base_C_Data_Type(&M_h_phi_restricted_to_Gamma, &V_h_phi_restricted_to_Gamma);
    Base_Matrix_D = new Base_D_Data_Type(&M_h_phi_restricted_to_Gamma, &Y_h_phi_restricted_to_Gamma);
    Base_Matrix_K = new Base_K_Data_Type(&Y_h_phi_restricted_to_Gamma, &Y_h_phi_restricted_to_Gamma);
    Base_Matrix_M = new Base_M_Data_Type(&V_h_phi_restricted_to_Omega, &V_h_phi_restricted_to_Omega);
    Base_Matrix_chi = new Base_chi_Data_Type(&M_h_phi_restricted_to_Gamma, &chi_Gamma_col_constant_phi);
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // pass pointers around
    Base_Matrix_A->geom_Gamma_embedded_in_Omega_restricted_to_Gamma = &geom_Gamma_embedded_in_Omega_restricted_to_Gamma;
    Base_Matrix_A->geom_Omega_embedded_in_Omega_restricted_to_Gamma = &geom_Omega_embedded_in_Omega_restricted_to_Gamma;
    Base_Matrix_A->geom_Omega_embedded_in_Omega_restricted_to_Omega = &geom_Omega_embedded_in_Omega_restricted_to_Omega;
    Base_Matrix_A->G_h_phi_restricted_to_Omega = &G_h_phi_restricted_to_Omega;
    Base_Matrix_A->M_h_phi_restricted_to_Gamma = &M_h_phi_restricted_to_Gamma;
    Base_Matrix_A->Q_h_phi_restricted_to_Omega = &Q_h_phi_restricted_to_Omega;
    Base_Matrix_A->V_h_phi_restricted_to_Gamma = &V_h_phi_restricted_to_Gamma;
    Base_Matrix_A->V_h_phi_restricted_to_Omega = &V_h_phi_restricted_to_Omega;
    Base_Matrix_A->Y_h_phi_restricted_to_Gamma = &Y_h_phi_restricted_to_Gamma;
    Base_Matrix_B->geom_Gamma_embedded_in_Omega_restricted_to_Gamma = &geom_Gamma_embedded_in_Omega_restricted_to_Gamma;
    Base_Matrix_B->geom_Omega_embedded_in_Omega_restricted_to_Gamma = &geom_Omega_embedded_in_Omega_restricted_to_Gamma;
    Base_Matrix_B->geom_Omega_embedded_in_Omega_restricted_to_Omega = &geom_Omega_embedded_in_Omega_restricted_to_Omega;
    Base_Matrix_B->G_h_phi_restricted_to_Omega = &G_h_phi_restricted_to_Omega;
    Base_Matrix_B->M_h_phi_restricted_to_Gamma = &M_h_phi_restricted_to_Gamma;
    Base_Matrix_B->Q_h_phi_restricted_to_Omega = &Q_h_phi_restricted_to_Omega;
    Base_Matrix_B->V_h_phi_restricted_to_Gamma = &V_h_phi_restricted_to_Gamma;
    Base_Matrix_B->V_h_phi_restricted_to_Omega = &V_h_phi_restricted_to_Omega;
    Base_Matrix_B->Y_h_phi_restricted_to_Gamma = &Y_h_phi_restricted_to_Gamma;
    Base_Matrix_C->geom_Gamma_embedded_in_Omega_restricted_to_Gamma = &geom_Gamma_embedded_in_Omega_restricted_to_Gamma;
    Base_Matrix_C->geom_Omega_embedded_in_Omega_restricted_to_Gamma = &geom_Omega_embedded_in_Omega_restricted_to_Gamma;
    Base_Matrix_C->geom_Omega_embedded_in_Omega_restricted_to_Omega = &geom_Omega_embedded_in_Omega_restricted_to_Omega;
    Base_Matrix_C->G_h_phi_restricted_to_Omega = &G_h_phi_restricted_to_Omega;
    Base_Matrix_C->M_h_phi_restricted_to_Gamma = &M_h_phi_restricted_to_Gamma;
    Base_Matrix_C->Q_h_phi_restricted_to_Omega = &Q_h_phi_restricted_to_Omega;
    Base_Matrix_C->V_h_phi_restricted_to_Gamma = &V_h_phi_restricted_to_Gamma;
    Base_Matrix_C->V_h_phi_restricted_to_Omega = &V_h_phi_restricted_to_Omega;
    Base_Matrix_C->Y_h_phi_restricted_to_Gamma = &Y_h_phi_restricted_to_Gamma;
    Base_Matrix_D->geom_Gamma_embedded_in_Omega_restricted_to_Gamma = &geom_Gamma_embedded_in_Omega_restricted_to_Gamma;
    Base_Matrix_D->geom_Omega_embedded_in_Omega_restricted_to_Gamma = &geom_Omega_embedded_in_Omega_restricted_to_Gamma;
    Base_Matrix_D->geom_Omega_embedded_in_Omega_restricted_to_Omega = &geom_Omega_embedded_in_Omega_restricted_to_Omega;
    Base_Matrix_D->G_h_phi_restricted_to_Omega = &G_h_phi_restricted_to_Omega;
    Base_Matrix_D->M_h_phi_restricted_to_Gamma = &M_h_phi_restricted_to_Gamma;
    Base_Matrix_D->Q_h_phi_restricted_to_Omega = &Q_h_phi_restricted_to_Omega;
    Base_Matrix_D->V_h_phi_restricted_to_Gamma = &V_h_phi_restricted_to_Gamma;
    Base_Matrix_D->V_h_phi_restricted_to_Omega = &V_h_phi_restricted_to_Omega;
    Base_Matrix_D->Y_h_phi_restricted_to_Gamma = &Y_h_phi_restricted_to_Gamma;
    Base_Matrix_K->geom_Gamma_embedded_in_Omega_restricted_to_Gamma = &geom_Gamma_embedded_in_Omega_restricted_to_Gamma;
    Base_Matrix_K->geom_Omega_embedded_in_Omega_restricted_to_Gamma = &geom_Omega_embedded_in_Omega_restricted_to_Gamma;
    Base_Matrix_K->geom_Omega_embedded_in_Omega_restricted_to_Omega = &geom_Omega_embedded_in_Omega_restricted_to_Omega;
    Base_Matrix_K->G_h_phi_restricted_to_Omega = &G_h_phi_restricted_to_Omega;
    Base_Matrix_K->M_h_phi_restricted_to_Gamma = &M_h_phi_restricted_to_Gamma;
    Base_Matrix_K->Q_h_phi_restricted_to_Omega = &Q_h_phi_restricted_to_Omega;
    Base_Matrix_K->V_h_phi_restricted_to_Gamma = &V_h_phi_restricted_to_Gamma;
    Base_Matrix_K->V_h_phi_restricted_to_Omega = &V_h_phi_restricted_to_Omega;
    Base_Matrix_K->Y_h_phi_restricted_to_Gamma = &Y_h_phi_restricted_to_Gamma;
    Base_Matrix_M->geom_Gamma_embedded_in_Omega_restricted_to_Gamma = &geom_Gamma_embedded_in_Omega_restricted_to_Gamma;
    Base_Matrix_M->geom_Omega_embedded_in_Omega_restricted_to_Gamma = &geom_Omega_embedded_in_Omega_restricted_to_Gamma;
    Base_Matrix_M->geom_Omega_embedded_in_Omega_restricted_to_Omega = &geom_Omega_embedded_in_Omega_restricted_to_Omega;
    Base_Matrix_M->G_h_phi_restricted_to_Omega = &G_h_phi_restricted_to_Omega;
    Base_Matrix_M->M_h_phi_restricted_to_Gamma = &M_h_phi_restricted_to_Gamma;
    Base_Matrix_M->Q_h_phi_restricted_to_Omega = &Q_h_phi_restricted_to_Omega;
    Base_Matrix_M->V_h_phi_restricted_to_Gamma = &V_h_phi_restricted_to_Gamma;
    Base_Matrix_M->V_h_phi_restricted_to_Omega = &V_h_phi_restricted_to_Omega;
    Base_Matrix_M->Y_h_phi_restricted_to_Gamma = &Y_h_phi_restricted_to_Gamma;
    Base_Matrix_chi->geom_Gamma_embedded_in_Omega_restricted_to_Gamma = &geom_Gamma_embedded_in_Omega_restricted_to_Gamma;
    Base_Matrix_chi->geom_Omega_embedded_in_Omega_restricted_to_Gamma = &geom_Omega_embedded_in_Omega_restricted_to_Gamma;
    Base_Matrix_chi->geom_Omega_embedded_in_Omega_restricted_to_Omega = &geom_Omega_embedded_in_Omega_restricted_to_Omega;
    Base_Matrix_chi->G_h_phi_restricted_to_Omega = &G_h_phi_restricted_to_Omega;
    Base_Matrix_chi->M_h_phi_restricted_to_Gamma = &M_h_phi_restricted_to_Gamma;
    Base_Matrix_chi->Q_h_phi_restricted_to_Omega = &Q_h_phi_restricted_to_Omega;
    Base_Matrix_chi->V_h_phi_restricted_to_Gamma = &V_h_phi_restricted_to_Gamma;
    Base_Matrix_chi->V_h_phi_restricted_to_Omega = &V_h_phi_restricted_to_Omega;
    Base_Matrix_chi->Y_h_phi_restricted_to_Gamma = &Y_h_phi_restricted_to_Gamma;
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // setup the block matrices to assemble
    Block_Assemble_Matrix_A = new Block_Assemble_A_Data_Type(&Sparse_Data_A, Base_Matrix_A);
    Block_Assemble_Matrix_B = new Block_Assemble_B_Data_Type(&Sparse_Data_B, Base_Matrix_B);
    Block_Assemble_Matrix_C = new Block_Assemble_C_Data_Type(&Sparse_Data_C, Base_Matrix_C);
    Block_Assemble_Matrix_D = new Block_Assemble_D_Data_Type(&Sparse_Data_D, Base_Matrix_D);
    Block_Assemble_Matrix_K = new Block_Assemble_K_Data_Type(&Sparse_Data_K, Base_Matrix_K);
    Block_Assemble_Matrix_M = new Block_Assemble_M_Data_Type(&Sparse_Data_M, Base_Matrix_M);
    Block_Assemble_Matrix_chi = new Block_Assemble_chi_Data_Type(&Sparse_Data_chi, Base_Matrix_chi);
    /*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/


/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* here is the main calling routine for assembling all FEM matrices */
void GFA::Assemble_Matrices ()
{
    // BEGIN: assemble matrices over the Integration Domain: Gamma

    if (Domain_Gamma_embedded_in_Omega_restricted_to_Gamma.Sub_Assem_List.empty())
        {
        mexPrintf("This integration domain is empty: Gamma\n");
        mexPrintf(" ... so no assembly necessary.\n");
        }

    // loop through each element
    for (std::vector<unsigned int>::iterator DoI_Ind = Domain_Gamma_embedded_in_Omega_restricted_to_Gamma.Sub_Assem_List.begin();
              DoI_Ind != Domain_Gamma_embedded_in_Omega_restricted_to_Gamma.Sub_Assem_List.end(); ++DoI_Ind)
        {
        Domain_Gamma_embedded_in_Omega_restricted_to_Gamma.Read_Embed_Data(*DoI_Ind);
        Domain_Omega_embedded_in_Omega_restricted_to_Gamma.Read_Embed_Data(*DoI_Ind);

        // get the local simplex transformation
        geom_Gamma_embedded_in_Omega_restricted_to_Gamma.Compute_Local_Transformation();
        geom_Omega_embedded_in_Omega_restricted_to_Gamma.Compute_Local_Transformation();

        // perform pre-computations with FE basis functions
        // NOTE: this must come before the external FE coefficient functions
        M_h_phi_restricted_to_Gamma.Transform_Basis_Functions();
        V_h_phi_restricted_to_Gamma.Transform_Basis_Functions();
        Y_h_phi_restricted_to_Gamma.Transform_Basis_Functions();


        // loop through the FE matrices to compute
        Base_Matrix_C->Tabulate_Tensor(geom_Gamma_embedded_in_Omega_restricted_to_Gamma);
        Base_Matrix_D->Tabulate_Tensor(geom_Gamma_embedded_in_Omega_restricted_to_Gamma);
        Base_Matrix_K->Tabulate_Tensor(geom_Gamma_embedded_in_Omega_restricted_to_Gamma);
        Base_Matrix_chi->Tabulate_Tensor(geom_Gamma_embedded_in_Omega_restricted_to_Gamma);

        // loop through the block FE matrices to assemble
        Block_Assemble_Matrix_C->Add_Entries_To_Global_Matrix_Gamma(Base_Matrix_C);
        Block_Assemble_Matrix_D->Add_Entries_To_Global_Matrix_Gamma(Base_Matrix_D);
        Block_Assemble_Matrix_K->Add_Entries_To_Global_Matrix_Gamma(Base_Matrix_K);
        Block_Assemble_Matrix_chi->Add_Entries_To_Global_Matrix_Gamma(Base_Matrix_chi);
        }
    // END: assemble matrices over the Integration Domain: Gamma

    // BEGIN: assemble matrices over the Integration Domain: Omega

    if (Domain_Omega_embedded_in_Omega_restricted_to_Omega.Sub_Assem_List.empty())
        {
        mexPrintf("This integration domain is empty: Omega\n");
        mexPrintf(" ... so no assembly necessary.\n");
        }

    // loop through each element
    for (std::vector<unsigned int>::iterator DoI_Ind = Domain_Omega_embedded_in_Omega_restricted_to_Omega.Sub_Assem_List.begin();
              DoI_Ind != Domain_Omega_embedded_in_Omega_restricted_to_Omega.Sub_Assem_List.end(); ++DoI_Ind)
        {
        Domain_Omega_embedded_in_Omega_restricted_to_Omega.Read_Embed_Data(*DoI_Ind);

        // get the local simplex transformation
        geom_Omega_embedded_in_Omega_restricted_to_Omega.Compute_Local_Transformation();

        // perform pre-computations with FE basis functions
        // NOTE: this must come before the external FE coefficient functions
        G_h_phi_restricted_to_Omega.Transform_Basis_Functions();
        Q_h_phi_restricted_to_Omega.Transform_Basis_Functions();
        V_h_phi_restricted_to_Omega.Transform_Basis_Functions();


        // loop through the FE matrices to compute
        Base_Matrix_A->Tabulate_Tensor(geom_Omega_embedded_in_Omega_restricted_to_Omega);
        Base_Matrix_B->Tabulate_Tensor(geom_Omega_embedded_in_Omega_restricted_to_Omega);
        Base_Matrix_M->Tabulate_Tensor(geom_Omega_embedded_in_Omega_restricted_to_Omega);

        // loop through the block FE matrices to assemble
        Block_Assemble_Matrix_A->Add_Entries_To_Global_Matrix_Omega(Base_Matrix_A);
        Block_Assemble_Matrix_B->Add_Entries_To_Global_Matrix_Omega(Base_Matrix_B);
        Block_Assemble_Matrix_M->Add_Entries_To_Global_Matrix_Omega(Base_Matrix_M);
        }
    // END: assemble matrices over the Integration Domain: Omega

}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* this outputs all FEM matrices as MATLAB sparse matrices */
void GFA::Output_Matrices (mxArray* plhs[])
{
    // declare internal matrix data storage pointer
    mxArray* Sparse_ptr;

    // create sparse MATLAB matrices and pass them back to MATLAB

    Sparse_ptr = Block_Assemble_Matrix_A->MAT->export_matrix();
    Output_Matrix(0, Sparse_ptr, mxCreateString(Block_Assemble_Matrix_A->Name), plhs[0]);
    delete(Block_Assemble_Matrix_A);

    Sparse_ptr = Block_Assemble_Matrix_B->MAT->export_matrix();
    Output_Matrix(1, Sparse_ptr, mxCreateString(Block_Assemble_Matrix_B->Name), plhs[0]);
    delete(Block_Assemble_Matrix_B);

    Sparse_ptr = Block_Assemble_Matrix_C->MAT->export_matrix();
    Output_Matrix(2, Sparse_ptr, mxCreateString(Block_Assemble_Matrix_C->Name), plhs[0]);
    delete(Block_Assemble_Matrix_C);

    Sparse_ptr = Block_Assemble_Matrix_D->MAT->export_matrix();
    Output_Matrix(3, Sparse_ptr, mxCreateString(Block_Assemble_Matrix_D->Name), plhs[0]);
    delete(Block_Assemble_Matrix_D);

    Sparse_ptr = Block_Assemble_Matrix_K->MAT->export_matrix();
    Output_Matrix(4, Sparse_ptr, mxCreateString(Block_Assemble_Matrix_K->Name), plhs[0]);
    delete(Block_Assemble_Matrix_K);

    Sparse_ptr = Block_Assemble_Matrix_M->MAT->export_matrix();
    Output_Matrix(5, Sparse_ptr, mxCreateString(Block_Assemble_Matrix_M->Name), plhs[0]);
    delete(Block_Assemble_Matrix_M);

    Sparse_ptr = Block_Assemble_Matrix_chi->MAT->export_matrix();
    Output_Matrix(6, Sparse_ptr, mxCreateString(Block_Assemble_Matrix_chi->Name), plhs[0]);
    delete(Block_Assemble_Matrix_chi);

}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/


/***************************************************************************************/
/* setup sparse FEM matrices to be output to MATLAB */
#define NUM_Matrix_Fieldnames (sizeof(Matrix_Fieldnames)/sizeof(*Matrix_Fieldnames))
void GFA::Init_Output_Matrices (mxArray* plhs[])               // output
{
    // // declare internal matrix data storage pointers
    // mxArray *Sparse_ptr;

    // declare constant arrays (see 'Generic_FEM_Assembly.h')
    const char *Matrix_Fieldnames[] = {OUT_MAT_str, OUT_FEM_NAME_str};

    // declare parameters for outputing structures to MATLAB
    mwSize Matrix_dims[2] = {1, 1}; // just initialize to a 1x1 struct

    // set the number of MATLAB structs to create
    Matrix_dims[1] = NUM_FEM_MAT;
    /*** setup LHS argument of MATLAB calling function (i.e. the output structs) ***/
    //                            2, 1xN struct,    X sub-fields,  field names
    plhs[0] = mxCreateStructArray(2, Matrix_dims, NUM_Matrix_Fieldnames, Matrix_Fieldnames);
}
/***************************************************************************************/


/***************************************************************************************/
/* output sparse FEM matrix to MATLAB */
void GFA::Output_Matrix(mwIndex index, mxArray* Sparse_ptr, mxArray* Matrix_Name,   // input
                        mxArray* mxOUT) // output
{
	// point the output MATLAB structure fields to the correct data blocks!
	mxSetFieldByNumber(mxOUT,index,mxGetFieldNumber(mxOUT,OUT_MAT_str), Sparse_ptr);
	mxSetFieldByNumber(mxOUT,index,mxGetFieldNumber(mxOUT,OUT_FEM_NAME_str), Matrix_Name);
}
/***************************************************************************************/


/***************************************************************************************/
/* this verifies the incoming FEM matrix index and name matches, and accesses the data
   (if appropriate) */
void GFA::Access_Previous_FEM_Matrix(
                          const mxArray* OLD_FEM, const char* Matrix_Name,  // inputs
						  const int& Array_Index,                           // inputs
						  PTR_TO_SPARSE&  Data)                             // outputs
{
	if (!mxIsEmpty(OLD_FEM))
		{
		const int Num_Prev_Matrices = (const int) mxGetNumberOfElements(OLD_FEM);
		if (Array_Index >= Num_Prev_Matrices)
			mexErrMsgTxt("Index exceeds the number of incoming FEM matrices!");

        // determine which index of the subdomain array is the one we want
		const mxArray* mxMAT_Name = mxGetField(OLD_FEM, (mwIndex)Array_Index, OUT_FEM_NAME_str);

		/* Copy the string data over... */
        mwSize name_len = mxGetNumberOfElements(mxMAT_Name) + 1;
        char* name_in   = (char*) mxCalloc(name_len, sizeof(char));
        if (mxGetString(mxMAT_Name, name_in, name_len) != 0)
            mexErrMsgIdAndTxt("MATLAB:explore:invalidStringArray","Could not convert Matrix_Name string data.");

        // if they match, then access the data
        const bool name_equal = (strcmp(Matrix_Name,name_in)==0);
        mxFree(name_in);
        if (name_equal)
            Read_Sparse_Ptr(OLD_FEM, Array_Index, Data);
		else // fail
			{
			mexPrintf("ERROR: The Matrix_Name: %s\n",Matrix_Name);
			mexPrintf("ERROR:     does not match: %s\n",name_in);
			mexPrintf("ERROR:     at OLD_FEM index: %d\n",Array_Index+1); // put into MATLAB style
			mexErrMsgTxt("Check your OLD_FEM data!");
			}
        }
	else
		Clear_Sparse_Ptr(Data);
}
/***************************************************************************************/


/***************************************************************************************/
/* read sparse pointer data, so we can reuse sparse data structure in assembly */
void GFA::Read_Sparse_Ptr (const mxArray* OLD_FEM, const int& Array_Index,   // inputs
						   PTR_TO_SPARSE&  Data)                             // outputs
{
	// indicate that there is a matrix
	Data.valid = true;

	// store the matrix name as a string
	const mxArray* String_ptr = mxGetField(OLD_FEM, (mwIndex)Array_Index, OUT_FEM_NAME_str);
	unsigned int buflen = ((unsigned int)mxGetN(String_ptr))*sizeof(mxChar) + 1;
	Data.name = (char*) mxMalloc((size_t)buflen);
	// copy name over
	const int status = mxGetString(String_ptr, Data.name, (mwSize)buflen);
	if (status==1) mexErrMsgTxt("FEM matrix string name not read in correctly!");

	// get pointer to sparse MATLAB matrix
	const mxArray* Sparse_ptr = mxGetField(OLD_FEM, (mwIndex)Array_Index, OUT_MAT_str);
	// store CSC matrix format info
	Data.m  = (int) mxGetM(Sparse_ptr);
	Data.n  = (int) mxGetN(Sparse_ptr);
	Data.jc = (mwIndex*) mxGetJc(Sparse_ptr);
	Data.ir = (mwIndex*) mxGetIr(Sparse_ptr);
	Data.pr =  (double*) mxGetPr(Sparse_ptr);
}
/***************************************************************************************/


/***************************************************************************************/
/* clear sparse pointer data */
void GFA::Clear_Sparse_Ptr (PTR_TO_SPARSE&  Data)
{
	// clear pointers
	Data.valid = false;
	Data.name  = NULL;
	Data.m     = 0;
	Data.n     = 0;
	Data.jc    = NULL;
	Data.ir    = NULL;
	Data.pr    = NULL;
}
/***************************************************************************************/

#undef GFA

/***/
