/*
============================================================================================
   This file contains an implementation of a derived C++ Class from the abstract base class
   in 'Block_Global_FE_Matrix.cc'.

   NOTE: portions of this code are automatically generated!

   Copyright (c) 06-14-2016,  Shawn W. Walker
============================================================================================
*/

/*------------ BEGIN: Auto Generate ------------*/
// Block Global matrix contains:
//
//Vector_Hess_Matrix
//
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
// define the name of the block FE matrix (should be the same as the filename of this file)
#define SpecificFEM        Block_Assemble_Vector_Hess_Matrix_Data_Type
#define SpecificFEM_str   "Vector_Hess_Matrix"

// the row function space is Type = CG, Name = "lagrange_deg2_dim1"

// set the number of blocks along the row dimension
#define ROW_Num_Block  1
// set the number of basis functions on each element
#define ROW_NB         3

// the col function space is Type = CG, Name = "lagrange_deg2_dim1"

// set the number of blocks along the col dimension
#define COL_Num_Block  1
// set the number of basis functions on each element
#define COL_NB  3
/*------------   END: Auto Generate ------------*/


/***************************************************************************************/
/* C++ (Specific) Block Global FE matrix class definition */
class SpecificFEM: public Block_Assemble_FE_MATRIX_Class // derive from base class
{
public:
    // data structure for sub-blocks of the block FE matrix:
    //      offsets for inserting sub-blocks into the global block matrix
    int     Block_Row_Shift[ROW_Num_Block];
    int     Block_Col_Shift[COL_Num_Block];

    /*------------ BEGIN: Auto Generate ------------*/
    // constructor
    SpecificFEM (const PTR_TO_SPARSE*, const Base_Vector_Hess_Matrix_Data_Type*);
    ~SpecificFEM (); // destructor
    void Init_Matrix_Assembler_Object(bool);
    void Add_Entries_To_Global_Matrix_Sigma (const Base_Vector_Hess_Matrix_Data_Type*);
    /*------------   END: Auto Generate ------------*/
private:
};


/***************************************************************************************/
/* constructor */
/*------------ BEGIN: Auto Generate ------------*/
SpecificFEM::SpecificFEM (const PTR_TO_SPARSE* Prev_Sparse_Data,
                          const Base_Vector_Hess_Matrix_Data_Type* Block_00
                          ) :
/*------------   END: Auto Generate ------------*/
Block_Assemble_FE_MATRIX_Class () // call the base class constructor
{
    // set the 'Name' of this Global matrix
    Name = (char*) SpecificFEM_str; // this should be similar to the Class identifier
    // record the number of block matrices (in the global matrix)
    Num_Blocks = 1;

    Sparse_Data = Prev_Sparse_Data;
    bool simple_assembler;
    /*------------ BEGIN: Auto Generate ------------*/
    simple_assembler = false; // use sparse matrix format

    // record the size of the block global matrix (only for ONE block)
    global_num_row = Block_00->get_global_num_row();
    global_num_col = Block_00->get_global_num_col();
    /*------------   END: Auto Generate ------------*/

    /*------------ BEGIN: Auto Generate ------------*/
    // input information for offsetting the sub-blocks
    Block_Row_Shift[0] = 0*Block_00->get_global_num_row();
    Block_Col_Shift[0] = 0*Block_00->get_global_num_col();
    /*------------   END: Auto Generate ------------*/

    Init_Matrix_Assembler_Object(simple_assembler);
}
/***************************************************************************************/


/***************************************************************************************/
/* destructor: should not usually need to be modified */
SpecificFEM::~SpecificFEM ()
{
    delete(MAT);
}
/***************************************************************************************/


/***************************************************************************************/
/* this initializes the object that handles matrix assembly */
void SpecificFEM::Init_Matrix_Assembler_Object(bool use_simple_assembler)
{
    if (use_simple_assembler)
        // create SimpleMatrixAssembler object
        MAT = new SimpleMatrixAssembler(global_num_row,global_num_col); // assemble from scratch
    else if (Sparse_Data->valid)
        {
        int str_diff = strcmp(Sparse_Data->name,Name);
        if (str_diff!=0)
            {
            mexPrintf("Matrix names do not match!\n");
            mexPrintf("Previously assembled matrix name is:  %s.\n", Sparse_Data->name);
            mexPrintf("The name SHOULD have been:            %s.\n", Name);
            mexErrMsgTxt("Check the previously assembled matrix structure!\n");
            }
        if (Sparse_Data->m!=global_num_row)
            {
            mexPrintf("Error with this matrix:   %s.\n", Name);
            mexErrMsgTxt("Number of rows in previous matrix does not match what the new matrix should be.\n");
            }
        if (Sparse_Data->n!=global_num_col)
            {
            mexPrintf("Error with this matrix:   %s.\n", Name);
            mexErrMsgTxt("Number of columns in previous matrix does not match what the new matrix should be.\n");
            }
        // create MatrixReassembler object (cf. David Bindel)
        MAT = new MatrixReassembler(Sparse_Data->jc,Sparse_Data->ir,Sparse_Data->pr,global_num_row,global_num_col);
        // use previous sparse structure
        }
    else // create MatrixAssembler object (cf. David Bindel)
        MAT = new MatrixAssembler(global_num_row,global_num_col); // assemble from scratch
}
/***************************************************************************************/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* assemble a local FE matrix on the domain Sigma */
void SpecificFEM::Add_Entries_To_Global_Matrix_Sigma(const Base_Vector_Hess_Matrix_Data_Type* Block_00)
{
    // get local to global index map for the current ROW element
    int  Row_Indices_0[ROW_NB];
    const int row_elem_index = Block_00->Vector_P2_phi_restricted_to_Sigma->Mesh->Domain->Sub_Cell_Index;
    Block_00->Vector_P2_phi_restricted_to_Sigma->Get_Local_to_Global_DoFmap(row_elem_index, Row_Indices_0);
    // shift Row_Indices_0 to account for the block matrix offset
    for (unsigned int ri = 0; ri < ROW_NB; ++ri)
        Row_Indices_0[ri] += Block_Row_Shift[0];

    // get local to global index map for the current COL element
    int  Col_Indices_0[COL_NB];
    const int col_elem_index = Block_00->Vector_P2_phi_restricted_to_Sigma->Mesh->Domain->Sub_Cell_Index;
    Block_00->Vector_P2_phi_restricted_to_Sigma->Get_Local_to_Global_DoFmap(col_elem_index, Col_Indices_0);
    // shift Col_Indices_0 to account for the block matrix offset
    for (unsigned int ci = 0; ci < COL_NB; ++ci)
        Col_Indices_0[ci] += Block_Col_Shift[0];

    // sort row indices (ascending order)
    int  Local_Row_Ind[ROW_NB] = {0, 1, 2};
    std::sort(Local_Row_Ind, Local_Row_Ind+ROW_NB,
              [&Row_Indices_0](int kk, int qq) { return (Row_Indices_0[kk] < Row_Indices_0[qq]); });

    // sort col indices (ascending order)
    int  Local_Col_Ind[COL_NB] = {0, 1, 2};
    std::sort(Local_Col_Ind, Local_Col_Ind+COL_NB,
              [&Col_Indices_0](int kk, int qq) { return (Col_Indices_0[kk] < Col_Indices_0[qq]); });

    // allocate (I,J,V) arrays to hold "big" local matrix
    int     COO_I[2*ROW_NB*COL_NB];
    int     COO_J[2*COL_NB];
    int     COO_J_IV_Range[2*COL_NB + 1]; // indicates what parts I,V correspond to J
    double  COO_V[2*ROW_NB*COL_NB];

    /* fill the (I,J,V) arrays (sorted) */

    // write column #0
    // write (I,J,V) data for Block_00->FE_Tensor_0, i.e. the (0,0) block
    // write the data directly
    COO_I[0] = Row_Indices_0[Local_Row_Ind[0]];
    COO_J[0] = Col_Indices_0[Local_Col_Ind[0]];
    COO_J_IV_Range[0] = 0;
    COO_V[0] = Block_00->FE_Tensor_0[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[0]];
    COO_I[1] = Row_Indices_0[Local_Row_Ind[1]];
    COO_V[1] = Block_00->FE_Tensor_0[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[1]];
    COO_I[2] = Row_Indices_0[Local_Row_Ind[2]];
    COO_V[2] = Block_00->FE_Tensor_0[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[2]];

    // write column #1
    // write (I,J,V) data for Block_00->FE_Tensor_0, i.e. the (0,0) block
    // write the data directly
    COO_I[3] = Row_Indices_0[Local_Row_Ind[0]];
    COO_J[1] = Col_Indices_0[Local_Col_Ind[1]];
    COO_J_IV_Range[1] = 3;
    COO_V[3] = Block_00->FE_Tensor_0[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[0]];
    COO_I[4] = Row_Indices_0[Local_Row_Ind[1]];
    COO_V[4] = Block_00->FE_Tensor_0[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[1]];
    COO_I[5] = Row_Indices_0[Local_Row_Ind[2]];
    COO_V[5] = Block_00->FE_Tensor_0[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[2]];

    // write column #2
    // write (I,J,V) data for Block_00->FE_Tensor_0, i.e. the (0,0) block
    // write the data directly
    COO_I[6] = Row_Indices_0[Local_Row_Ind[0]];
    COO_J[2] = Col_Indices_0[Local_Col_Ind[2]];
    COO_J_IV_Range[2] = 6;
    COO_V[6] = Block_00->FE_Tensor_0[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[0]];
    COO_I[7] = Row_Indices_0[Local_Row_Ind[1]];
    COO_V[7] = Block_00->FE_Tensor_0[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[1]];
    COO_I[8] = Row_Indices_0[Local_Row_Ind[2]];
    COO_V[8] = Block_00->FE_Tensor_0[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[2]];

    // write column #0
    // write (I,J,V) data for Block_00->FE_Tensor_1, i.e. the (1,1) block
    // actually:  copy Block_00->FE_Tensor_0 to Block_00->FE_Tensor_1
    // write the data directly
    COO_I[9] = Row_Indices_0[Local_Row_Ind[0]] + Block_00->Row_Shift[1];
    COO_J[3] = Col_Indices_0[Local_Col_Ind[0]] + Block_00->Col_Shift[1];
    COO_J_IV_Range[3] = 9;
    COO_V[9] = Block_00->FE_Tensor_0[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[0]];
    COO_I[10] = Row_Indices_0[Local_Row_Ind[1]] + Block_00->Row_Shift[1];
    COO_V[10] = Block_00->FE_Tensor_0[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[1]];
    COO_I[11] = Row_Indices_0[Local_Row_Ind[2]] + Block_00->Row_Shift[1];
    COO_V[11] = Block_00->FE_Tensor_0[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[2]];

    // write column #1
    // write (I,J,V) data for Block_00->FE_Tensor_1, i.e. the (1,1) block
    // actually:  copy Block_00->FE_Tensor_0 to Block_00->FE_Tensor_1
    // write the data directly
    COO_I[12] = Row_Indices_0[Local_Row_Ind[0]] + Block_00->Row_Shift[1];
    COO_J[4] = Col_Indices_0[Local_Col_Ind[1]] + Block_00->Col_Shift[1];
    COO_J_IV_Range[4] = 12;
    COO_V[12] = Block_00->FE_Tensor_0[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[0]];
    COO_I[13] = Row_Indices_0[Local_Row_Ind[1]] + Block_00->Row_Shift[1];
    COO_V[13] = Block_00->FE_Tensor_0[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[1]];
    COO_I[14] = Row_Indices_0[Local_Row_Ind[2]] + Block_00->Row_Shift[1];
    COO_V[14] = Block_00->FE_Tensor_0[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[2]];

    // write column #2
    // write (I,J,V) data for Block_00->FE_Tensor_1, i.e. the (1,1) block
    // actually:  copy Block_00->FE_Tensor_0 to Block_00->FE_Tensor_1
    // write the data directly
    COO_I[15] = Row_Indices_0[Local_Row_Ind[0]] + Block_00->Row_Shift[1];
    COO_J[5] = Col_Indices_0[Local_Col_Ind[2]] + Block_00->Col_Shift[1];
    COO_J_IV_Range[5] = 15;
    COO_V[15] = Block_00->FE_Tensor_0[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[0]];
    COO_I[16] = Row_Indices_0[Local_Row_Ind[1]] + Block_00->Row_Shift[1];
    COO_V[16] = Block_00->FE_Tensor_0[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[1]];
    COO_I[17] = Row_Indices_0[Local_Row_Ind[2]] + Block_00->Row_Shift[1];
    COO_V[17] = Block_00->FE_Tensor_0[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[2]];

    COO_J_IV_Range[6] = 18; // end of range
    // now insert into the matrix!
    MAT->add_entries(COO_I, COO_J, COO_J_IV_Range, COO_V, 2*COL_NB);

}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

// remove those macros!
#undef SpecificFEM
#undef SpecificFEM_str

#undef ROW_Num_Block
#undef ROW_NB
#undef COL_Num_Block
#undef COL_NB

/***/

