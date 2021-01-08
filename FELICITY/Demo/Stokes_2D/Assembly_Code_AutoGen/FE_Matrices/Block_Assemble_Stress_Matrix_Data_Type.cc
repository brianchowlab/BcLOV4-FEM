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
//Stress_Matrix
//
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
// define the name of the block FE matrix (should be the same as the filename of this file)
#define SpecificFEM        Block_Assemble_Stress_Matrix_Data_Type
#define SpecificFEM_str   "Stress_Matrix"

// the row function space is Type = CG, Name = "lagrange_deg2_dim2"

// set the number of blocks along the row dimension
#define ROW_Num_Block  1
// set the number of basis functions on each element
#define ROW_NB         6

// the col function space is Type = CG, Name = "lagrange_deg2_dim2"

// set the number of blocks along the col dimension
#define COL_Num_Block  1
// set the number of basis functions on each element
#define COL_NB  6
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
    SpecificFEM (const PTR_TO_SPARSE*, const Base_Stress_Matrix_Data_Type*);
    ~SpecificFEM (); // destructor
    void Init_Matrix_Assembler_Object(bool);
    void Add_Entries_To_Global_Matrix_Omega (const Base_Stress_Matrix_Data_Type*);
    /*------------   END: Auto Generate ------------*/
private:
};


/***************************************************************************************/
/* constructor */
/*------------ BEGIN: Auto Generate ------------*/
SpecificFEM::SpecificFEM (const PTR_TO_SPARSE* Prev_Sparse_Data,
                          const Base_Stress_Matrix_Data_Type* Block_00
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
/* assemble a local FE matrix on the domain Omega */
void SpecificFEM::Add_Entries_To_Global_Matrix_Omega(const Base_Stress_Matrix_Data_Type* Block_00)
{
    // get local to global index map for the current ROW element
    int  Row_Indices_0[ROW_NB];
    const int row_elem_index = Block_00->Vector_P2_phi_restricted_to_Omega->Mesh->Domain->Sub_Cell_Index;
    Block_00->Vector_P2_phi_restricted_to_Omega->Get_Local_to_Global_DoFmap(row_elem_index, Row_Indices_0);
    // shift Row_Indices_0 to account for the block matrix offset
    for (unsigned int ri = 0; ri < ROW_NB; ++ri)
        Row_Indices_0[ri] += Block_Row_Shift[0];

    // get local to global index map for the current COL element
    int  Col_Indices_0[COL_NB];
    const int col_elem_index = Block_00->Vector_P2_phi_restricted_to_Omega->Mesh->Domain->Sub_Cell_Index;
    Block_00->Vector_P2_phi_restricted_to_Omega->Get_Local_to_Global_DoFmap(col_elem_index, Col_Indices_0);
    // shift Col_Indices_0 to account for the block matrix offset
    for (unsigned int ci = 0; ci < COL_NB; ++ci)
        Col_Indices_0[ci] += Block_Col_Shift[0];

    // sort row indices (ascending order)
    int  Local_Row_Ind[ROW_NB] = {0, 1, 2, 3, 4, 5};
    std::sort(Local_Row_Ind, Local_Row_Ind+ROW_NB,
              [&Row_Indices_0](int kk, int qq) { return (Row_Indices_0[kk] < Row_Indices_0[qq]); });

    // sort col indices (ascending order)
    int  Local_Col_Ind[COL_NB] = {0, 1, 2, 3, 4, 5};
    std::sort(Local_Col_Ind, Local_Col_Ind+COL_NB,
              [&Col_Indices_0](int kk, int qq) { return (Col_Indices_0[kk] < Col_Indices_0[qq]); });

    // allocate (I,J,V) arrays to hold "big" local matrix
    int     COO_I[4*ROW_NB*COL_NB];
    int     COO_J[2*COL_NB];
    int     COO_J_IV_Range[2*COL_NB + 1]; // indicates what parts I,V correspond to J
    double  COO_V[4*ROW_NB*COL_NB];

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
    COO_I[3] = Row_Indices_0[Local_Row_Ind[3]];
    COO_V[3] = Block_00->FE_Tensor_0[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[3]];
    COO_I[4] = Row_Indices_0[Local_Row_Ind[4]];
    COO_V[4] = Block_00->FE_Tensor_0[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[4]];
    COO_I[5] = Row_Indices_0[Local_Row_Ind[5]];
    COO_V[5] = Block_00->FE_Tensor_0[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[5]];

    // write (I,J,V) data for Block_00->FE_Tensor_1, i.e. the (1,0) block
    // write the data directly
    COO_I[6] = Row_Indices_0[Local_Row_Ind[0]] + Block_00->Row_Shift[1];
    COO_V[6] = Block_00->FE_Tensor_1[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[0]];
    COO_I[7] = Row_Indices_0[Local_Row_Ind[1]] + Block_00->Row_Shift[1];
    COO_V[7] = Block_00->FE_Tensor_1[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[1]];
    COO_I[8] = Row_Indices_0[Local_Row_Ind[2]] + Block_00->Row_Shift[1];
    COO_V[8] = Block_00->FE_Tensor_1[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[2]];
    COO_I[9] = Row_Indices_0[Local_Row_Ind[3]] + Block_00->Row_Shift[1];
    COO_V[9] = Block_00->FE_Tensor_1[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[3]];
    COO_I[10] = Row_Indices_0[Local_Row_Ind[4]] + Block_00->Row_Shift[1];
    COO_V[10] = Block_00->FE_Tensor_1[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[4]];
    COO_I[11] = Row_Indices_0[Local_Row_Ind[5]] + Block_00->Row_Shift[1];
    COO_V[11] = Block_00->FE_Tensor_1[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[5]];

    // write column #1
    // write (I,J,V) data for Block_00->FE_Tensor_0, i.e. the (0,0) block
    // write the data directly
    COO_I[12] = Row_Indices_0[Local_Row_Ind[0]];
    COO_J[1] = Col_Indices_0[Local_Col_Ind[1]];
    COO_J_IV_Range[1] = 12;
    COO_V[12] = Block_00->FE_Tensor_0[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[0]];
    COO_I[13] = Row_Indices_0[Local_Row_Ind[1]];
    COO_V[13] = Block_00->FE_Tensor_0[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[1]];
    COO_I[14] = Row_Indices_0[Local_Row_Ind[2]];
    COO_V[14] = Block_00->FE_Tensor_0[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[2]];
    COO_I[15] = Row_Indices_0[Local_Row_Ind[3]];
    COO_V[15] = Block_00->FE_Tensor_0[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[3]];
    COO_I[16] = Row_Indices_0[Local_Row_Ind[4]];
    COO_V[16] = Block_00->FE_Tensor_0[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[4]];
    COO_I[17] = Row_Indices_0[Local_Row_Ind[5]];
    COO_V[17] = Block_00->FE_Tensor_0[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[5]];

    // write (I,J,V) data for Block_00->FE_Tensor_1, i.e. the (1,0) block
    // write the data directly
    COO_I[18] = Row_Indices_0[Local_Row_Ind[0]] + Block_00->Row_Shift[1];
    COO_V[18] = Block_00->FE_Tensor_1[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[0]];
    COO_I[19] = Row_Indices_0[Local_Row_Ind[1]] + Block_00->Row_Shift[1];
    COO_V[19] = Block_00->FE_Tensor_1[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[1]];
    COO_I[20] = Row_Indices_0[Local_Row_Ind[2]] + Block_00->Row_Shift[1];
    COO_V[20] = Block_00->FE_Tensor_1[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[2]];
    COO_I[21] = Row_Indices_0[Local_Row_Ind[3]] + Block_00->Row_Shift[1];
    COO_V[21] = Block_00->FE_Tensor_1[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[3]];
    COO_I[22] = Row_Indices_0[Local_Row_Ind[4]] + Block_00->Row_Shift[1];
    COO_V[22] = Block_00->FE_Tensor_1[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[4]];
    COO_I[23] = Row_Indices_0[Local_Row_Ind[5]] + Block_00->Row_Shift[1];
    COO_V[23] = Block_00->FE_Tensor_1[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[5]];

    // write column #2
    // write (I,J,V) data for Block_00->FE_Tensor_0, i.e. the (0,0) block
    // write the data directly
    COO_I[24] = Row_Indices_0[Local_Row_Ind[0]];
    COO_J[2] = Col_Indices_0[Local_Col_Ind[2]];
    COO_J_IV_Range[2] = 24;
    COO_V[24] = Block_00->FE_Tensor_0[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[0]];
    COO_I[25] = Row_Indices_0[Local_Row_Ind[1]];
    COO_V[25] = Block_00->FE_Tensor_0[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[1]];
    COO_I[26] = Row_Indices_0[Local_Row_Ind[2]];
    COO_V[26] = Block_00->FE_Tensor_0[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[2]];
    COO_I[27] = Row_Indices_0[Local_Row_Ind[3]];
    COO_V[27] = Block_00->FE_Tensor_0[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[3]];
    COO_I[28] = Row_Indices_0[Local_Row_Ind[4]];
    COO_V[28] = Block_00->FE_Tensor_0[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[4]];
    COO_I[29] = Row_Indices_0[Local_Row_Ind[5]];
    COO_V[29] = Block_00->FE_Tensor_0[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[5]];

    // write (I,J,V) data for Block_00->FE_Tensor_1, i.e. the (1,0) block
    // write the data directly
    COO_I[30] = Row_Indices_0[Local_Row_Ind[0]] + Block_00->Row_Shift[1];
    COO_V[30] = Block_00->FE_Tensor_1[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[0]];
    COO_I[31] = Row_Indices_0[Local_Row_Ind[1]] + Block_00->Row_Shift[1];
    COO_V[31] = Block_00->FE_Tensor_1[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[1]];
    COO_I[32] = Row_Indices_0[Local_Row_Ind[2]] + Block_00->Row_Shift[1];
    COO_V[32] = Block_00->FE_Tensor_1[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[2]];
    COO_I[33] = Row_Indices_0[Local_Row_Ind[3]] + Block_00->Row_Shift[1];
    COO_V[33] = Block_00->FE_Tensor_1[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[3]];
    COO_I[34] = Row_Indices_0[Local_Row_Ind[4]] + Block_00->Row_Shift[1];
    COO_V[34] = Block_00->FE_Tensor_1[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[4]];
    COO_I[35] = Row_Indices_0[Local_Row_Ind[5]] + Block_00->Row_Shift[1];
    COO_V[35] = Block_00->FE_Tensor_1[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[5]];

    // write column #3
    // write (I,J,V) data for Block_00->FE_Tensor_0, i.e. the (0,0) block
    // write the data directly
    COO_I[36] = Row_Indices_0[Local_Row_Ind[0]];
    COO_J[3] = Col_Indices_0[Local_Col_Ind[3]];
    COO_J_IV_Range[3] = 36;
    COO_V[36] = Block_00->FE_Tensor_0[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[0]];
    COO_I[37] = Row_Indices_0[Local_Row_Ind[1]];
    COO_V[37] = Block_00->FE_Tensor_0[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[1]];
    COO_I[38] = Row_Indices_0[Local_Row_Ind[2]];
    COO_V[38] = Block_00->FE_Tensor_0[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[2]];
    COO_I[39] = Row_Indices_0[Local_Row_Ind[3]];
    COO_V[39] = Block_00->FE_Tensor_0[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[3]];
    COO_I[40] = Row_Indices_0[Local_Row_Ind[4]];
    COO_V[40] = Block_00->FE_Tensor_0[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[4]];
    COO_I[41] = Row_Indices_0[Local_Row_Ind[5]];
    COO_V[41] = Block_00->FE_Tensor_0[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[5]];

    // write (I,J,V) data for Block_00->FE_Tensor_1, i.e. the (1,0) block
    // write the data directly
    COO_I[42] = Row_Indices_0[Local_Row_Ind[0]] + Block_00->Row_Shift[1];
    COO_V[42] = Block_00->FE_Tensor_1[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[0]];
    COO_I[43] = Row_Indices_0[Local_Row_Ind[1]] + Block_00->Row_Shift[1];
    COO_V[43] = Block_00->FE_Tensor_1[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[1]];
    COO_I[44] = Row_Indices_0[Local_Row_Ind[2]] + Block_00->Row_Shift[1];
    COO_V[44] = Block_00->FE_Tensor_1[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[2]];
    COO_I[45] = Row_Indices_0[Local_Row_Ind[3]] + Block_00->Row_Shift[1];
    COO_V[45] = Block_00->FE_Tensor_1[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[3]];
    COO_I[46] = Row_Indices_0[Local_Row_Ind[4]] + Block_00->Row_Shift[1];
    COO_V[46] = Block_00->FE_Tensor_1[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[4]];
    COO_I[47] = Row_Indices_0[Local_Row_Ind[5]] + Block_00->Row_Shift[1];
    COO_V[47] = Block_00->FE_Tensor_1[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[5]];

    // write column #4
    // write (I,J,V) data for Block_00->FE_Tensor_0, i.e. the (0,0) block
    // write the data directly
    COO_I[48] = Row_Indices_0[Local_Row_Ind[0]];
    COO_J[4] = Col_Indices_0[Local_Col_Ind[4]];
    COO_J_IV_Range[4] = 48;
    COO_V[48] = Block_00->FE_Tensor_0[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[0]];
    COO_I[49] = Row_Indices_0[Local_Row_Ind[1]];
    COO_V[49] = Block_00->FE_Tensor_0[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[1]];
    COO_I[50] = Row_Indices_0[Local_Row_Ind[2]];
    COO_V[50] = Block_00->FE_Tensor_0[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[2]];
    COO_I[51] = Row_Indices_0[Local_Row_Ind[3]];
    COO_V[51] = Block_00->FE_Tensor_0[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[3]];
    COO_I[52] = Row_Indices_0[Local_Row_Ind[4]];
    COO_V[52] = Block_00->FE_Tensor_0[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[4]];
    COO_I[53] = Row_Indices_0[Local_Row_Ind[5]];
    COO_V[53] = Block_00->FE_Tensor_0[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[5]];

    // write (I,J,V) data for Block_00->FE_Tensor_1, i.e. the (1,0) block
    // write the data directly
    COO_I[54] = Row_Indices_0[Local_Row_Ind[0]] + Block_00->Row_Shift[1];
    COO_V[54] = Block_00->FE_Tensor_1[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[0]];
    COO_I[55] = Row_Indices_0[Local_Row_Ind[1]] + Block_00->Row_Shift[1];
    COO_V[55] = Block_00->FE_Tensor_1[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[1]];
    COO_I[56] = Row_Indices_0[Local_Row_Ind[2]] + Block_00->Row_Shift[1];
    COO_V[56] = Block_00->FE_Tensor_1[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[2]];
    COO_I[57] = Row_Indices_0[Local_Row_Ind[3]] + Block_00->Row_Shift[1];
    COO_V[57] = Block_00->FE_Tensor_1[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[3]];
    COO_I[58] = Row_Indices_0[Local_Row_Ind[4]] + Block_00->Row_Shift[1];
    COO_V[58] = Block_00->FE_Tensor_1[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[4]];
    COO_I[59] = Row_Indices_0[Local_Row_Ind[5]] + Block_00->Row_Shift[1];
    COO_V[59] = Block_00->FE_Tensor_1[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[5]];

    // write column #5
    // write (I,J,V) data for Block_00->FE_Tensor_0, i.e. the (0,0) block
    // write the data directly
    COO_I[60] = Row_Indices_0[Local_Row_Ind[0]];
    COO_J[5] = Col_Indices_0[Local_Col_Ind[5]];
    COO_J_IV_Range[5] = 60;
    COO_V[60] = Block_00->FE_Tensor_0[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[0]];
    COO_I[61] = Row_Indices_0[Local_Row_Ind[1]];
    COO_V[61] = Block_00->FE_Tensor_0[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[1]];
    COO_I[62] = Row_Indices_0[Local_Row_Ind[2]];
    COO_V[62] = Block_00->FE_Tensor_0[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[2]];
    COO_I[63] = Row_Indices_0[Local_Row_Ind[3]];
    COO_V[63] = Block_00->FE_Tensor_0[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[3]];
    COO_I[64] = Row_Indices_0[Local_Row_Ind[4]];
    COO_V[64] = Block_00->FE_Tensor_0[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[4]];
    COO_I[65] = Row_Indices_0[Local_Row_Ind[5]];
    COO_V[65] = Block_00->FE_Tensor_0[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[5]];

    // write (I,J,V) data for Block_00->FE_Tensor_1, i.e. the (1,0) block
    // write the data directly
    COO_I[66] = Row_Indices_0[Local_Row_Ind[0]] + Block_00->Row_Shift[1];
    COO_V[66] = Block_00->FE_Tensor_1[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[0]];
    COO_I[67] = Row_Indices_0[Local_Row_Ind[1]] + Block_00->Row_Shift[1];
    COO_V[67] = Block_00->FE_Tensor_1[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[1]];
    COO_I[68] = Row_Indices_0[Local_Row_Ind[2]] + Block_00->Row_Shift[1];
    COO_V[68] = Block_00->FE_Tensor_1[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[2]];
    COO_I[69] = Row_Indices_0[Local_Row_Ind[3]] + Block_00->Row_Shift[1];
    COO_V[69] = Block_00->FE_Tensor_1[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[3]];
    COO_I[70] = Row_Indices_0[Local_Row_Ind[4]] + Block_00->Row_Shift[1];
    COO_V[70] = Block_00->FE_Tensor_1[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[4]];
    COO_I[71] = Row_Indices_0[Local_Row_Ind[5]] + Block_00->Row_Shift[1];
    COO_V[71] = Block_00->FE_Tensor_1[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[5]];

    // write column #0
    // write (I,J,V) data for Block_00->FE_Tensor_2, i.e. the (0,1) block
    // actually:  copy Block_00->FE_Tensor_1 to Block_00->FE_Tensor_2
    // write the data directly (but take the transpose!)
    COO_I[72] = Row_Indices_0[Local_Row_Ind[0]];
    COO_J[6] = Col_Indices_0[Local_Col_Ind[0]] + Block_00->Col_Shift[1];
    COO_J_IV_Range[6] = 72;
    COO_V[72] = Block_00->FE_Tensor_1[Local_Row_Ind[0]*COL_NB + Local_Col_Ind[0]];
    COO_I[73] = Row_Indices_0[Local_Row_Ind[1]];
    COO_V[73] = Block_00->FE_Tensor_1[Local_Row_Ind[1]*COL_NB + Local_Col_Ind[0]];
    COO_I[74] = Row_Indices_0[Local_Row_Ind[2]];
    COO_V[74] = Block_00->FE_Tensor_1[Local_Row_Ind[2]*COL_NB + Local_Col_Ind[0]];
    COO_I[75] = Row_Indices_0[Local_Row_Ind[3]];
    COO_V[75] = Block_00->FE_Tensor_1[Local_Row_Ind[3]*COL_NB + Local_Col_Ind[0]];
    COO_I[76] = Row_Indices_0[Local_Row_Ind[4]];
    COO_V[76] = Block_00->FE_Tensor_1[Local_Row_Ind[4]*COL_NB + Local_Col_Ind[0]];
    COO_I[77] = Row_Indices_0[Local_Row_Ind[5]];
    COO_V[77] = Block_00->FE_Tensor_1[Local_Row_Ind[5]*COL_NB + Local_Col_Ind[0]];

    // write (I,J,V) data for Block_00->FE_Tensor_3, i.e. the (1,1) block
    // write the data directly
    COO_I[78] = Row_Indices_0[Local_Row_Ind[0]] + Block_00->Row_Shift[1];
    COO_V[78] = Block_00->FE_Tensor_3[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[0]];
    COO_I[79] = Row_Indices_0[Local_Row_Ind[1]] + Block_00->Row_Shift[1];
    COO_V[79] = Block_00->FE_Tensor_3[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[1]];
    COO_I[80] = Row_Indices_0[Local_Row_Ind[2]] + Block_00->Row_Shift[1];
    COO_V[80] = Block_00->FE_Tensor_3[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[2]];
    COO_I[81] = Row_Indices_0[Local_Row_Ind[3]] + Block_00->Row_Shift[1];
    COO_V[81] = Block_00->FE_Tensor_3[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[3]];
    COO_I[82] = Row_Indices_0[Local_Row_Ind[4]] + Block_00->Row_Shift[1];
    COO_V[82] = Block_00->FE_Tensor_3[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[4]];
    COO_I[83] = Row_Indices_0[Local_Row_Ind[5]] + Block_00->Row_Shift[1];
    COO_V[83] = Block_00->FE_Tensor_3[Local_Col_Ind[0]*ROW_NB + Local_Row_Ind[5]];

    // write column #1
    // write (I,J,V) data for Block_00->FE_Tensor_2, i.e. the (0,1) block
    // actually:  copy Block_00->FE_Tensor_1 to Block_00->FE_Tensor_2
    // write the data directly (but take the transpose!)
    COO_I[84] = Row_Indices_0[Local_Row_Ind[0]];
    COO_J[7] = Col_Indices_0[Local_Col_Ind[1]] + Block_00->Col_Shift[1];
    COO_J_IV_Range[7] = 84;
    COO_V[84] = Block_00->FE_Tensor_1[Local_Row_Ind[0]*COL_NB + Local_Col_Ind[1]];
    COO_I[85] = Row_Indices_0[Local_Row_Ind[1]];
    COO_V[85] = Block_00->FE_Tensor_1[Local_Row_Ind[1]*COL_NB + Local_Col_Ind[1]];
    COO_I[86] = Row_Indices_0[Local_Row_Ind[2]];
    COO_V[86] = Block_00->FE_Tensor_1[Local_Row_Ind[2]*COL_NB + Local_Col_Ind[1]];
    COO_I[87] = Row_Indices_0[Local_Row_Ind[3]];
    COO_V[87] = Block_00->FE_Tensor_1[Local_Row_Ind[3]*COL_NB + Local_Col_Ind[1]];
    COO_I[88] = Row_Indices_0[Local_Row_Ind[4]];
    COO_V[88] = Block_00->FE_Tensor_1[Local_Row_Ind[4]*COL_NB + Local_Col_Ind[1]];
    COO_I[89] = Row_Indices_0[Local_Row_Ind[5]];
    COO_V[89] = Block_00->FE_Tensor_1[Local_Row_Ind[5]*COL_NB + Local_Col_Ind[1]];

    // write (I,J,V) data for Block_00->FE_Tensor_3, i.e. the (1,1) block
    // write the data directly
    COO_I[90] = Row_Indices_0[Local_Row_Ind[0]] + Block_00->Row_Shift[1];
    COO_V[90] = Block_00->FE_Tensor_3[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[0]];
    COO_I[91] = Row_Indices_0[Local_Row_Ind[1]] + Block_00->Row_Shift[1];
    COO_V[91] = Block_00->FE_Tensor_3[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[1]];
    COO_I[92] = Row_Indices_0[Local_Row_Ind[2]] + Block_00->Row_Shift[1];
    COO_V[92] = Block_00->FE_Tensor_3[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[2]];
    COO_I[93] = Row_Indices_0[Local_Row_Ind[3]] + Block_00->Row_Shift[1];
    COO_V[93] = Block_00->FE_Tensor_3[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[3]];
    COO_I[94] = Row_Indices_0[Local_Row_Ind[4]] + Block_00->Row_Shift[1];
    COO_V[94] = Block_00->FE_Tensor_3[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[4]];
    COO_I[95] = Row_Indices_0[Local_Row_Ind[5]] + Block_00->Row_Shift[1];
    COO_V[95] = Block_00->FE_Tensor_3[Local_Col_Ind[1]*ROW_NB + Local_Row_Ind[5]];

    // write column #2
    // write (I,J,V) data for Block_00->FE_Tensor_2, i.e. the (0,1) block
    // actually:  copy Block_00->FE_Tensor_1 to Block_00->FE_Tensor_2
    // write the data directly (but take the transpose!)
    COO_I[96] = Row_Indices_0[Local_Row_Ind[0]];
    COO_J[8] = Col_Indices_0[Local_Col_Ind[2]] + Block_00->Col_Shift[1];
    COO_J_IV_Range[8] = 96;
    COO_V[96] = Block_00->FE_Tensor_1[Local_Row_Ind[0]*COL_NB + Local_Col_Ind[2]];
    COO_I[97] = Row_Indices_0[Local_Row_Ind[1]];
    COO_V[97] = Block_00->FE_Tensor_1[Local_Row_Ind[1]*COL_NB + Local_Col_Ind[2]];
    COO_I[98] = Row_Indices_0[Local_Row_Ind[2]];
    COO_V[98] = Block_00->FE_Tensor_1[Local_Row_Ind[2]*COL_NB + Local_Col_Ind[2]];
    COO_I[99] = Row_Indices_0[Local_Row_Ind[3]];
    COO_V[99] = Block_00->FE_Tensor_1[Local_Row_Ind[3]*COL_NB + Local_Col_Ind[2]];
    COO_I[100] = Row_Indices_0[Local_Row_Ind[4]];
    COO_V[100] = Block_00->FE_Tensor_1[Local_Row_Ind[4]*COL_NB + Local_Col_Ind[2]];
    COO_I[101] = Row_Indices_0[Local_Row_Ind[5]];
    COO_V[101] = Block_00->FE_Tensor_1[Local_Row_Ind[5]*COL_NB + Local_Col_Ind[2]];

    // write (I,J,V) data for Block_00->FE_Tensor_3, i.e. the (1,1) block
    // write the data directly
    COO_I[102] = Row_Indices_0[Local_Row_Ind[0]] + Block_00->Row_Shift[1];
    COO_V[102] = Block_00->FE_Tensor_3[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[0]];
    COO_I[103] = Row_Indices_0[Local_Row_Ind[1]] + Block_00->Row_Shift[1];
    COO_V[103] = Block_00->FE_Tensor_3[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[1]];
    COO_I[104] = Row_Indices_0[Local_Row_Ind[2]] + Block_00->Row_Shift[1];
    COO_V[104] = Block_00->FE_Tensor_3[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[2]];
    COO_I[105] = Row_Indices_0[Local_Row_Ind[3]] + Block_00->Row_Shift[1];
    COO_V[105] = Block_00->FE_Tensor_3[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[3]];
    COO_I[106] = Row_Indices_0[Local_Row_Ind[4]] + Block_00->Row_Shift[1];
    COO_V[106] = Block_00->FE_Tensor_3[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[4]];
    COO_I[107] = Row_Indices_0[Local_Row_Ind[5]] + Block_00->Row_Shift[1];
    COO_V[107] = Block_00->FE_Tensor_3[Local_Col_Ind[2]*ROW_NB + Local_Row_Ind[5]];

    // write column #3
    // write (I,J,V) data for Block_00->FE_Tensor_2, i.e. the (0,1) block
    // actually:  copy Block_00->FE_Tensor_1 to Block_00->FE_Tensor_2
    // write the data directly (but take the transpose!)
    COO_I[108] = Row_Indices_0[Local_Row_Ind[0]];
    COO_J[9] = Col_Indices_0[Local_Col_Ind[3]] + Block_00->Col_Shift[1];
    COO_J_IV_Range[9] = 108;
    COO_V[108] = Block_00->FE_Tensor_1[Local_Row_Ind[0]*COL_NB + Local_Col_Ind[3]];
    COO_I[109] = Row_Indices_0[Local_Row_Ind[1]];
    COO_V[109] = Block_00->FE_Tensor_1[Local_Row_Ind[1]*COL_NB + Local_Col_Ind[3]];
    COO_I[110] = Row_Indices_0[Local_Row_Ind[2]];
    COO_V[110] = Block_00->FE_Tensor_1[Local_Row_Ind[2]*COL_NB + Local_Col_Ind[3]];
    COO_I[111] = Row_Indices_0[Local_Row_Ind[3]];
    COO_V[111] = Block_00->FE_Tensor_1[Local_Row_Ind[3]*COL_NB + Local_Col_Ind[3]];
    COO_I[112] = Row_Indices_0[Local_Row_Ind[4]];
    COO_V[112] = Block_00->FE_Tensor_1[Local_Row_Ind[4]*COL_NB + Local_Col_Ind[3]];
    COO_I[113] = Row_Indices_0[Local_Row_Ind[5]];
    COO_V[113] = Block_00->FE_Tensor_1[Local_Row_Ind[5]*COL_NB + Local_Col_Ind[3]];

    // write (I,J,V) data for Block_00->FE_Tensor_3, i.e. the (1,1) block
    // write the data directly
    COO_I[114] = Row_Indices_0[Local_Row_Ind[0]] + Block_00->Row_Shift[1];
    COO_V[114] = Block_00->FE_Tensor_3[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[0]];
    COO_I[115] = Row_Indices_0[Local_Row_Ind[1]] + Block_00->Row_Shift[1];
    COO_V[115] = Block_00->FE_Tensor_3[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[1]];
    COO_I[116] = Row_Indices_0[Local_Row_Ind[2]] + Block_00->Row_Shift[1];
    COO_V[116] = Block_00->FE_Tensor_3[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[2]];
    COO_I[117] = Row_Indices_0[Local_Row_Ind[3]] + Block_00->Row_Shift[1];
    COO_V[117] = Block_00->FE_Tensor_3[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[3]];
    COO_I[118] = Row_Indices_0[Local_Row_Ind[4]] + Block_00->Row_Shift[1];
    COO_V[118] = Block_00->FE_Tensor_3[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[4]];
    COO_I[119] = Row_Indices_0[Local_Row_Ind[5]] + Block_00->Row_Shift[1];
    COO_V[119] = Block_00->FE_Tensor_3[Local_Col_Ind[3]*ROW_NB + Local_Row_Ind[5]];

    // write column #4
    // write (I,J,V) data for Block_00->FE_Tensor_2, i.e. the (0,1) block
    // actually:  copy Block_00->FE_Tensor_1 to Block_00->FE_Tensor_2
    // write the data directly (but take the transpose!)
    COO_I[120] = Row_Indices_0[Local_Row_Ind[0]];
    COO_J[10] = Col_Indices_0[Local_Col_Ind[4]] + Block_00->Col_Shift[1];
    COO_J_IV_Range[10] = 120;
    COO_V[120] = Block_00->FE_Tensor_1[Local_Row_Ind[0]*COL_NB + Local_Col_Ind[4]];
    COO_I[121] = Row_Indices_0[Local_Row_Ind[1]];
    COO_V[121] = Block_00->FE_Tensor_1[Local_Row_Ind[1]*COL_NB + Local_Col_Ind[4]];
    COO_I[122] = Row_Indices_0[Local_Row_Ind[2]];
    COO_V[122] = Block_00->FE_Tensor_1[Local_Row_Ind[2]*COL_NB + Local_Col_Ind[4]];
    COO_I[123] = Row_Indices_0[Local_Row_Ind[3]];
    COO_V[123] = Block_00->FE_Tensor_1[Local_Row_Ind[3]*COL_NB + Local_Col_Ind[4]];
    COO_I[124] = Row_Indices_0[Local_Row_Ind[4]];
    COO_V[124] = Block_00->FE_Tensor_1[Local_Row_Ind[4]*COL_NB + Local_Col_Ind[4]];
    COO_I[125] = Row_Indices_0[Local_Row_Ind[5]];
    COO_V[125] = Block_00->FE_Tensor_1[Local_Row_Ind[5]*COL_NB + Local_Col_Ind[4]];

    // write (I,J,V) data for Block_00->FE_Tensor_3, i.e. the (1,1) block
    // write the data directly
    COO_I[126] = Row_Indices_0[Local_Row_Ind[0]] + Block_00->Row_Shift[1];
    COO_V[126] = Block_00->FE_Tensor_3[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[0]];
    COO_I[127] = Row_Indices_0[Local_Row_Ind[1]] + Block_00->Row_Shift[1];
    COO_V[127] = Block_00->FE_Tensor_3[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[1]];
    COO_I[128] = Row_Indices_0[Local_Row_Ind[2]] + Block_00->Row_Shift[1];
    COO_V[128] = Block_00->FE_Tensor_3[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[2]];
    COO_I[129] = Row_Indices_0[Local_Row_Ind[3]] + Block_00->Row_Shift[1];
    COO_V[129] = Block_00->FE_Tensor_3[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[3]];
    COO_I[130] = Row_Indices_0[Local_Row_Ind[4]] + Block_00->Row_Shift[1];
    COO_V[130] = Block_00->FE_Tensor_3[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[4]];
    COO_I[131] = Row_Indices_0[Local_Row_Ind[5]] + Block_00->Row_Shift[1];
    COO_V[131] = Block_00->FE_Tensor_3[Local_Col_Ind[4]*ROW_NB + Local_Row_Ind[5]];

    // write column #5
    // write (I,J,V) data for Block_00->FE_Tensor_2, i.e. the (0,1) block
    // actually:  copy Block_00->FE_Tensor_1 to Block_00->FE_Tensor_2
    // write the data directly (but take the transpose!)
    COO_I[132] = Row_Indices_0[Local_Row_Ind[0]];
    COO_J[11] = Col_Indices_0[Local_Col_Ind[5]] + Block_00->Col_Shift[1];
    COO_J_IV_Range[11] = 132;
    COO_V[132] = Block_00->FE_Tensor_1[Local_Row_Ind[0]*COL_NB + Local_Col_Ind[5]];
    COO_I[133] = Row_Indices_0[Local_Row_Ind[1]];
    COO_V[133] = Block_00->FE_Tensor_1[Local_Row_Ind[1]*COL_NB + Local_Col_Ind[5]];
    COO_I[134] = Row_Indices_0[Local_Row_Ind[2]];
    COO_V[134] = Block_00->FE_Tensor_1[Local_Row_Ind[2]*COL_NB + Local_Col_Ind[5]];
    COO_I[135] = Row_Indices_0[Local_Row_Ind[3]];
    COO_V[135] = Block_00->FE_Tensor_1[Local_Row_Ind[3]*COL_NB + Local_Col_Ind[5]];
    COO_I[136] = Row_Indices_0[Local_Row_Ind[4]];
    COO_V[136] = Block_00->FE_Tensor_1[Local_Row_Ind[4]*COL_NB + Local_Col_Ind[5]];
    COO_I[137] = Row_Indices_0[Local_Row_Ind[5]];
    COO_V[137] = Block_00->FE_Tensor_1[Local_Row_Ind[5]*COL_NB + Local_Col_Ind[5]];

    // write (I,J,V) data for Block_00->FE_Tensor_3, i.e. the (1,1) block
    // write the data directly
    COO_I[138] = Row_Indices_0[Local_Row_Ind[0]] + Block_00->Row_Shift[1];
    COO_V[138] = Block_00->FE_Tensor_3[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[0]];
    COO_I[139] = Row_Indices_0[Local_Row_Ind[1]] + Block_00->Row_Shift[1];
    COO_V[139] = Block_00->FE_Tensor_3[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[1]];
    COO_I[140] = Row_Indices_0[Local_Row_Ind[2]] + Block_00->Row_Shift[1];
    COO_V[140] = Block_00->FE_Tensor_3[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[2]];
    COO_I[141] = Row_Indices_0[Local_Row_Ind[3]] + Block_00->Row_Shift[1];
    COO_V[141] = Block_00->FE_Tensor_3[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[3]];
    COO_I[142] = Row_Indices_0[Local_Row_Ind[4]] + Block_00->Row_Shift[1];
    COO_V[142] = Block_00->FE_Tensor_3[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[4]];
    COO_I[143] = Row_Indices_0[Local_Row_Ind[5]] + Block_00->Row_Shift[1];
    COO_V[143] = Block_00->FE_Tensor_3[Local_Col_Ind[5]*ROW_NB + Local_Row_Ind[5]];

    COO_J_IV_Range[12] = 144; // end of range
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

