/*
============================================================================================
   This class points to outside MATLAB data which contains info about the
   finite element data.  It will also create a DoF numbering for that element.
   
   NOTE: portions of this code are automatically generated!
   
   Copyright (c) 10-16-2016,  Shawn W. Walker
============================================================================================
*/

/*------------ BEGIN: Auto Generate ------------*/
// set the type of element
#define ELEM_Name        "Elem4_Test"
// set allocator class name
#define EDA               Elem4_Test_DoF_Allocator
// set nodal topology data type name
#define NODAL_TOPOLOGY    Elem4_Test_nodal_top
// set permutation map data type names
#define SINGLE_EDGE_PERM_MAP    Elem4_Test_single_edge_perm_map
#define EDGE_DOF_PERMUTATION    Elem4_Test_edge_dof_perm
#define SINGLE_FACE_PERM_MAP    Elem4_Test_single_face_perm_map
#define FACE_DOF_PERMUTATION    Elem4_Test_face_dof_perm

// set the number of DoF sets
#define Num_Vtx_Sets        1
#define Num_Edge_Sets       0
#define Num_Face_Sets       1
#define Num_Tet_Sets        0

// set the max number (over all sets) of DoFs per entity per set
#define Max_DoF_Per_Vtx     1
#define Max_DoF_Per_Edge    0
#define Max_DoF_Per_Face    1
#define Max_DoF_Per_Tet     0

// set the total number of DoFs per entity
#define Num_DoF_Per_Vtx     1
#define Num_DoF_Per_Edge    0
#define Num_DoF_Per_Face    1
#define Num_DoF_Per_Tet     0

// set the TOTAL number of DoFs per cell (element)
#define Total_DoF_Per_Cell  4
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
// set (intrinsic) dimension and domain
#define Dimension    2
#define Domain_str  "triangle"

// set the number of vertices
#define NUM_VTX    3
// set the number of edges
#define NUM_EDGE   3
// set the number of faces
#define NUM_FACE   1
// set the number of tetrahedrons
#define NUM_TET    0
/*------------   END: Auto Generate ------------*/

/***************************************************************************************/
class EDA
{
public:
    char*    Name;              // name of finite element space
    int      Dim;               // intrinsic dimension
    char*    Domain;            // type of domain: "interval", "triangle", "tetrahedron"

    EDA ();  // constructor
    ~EDA (); // DE-structor
    mxArray* Init_DoF_Map(int);
    void  Fill_DoF_Map  (TRIANGLE_EDGE_SEARCH*);
    void  Error_Check_DoF_Map(const int&);

private:
    int*  cell_dof[Total_DoF_Per_Cell+1];
};

/***************************************************************************************/
/* constructor */
EDA::EDA ()
{
    // setup
    Name              = (char *) ELEM_Name;
    Dim               = Dimension;
    Domain            = (char *) Domain_str;
}
/***************************************************************************************/


/***************************************************************************************/
/* destructor */
EDA::~EDA ()
{
}
/***************************************************************************************/


/***************************************************************************************/
/* initialize FEM DoF map to all zeros */
mxArray* EDA::Init_DoF_Map(int  Num_Cell)     // output to MATLAB
{
    mxArray* DoF_Map;

    // BEGIN: allocate and access output data
    DoF_Map = mxCreateNumericMatrix((mwSize)Num_Cell, (mwSize)Total_DoF_Per_Cell, mxUINT32_CLASS, mxREAL);
    // access the data
    // split up the columns of the data
    cell_dof[0] = NULL; // not used!
    cell_dof[1] = (int *) mxGetPr(DoF_Map);
    for (int i = 2; (i <= Total_DoF_Per_Cell); i++) // note: off by one because of C-style indexing!
        cell_dof[i] = cell_dof[i-1] + Num_Cell;
    // END: allocate and access output data

    return DoF_Map;
}
/***************************************************************************************/


/***************************************************************************************/
/* assign DoFs (DG-style) */
void EDA::Fill_DoF_Map(TRIANGLE_EDGE_SEARCH*  Tri_Edge_Search)
{
    // output what we are doing
    printf("Generating DoFs on ");
    printf(Domain);
    printf("s for the Finite Element: ");
    printf(Name);
    printf(".\n");

    // write the local to global DoF mapping
    int DoF_Counter = 0;
    for (int ind=0; ind < Tri_Edge_Search->Num_Tri; ind++)
        {
        for (int i = 1; (i <= Total_DoF_Per_Cell); i++) // note: off by one because of C-style indexing!
            {
            DoF_Counter++;
            cell_dof[i][ind] = DoF_Counter;
            }
        }
    Error_Check_DoF_Map(Tri_Edge_Search->Num_Tri);
}
/***************************************************************************************/


/***************************************************************************************/
/* perform error checks on Finite Element Degree-of-Freedom (DoF) map */
void EDA::Error_Check_DoF_Map(const int&  Num_Cell)
{
    // BEGIN: error check
    const unsigned int LENGTH = Total_DoF_Per_Cell*Num_Cell;
    const unsigned int MIN_DoF = *min_element(cell_dof[1], cell_dof[1]+LENGTH);
    if (MIN_DoF < 1)
        {
        mexPrintf("ERROR: DoFmap references DoF with index 0!\n");
        mexPrintf("ERROR: Some DoFs were never allocated!\n");
        mexErrMsgTxt("ERROR: Make sure your mesh describes a domain that is a manifold!");
        }
    // END: error check
}
/***************************************************************************************/


#undef ELEM_Name
#undef EDA
#undef NODAL_TOPOLOGY
#undef SINGLE_EDGE_PERM_MAP
#undef EDGE_DOF_PERMUTATION
#undef SINGLE_FACE_PERM_MAP
#undef FACE_DOF_PERMUTATION
#undef Dimension
#undef Domain_str

#undef NUM_VTX
#undef NUM_EDGE
#undef NUM_FACE
#undef NUM_TET

#undef Num_Vtx_Sets
#undef Num_Edge_Sets
#undef Num_Face_Sets
#undef Num_Tet_Sets

#undef Max_DoF_Per_Vtx
#undef Max_DoF_Per_Edge
#undef Max_DoF_Per_Face
#undef Max_DoF_Per_Tet

#undef Num_DoF_Per_Vtx
#undef Num_DoF_Per_Edge
#undef Num_DoF_Per_Face
#undef Num_DoF_Per_Tet

#undef Total_DoF_Per_Cell

/***/