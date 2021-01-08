/*
============================================================================================
   This class contains data about a given FE function space, and methods for computing
   transformations of the local basis functions.

   This code references the header files:

   matrix_vector_defn.h
   matrix_vector_ops.h
   geometric_computations.h
   basis_function_computations.h


   NOTE: portions of this code are automatically generated!

   Copyright (c) 01-15-2018,  Shawn W. Walker
============================================================================================
*/

/*------------ BEGIN: Auto Generate ------------*/
// define the name of the FE basis function (should be the same as the filename of this file)
#define SpecificFUNC        Data_Type_Ned2_phi_restricted_to_Omega
#define SpecificFUNC_str   "Data_Type_Ned2_phi_restricted_to_Omega"

// set the type of function space
#define SPACE_type  "CG - nedelec_1stkind_deg2_dim3"
// set the name of function space
#define SPACE_name  "Ned2"

// set the Subdomain topological dimension
#define SUB_TD  3
// set the Domain of Integration (DoI) topological dimension
#define DOI_TD  3
// set the geometric dimension
#define GD  3
// set the number of cartesian tuple components (m*n) = 1 * 1
#define NC  1
// NOTE: the (i,j) tuple component is accessed by the linear index k = i + (j-1)*m
// set the number of quad points
#define NQ  10
// set the number of basis functions
#define NB  20
/*------------   END: Auto Generate ------------*/

/* C++ (Specific) FE Function class definition */
class SpecificFUNC: public ABSTRACT_FEM_Function_Class // derive from base class
{
public:
    int*     Elem_DoF[NB];    // element DoF list

    // data structure containing information on the function evaluations.
    // Note: this data is evaluated at several quadrature points!
    // vector valued H(curl) basis functions
    VEC_3x1 Func_vv_Value[NB][NQ];
    // curl of vector valued H(curl) basis functions (curl is a vector in 3-D)
    VEC_3x1 Func_vv_Curl[NB][NQ];

    // constructor
    SpecificFUNC ();
    ~SpecificFUNC (); // destructor
    void Setup_Function_Space(const mxArray*);
    void Get_Local_to_Global_DoFmap(const int&, int*) const;
                   // need the "const" to ENSURE that nothing in this object will change!
    void Transform_Basis_Functions();
    const CLASS_geom_Omega_embedded_in_Omega_restricted_to_Omega*  Mesh;

private:
    void Determine_Element_Order(const int Mesh_Vertex[SUB_TD+1], bool& Std_Elem_Order) const;
    void Map_Basis_p1();
};

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* constructor */
SpecificFUNC::SpecificFUNC () :
ABSTRACT_FEM_Function_Class () // call the base class constructor
{
    Name       = (char*) SpecificFUNC_str;
    Type       = (char*) SPACE_type;
    Space_Name = (char*) SPACE_name;
    Sub_TopDim = SUB_TD;
    DoI_TopDim = DOI_TD;
    GeoDim     = GD;
    Num_Basis  = NB;
    Num_Comp   = NC;
    Num_QP     = NQ;
    Mesh       = NULL;

    // init DoF information to NULL
    for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
        Elem_DoF[basis_i] = NULL;

    // init everything to zero
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
    for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
        Func_vv_Value[basis_i][qp_i].Set_To_Zero();
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
    for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
        Func_vv_Curl[basis_i][qp_i].Set_To_Zero();
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/


/***************************************************************************************/
/* DE-structor */
SpecificFUNC::~SpecificFUNC ()
{
}
/***************************************************************************************/


/***************************************************************************************/
/* put incoming function data from MATLAB into a nice struct  */
void SpecificFUNC::Setup_Function_Space(const mxArray* Elem)          // inputs
{
    Init_Function_Space(Elem);

    // split up the columns of the element data
    Elem_DoF[0] = (int *) mxGetPr(Elem);
    for (int basis_i = 1; (basis_i < Num_Basis); basis_i++)
        Elem_DoF[basis_i] = Elem_DoF[basis_i-1] + Num_Elem;
}
/***************************************************************************************/


/***************************************************************************************/
/* get the local DoFs on the given element.
   Note: elem_index is in the   C-style (i.e. 0 <= elem_index <= Num_Elem - 1),
         Indices is in the MATLAB-style (i.e. 1 <= Indices[:] <= max(Elem_DoF)). */
void SpecificFUNC::Get_Local_to_Global_DoFmap(const int& elem_index, int* Indices) const  // inputs
{
    /* error check: */
    if (elem_index < 0)
        {
        mexPrintf("ERROR: Given element index #%d is not positive. It must be > 0!\n",elem_index+1);
        mexPrintf("ERROR: There is an issue with the Finite Element Space = %s!\n",Space_Name);
        mexErrMsgTxt("ERROR: Make sure your inputs are valid!");
        }
    else if (elem_index >= Num_Elem)
        {
        mexPrintf("ERROR: Given element index #%d exceeds the number of elements in the finite element (FE) space.\n",elem_index+1);
		mexPrintf("It must be <= %d!  OR  Your FE space DoFmap is not defined correctly!\n",Num_Elem);
		mexPrintf("   For example, the number of rows in DoFmap should *equal*\n");
		mexPrintf("       the number of mesh elements in the (sub)-domain.\n");
        mexPrintf("ERROR: There is an issue with the Finite Element Space = %s!\n",Space_Name);
        mexErrMsgTxt("ERROR: Make sure your inputs are valid!");
        }

    // get local to global index map for the current element
    for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
        {
        int DoF_index = Elem_DoF[basis_i][elem_index] - 1; // shifted for C - style indexing
        Indices[basis_i] = DoF_index;
        }
}
/***************************************************************************************/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* determine ordering of the current mesh element, so we can choose the correct
   reference element (and associated basis functions).
   Let [V_1, V_2, V_3, V_4] be the global vertex indices of the current mesh element.
        Only two orderings are allowed:
        V_1 < V_2 < V_3 < V_4, (ascending order);
        V_1 < V_3 < V_2 < V_4, (mirror reflection ascending order).
        Any other order outputs an error message.

        This is done because we use the "Paul Wesson Trick" to choose the orientation
        of tetrahedral edges, as well as the tangent basis of each tetrahedral face.
        For more information, see pages 142-143 in
                "Finite Element Methods for Maxwell's Equations", by Peter Monk.  */
void SpecificFUNC::Determine_Element_Order(const int Mesh_Vertex[SUB_TD+1], bool& Std_Elem_Order) const
{
    // determine ascending order of local mesh vertices
    int Sorted_Indices[SUB_TD+1];

    Sort_Four_Ints(Mesh_Vertex[0], Mesh_Vertex[1], Mesh_Vertex[2], Mesh_Vertex[3], Sorted_Indices);

    // determine which vertex order case we have
    const bool ascending = (Sorted_Indices[0]==0) & 
                           (Sorted_Indices[1]==1) & 
                           (Sorted_Indices[2]==2) & 
                           (Sorted_Indices[3]==3);
    const bool mirror    = (Sorted_Indices[0]==0) & 
                           (Sorted_Indices[1]==2) & 
                           (Sorted_Indices[2]==1) & 
                           (Sorted_Indices[3]==3);

    if (ascending)
        Std_Elem_Order = true;
    else if (mirror)
        Std_Elem_Order = false;
    else
        {
        mexPrintf("FELICITY ERROR\n");
        mexPrintf("-----------------------------------------------------------------------------\n");
        mexPrintf("Given mesh element [V_1, V_2, V_3, V_4] does not satisfy either ordering:\n");
        mexPrintf("        V_1 < V_2 < V_3 < V_4, (ascending order);\n");
        mexPrintf("        V_1 < V_3 < V_2 < V_4, (mirror reflection ascending order).\n");
        mexPrintf("The actual order is:\n");
        mexPrintf("        V_%d < V_%d < V_%d < V_%d.\n",Sorted_Indices[0]+1,Sorted_Indices[1]+1,Sorted_Indices[2]+1,Sorted_Indices[3]+1);
        mexPrintf("\n");
        mexPrintf("You must sort your triangulation!\n");
        mexPrintf("    I.e. use the MeshTetrahedron method 'Order_Cell_Vertices_For_Hcurl'.\n");
        mexPrintf("-----------------------------------------------------------------------------\n");
        mexErrMsgTxt("Triangulation is not sorted correctly!\n");
        }
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* compute the correct local transformation */
void SpecificFUNC::Transform_Basis_Functions()
{
    Map_Basis_p1();
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/


/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* map basis functions from the standard reference element
       to an actual element in the Domain.     */
void SpecificFUNC::Map_Basis_p1()
{

/*------------ BEGIN: Auto Generate ------------*/
    /* call sub-routine to determine current mesh element ordering */
    /*      (This is related to the Paul Wesson trick!)  */
    bool Std_Elem_Order;
    // retrieve global vertex indices on the current cell
    int Vtx_Ind[SUB_TD+1];
    Mesh->Get_Current_Cell_Vertex_Indices(Vtx_Ind);
    Determine_Element_Order(Vtx_Ind, Std_Elem_Order);
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
// Local Element defined on Subdomain: CG, nedelec_1stkind_deg2_dim3
// the Subdomain             has topological dimension = 3
// the Domain of Integration has topological dimension = 3
// geometric dimension = 3
// Number of Quadrature Points = 10

    // get "Val" of basis functions
    VEC_3x1 phi_Val[NQ][NB];
    // get "Curl" of basis functions
    VEC_3x1 phi_Curl[NQ][NB];

    // choose basis evaluations based on element order
    if (Std_Elem_Order) // V_1 < V_2 < V_3 < V_4
        {
        phi_Val[0][0].Set_Equal_To(-8.92550232934685051E-01, -4.82681821336412742E-01, -4.82681821336412742E-01);
        phi_Val[0][1].Set_Equal_To(6.77650698804054930E-01, 3.11184295615549111E-01, 3.11184295615549111E-01);
        phi_Val[0][2].Set_Equal_To(-1.22155467729501976E-01, -5.32023879327774285E-01, -1.22155467729501976E-01);
        phi_Val[0][3].Set_Equal_To(-4.09868411598272309E-01, -5.32023879327774285E-01, -4.09868411598272309E-01);
        phi_Val[0][4].Set_Equal_To(-1.22155467729501976E-01, -1.22155467729501976E-01, -5.32023879327774285E-01);
        phi_Val[0][5].Set_Equal_To(-4.09868411598272309E-01, -4.09868411598272309E-01, -5.32023879327774285E-01);
        phi_Val[0][6].Set_Equal_To(-3.66466403188505874E-01, 3.11184295615549111E-01, 0.00000000000000000E+00);
        phi_Val[0][7].Set_Equal_To(4.09868411598272309E-01, -4.82681821336412631E-01, 0.00000000000000000E+00);
        phi_Val[0][8].Set_Equal_To(0.00000000000000000E+00, 1.22155467729501949E-01, -4.09868411598272309E-01);
        phi_Val[0][9].Set_Equal_To(0.00000000000000000E+00, 4.09868411598272309E-01, -1.22155467729501949E-01);
        phi_Val[0][10].Set_Equal_To(4.09868411598272309E-01, 0.00000000000000000E+00, -4.82681821336412631E-01);
        phi_Val[0][11].Set_Equal_To(-3.66466403188505874E-01, 0.00000000000000000E+00, 3.11184295615549111E-01);
        phi_Val[0][12].Set_Equal_To(-8.27787380696342195E-02, 6.54179347057276317E-01, -3.27089673528638158E-01);
        phi_Val[0][13].Set_Equal_To(8.27787380696342195E-02, 2.48336214208902617E-01, 1.59712767382897074E-17);
        phi_Val[0][14].Set_Equal_To(4.92647149667906514E-01, 3.27089673528638158E-01, 6.31084722476465350E-17);
        phi_Val[0][15].Set_Equal_To(4.92647149667906514E-01, 6.31084722476465350E-17, 3.27089673528638158E-01);
        phi_Val[0][16].Set_Equal_To(-8.27787380696342195E-02, -3.27089673528638158E-01, 6.54179347057276317E-01);
        phi_Val[0][17].Set_Equal_To(8.27787380696342195E-02, 1.59712767382897074E-17, 2.48336214208902617E-01);
        phi_Val[0][18].Set_Equal_To(2.44310935459003925E-01, 3.27089673528638158E-01, 9.81269020585914253E-01);
        phi_Val[0][19].Set_Equal_To(2.44310935459003925E-01, 9.81269020585914253E-01, 3.27089673528638158E-01);
        phi_Val[1][0].Set_Equal_To(6.77650698804054819E-01, 3.66466403188505818E-01, 3.66466403188505818E-01);
        phi_Val[1][1].Set_Equal_To(-8.92550232934684940E-01, -4.09868411598272309E-01, -4.09868411598272309E-01);
        phi_Val[1][2].Set_Equal_To(3.66466403188505818E-01, 6.77650698804054819E-01, 3.66466403188505818E-01);
        phi_Val[1][3].Set_Equal_To(-4.09868411598272309E-01, -8.92550232934684940E-01, -4.09868411598272309E-01);
        phi_Val[1][4].Set_Equal_To(3.66466403188505818E-01, 3.66466403188505818E-01, 6.77650698804054819E-01);
        phi_Val[1][5].Set_Equal_To(-4.09868411598272309E-01, -4.09868411598272309E-01, -8.92550232934684940E-01);
        phi_Val[1][6].Set_Equal_To(1.22155467729501949E-01, -4.09868411598272309E-01, 0.00000000000000000E+00);
        phi_Val[1][7].Set_Equal_To(4.09868411598272309E-01, -1.22155467729501949E-01, 0.00000000000000000E+00);
        phi_Val[1][8].Set_Equal_To(0.00000000000000000E+00, 1.22155467729501949E-01, -4.09868411598272309E-01);
        phi_Val[1][9].Set_Equal_To(0.00000000000000000E+00, 4.09868411598272309E-01, -1.22155467729501949E-01);
        phi_Val[1][10].Set_Equal_To(4.09868411598272309E-01, 0.00000000000000000E+00, -1.22155467729501949E-01);
        phi_Val[1][11].Set_Equal_To(1.22155467729501949E-01, 0.00000000000000000E+00, -4.09868411598272309E-01);
        phi_Val[1][12].Set_Equal_To(-8.27787380696342195E-02, 1.65557476139268439E-01, -8.27787380696342195E-02);
        phi_Val[1][13].Set_Equal_To(8.27787380696342195E-02, 7.36958085126910412E-01, -2.44310935459003897E-01);
        phi_Val[1][14].Set_Equal_To(7.36958085126910412E-01, 8.27787380696342195E-02, -2.44310935459003897E-01);
        phi_Val[1][15].Set_Equal_To(7.36958085126910412E-01, -2.44310935459003897E-01, 8.27787380696342195E-02);
        phi_Val[1][16].Set_Equal_To(-8.27787380696342195E-02, -8.27787380696342195E-02, 1.65557476139268439E-01);
        phi_Val[1][17].Set_Equal_To(8.27787380696342195E-02, -2.44310935459003897E-01, 7.36958085126910412E-01);
        phi_Val[1][18].Set_Equal_To(-2.44310935459003897E-01, 8.27787380696342195E-02, 7.36958085126910412E-01);
        phi_Val[1][19].Set_Equal_To(-2.44310935459003897E-01, 7.36958085126910412E-01, 8.27787380696342195E-02);
        phi_Val[2][0].Set_Equal_To(-5.32023879327774285E-01, -1.22155467729501976E-01, -1.22155467729501976E-01);
        phi_Val[2][1].Set_Equal_To(-5.32023879327774285E-01, -4.09868411598272309E-01, -4.09868411598272309E-01);
        phi_Val[2][2].Set_Equal_To(-1.22155467729501976E-01, -5.32023879327774285E-01, -1.22155467729501976E-01);
        phi_Val[2][3].Set_Equal_To(-4.09868411598272309E-01, -5.32023879327774285E-01, -4.09868411598272309E-01);
        phi_Val[2][4].Set_Equal_To(-4.82681821336412742E-01, -4.82681821336412742E-01, -8.92550232934685051E-01);
        phi_Val[2][5].Set_Equal_To(3.11184295615549111E-01, 3.11184295615549111E-01, 6.77650698804054930E-01);
        phi_Val[2][6].Set_Equal_To(1.22155467729501949E-01, -4.09868411598272309E-01, 0.00000000000000000E+00);
        phi_Val[2][7].Set_Equal_To(4.09868411598272309E-01, -1.22155467729501949E-01, 0.00000000000000000E+00);
        phi_Val[2][8].Set_Equal_To(0.00000000000000000E+00, 4.82681821336412631E-01, -4.09868411598272309E-01);
        phi_Val[2][9].Set_Equal_To(0.00000000000000000E+00, -3.11184295615549111E-01, 3.66466403188505874E-01);
        phi_Val[2][10].Set_Equal_To(-3.11184295615549111E-01, 0.00000000000000000E+00, 3.66466403188505874E-01);
        phi_Val[2][11].Set_Equal_To(4.82681821336412631E-01, 0.00000000000000000E+00, -4.09868411598272309E-01);
        phi_Val[2][12].Set_Equal_To(-3.27089673528638158E-01, 6.54179347057276317E-01, -8.27787380696342195E-02);
        phi_Val[2][13].Set_Equal_To(3.27089673528638158E-01, 9.81269020585914253E-01, 2.44310935459003925E-01);
        phi_Val[2][14].Set_Equal_To(9.81269020585914253E-01, 3.27089673528638158E-01, 2.44310935459003925E-01);
        phi_Val[2][15].Set_Equal_To(2.48336214208902617E-01, 1.59712767382897074E-17, 8.27787380696342195E-02);
        phi_Val[2][16].Set_Equal_To(-3.27089673528638158E-01, -3.27089673528638158E-01, 1.65557476139268439E-01);
        phi_Val[2][17].Set_Equal_To(3.27089673528638158E-01, 6.31084722476465350E-17, 4.92647149667906514E-01);
        phi_Val[2][18].Set_Equal_To(6.31084722476465350E-17, 3.27089673528638158E-01, 4.92647149667906514E-01);
        phi_Val[2][19].Set_Equal_To(1.59712767382897074E-17, 2.48336214208902617E-01, 8.27787380696342195E-02);
        phi_Val[3][0].Set_Equal_To(-5.32023879327774285E-01, -1.22155467729501976E-01, -1.22155467729501976E-01);
        phi_Val[3][1].Set_Equal_To(-5.32023879327774285E-01, -4.09868411598272309E-01, -4.09868411598272309E-01);
        phi_Val[3][2].Set_Equal_To(-4.82681821336412742E-01, -8.92550232934685051E-01, -4.82681821336412742E-01);
        phi_Val[3][3].Set_Equal_To(3.11184295615549111E-01, 6.77650698804054930E-01, 3.11184295615549111E-01);
        phi_Val[3][4].Set_Equal_To(-1.22155467729501976E-01, -1.22155467729501976E-01, -5.32023879327774285E-01);
        phi_Val[3][5].Set_Equal_To(-4.09868411598272309E-01, -4.09868411598272309E-01, -5.32023879327774285E-01);
        phi_Val[3][6].Set_Equal_To(4.82681821336412631E-01, -4.09868411598272309E-01, 0.00000000000000000E+00);
        phi_Val[3][7].Set_Equal_To(-3.11184295615549111E-01, 3.66466403188505874E-01, 0.00000000000000000E+00);
        phi_Val[3][8].Set_Equal_To(0.00000000000000000E+00, -3.66466403188505874E-01, 3.11184295615549111E-01);
        phi_Val[3][9].Set_Equal_To(0.00000000000000000E+00, 4.09868411598272309E-01, -4.82681821336412631E-01);
        phi_Val[3][10].Set_Equal_To(4.09868411598272309E-01, 0.00000000000000000E+00, -1.22155467729501949E-01);
        phi_Val[3][11].Set_Equal_To(1.22155467729501949E-01, 0.00000000000000000E+00, -4.09868411598272309E-01);
        phi_Val[3][12].Set_Equal_To(-3.27089673528638158E-01, 1.65557476139268439E-01, -3.27089673528638158E-01);
        phi_Val[3][13].Set_Equal_To(3.27089673528638158E-01, 4.92647149667906514E-01, 6.31084722476465350E-17);
        phi_Val[3][14].Set_Equal_To(2.48336214208902617E-01, 8.27787380696342195E-02, 1.59712767382897074E-17);
        phi_Val[3][15].Set_Equal_To(9.81269020585914253E-01, 2.44310935459003925E-01, 3.27089673528638158E-01);
        phi_Val[3][16].Set_Equal_To(-3.27089673528638158E-01, -8.27787380696342195E-02, 6.54179347057276317E-01);
        phi_Val[3][17].Set_Equal_To(3.27089673528638158E-01, 2.44310935459003925E-01, 9.81269020585914253E-01);
        phi_Val[3][18].Set_Equal_To(1.59712767382897074E-17, 8.27787380696342195E-02, 2.48336214208902617E-01);
        phi_Val[3][19].Set_Equal_To(6.31084722476465350E-17, 4.92647149667906514E-01, 3.27089673528638158E-01);
        phi_Val[4][0].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][1].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][2].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[4][3].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][4].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[4][5].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][6].Set_Equal_To(1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][8].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][9].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[4][10].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][11].Set_Equal_To(1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][12].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][13].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[4][14].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][15].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][16].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][17].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[4][18].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][19].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][0].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[5][1].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][2].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][3].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][4].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[5][5].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][7].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][8].Set_Equal_To(0.00000000000000000E+00, 1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][9].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][10].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[5][11].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][12].Set_Equal_To(0.00000000000000000E+00, 2.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][13].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][14].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[5][15].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][16].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][17].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][18].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[5][19].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][0].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[6][1].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][2].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[6][3].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][4].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][5].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][6].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][7].Set_Equal_To(0.00000000000000000E+00, 1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][8].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][9].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[6][10].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[6][11].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[6][13].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][14].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][15].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[6][16].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 2.00000000000000000E+00);
        phi_Val[6][17].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][18].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][19].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[7][0].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[7][1].Set_Equal_To(1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][2].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][3].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][4].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][5].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[7][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][7].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][8].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][9].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][10].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[7][11].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][13].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][14].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[7][15].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][16].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][17].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][18].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 2.00000000000000000E+00);
        phi_Val[7][19].Set_Equal_To(0.00000000000000000E+00, 2.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][0].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][1].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][2].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[8][3].Set_Equal_To(0.00000000000000000E+00, 1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][4].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][5].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[8][6].Set_Equal_To(1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][8].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][9].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[8][10].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][11].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][13].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[8][14].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][15].Set_Equal_To(2.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][16].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][17].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 2.00000000000000000E+00);
        phi_Val[8][18].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][19].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][0].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][1].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][2].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][3].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][4].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[9][5].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[9][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][8].Set_Equal_To(0.00000000000000000E+00, 1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][9].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][10].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][11].Set_Equal_To(1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][13].Set_Equal_To(0.00000000000000000E+00, 2.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][14].Set_Equal_To(2.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][15].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][16].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][17].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][18].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][19].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);

        phi_Curl[0][0].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475601E+00, -2.54744467357475601E+00);
        phi_Curl[0][1].Set_Equal_To(0.00000000000000000E+00, -7.64233402072426671E+00, 7.64233402072426671E+00);
        phi_Curl[0][2].Set_Equal_To(-2.54744467357475601E+00, 0.00000000000000000E+00, 2.54744467357475601E+00);
        phi_Curl[0][3].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 2.54744467357475557E+00);
        phi_Curl[0][4].Set_Equal_To(2.54744467357475601E+00, -2.54744467357475601E+00, 0.00000000000000000E+00);
        phi_Curl[0][5].Set_Equal_To(2.54744467357475557E+00, -2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[0][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 7.64233402072426671E+00);
        phi_Curl[0][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -2.54744467357475557E+00);
        phi_Curl[0][8].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[0][9].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[0][10].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[0][11].Set_Equal_To(0.00000000000000000E+00, -7.64233402072426671E+00, 0.00000000000000000E+00);
        phi_Curl[0][12].Set_Equal_To(-6.82116701036213335E+00, 0.00000000000000000E+00, 1.72627766321262222E+00);
        phi_Curl[0][13].Set_Equal_To(3.33066907387546962E-16, 0.00000000000000000E+00, -1.72627766321262222E+00);
        phi_Curl[0][14].Set_Equal_To(0.00000000000000000E+00, -3.33066907387546962E-16, 1.72627766321262222E+00);
        phi_Curl[0][15].Set_Equal_To(0.00000000000000000E+00, -1.72627766321262222E+00, 3.33066907387546962E-16);
        phi_Curl[0][16].Set_Equal_To(6.82116701036213335E+00, -1.72627766321262222E+00, 0.00000000000000000E+00);
        phi_Curl[0][17].Set_Equal_To(-3.33066907387546962E-16, 1.72627766321262222E+00, 0.00000000000000000E+00);
        phi_Curl[0][18].Set_Equal_To(-6.82116701036213335E+00, 5.09488934714951114E+00, 0.00000000000000000E+00);
        phi_Curl[0][19].Set_Equal_To(6.82116701036213335E+00, 0.00000000000000000E+00, -5.09488934714951114E+00);
        phi_Curl[1][0].Set_Equal_To(0.00000000000000000E+00, -7.64233402072426582E+00, 7.64233402072426582E+00);
        phi_Curl[1][1].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, -2.54744467357475557E+00);
        phi_Curl[1][2].Set_Equal_To(7.64233402072426582E+00, 0.00000000000000000E+00, -7.64233402072426582E+00);
        phi_Curl[1][3].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 2.54744467357475557E+00);
        phi_Curl[1][4].Set_Equal_To(-7.64233402072426582E+00, 7.64233402072426582E+00, 0.00000000000000000E+00);
        phi_Curl[1][5].Set_Equal_To(2.54744467357475557E+00, -2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[1][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -2.54744467357475557E+00);
        phi_Curl[1][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -2.54744467357475557E+00);
        phi_Curl[1][8].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[1][9].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[1][10].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[1][11].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[1][12].Set_Equal_To(-1.72627766321262222E+00, 0.00000000000000000E+00, 1.72627766321262222E+00);
        phi_Curl[1][13].Set_Equal_To(-5.09488934714951114E+00, 0.00000000000000000E+00, -1.72627766321262222E+00);
        phi_Curl[1][14].Set_Equal_To(0.00000000000000000E+00, 5.09488934714951114E+00, 1.72627766321262222E+00);
        phi_Curl[1][15].Set_Equal_To(0.00000000000000000E+00, -1.72627766321262222E+00, -5.09488934714951114E+00);
        phi_Curl[1][16].Set_Equal_To(1.72627766321262222E+00, -1.72627766321262222E+00, 0.00000000000000000E+00);
        phi_Curl[1][17].Set_Equal_To(5.09488934714951114E+00, 1.72627766321262222E+00, 0.00000000000000000E+00);
        phi_Curl[1][18].Set_Equal_To(-1.72627766321262222E+00, -5.09488934714951114E+00, 0.00000000000000000E+00);
        phi_Curl[1][19].Set_Equal_To(1.72627766321262222E+00, 0.00000000000000000E+00, 5.09488934714951114E+00);
        phi_Curl[2][0].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475601E+00, -2.54744467357475601E+00);
        phi_Curl[2][1].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, -2.54744467357475557E+00);
        phi_Curl[2][2].Set_Equal_To(-2.54744467357475601E+00, 0.00000000000000000E+00, 2.54744467357475601E+00);
        phi_Curl[2][3].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 2.54744467357475557E+00);
        phi_Curl[2][4].Set_Equal_To(2.54744467357475601E+00, -2.54744467357475601E+00, 0.00000000000000000E+00);
        phi_Curl[2][5].Set_Equal_To(-7.64233402072426671E+00, 7.64233402072426671E+00, 0.00000000000000000E+00);
        phi_Curl[2][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -2.54744467357475557E+00);
        phi_Curl[2][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -2.54744467357475557E+00);
        phi_Curl[2][8].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[2][9].Set_Equal_To(7.64233402072426671E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[2][10].Set_Equal_To(0.00000000000000000E+00, -7.64233402072426671E+00, 0.00000000000000000E+00);
        phi_Curl[2][11].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[2][12].Set_Equal_To(-1.72627766321262222E+00, 0.00000000000000000E+00, 6.82116701036213335E+00);
        phi_Curl[2][13].Set_Equal_To(5.09488934714951114E+00, 0.00000000000000000E+00, -6.82116701036213335E+00);
        phi_Curl[2][14].Set_Equal_To(0.00000000000000000E+00, -5.09488934714951114E+00, 6.82116701036213335E+00);
        phi_Curl[2][15].Set_Equal_To(0.00000000000000000E+00, -1.72627766321262222E+00, 3.33066907387546962E-16);
        phi_Curl[2][16].Set_Equal_To(1.72627766321262222E+00, -1.72627766321262222E+00, 0.00000000000000000E+00);
        phi_Curl[2][17].Set_Equal_To(-3.33066907387546962E-16, 1.72627766321262222E+00, 0.00000000000000000E+00);
        phi_Curl[2][18].Set_Equal_To(-1.72627766321262222E+00, 3.33066907387546962E-16, 0.00000000000000000E+00);
        phi_Curl[2][19].Set_Equal_To(1.72627766321262222E+00, 0.00000000000000000E+00, -3.33066907387546962E-16);
        phi_Curl[3][0].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475601E+00, -2.54744467357475601E+00);
        phi_Curl[3][1].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, -2.54744467357475557E+00);
        phi_Curl[3][2].Set_Equal_To(-2.54744467357475601E+00, 0.00000000000000000E+00, 2.54744467357475601E+00);
        phi_Curl[3][3].Set_Equal_To(7.64233402072426671E+00, 0.00000000000000000E+00, -7.64233402072426671E+00);
        phi_Curl[3][4].Set_Equal_To(2.54744467357475601E+00, -2.54744467357475601E+00, 0.00000000000000000E+00);
        phi_Curl[3][5].Set_Equal_To(2.54744467357475557E+00, -2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[3][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -2.54744467357475557E+00);
        phi_Curl[3][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 7.64233402072426671E+00);
        phi_Curl[3][8].Set_Equal_To(7.64233402072426671E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[3][9].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[3][10].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[3][11].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[3][12].Set_Equal_To(-1.72627766321262222E+00, 0.00000000000000000E+00, 1.72627766321262222E+00);
        phi_Curl[3][13].Set_Equal_To(3.33066907387546962E-16, 0.00000000000000000E+00, -1.72627766321262222E+00);
        phi_Curl[3][14].Set_Equal_To(0.00000000000000000E+00, -3.33066907387546962E-16, 1.72627766321262222E+00);
        phi_Curl[3][15].Set_Equal_To(0.00000000000000000E+00, -6.82116701036213335E+00, 5.09488934714951114E+00);
        phi_Curl[3][16].Set_Equal_To(1.72627766321262222E+00, -6.82116701036213335E+00, 0.00000000000000000E+00);
        phi_Curl[3][17].Set_Equal_To(-5.09488934714951114E+00, 6.82116701036213335E+00, 0.00000000000000000E+00);
        phi_Curl[3][18].Set_Equal_To(-1.72627766321262222E+00, 3.33066907387546962E-16, 0.00000000000000000E+00);
        phi_Curl[3][19].Set_Equal_To(1.72627766321262222E+00, 0.00000000000000000E+00, -3.33066907387546962E-16);
        phi_Curl[4][0].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[4][1].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[4][2].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[4][3].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[4][4].Set_Equal_To(6.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][5].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[4][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[4][8].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][9].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][10].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][11].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[4][13].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[4][14].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[4][15].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[4][16].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][17].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][18].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][19].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][0].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[5][1].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[5][2].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[5][3].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[5][4].Set_Equal_To(6.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][5].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[5][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[5][8].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][9].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][10].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][11].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][12].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[5][13].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[5][14].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[5][15].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][16].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][17].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][18].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][19].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[6][0].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[6][1].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[6][2].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[6][3].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[6][4].Set_Equal_To(6.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][5].Set_Equal_To(6.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[6][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[6][8].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][9].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][10].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][11].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][12].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][13].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][14].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][15].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[6][16].Set_Equal_To(6.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][17].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][18].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][19].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[7][0].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[7][1].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[7][2].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[7][3].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[7][4].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][5].Set_Equal_To(6.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[7][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[7][8].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][9].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][10].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][11].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][12].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][13].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][14].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][15].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[7][16].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][17].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][18].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][19].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][0].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[8][1].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[8][2].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[8][3].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[8][4].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][5].Set_Equal_To(6.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[8][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[8][8].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][9].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][10].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][11].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][13].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][14].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][15].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][16].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][17].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][18].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][19].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[9][0].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[9][1].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[9][2].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[9][3].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[9][4].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][5].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[9][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[9][8].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][9].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][10].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][11].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[9][13].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[9][14].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[9][15].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[9][16].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][17].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][18].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][19].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);


        }
    else // V_1 < V_3 < V_2 < V_4
        {
        phi_Val[0][0].Set_Equal_To(-8.92550232934685051E-01, -4.82681821336412742E-01, -4.82681821336412742E-01);
        phi_Val[0][1].Set_Equal_To(6.77650698804054930E-01, 3.11184295615549111E-01, 3.11184295615549111E-01);
        phi_Val[0][2].Set_Equal_To(-1.22155467729501976E-01, -5.32023879327774285E-01, -1.22155467729501976E-01);
        phi_Val[0][3].Set_Equal_To(-4.09868411598272309E-01, -5.32023879327774285E-01, -4.09868411598272309E-01);
        phi_Val[0][4].Set_Equal_To(-1.22155467729501976E-01, -1.22155467729501976E-01, -5.32023879327774285E-01);
        phi_Val[0][5].Set_Equal_To(-4.09868411598272309E-01, -4.09868411598272309E-01, -5.32023879327774285E-01);
        phi_Val[0][6].Set_Equal_To(3.66466403188505874E-01, -3.11184295615549111E-01, 0.00000000000000000E+00);
        phi_Val[0][7].Set_Equal_To(-4.09868411598272309E-01, 4.82681821336412631E-01, 0.00000000000000000E+00);
        phi_Val[0][8].Set_Equal_To(0.00000000000000000E+00, 1.22155467729501949E-01, -4.09868411598272309E-01);
        phi_Val[0][9].Set_Equal_To(0.00000000000000000E+00, 4.09868411598272309E-01, -1.22155467729501949E-01);
        phi_Val[0][10].Set_Equal_To(4.09868411598272309E-01, 0.00000000000000000E+00, -4.82681821336412631E-01);
        phi_Val[0][11].Set_Equal_To(-3.66466403188505874E-01, 0.00000000000000000E+00, 3.11184295615549111E-01);
        phi_Val[0][12].Set_Equal_To(1.65557476139268439E-01, -3.27089673528638158E-01, -3.27089673528638158E-01);
        phi_Val[0][13].Set_Equal_To(8.27787380696342195E-02, 2.48336214208902617E-01, 1.59712767382897074E-17);
        phi_Val[0][14].Set_Equal_To(4.92647149667906514E-01, 3.27089673528638158E-01, 6.31084722476465350E-17);
        phi_Val[0][15].Set_Equal_To(2.44310935459003925E-01, 9.81269020585914253E-01, 3.27089673528638158E-01);
        phi_Val[0][16].Set_Equal_To(-8.27787380696342195E-02, -3.27089673528638158E-01, 6.54179347057276317E-01);
        phi_Val[0][17].Set_Equal_To(8.27787380696342195E-02, 1.59712767382897074E-17, 2.48336214208902617E-01);
        phi_Val[0][18].Set_Equal_To(2.44310935459003925E-01, 3.27089673528638158E-01, 9.81269020585914253E-01);
        phi_Val[0][19].Set_Equal_To(4.92647149667906514E-01, 6.31084722476465350E-17, 3.27089673528638158E-01);
        phi_Val[1][0].Set_Equal_To(6.77650698804054819E-01, 3.66466403188505818E-01, 3.66466403188505818E-01);
        phi_Val[1][1].Set_Equal_To(-8.92550232934684940E-01, -4.09868411598272309E-01, -4.09868411598272309E-01);
        phi_Val[1][2].Set_Equal_To(3.66466403188505818E-01, 6.77650698804054819E-01, 3.66466403188505818E-01);
        phi_Val[1][3].Set_Equal_To(-4.09868411598272309E-01, -8.92550232934684940E-01, -4.09868411598272309E-01);
        phi_Val[1][4].Set_Equal_To(3.66466403188505818E-01, 3.66466403188505818E-01, 6.77650698804054819E-01);
        phi_Val[1][5].Set_Equal_To(-4.09868411598272309E-01, -4.09868411598272309E-01, -8.92550232934684940E-01);
        phi_Val[1][6].Set_Equal_To(-1.22155467729501949E-01, 4.09868411598272309E-01, 0.00000000000000000E+00);
        phi_Val[1][7].Set_Equal_To(-4.09868411598272309E-01, 1.22155467729501949E-01, 0.00000000000000000E+00);
        phi_Val[1][8].Set_Equal_To(0.00000000000000000E+00, 1.22155467729501949E-01, -4.09868411598272309E-01);
        phi_Val[1][9].Set_Equal_To(0.00000000000000000E+00, 4.09868411598272309E-01, -1.22155467729501949E-01);
        phi_Val[1][10].Set_Equal_To(4.09868411598272309E-01, 0.00000000000000000E+00, -1.22155467729501949E-01);
        phi_Val[1][11].Set_Equal_To(1.22155467729501949E-01, 0.00000000000000000E+00, -4.09868411598272309E-01);
        phi_Val[1][12].Set_Equal_To(1.65557476139268439E-01, -8.27787380696342195E-02, -8.27787380696342195E-02);
        phi_Val[1][13].Set_Equal_To(8.27787380696342195E-02, 7.36958085126910412E-01, -2.44310935459003897E-01);
        phi_Val[1][14].Set_Equal_To(7.36958085126910412E-01, 8.27787380696342195E-02, -2.44310935459003897E-01);
        phi_Val[1][15].Set_Equal_To(-2.44310935459003897E-01, 7.36958085126910412E-01, 8.27787380696342195E-02);
        phi_Val[1][16].Set_Equal_To(-8.27787380696342195E-02, -8.27787380696342195E-02, 1.65557476139268439E-01);
        phi_Val[1][17].Set_Equal_To(8.27787380696342195E-02, -2.44310935459003897E-01, 7.36958085126910412E-01);
        phi_Val[1][18].Set_Equal_To(-2.44310935459003897E-01, 8.27787380696342195E-02, 7.36958085126910412E-01);
        phi_Val[1][19].Set_Equal_To(7.36958085126910412E-01, -2.44310935459003897E-01, 8.27787380696342195E-02);
        phi_Val[2][0].Set_Equal_To(-5.32023879327774285E-01, -1.22155467729501976E-01, -1.22155467729501976E-01);
        phi_Val[2][1].Set_Equal_To(-5.32023879327774285E-01, -4.09868411598272309E-01, -4.09868411598272309E-01);
        phi_Val[2][2].Set_Equal_To(-1.22155467729501976E-01, -5.32023879327774285E-01, -1.22155467729501976E-01);
        phi_Val[2][3].Set_Equal_To(-4.09868411598272309E-01, -5.32023879327774285E-01, -4.09868411598272309E-01);
        phi_Val[2][4].Set_Equal_To(-4.82681821336412742E-01, -4.82681821336412742E-01, -8.92550232934685051E-01);
        phi_Val[2][5].Set_Equal_To(3.11184295615549111E-01, 3.11184295615549111E-01, 6.77650698804054930E-01);
        phi_Val[2][6].Set_Equal_To(-1.22155467729501949E-01, 4.09868411598272309E-01, 0.00000000000000000E+00);
        phi_Val[2][7].Set_Equal_To(-4.09868411598272309E-01, 1.22155467729501949E-01, 0.00000000000000000E+00);
        phi_Val[2][8].Set_Equal_To(0.00000000000000000E+00, 4.82681821336412631E-01, -4.09868411598272309E-01);
        phi_Val[2][9].Set_Equal_To(0.00000000000000000E+00, -3.11184295615549111E-01, 3.66466403188505874E-01);
        phi_Val[2][10].Set_Equal_To(-3.11184295615549111E-01, 0.00000000000000000E+00, 3.66466403188505874E-01);
        phi_Val[2][11].Set_Equal_To(4.82681821336412631E-01, 0.00000000000000000E+00, -4.09868411598272309E-01);
        phi_Val[2][12].Set_Equal_To(6.54179347057276317E-01, -3.27089673528638158E-01, -8.27787380696342195E-02);
        phi_Val[2][13].Set_Equal_To(3.27089673528638158E-01, 9.81269020585914253E-01, 2.44310935459003925E-01);
        phi_Val[2][14].Set_Equal_To(9.81269020585914253E-01, 3.27089673528638158E-01, 2.44310935459003925E-01);
        phi_Val[2][15].Set_Equal_To(1.59712767382897074E-17, 2.48336214208902617E-01, 8.27787380696342195E-02);
        phi_Val[2][16].Set_Equal_To(-3.27089673528638158E-01, -3.27089673528638158E-01, 1.65557476139268439E-01);
        phi_Val[2][17].Set_Equal_To(3.27089673528638158E-01, 6.31084722476465350E-17, 4.92647149667906514E-01);
        phi_Val[2][18].Set_Equal_To(6.31084722476465350E-17, 3.27089673528638158E-01, 4.92647149667906514E-01);
        phi_Val[2][19].Set_Equal_To(2.48336214208902617E-01, 1.59712767382897074E-17, 8.27787380696342195E-02);
        phi_Val[3][0].Set_Equal_To(-5.32023879327774285E-01, -1.22155467729501976E-01, -1.22155467729501976E-01);
        phi_Val[3][1].Set_Equal_To(-5.32023879327774285E-01, -4.09868411598272309E-01, -4.09868411598272309E-01);
        phi_Val[3][2].Set_Equal_To(-4.82681821336412742E-01, -8.92550232934685051E-01, -4.82681821336412742E-01);
        phi_Val[3][3].Set_Equal_To(3.11184295615549111E-01, 6.77650698804054930E-01, 3.11184295615549111E-01);
        phi_Val[3][4].Set_Equal_To(-1.22155467729501976E-01, -1.22155467729501976E-01, -5.32023879327774285E-01);
        phi_Val[3][5].Set_Equal_To(-4.09868411598272309E-01, -4.09868411598272309E-01, -5.32023879327774285E-01);
        phi_Val[3][6].Set_Equal_To(-4.82681821336412631E-01, 4.09868411598272309E-01, 0.00000000000000000E+00);
        phi_Val[3][7].Set_Equal_To(3.11184295615549111E-01, -3.66466403188505874E-01, 0.00000000000000000E+00);
        phi_Val[3][8].Set_Equal_To(0.00000000000000000E+00, -3.66466403188505874E-01, 3.11184295615549111E-01);
        phi_Val[3][9].Set_Equal_To(0.00000000000000000E+00, 4.09868411598272309E-01, -4.82681821336412631E-01);
        phi_Val[3][10].Set_Equal_To(4.09868411598272309E-01, 0.00000000000000000E+00, -1.22155467729501949E-01);
        phi_Val[3][11].Set_Equal_To(1.22155467729501949E-01, 0.00000000000000000E+00, -4.09868411598272309E-01);
        phi_Val[3][12].Set_Equal_To(6.54179347057276317E-01, -8.27787380696342195E-02, -3.27089673528638158E-01);
        phi_Val[3][13].Set_Equal_To(3.27089673528638158E-01, 4.92647149667906514E-01, 6.31084722476465350E-17);
        phi_Val[3][14].Set_Equal_To(2.48336214208902617E-01, 8.27787380696342195E-02, 1.59712767382897074E-17);
        phi_Val[3][15].Set_Equal_To(6.31084722476465350E-17, 4.92647149667906514E-01, 3.27089673528638158E-01);
        phi_Val[3][16].Set_Equal_To(-3.27089673528638158E-01, -8.27787380696342195E-02, 6.54179347057276317E-01);
        phi_Val[3][17].Set_Equal_To(3.27089673528638158E-01, 2.44310935459003925E-01, 9.81269020585914253E-01);
        phi_Val[3][18].Set_Equal_To(1.59712767382897074E-17, 8.27787380696342195E-02, 2.48336214208902617E-01);
        phi_Val[3][19].Set_Equal_To(9.81269020585914253E-01, 2.44310935459003925E-01, 3.27089673528638158E-01);
        phi_Val[4][0].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][1].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][2].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[4][3].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][4].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[4][5].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][6].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][8].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][9].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[4][10].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][11].Set_Equal_To(1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][12].Set_Equal_To(2.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][13].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[4][14].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][15].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][16].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][17].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[4][18].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[4][19].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][0].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[5][1].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][2].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][3].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][4].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[5][5].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][7].Set_Equal_To(0.00000000000000000E+00, 1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][8].Set_Equal_To(0.00000000000000000E+00, 1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][9].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][10].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[5][11].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][12].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][13].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][14].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[5][15].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][16].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][17].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[5][18].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[5][19].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][0].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[6][1].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][2].Set_Equal_To(-1.00000000000000000E+00, -1.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[6][3].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][4].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][5].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][6].Set_Equal_To(1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][7].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][8].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][9].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[6][10].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[6][11].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[6][13].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][14].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][15].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[6][16].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 2.00000000000000000E+00);
        phi_Val[6][17].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][18].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[6][19].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[7][0].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[7][1].Set_Equal_To(1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][2].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][3].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][4].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][5].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[7][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][7].Set_Equal_To(0.00000000000000000E+00, 1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][8].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][9].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][10].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[7][11].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][13].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][14].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[7][15].Set_Equal_To(0.00000000000000000E+00, 2.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][16].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][17].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[7][18].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 2.00000000000000000E+00);
        phi_Val[7][19].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][0].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][1].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][2].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[8][3].Set_Equal_To(0.00000000000000000E+00, 1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][4].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][5].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[8][6].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][8].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][9].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[8][10].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][11].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][13].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -1.00000000000000000E+00);
        phi_Val[8][14].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][15].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][16].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][17].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 2.00000000000000000E+00);
        phi_Val[8][18].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[8][19].Set_Equal_To(2.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][0].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][1].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][2].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][3].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][4].Set_Equal_To(1.00000000000000000E+00, 1.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[9][5].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 1.00000000000000000E+00);
        phi_Val[9][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][8].Set_Equal_To(0.00000000000000000E+00, 1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][9].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][10].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][11].Set_Equal_To(1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][13].Set_Equal_To(0.00000000000000000E+00, 2.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][14].Set_Equal_To(2.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][15].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][16].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][17].Set_Equal_To(0.00000000000000000E+00, -1.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][18].Set_Equal_To(-1.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Val[9][19].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);

        phi_Curl[0][0].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475601E+00, -2.54744467357475601E+00);
        phi_Curl[0][1].Set_Equal_To(0.00000000000000000E+00, -7.64233402072426671E+00, 7.64233402072426671E+00);
        phi_Curl[0][2].Set_Equal_To(-2.54744467357475601E+00, 0.00000000000000000E+00, 2.54744467357475601E+00);
        phi_Curl[0][3].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 2.54744467357475557E+00);
        phi_Curl[0][4].Set_Equal_To(2.54744467357475601E+00, -2.54744467357475601E+00, 0.00000000000000000E+00);
        phi_Curl[0][5].Set_Equal_To(2.54744467357475557E+00, -2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[0][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -7.64233402072426671E+00);
        phi_Curl[0][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 2.54744467357475557E+00);
        phi_Curl[0][8].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[0][9].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[0][10].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[0][11].Set_Equal_To(0.00000000000000000E+00, -7.64233402072426671E+00, 0.00000000000000000E+00);
        phi_Curl[0][12].Set_Equal_To(0.00000000000000000E+00, 1.72627766321262222E+00, -1.72627766321262222E+00);
        phi_Curl[0][13].Set_Equal_To(3.33066907387546962E-16, 0.00000000000000000E+00, -1.72627766321262222E+00);
        phi_Curl[0][14].Set_Equal_To(0.00000000000000000E+00, -3.33066907387546962E-16, 1.72627766321262222E+00);
        phi_Curl[0][15].Set_Equal_To(6.82116701036213335E+00, 0.00000000000000000E+00, -5.09488934714951114E+00);
        phi_Curl[0][16].Set_Equal_To(6.82116701036213335E+00, -1.72627766321262222E+00, 0.00000000000000000E+00);
        phi_Curl[0][17].Set_Equal_To(-3.33066907387546962E-16, 1.72627766321262222E+00, 0.00000000000000000E+00);
        phi_Curl[0][18].Set_Equal_To(-6.82116701036213335E+00, 5.09488934714951114E+00, 0.00000000000000000E+00);
        phi_Curl[0][19].Set_Equal_To(0.00000000000000000E+00, -1.72627766321262222E+00, 3.33066907387546962E-16);
        phi_Curl[1][0].Set_Equal_To(0.00000000000000000E+00, -7.64233402072426582E+00, 7.64233402072426582E+00);
        phi_Curl[1][1].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, -2.54744467357475557E+00);
        phi_Curl[1][2].Set_Equal_To(7.64233402072426582E+00, 0.00000000000000000E+00, -7.64233402072426582E+00);
        phi_Curl[1][3].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 2.54744467357475557E+00);
        phi_Curl[1][4].Set_Equal_To(-7.64233402072426582E+00, 7.64233402072426582E+00, 0.00000000000000000E+00);
        phi_Curl[1][5].Set_Equal_To(2.54744467357475557E+00, -2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[1][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 2.54744467357475557E+00);
        phi_Curl[1][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 2.54744467357475557E+00);
        phi_Curl[1][8].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[1][9].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[1][10].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[1][11].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[1][12].Set_Equal_To(0.00000000000000000E+00, 1.72627766321262222E+00, -1.72627766321262222E+00);
        phi_Curl[1][13].Set_Equal_To(-5.09488934714951114E+00, 0.00000000000000000E+00, -1.72627766321262222E+00);
        phi_Curl[1][14].Set_Equal_To(0.00000000000000000E+00, 5.09488934714951114E+00, 1.72627766321262222E+00);
        phi_Curl[1][15].Set_Equal_To(1.72627766321262222E+00, 0.00000000000000000E+00, 5.09488934714951114E+00);
        phi_Curl[1][16].Set_Equal_To(1.72627766321262222E+00, -1.72627766321262222E+00, 0.00000000000000000E+00);
        phi_Curl[1][17].Set_Equal_To(5.09488934714951114E+00, 1.72627766321262222E+00, 0.00000000000000000E+00);
        phi_Curl[1][18].Set_Equal_To(-1.72627766321262222E+00, -5.09488934714951114E+00, 0.00000000000000000E+00);
        phi_Curl[1][19].Set_Equal_To(0.00000000000000000E+00, -1.72627766321262222E+00, -5.09488934714951114E+00);
        phi_Curl[2][0].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475601E+00, -2.54744467357475601E+00);
        phi_Curl[2][1].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, -2.54744467357475557E+00);
        phi_Curl[2][2].Set_Equal_To(-2.54744467357475601E+00, 0.00000000000000000E+00, 2.54744467357475601E+00);
        phi_Curl[2][3].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 2.54744467357475557E+00);
        phi_Curl[2][4].Set_Equal_To(2.54744467357475601E+00, -2.54744467357475601E+00, 0.00000000000000000E+00);
        phi_Curl[2][5].Set_Equal_To(-7.64233402072426671E+00, 7.64233402072426671E+00, 0.00000000000000000E+00);
        phi_Curl[2][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 2.54744467357475557E+00);
        phi_Curl[2][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 2.54744467357475557E+00);
        phi_Curl[2][8].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[2][9].Set_Equal_To(7.64233402072426671E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[2][10].Set_Equal_To(0.00000000000000000E+00, -7.64233402072426671E+00, 0.00000000000000000E+00);
        phi_Curl[2][11].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[2][12].Set_Equal_To(0.00000000000000000E+00, 1.72627766321262222E+00, -6.82116701036213335E+00);
        phi_Curl[2][13].Set_Equal_To(5.09488934714951114E+00, 0.00000000000000000E+00, -6.82116701036213335E+00);
        phi_Curl[2][14].Set_Equal_To(0.00000000000000000E+00, -5.09488934714951114E+00, 6.82116701036213335E+00);
        phi_Curl[2][15].Set_Equal_To(1.72627766321262222E+00, 0.00000000000000000E+00, -3.33066907387546962E-16);
        phi_Curl[2][16].Set_Equal_To(1.72627766321262222E+00, -1.72627766321262222E+00, 0.00000000000000000E+00);
        phi_Curl[2][17].Set_Equal_To(-3.33066907387546962E-16, 1.72627766321262222E+00, 0.00000000000000000E+00);
        phi_Curl[2][18].Set_Equal_To(-1.72627766321262222E+00, 3.33066907387546962E-16, 0.00000000000000000E+00);
        phi_Curl[2][19].Set_Equal_To(0.00000000000000000E+00, -1.72627766321262222E+00, 3.33066907387546962E-16);
        phi_Curl[3][0].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475601E+00, -2.54744467357475601E+00);
        phi_Curl[3][1].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, -2.54744467357475557E+00);
        phi_Curl[3][2].Set_Equal_To(-2.54744467357475601E+00, 0.00000000000000000E+00, 2.54744467357475601E+00);
        phi_Curl[3][3].Set_Equal_To(7.64233402072426671E+00, 0.00000000000000000E+00, -7.64233402072426671E+00);
        phi_Curl[3][4].Set_Equal_To(2.54744467357475601E+00, -2.54744467357475601E+00, 0.00000000000000000E+00);
        phi_Curl[3][5].Set_Equal_To(2.54744467357475557E+00, -2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[3][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 2.54744467357475557E+00);
        phi_Curl[3][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -7.64233402072426671E+00);
        phi_Curl[3][8].Set_Equal_To(7.64233402072426671E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[3][9].Set_Equal_To(-2.54744467357475557E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[3][10].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[3][11].Set_Equal_To(0.00000000000000000E+00, 2.54744467357475557E+00, 0.00000000000000000E+00);
        phi_Curl[3][12].Set_Equal_To(0.00000000000000000E+00, 6.82116701036213335E+00, -1.72627766321262222E+00);
        phi_Curl[3][13].Set_Equal_To(3.33066907387546962E-16, 0.00000000000000000E+00, -1.72627766321262222E+00);
        phi_Curl[3][14].Set_Equal_To(0.00000000000000000E+00, -3.33066907387546962E-16, 1.72627766321262222E+00);
        phi_Curl[3][15].Set_Equal_To(1.72627766321262222E+00, 0.00000000000000000E+00, -3.33066907387546962E-16);
        phi_Curl[3][16].Set_Equal_To(1.72627766321262222E+00, -6.82116701036213335E+00, 0.00000000000000000E+00);
        phi_Curl[3][17].Set_Equal_To(-5.09488934714951114E+00, 6.82116701036213335E+00, 0.00000000000000000E+00);
        phi_Curl[3][18].Set_Equal_To(-1.72627766321262222E+00, 3.33066907387546962E-16, 0.00000000000000000E+00);
        phi_Curl[3][19].Set_Equal_To(0.00000000000000000E+00, -6.82116701036213335E+00, 5.09488934714951114E+00);
        phi_Curl[4][0].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[4][1].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[4][2].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[4][3].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[4][4].Set_Equal_To(6.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][5].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[4][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[4][8].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][9].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][10].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][11].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][12].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[4][13].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[4][14].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[4][15].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][16].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][17].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][18].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[4][19].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[5][0].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[5][1].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[5][2].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[5][3].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[5][4].Set_Equal_To(6.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][5].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[5][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[5][8].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][9].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][10].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][11].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[5][13].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[5][14].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[5][15].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[5][16].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][17].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][18].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[5][19].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][0].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[6][1].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[6][2].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[6][3].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[6][4].Set_Equal_To(6.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][5].Set_Equal_To(6.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[6][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[6][8].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][9].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][10].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][11].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][12].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][13].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][14].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][15].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[6][16].Set_Equal_To(6.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][17].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][18].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[6][19].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[7][0].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[7][1].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[7][2].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[7][3].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[7][4].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][5].Set_Equal_To(6.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[7][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[7][8].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][9].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][10].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][11].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][13].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][14].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][15].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][16].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][17].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][18].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[7][19].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[8][0].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[8][1].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[8][2].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[8][3].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[8][4].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][5].Set_Equal_To(6.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[8][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[8][8].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][9].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][10].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][11].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][12].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][13].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][14].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][15].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[8][16].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][17].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][18].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[8][19].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][0].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[9][1].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[9][2].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[9][3].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[9][4].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][5].Set_Equal_To(-6.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][6].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[9][7].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[9][8].Set_Equal_To(-6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][9].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][10].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][11].Set_Equal_To(0.00000000000000000E+00, 6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][12].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[9][13].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);
        phi_Curl[9][14].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[9][15].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 6.00000000000000000E+00);
        phi_Curl[9][16].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][17].Set_Equal_To(6.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][18].Set_Equal_To(0.00000000000000000E+00, -6.00000000000000000E+00, 0.00000000000000000E+00);
        phi_Curl[9][19].Set_Equal_To(0.00000000000000000E+00, 0.00000000000000000E+00, -6.00000000000000000E+00);


        }
    // END: choose basis evaluations based on element order
// // set of quadrature points
// static const double Quad_Points[NQ][SUB_TD] = { \
//     {5.68430584196844446E-01, 1.43856471934385194E-01, 1.43856471934385194E-01}, \
//     {1.43856471934385194E-01, 1.43856471934385194E-01, 1.43856471934385194E-01}, \
//     {1.43856471934385194E-01, 1.43856471934385194E-01, 5.68430584196844446E-01}, \
//     {1.43856471934385194E-01, 5.68430584196844446E-01, 1.43856471934385194E-01}, \
//     {0.00000000000000000E+00, 5.00000000000000000E-01, 5.00000000000000000E-01}, \
//     {5.00000000000000000E-01, 0.00000000000000000E+00, 5.00000000000000000E-01}, \
//     {5.00000000000000000E-01, 5.00000000000000000E-01, 0.00000000000000000E+00}, \
//     {5.00000000000000000E-01, 0.00000000000000000E+00, 0.00000000000000000E+00}, \
//     {0.00000000000000000E+00, 5.00000000000000000E-01, 0.00000000000000000E+00}, \
//     {0.00000000000000000E+00, 0.00000000000000000E+00, 5.00000000000000000E-01}  \
//     };

// // set of quadrature weights
// static const double Quad_Weights[NQ] = { \
//     3.62941783134008988E-02, \
//     3.62941783134008988E-02, \
//     3.62941783134008988E-02, \
//     3.62941783134008988E-02, \
//     3.58165890217718337E-03, \
//     3.58165890217718337E-03, \
//     3.58165890217718337E-03, \
//     3.58165890217718337E-03, \
//     3.58165890217718337E-03, \
//     3.58165890217718337E-03  \
//     };
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
    /*** compute basis function quantities ***/
    // map basis vectors over (indexing is in the C style)
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        // evaluate for each basis function
        for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
            {
            // multiply by inverse Jacobian matrix transpose
            Mat_Transpose_Vec(Mesh->Map_PHI_Inv_Grad[0], phi_Val[qp_i][basis_i], Func_vv_Value[basis_i][qp_i]);
            }
        }
    // map curl of basis vectors over (indexing is in the C style)
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        // evaluate for each basis function
        for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
            {
            // compute curl of basis vector
            VEC_3x1 curl_temp;
            // pre-multiply by 1/det(Jac)
            Scalar_Mult_Vector(phi_Curl[qp_i][basis_i], Mesh->Map_Inv_Det_Jac[0].a, curl_temp);
            // multiply by Jacobian matrix
            Mat_Vec(Mesh->Map_PHI_Grad[0], curl_temp, Func_vv_Curl[basis_i][qp_i]);
            }
        }
/*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

// remove those macros!
#undef SpecificFUNC
#undef SpecificFUNC_str

#undef SPACE_type
#undef SPACE_name
#undef SUB_TD
#undef DOI_TD
#undef GD
#undef NC
#undef NB
#undef NQ

/***/
