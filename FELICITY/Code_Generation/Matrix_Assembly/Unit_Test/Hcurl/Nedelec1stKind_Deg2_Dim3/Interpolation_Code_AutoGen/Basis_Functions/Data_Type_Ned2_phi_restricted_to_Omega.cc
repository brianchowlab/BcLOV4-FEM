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
#define NQ  1
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

    // get "Val" of basis functions
    VEC_3x1 phi_Val[NB];

    // convenience variables
    const double& x_hat = local_coord[0];
    const double& y_hat = local_coord[1];
    const double& z_hat = local_coord[2];

    // choose basis evaluations based on element order
    if (Std_Elem_Order) // V_1 < V_2 < V_3 < V_4
        {
        phi_Val[0].Set_Equal_To(x_hat*-6.0-y_hat*1.2E+1-z_hat*1.2E+1+x_hat*y_hat*8.0+x_hat*z_hat*8.0+y_hat*z_hat*1.6E+1+(y_hat*y_hat)*8.0+(z_hat*z_hat)*8.0+4.0, x_hat*(x_hat*4.0+y_hat*4.0+z_hat*4.0-3.0)*-2.0, x_hat*(x_hat*4.0+y_hat*4.0+z_hat*4.0-3.0)*-2.0);
        phi_Val[1].Set_Equal_To(x_hat*6.0+y_hat*2.0+z_hat*2.0-x_hat*y_hat*8.0-x_hat*z_hat*8.0-2.0, x_hat*(x_hat*2.0-1.0)*4.0, x_hat*(x_hat*2.0-1.0)*4.0);
        phi_Val[2].Set_Equal_To(y_hat*(x_hat*4.0+y_hat*4.0+z_hat*4.0-3.0)*-2.0, x_hat*-1.2E+1-y_hat*6.0-z_hat*1.2E+1+x_hat*y_hat*8.0+x_hat*z_hat*1.6E+1+y_hat*z_hat*8.0+(x_hat*x_hat)*8.0+(z_hat*z_hat)*8.0+4.0, y_hat*(x_hat*4.0+y_hat*4.0+z_hat*4.0-3.0)*-2.0);
        phi_Val[3].Set_Equal_To(y_hat*(y_hat*2.0-1.0)*4.0, x_hat*2.0+y_hat*6.0+z_hat*2.0-x_hat*y_hat*8.0-y_hat*z_hat*8.0-2.0, y_hat*(y_hat*2.0-1.0)*4.0);
        phi_Val[4].Set_Equal_To(z_hat*(x_hat*4.0+y_hat*4.0+z_hat*4.0-3.0)*-2.0, z_hat*(x_hat*4.0+y_hat*4.0+z_hat*4.0-3.0)*-2.0, x_hat*-1.2E+1-y_hat*1.2E+1-z_hat*6.0+x_hat*y_hat*1.6E+1+x_hat*z_hat*8.0+y_hat*z_hat*8.0+(x_hat*x_hat)*8.0+(y_hat*y_hat)*8.0+4.0);
        phi_Val[5].Set_Equal_To(z_hat*(z_hat*2.0-1.0)*4.0, z_hat*(z_hat*2.0-1.0)*4.0, x_hat*2.0+y_hat*2.0+z_hat*6.0-x_hat*z_hat*8.0-y_hat*z_hat*8.0-2.0);
        phi_Val[6].Set_Equal_To(y_hat*(x_hat*4.0-1.0)*-2.0, x_hat*(x_hat*2.0-1.0)*4.0, 0.0);
        phi_Val[7].Set_Equal_To(y_hat*(y_hat*2.0-1.0)*-4.0, x_hat*(y_hat*4.0-1.0)*2.0, 0.0);
        phi_Val[8].Set_Equal_To(0.0, z_hat*(y_hat*4.0-1.0)*-2.0, y_hat*(y_hat*2.0-1.0)*4.0);
        phi_Val[9].Set_Equal_To(0.0, z_hat*(z_hat*2.0-1.0)*-4.0, y_hat*(z_hat*4.0-1.0)*2.0);
        phi_Val[10].Set_Equal_To(z_hat*(z_hat*2.0-1.0)*-4.0, 0.0, x_hat*(z_hat*4.0-1.0)*2.0);
        phi_Val[11].Set_Equal_To(z_hat*(x_hat*4.0-1.0)*-2.0, 0.0, x_hat*(x_hat*2.0-1.0)*4.0);
        phi_Val[12].Set_Equal_To(y_hat*z_hat*-4.0, x_hat*z_hat*8.0, x_hat*y_hat*-4.0);
        phi_Val[13].Set_Equal_To(y_hat*z_hat*4.0, z_hat*(x_hat*2.0+y_hat+z_hat*2.0-2.0)*-4.0, y_hat*(x_hat+y_hat+z_hat*2.0-1.0)*4.0);
        phi_Val[14].Set_Equal_To(z_hat*(x_hat+y_hat*2.0+z_hat*2.0-2.0)*-4.0, x_hat*z_hat*4.0, x_hat*(x_hat+y_hat+z_hat*2.0-1.0)*4.0);
        phi_Val[15].Set_Equal_To(y_hat*(x_hat+y_hat*2.0+z_hat*2.0-2.0)*-4.0, x_hat*(x_hat+y_hat*2.0+z_hat-1.0)*4.0, x_hat*y_hat*4.0);
        phi_Val[16].Set_Equal_To(y_hat*z_hat*-4.0, x_hat*z_hat*-4.0, x_hat*y_hat*8.0);
        phi_Val[17].Set_Equal_To(y_hat*z_hat*4.0, z_hat*(x_hat+y_hat*2.0+z_hat-1.0)*4.0, y_hat*(x_hat*2.0+y_hat*2.0+z_hat-2.0)*-4.0);
        phi_Val[18].Set_Equal_To(z_hat*(x_hat*2.0+y_hat+z_hat-1.0)*4.0, x_hat*z_hat*4.0, x_hat*(x_hat*2.0+y_hat*2.0+z_hat-2.0)*-4.0);
        phi_Val[19].Set_Equal_To(y_hat*(x_hat*2.0+y_hat+z_hat-1.0)*4.0, x_hat*(x_hat*2.0+y_hat+z_hat*2.0-2.0)*-4.0, x_hat*y_hat*4.0);

        }
    else // V_1 < V_3 < V_2 < V_4
        {
        phi_Val[0].Set_Equal_To(x_hat*-6.0-y_hat*1.2E+1-z_hat*1.2E+1+x_hat*y_hat*8.0+x_hat*z_hat*8.0+y_hat*z_hat*1.6E+1+(y_hat*y_hat)*8.0+(z_hat*z_hat)*8.0+4.0, x_hat*(x_hat*4.0+y_hat*4.0+z_hat*4.0-3.0)*-2.0, x_hat*(x_hat*4.0+y_hat*4.0+z_hat*4.0-3.0)*-2.0);
        phi_Val[1].Set_Equal_To(x_hat*6.0+y_hat*2.0+z_hat*2.0-x_hat*y_hat*8.0-x_hat*z_hat*8.0-2.0, x_hat*(x_hat*2.0-1.0)*4.0, x_hat*(x_hat*2.0-1.0)*4.0);
        phi_Val[2].Set_Equal_To(y_hat*(x_hat*4.0+y_hat*4.0+z_hat*4.0-3.0)*-2.0, x_hat*-1.2E+1-y_hat*6.0-z_hat*1.2E+1+x_hat*y_hat*8.0+x_hat*z_hat*1.6E+1+y_hat*z_hat*8.0+(x_hat*x_hat)*8.0+(z_hat*z_hat)*8.0+4.0, y_hat*(x_hat*4.0+y_hat*4.0+z_hat*4.0-3.0)*-2.0);
        phi_Val[3].Set_Equal_To(y_hat*(y_hat*2.0-1.0)*4.0, x_hat*2.0+y_hat*6.0+z_hat*2.0-x_hat*y_hat*8.0-y_hat*z_hat*8.0-2.0, y_hat*(y_hat*2.0-1.0)*4.0);
        phi_Val[4].Set_Equal_To(z_hat*(x_hat*4.0+y_hat*4.0+z_hat*4.0-3.0)*-2.0, z_hat*(x_hat*4.0+y_hat*4.0+z_hat*4.0-3.0)*-2.0, x_hat*-1.2E+1-y_hat*1.2E+1-z_hat*6.0+x_hat*y_hat*1.6E+1+x_hat*z_hat*8.0+y_hat*z_hat*8.0+(x_hat*x_hat)*8.0+(y_hat*y_hat)*8.0+4.0);
        phi_Val[5].Set_Equal_To(z_hat*(z_hat*2.0-1.0)*4.0, z_hat*(z_hat*2.0-1.0)*4.0, x_hat*2.0+y_hat*2.0+z_hat*6.0-x_hat*z_hat*8.0-y_hat*z_hat*8.0-2.0);
        phi_Val[6].Set_Equal_To(y_hat*(x_hat*4.0-1.0)*2.0, x_hat*(x_hat*2.0-1.0)*-4.0, 0.0);
        phi_Val[7].Set_Equal_To(y_hat*(y_hat*2.0-1.0)*4.0, x_hat*(y_hat*4.0-1.0)*-2.0, 0.0);
        phi_Val[8].Set_Equal_To(0.0, z_hat*(y_hat*4.0-1.0)*-2.0, y_hat*(y_hat*2.0-1.0)*4.0);
        phi_Val[9].Set_Equal_To(0.0, z_hat*(z_hat*2.0-1.0)*-4.0, y_hat*(z_hat*4.0-1.0)*2.0);
        phi_Val[10].Set_Equal_To(z_hat*(z_hat*2.0-1.0)*-4.0, 0.0, x_hat*(z_hat*4.0-1.0)*2.0);
        phi_Val[11].Set_Equal_To(z_hat*(x_hat*4.0-1.0)*-2.0, 0.0, x_hat*(x_hat*2.0-1.0)*4.0);
        phi_Val[12].Set_Equal_To(y_hat*z_hat*8.0, x_hat*z_hat*-4.0, x_hat*y_hat*-4.0);
        phi_Val[13].Set_Equal_To(y_hat*z_hat*4.0, z_hat*(x_hat*2.0+y_hat+z_hat*2.0-2.0)*-4.0, y_hat*(x_hat+y_hat+z_hat*2.0-1.0)*4.0);
        phi_Val[14].Set_Equal_To(z_hat*(x_hat+y_hat*2.0+z_hat*2.0-2.0)*-4.0, x_hat*z_hat*4.0, x_hat*(x_hat+y_hat+z_hat*2.0-1.0)*4.0);
        phi_Val[15].Set_Equal_To(y_hat*(x_hat*2.0+y_hat+z_hat-1.0)*4.0, x_hat*(x_hat*2.0+y_hat+z_hat*2.0-2.0)*-4.0, x_hat*y_hat*4.0);
        phi_Val[16].Set_Equal_To(y_hat*z_hat*-4.0, x_hat*z_hat*-4.0, x_hat*y_hat*8.0);
        phi_Val[17].Set_Equal_To(y_hat*z_hat*4.0, z_hat*(x_hat+y_hat*2.0+z_hat-1.0)*4.0, y_hat*(x_hat*2.0+y_hat*2.0+z_hat-2.0)*-4.0);
        phi_Val[18].Set_Equal_To(z_hat*(x_hat*2.0+y_hat+z_hat-1.0)*4.0, x_hat*z_hat*4.0, x_hat*(x_hat*2.0+y_hat*2.0+z_hat-2.0)*-4.0);
        phi_Val[19].Set_Equal_To(y_hat*(x_hat+y_hat*2.0+z_hat*2.0-2.0)*-4.0, x_hat*(x_hat+y_hat*2.0+z_hat-1.0)*4.0, x_hat*y_hat*4.0);

        }
    // END: choose basis evaluations based on element order
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
            Mat_Transpose_Vec(Mesh->Map_PHI_Inv_Grad[0], phi_Val[basis_i], Func_vv_Value[basis_i][qp_i]);
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
