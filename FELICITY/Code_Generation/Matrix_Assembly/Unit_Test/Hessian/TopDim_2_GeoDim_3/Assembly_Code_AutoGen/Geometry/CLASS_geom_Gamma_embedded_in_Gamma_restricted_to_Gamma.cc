/*
============================================================================================
   This class accesses the mesh geometry data and computes the local transformation from
   the reference element to a `general' element in the mesh.
   
   Several things are computed, such as the gradient of the transformation and Jacobian,
   as well as many other ``custom items''.
   
   This code references the header files:
   
   matrix_vector_defn.h
   matrix_vector_ops.h
   geometric_computations.h
   

   NOTE: portions of this code are automatically generated!

   Copyright (c) 06-06-2012,  Shawn W. Walker
============================================================================================
*/

/*------------ BEGIN: Auto Generate ------------*/
// set mesh geometry data type name
#define MGC        CLASS_geom_Gamma_embedded_in_Gamma_restricted_to_Gamma
// set the type of map
#define MAP_type  "CG - lagrange_deg2_dim2"

// set the Global mesh topological dimension
#define GLOBAL_TD  2
// set the Subdomain topological dimension
#define SUB_TD  2
// set the Domain of Integration (DoI) topological dimension
#define DOI_TD  2
// set the (ambient) geometric dimension
#define GD  3
// set the number of quad points
#define NQ  9
// set the number of basis functions
#define NB  6
// set whether to access the local simplex (facet) orientation
#define ORIENT  false
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
/* C++ (Specific) Mesh class definition */
class MGC: public Abstract_MESH_GEOMETRY_Class // derive from base class
{
public:
    /***************************************************************************************/
    // data structures containing information about the local mapping from a
    // reference simplex to the actual element in the mesh.
    // Note: this data is evaluated at several quadrature points.
    // gradient of the transformation (matrix)
    MAT_3x2 Map_PHI_Grad[NQ];
    // metric tensor of the map
    MAT_2x2 Map_PHI_Metric[NQ];
    // determinant of the metric matrix
    SCALAR Map_Det_Metric[NQ];
    // inverse of determinant of Metric
    SCALAR Map_Inv_Det_Metric[NQ];
    // inverse of the Metric tensor
    MAT_2x2 Map_PHI_Inv_Metric[NQ];
    // determinant of the transformation (Jacobian)
    SCALAR Map_Det_Jac[NQ];
    // determinant of Jacobian multiplied by quadrature weight
    SCALAR Map_Det_Jac_w_Weight[NQ];
    // hessian of the transformation
    MAT_3x2x2 Map_PHI_Hess[NQ];
    // gradient (in local coordinates) of the metric tensor of the map
    MAT_2x2x2 Map_PHI_Grad_Metric[NQ];
    // Christoffel symbols of the 2nd kind \Gamma^k_{i,j} (in local coordinates)
    MAT_2x2x2 Map_PHI_Christoffel_2nd_Kind[NQ];

    MGC (); // constructor
    ~MGC ();   // DE-structor
    void Setup_Mesh_Geometry(const mxArray*, const mxArray*, const mxArray*);
    void Compute_Local_Transformation();
    void Get_Current_Cell_Vertex_Indices(int Vtx_Indices[SUB_TD+1]) const
        {
        // transfer global vertex indices from kc[:]
        // Note: this code is custom made for this geometry class
        Vtx_Indices[0] = kc[0];
        Vtx_Indices[1] = kc[1];
        Vtx_Indices[2] = kc[2];
        }
    double Orientation[SUB_TD+1]; // mesh "facet" orientation direction for the current subdomain element
    // pointer to Domain class
    const CLASS_Domain_Gamma_embedded_in_Gamma_restricted_to_Gamma*   Domain;

private:
    double*  Node_Value[GD];    // mesh node values
    int*     Elem_DoF[NB];      // element DoF list
    bool*    Elem_Orient[SUB_TD+1]; // element facet orientation
                                    // true  = face has an outward normal vector (+1)
                                    // false = face has an  inward normal vector (-1)
    int      kc[NB];            // for storing the local mesh element DoFmap

    void Get_Local_to_Global_DoFmap(const int&, int*);
    void Get_Local_Orientation(const int&);
    void Compute_Map_p1(const int&);
};
/*------------   END: Auto Generate ------------*/

/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* constructor */
MGC::MGC () : Abstract_MESH_GEOMETRY_Class ()
{
    Type      = (char*) MAP_type;
    Global_TopDim = GLOBAL_TD;
    Sub_TopDim    = SUB_TD;
    DoI_TopDim    = DOI_TD;
    GeoDim    = GD;
    Num_QP    = NQ;
    Num_Basis = NB;

    // init mesh information to NULL
    Num_Nodes = 0;
    Num_Elem  = 0;
    for (int gd_i = 0; (gd_i < GeoDim); gd_i++)
        Node_Value[gd_i] = NULL;
    for (int basis_i = 0; (basis_i < Num_Basis); basis_i++)
        {
        Elem_DoF[basis_i] = NULL;
        kc[basis_i]       = -1; // a NULL value
        }
    for (int o_i = 0; (o_i < (Sub_TopDim+1)); o_i++)
        {
        Orientation[o_i] = +1.0;
        Elem_Orient[o_i] = NULL;
        }

    // init everything to zero
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        Map_PHI_Grad[qp_i].Set_To_Zero();
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        Map_PHI_Metric[qp_i].Set_To_Zero();
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        Map_Det_Metric[qp_i].Set_To_Zero();
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        Map_Inv_Det_Metric[qp_i].Set_To_Zero();
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        Map_PHI_Inv_Metric[qp_i].Set_To_Zero();
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        Map_Det_Jac[qp_i].Set_To_Zero();
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        Map_Det_Jac_w_Weight[qp_i].Set_To_Zero();
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        Map_PHI_Hess[qp_i].Set_To_Zero();
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        Map_PHI_Grad_Metric[qp_i].Set_To_Zero();
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        Map_PHI_Christoffel_2nd_Kind[qp_i].Set_To_Zero();
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/


/***************************************************************************************/
/* DE-structor */
MGC::~MGC ()
{
}
/***************************************************************************************/


/***************************************************************************************/
/* put incoming mesh data from MATLAB into a nice struct  */
void MGC::Setup_Mesh_Geometry(const mxArray *Vtx,       // inputs
                              const mxArray *Elem,      // inputs
                              const mxArray *Orient)    // inputs
{
    // get the ambient geometric dimension
    int CHK_GeoDim = (int) mxGetN(Vtx);
    // get the number of vertices
    Num_Nodes = (int) mxGetM(Vtx);
    // get the number of elements
    Num_Elem  = (int) mxGetM(Elem);
    // get the number of basis functions for each element
    int CHK_Num_Basis = (int) mxGetN(Elem);
    // get the number of rows, columns in the orientation data
    int CHK_Num_Row_Orient = (int) mxGetM(Orient);
    int CHK_Num_Col_Orient = (int) mxGetN(Orient);

    /* BEGIN: Simple Error Checking */
    if (mxGetClassID(Elem)!=mxUINT32_CLASS) mexErrMsgTxt("ERROR: Geometry DoFmap must be uint32!");
    if (CHK_GeoDim != GeoDim)
        {
        mexPrintf("ERROR: Vertex Coordinate List has %d columns; expected %d columns.\n", CHK_GeoDim, GeoDim);
        mexErrMsgTxt("ERROR: ambient geometric dimension must match!");
        }
    if (CHK_Num_Basis!=Num_Basis)
        {
        mexPrintf("ERROR: Mesh DoFmap has %d columns; expected %d columns.\n", CHK_Num_Basis, Num_Basis);
        mexPrintf("ERROR: A common reason for this error is you are using a finite element space\n");
        mexPrintf("ERROR:     to represent the mesh that is higher order than piecewise linear\n");
        mexPrintf("ERROR:     and you forgot to create a distinct DoFmap for that space.\n");
        mexPrintf("ERROR: You *cannot* just use the plain triangulation data!\n");
        mexPrintf("ERROR:     That only works for linear elements.\n");
        mexErrMsgTxt("ERROR: number of basis functions describing geometry must match!");
        }
    if (ORIENT) // if we should access orientation data, then make some checks
        {
        if (mxGetClassID(Orient)!=mxLOGICAL_CLASS) mexErrMsgTxt("ERROR: Mesh Orientation must be logical!");
        if (CHK_Num_Row_Orient!=Num_Elem)
            {
            mexPrintf("ERROR: Mesh Orientation has %d rows; expected %d rows.\n", CHK_Num_Row_Orient, Num_Elem);
            mexErrMsgTxt("ERROR: Orientation rows should match the Mesh DoFmap rows!");
            }
        if (CHK_Num_Col_Orient!=(Sub_TopDim+1))
            {
            mexPrintf("ERROR: Mesh Orientation has %d columns; expected %d columns.\n", CHK_Num_Col_Orient, (Sub_TopDim+1));
            mexErrMsgTxt("ERROR: Orientation cols should match the Mesh topological dimension + 1!");
            }
        }
    /* END: Simple Error Checking */


    // split up the columns of the node data
    Node_Value[0] = mxGetPr(Vtx);
    for (int gd_i = 1; (gd_i < GeoDim); gd_i++)
        Node_Value[gd_i] = Node_Value[gd_i-1] + Num_Nodes;

    // split up the columns of the element data
    Elem_DoF[0] = (int *) mxGetPr(Elem);
    for (int basis_i = 1; (basis_i < Num_Basis); basis_i++)
        Elem_DoF[basis_i] = Elem_DoF[basis_i-1] + Num_Elem;

    // split up the columns of the element (facet) orientation data
    if (ORIENT)
        {
        Elem_Orient[0] = (bool *) mxGetPr(Orient);
        for (int o_i = 1; (o_i < (Sub_TopDim+1)); o_i++)
            Elem_Orient[o_i] = Elem_Orient[o_i-1] + Num_Elem;
        }

    // get maximum DoF present in Elem
    int Elem_Num_Nodes  = *std::max_element(Elem_DoF[0],Elem_DoF[0] + (Num_Elem*Num_Basis));
    int Min_DoF         = *std::min_element(Elem_DoF[0],Elem_DoF[0] + (Num_Elem*Num_Basis));
    if ((Min_DoF < 1) || (Elem_Num_Nodes < 1))
        {
        mexPrintf("ERROR: There are Mesh DoFs that have indices < 1!\n");
        mexPrintf("ERROR: There are problems with this Mesh Type = %s!\n",Type);
        mexErrMsgTxt("ERROR: Fix your Mesh DoFmap!");
        }
    if (Elem_Num_Nodes > Num_Nodes)
        {
        mexPrintf("ERROR: There are Mesh DoFs that have indices > number of Mesh Values!\n");
        mexPrintf("ERROR: There are problems with this Mesh Type = %s!\n",Type);
        mexErrMsgTxt("ERROR: Fix your Mesh Values or DoFmap!");
        }
}
/***************************************************************************************/


/***************************************************************************************/
/* get the local DoFs on the given cell (element).
   Note: elem_index is in the   C-style (i.e. 0 <= elem_index <= Num_Elem - 1),
         Indices is in the MATLAB-style (i.e. 1 <= Indices[:] <= max(Elem_DoF)). */
void MGC::Get_Local_to_Global_DoFmap(const int& elem_index, int* Indices)  // inputs
{
    /* error check: */
    if (elem_index < 0)
        {
        mexPrintf("ERROR: Given cell index #%d is not positive. It must be > 0!\n",elem_index+1);
        mexPrintf("ERROR: There is an issue with a mesh of Type = %s!\n",Type);
        mexErrMsgTxt("ERROR: Make sure your inputs are valid!");
        }
    else if (elem_index >= Num_Elem)
        {
        mexPrintf("ERROR: Given cell index #%d exceeds the number of mesh cells. It must be <= %d!\n",elem_index+1,Num_Elem);
        mexPrintf("ERROR: There is an issue with a mesh of Type = %s!\n",Type);
        mexErrMsgTxt("ERROR: Make sure your inputs are valid!");
        }

    // get local to global indexing for geometry DoFmap
    for (int basis_i = 0; (basis_i < NB); basis_i++)
        {
        int DoF_index = Elem_DoF[basis_i][elem_index] - 1; // shifted for C - style indexing
        Indices[basis_i] = DoF_index;
        }
}
/***************************************************************************************/


/***************************************************************************************/
/* get the local orientation of the given element */
void MGC::Get_Local_Orientation(const int& elem_index)  // inputs
{
    // translate logical info to +/- 1.0
    for (int o_i = 0; (o_i < (SUB_TD+1)); o_i++)
        {
        const bool Orient_TF = Elem_Orient[o_i][elem_index];
        if (Orient_TF) Orientation[o_i] = +1.0;
        else           Orientation[o_i] = -1.0;
        }
}
/***************************************************************************************/


/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* compute the local transformation based on the given mesh entity */
void MGC::Compute_Local_Transformation()
{
    // read in the embedding info
    const int Global_Cell_Index = Domain->Global_Cell_Index;

    // compute "facet" orientation directions of current simplex
    if (ORIENT)  Get_Local_Orientation(Global_Cell_Index);

    /* compute local map */

    Compute_Map_p1(Global_Cell_Index);
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/


/*------------ BEGIN: Auto Generate ------------*/
/***************************************************************************************/
/* compute the local transformation from the standard reference element
          to an element in the Mesh    */
void MGC::Compute_Map_p1(const int& elem_index)           // current mesh element index
{
Get_Local_to_Global_DoFmap(elem_index, kc);

/*------------ BEGIN: Auto Generate ------------*/
// Intrinsic Map
// Local Map Element on Global mesh: CG, lagrange_deg2_dim2
// the Global mesh           has topological dimension = 2
// the Subdomain             has topological dimension = 2
// the Domain of Integration has topological dimension = 2
// Map has geometric dimension = 3
// Number of Quadrature Points = 9

    // Value of basis function, derivatives = [0  0  0], at quadrature points
    static const double Geo_Basis_Val_0_0_0[NQ][NB] = { \
    {-5.46685624375002968E-02, -9.37247465167690974E-02, -5.46685624375002552E-02, 2.18674249750001104E-01, 7.65713371891767469E-01, 2.18674249750001076E-01}, \
    {-5.46685624375002968E-02, -5.46685624375002552E-02, -9.37247465167690974E-02, 2.18674249750001104E-01, 2.18674249750001076E-01, 7.65713371891767469E-01}, \
    {-9.37247465167690696E-02, -5.46685624375002552E-02, -5.46685624375002552E-02, 7.65713371891767580E-01, 2.18674249750001021E-01, 2.18674249750001021E-01}, \
    {-3.46683066179297728E-02, 4.73664507650718436E-01, -1.10689039231616063E-01, 5.27401383462791751E-01, 2.47965497801223254E-02, 1.19494904955913278E-01}, \
    {-1.10689039231616063E-01, 4.73664507650718436E-01, -3.46683066179297519E-02, 1.19494904955913195E-01, 2.47965497801223116E-02, 5.27401383462791862E-01}, \
    {-3.46683066179297728E-02, -1.10689039231616063E-01, 4.73664507650718436E-01, 5.27401383462791751E-01, 1.19494904955913278E-01, 2.47965497801223254E-02}, \
    {4.73664507650718491E-01, -1.10689039231616063E-01, -3.46683066179297519E-02, 2.47965497801223081E-02, 1.19494904955913195E-01, 5.27401383462791862E-01}, \
    {-1.10689039231616063E-01, -3.46683066179297519E-02, 4.73664507650718436E-01, 1.19494904955913195E-01, 5.27401383462791862E-01, 2.47965497801223116E-02}, \
    {4.73664507650718491E-01, -3.46683066179297519E-02, -1.10689039231616063E-01, 2.47965497801223081E-02, 5.27401383462791862E-01, 1.19494904955913195E-01}  \
    };

    // Value of basis function, derivatives = [0  1  0], at quadrature points
    static const double Geo_Basis_Val_0_1_0[NQ][NB] = { \
    {-7.50100993533535876E-01, 0.00000000000000000E+00, 7.50100993533536098E-01, 4.99798012932928026E-01, -2.22044604925031308E-16, -4.99798012932928026E-01}, \
    {-7.50100993533535876E-01, 0.00000000000000000E+00, -5.00201987067071974E-01, 1.75010099353353610E+00, 1.25030298060060785E+00, -1.75010099353353610E+00}, \
    {5.00201987067072196E-01, 0.00000000000000000E+00, 7.50100993533536098E-01, 1.75010099353353610E+00, -1.25030298060060829E+00, -1.75010099353353610E+00}, \
    {8.50090316999647877E-01, 0.00000000000000000E+00, -3.38360290440635980E-01, 3.18845060744028386E+00, -5.11730026559011897E-01, -3.18845060744028386E+00}, \
    {3.38360290440635869E-01, 0.00000000000000000E+00, -8.50090316999647988E-01, 3.18845060744028386E+00, 5.11730026559012119E-01, -3.18845060744028386E+00}, \
    {8.50090316999647877E-01, 0.00000000000000000E+00, 2.18845060744028386E+00, 6.61639709559364020E-01, -3.03854092443993196E+00, -6.61639709559364020E-01}, \
    {-2.18845060744028386E+00, 0.00000000000000000E+00, -8.50090316999647988E-01, 6.61639709559364020E-01, 3.03854092443993196E+00, -6.61639709559364020E-01}, \
    {3.38360290440635869E-01, 0.00000000000000000E+00, 2.18845060744028386E+00, 1.49909683000352012E-01, -2.52681089788091962E+00, -1.49909683000352012E-01}, \
    {-2.18845060744028386E+00, 0.00000000000000000E+00, -3.38360290440635980E-01, 1.49909683000352012E-01, 2.52681089788092006E+00, -1.49909683000352012E-01}  \
    };

    // Value of basis function, derivatives = [0  2  0], at quadrature points
    static const double Geo_Basis_Val_0_2_0[NQ][NB] = { \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00, 0.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00, 0.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00, 0.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00, 0.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00, 0.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00, 0.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00, 0.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00, 0.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00, 0.00000000000000000E+00}  \
    };

    // Value of basis function, derivatives = [1  0  0], at quadrature points
    static const double Geo_Basis_Val_1_0_0[NQ][NB] = { \
    {-7.50100993533535876E-01, -5.00201987067071974E-01, 0.00000000000000000E+00, 1.75010099353353610E+00, -1.75010099353353610E+00, 1.25030298060060785E+00}, \
    {-7.50100993533535876E-01, 7.50100993533536098E-01, 0.00000000000000000E+00, 4.99798012932928026E-01, -4.99798012932928026E-01, -2.22044604925031308E-16}, \
    {5.00201987067072196E-01, 7.50100993533536098E-01, 0.00000000000000000E+00, 1.75010099353353610E+00, -1.75010099353353610E+00, -1.25030298060060829E+00}, \
    {8.50090316999647877E-01, 2.18845060744028386E+00, 0.00000000000000000E+00, 6.61639709559364020E-01, -6.61639709559364020E-01, -3.03854092443993196E+00}, \
    {3.38360290440635869E-01, 2.18845060744028386E+00, 0.00000000000000000E+00, 1.49909683000352012E-01, -1.49909683000352012E-01, -2.52681089788091962E+00}, \
    {8.50090316999647877E-01, -3.38360290440635980E-01, 0.00000000000000000E+00, 3.18845060744028386E+00, -3.18845060744028386E+00, -5.11730026559011897E-01}, \
    {-2.18845060744028386E+00, -3.38360290440635980E-01, 0.00000000000000000E+00, 1.49909683000352012E-01, -1.49909683000352012E-01, 2.52681089788092006E+00}, \
    {3.38360290440635869E-01, -8.50090316999647988E-01, 0.00000000000000000E+00, 3.18845060744028386E+00, -3.18845060744028386E+00, 5.11730026559012119E-01}, \
    {-2.18845060744028386E+00, -8.50090316999647988E-01, 0.00000000000000000E+00, 6.61639709559364020E-01, -6.61639709559364020E-01, 3.03854092443993196E+00}  \
    };

    // Value of basis function, derivatives = [1  1  0], at quadrature points
    static const double Geo_Basis_Val_1_1_0[NQ][NB] = { \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, -4.00000000000000000E+00, -4.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, -4.00000000000000000E+00, -4.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, -4.00000000000000000E+00, -4.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, -4.00000000000000000E+00, -4.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, -4.00000000000000000E+00, -4.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, -4.00000000000000000E+00, -4.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, -4.00000000000000000E+00, -4.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, -4.00000000000000000E+00, -4.00000000000000000E+00}, \
    {4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 4.00000000000000000E+00, -4.00000000000000000E+00, -4.00000000000000000E+00}  \
    };

    // Value of basis function, derivatives = [2  0  0], at quadrature points
    static const double Geo_Basis_Val_2_0_0[NQ][NB] = { \
    {4.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00}, \
    {4.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00}, \
    {4.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00}, \
    {4.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00}, \
    {4.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00}, \
    {4.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00}, \
    {4.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00}, \
    {4.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00}, \
    {4.00000000000000000E+00, 4.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, 0.00000000000000000E+00, -8.00000000000000000E+00}  \
    };

// // set of quadrature points
// static const double Quad_Points[NQ][GLOBAL_TD] = { \
//     {1.24949503233232007E-01, 4.37525248383384024E-01}, \
//     {4.37525248383384024E-01, 1.24949503233232007E-01}, \
//     {4.37525248383384024E-01, 4.37525248383384024E-01}, \
//     {7.97112651860070964E-01, 1.65409927389841005E-01}, \
//     {7.97112651860070964E-01, 3.74774207500880030E-02}, \
//     {1.65409927389841005E-01, 7.97112651860070964E-01}, \
//     {1.65409927389841005E-01, 3.74774207500880030E-02}, \
//     {3.74774207500880030E-02, 7.97112651860070964E-01}, \
//     {3.74774207500880030E-02, 1.65409927389841005E-01}  \
//     };
// set of quadrature weights
static const double Quad_Weights[NQ] = { \
    1.02975252380443499E-01, \
    1.02975252380443499E-01, \
    1.02975252380443499E-01, \
    3.18457071431115027E-02, \
    3.18457071431115027E-02, \
    3.18457071431115027E-02, \
    3.18457071431115027E-02, \
    3.18457071431115027E-02, \
    3.18457071431115027E-02  \
    };
/*------------   END: Auto Generate ------------*/
/*------------ BEGIN: Auto Generate ------------*/
    /*** compute geometric quantities ***/
    // compute the gradient of the local map
    // note: indexing is in the C style
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        // sum over basis functions
        MAT_3x2& Map_PHI_Grad_qp_i = Map_PHI_Grad[qp_i];
        Map_PHI_Grad_qp_i.m[0][0] = Geo_Basis_Val_1_0_0[qp_i][0]*Node_Value[0][kc[0]]+Geo_Basis_Val_1_0_0[qp_i][1]*Node_Value[0][kc[1]]+Geo_Basis_Val_1_0_0[qp_i][2]*Node_Value[0][kc[2]]+Geo_Basis_Val_1_0_0[qp_i][3]*Node_Value[0][kc[3]]+Geo_Basis_Val_1_0_0[qp_i][4]*Node_Value[0][kc[4]]+Geo_Basis_Val_1_0_0[qp_i][5]*Node_Value[0][kc[5]];
        Map_PHI_Grad_qp_i.m[0][1] = Geo_Basis_Val_0_1_0[qp_i][0]*Node_Value[0][kc[0]]+Geo_Basis_Val_0_1_0[qp_i][1]*Node_Value[0][kc[1]]+Geo_Basis_Val_0_1_0[qp_i][2]*Node_Value[0][kc[2]]+Geo_Basis_Val_0_1_0[qp_i][3]*Node_Value[0][kc[3]]+Geo_Basis_Val_0_1_0[qp_i][4]*Node_Value[0][kc[4]]+Geo_Basis_Val_0_1_0[qp_i][5]*Node_Value[0][kc[5]];
        Map_PHI_Grad_qp_i.m[1][0] = Geo_Basis_Val_1_0_0[qp_i][0]*Node_Value[1][kc[0]]+Geo_Basis_Val_1_0_0[qp_i][1]*Node_Value[1][kc[1]]+Geo_Basis_Val_1_0_0[qp_i][2]*Node_Value[1][kc[2]]+Geo_Basis_Val_1_0_0[qp_i][3]*Node_Value[1][kc[3]]+Geo_Basis_Val_1_0_0[qp_i][4]*Node_Value[1][kc[4]]+Geo_Basis_Val_1_0_0[qp_i][5]*Node_Value[1][kc[5]];
        Map_PHI_Grad_qp_i.m[1][1] = Geo_Basis_Val_0_1_0[qp_i][0]*Node_Value[1][kc[0]]+Geo_Basis_Val_0_1_0[qp_i][1]*Node_Value[1][kc[1]]+Geo_Basis_Val_0_1_0[qp_i][2]*Node_Value[1][kc[2]]+Geo_Basis_Val_0_1_0[qp_i][3]*Node_Value[1][kc[3]]+Geo_Basis_Val_0_1_0[qp_i][4]*Node_Value[1][kc[4]]+Geo_Basis_Val_0_1_0[qp_i][5]*Node_Value[1][kc[5]];
        Map_PHI_Grad_qp_i.m[2][0] = Geo_Basis_Val_1_0_0[qp_i][0]*Node_Value[2][kc[0]]+Geo_Basis_Val_1_0_0[qp_i][1]*Node_Value[2][kc[1]]+Geo_Basis_Val_1_0_0[qp_i][2]*Node_Value[2][kc[2]]+Geo_Basis_Val_1_0_0[qp_i][3]*Node_Value[2][kc[3]]+Geo_Basis_Val_1_0_0[qp_i][4]*Node_Value[2][kc[4]]+Geo_Basis_Val_1_0_0[qp_i][5]*Node_Value[2][kc[5]];
        Map_PHI_Grad_qp_i.m[2][1] = Geo_Basis_Val_0_1_0[qp_i][0]*Node_Value[2][kc[0]]+Geo_Basis_Val_0_1_0[qp_i][1]*Node_Value[2][kc[1]]+Geo_Basis_Val_0_1_0[qp_i][2]*Node_Value[2][kc[2]]+Geo_Basis_Val_0_1_0[qp_i][3]*Node_Value[2][kc[3]]+Geo_Basis_Val_0_1_0[qp_i][4]*Node_Value[2][kc[4]]+Geo_Basis_Val_0_1_0[qp_i][5]*Node_Value[2][kc[5]];
        }
    // compute metric tensor from jacobian matrix
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        //
        MAT_2x2& Map_PHI_Metric_qp_i = Map_PHI_Metric[qp_i];
        Mat_Transpose_Mat_Self(Map_PHI_Grad[qp_i], Map_PHI_Metric_qp_i);
        }
    // compute determinant of Metric: det(PHI_Metric)
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        //
        SCALAR& Map_Det_Metric_qp_i = Map_Det_Metric[qp_i];
        Map_Det_Metric_qp_i.a = Determinant(Map_PHI_Metric[qp_i]);
        }
    // compute 1 / det(Metric)
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        //
        SCALAR& Map_Inv_Det_Metric_qp_i = Map_Inv_Det_Metric[qp_i];
        Map_Inv_Det_Metric_qp_i.a = 1.0 / Map_Det_Metric[qp_i].a;
        }
    // compute inverse of the Metric tensor
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        //
        MAT_2x2& Map_PHI_Inv_Metric_qp_i = Map_PHI_Inv_Metric[qp_i];
        Matrix_Inverse(Map_PHI_Metric[qp_i], Map_Inv_Det_Metric[qp_i], Map_PHI_Inv_Metric_qp_i);
        }
    // compute determinant of Jacobian
    // note: det(Jac) = sqrt(det(PHI_Metric))
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        //
        SCALAR& Map_Det_Jac_qp_i = Map_Det_Jac[qp_i];
        Map_Det_Jac_qp_i.a = sqrt(Map_Det_Metric[qp_i].a);
        }
    // multiply det(jacobian) by quadrature weight
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        //
        SCALAR& Map_Det_Jac_w_Weight_qp_i = Map_Det_Jac_w_Weight[qp_i];
        Map_Det_Jac_w_Weight_qp_i.a = Map_Det_Jac[qp_i].a * Quad_Weights[qp_i];
        }
    // compute the hessian of the local map
    // note: indexing is in the C style
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        // sum over basis functions
        MAT_3x2x2& Map_PHI_Hess_qp_i = Map_PHI_Hess[qp_i];
        Map_PHI_Hess_qp_i.m[0][0][0] = Geo_Basis_Val_2_0_0[qp_i][0] * Node_Value[0][kc[0]] + Geo_Basis_Val_2_0_0[qp_i][1] * Node_Value[0][kc[1]] + Geo_Basis_Val_2_0_0[qp_i][2] * Node_Value[0][kc[2]] + Geo_Basis_Val_2_0_0[qp_i][3] * Node_Value[0][kc[3]] + Geo_Basis_Val_2_0_0[qp_i][4] * Node_Value[0][kc[4]] + Geo_Basis_Val_2_0_0[qp_i][5] * Node_Value[0][kc[5]];
        Map_PHI_Hess_qp_i.m[0][0][1] = Geo_Basis_Val_1_1_0[qp_i][0] * Node_Value[0][kc[0]] + Geo_Basis_Val_1_1_0[qp_i][1] * Node_Value[0][kc[1]] + Geo_Basis_Val_1_1_0[qp_i][2] * Node_Value[0][kc[2]] + Geo_Basis_Val_1_1_0[qp_i][3] * Node_Value[0][kc[3]] + Geo_Basis_Val_1_1_0[qp_i][4] * Node_Value[0][kc[4]] + Geo_Basis_Val_1_1_0[qp_i][5] * Node_Value[0][kc[5]];
        Map_PHI_Hess_qp_i.m[0][1][0] = Map_PHI_Hess_qp_i.m[0][0][1]; // symmetry!
        Map_PHI_Hess_qp_i.m[0][1][1] = Geo_Basis_Val_0_2_0[qp_i][0] * Node_Value[0][kc[0]] + Geo_Basis_Val_0_2_0[qp_i][1] * Node_Value[0][kc[1]] + Geo_Basis_Val_0_2_0[qp_i][2] * Node_Value[0][kc[2]] + Geo_Basis_Val_0_2_0[qp_i][3] * Node_Value[0][kc[3]] + Geo_Basis_Val_0_2_0[qp_i][4] * Node_Value[0][kc[4]] + Geo_Basis_Val_0_2_0[qp_i][5] * Node_Value[0][kc[5]];
        Map_PHI_Hess_qp_i.m[1][0][0] = Geo_Basis_Val_2_0_0[qp_i][0] * Node_Value[1][kc[0]] + Geo_Basis_Val_2_0_0[qp_i][1] * Node_Value[1][kc[1]] + Geo_Basis_Val_2_0_0[qp_i][2] * Node_Value[1][kc[2]] + Geo_Basis_Val_2_0_0[qp_i][3] * Node_Value[1][kc[3]] + Geo_Basis_Val_2_0_0[qp_i][4] * Node_Value[1][kc[4]] + Geo_Basis_Val_2_0_0[qp_i][5] * Node_Value[1][kc[5]];
        Map_PHI_Hess_qp_i.m[1][0][1] = Geo_Basis_Val_1_1_0[qp_i][0] * Node_Value[1][kc[0]] + Geo_Basis_Val_1_1_0[qp_i][1] * Node_Value[1][kc[1]] + Geo_Basis_Val_1_1_0[qp_i][2] * Node_Value[1][kc[2]] + Geo_Basis_Val_1_1_0[qp_i][3] * Node_Value[1][kc[3]] + Geo_Basis_Val_1_1_0[qp_i][4] * Node_Value[1][kc[4]] + Geo_Basis_Val_1_1_0[qp_i][5] * Node_Value[1][kc[5]];
        Map_PHI_Hess_qp_i.m[1][1][0] = Map_PHI_Hess_qp_i.m[1][0][1]; // symmetry!
        Map_PHI_Hess_qp_i.m[1][1][1] = Geo_Basis_Val_0_2_0[qp_i][0] * Node_Value[1][kc[0]] + Geo_Basis_Val_0_2_0[qp_i][1] * Node_Value[1][kc[1]] + Geo_Basis_Val_0_2_0[qp_i][2] * Node_Value[1][kc[2]] + Geo_Basis_Val_0_2_0[qp_i][3] * Node_Value[1][kc[3]] + Geo_Basis_Val_0_2_0[qp_i][4] * Node_Value[1][kc[4]] + Geo_Basis_Val_0_2_0[qp_i][5] * Node_Value[1][kc[5]];
        Map_PHI_Hess_qp_i.m[2][0][0] = Geo_Basis_Val_2_0_0[qp_i][0] * Node_Value[2][kc[0]] + Geo_Basis_Val_2_0_0[qp_i][1] * Node_Value[2][kc[1]] + Geo_Basis_Val_2_0_0[qp_i][2] * Node_Value[2][kc[2]] + Geo_Basis_Val_2_0_0[qp_i][3] * Node_Value[2][kc[3]] + Geo_Basis_Val_2_0_0[qp_i][4] * Node_Value[2][kc[4]] + Geo_Basis_Val_2_0_0[qp_i][5] * Node_Value[2][kc[5]];
        Map_PHI_Hess_qp_i.m[2][0][1] = Geo_Basis_Val_1_1_0[qp_i][0] * Node_Value[2][kc[0]] + Geo_Basis_Val_1_1_0[qp_i][1] * Node_Value[2][kc[1]] + Geo_Basis_Val_1_1_0[qp_i][2] * Node_Value[2][kc[2]] + Geo_Basis_Val_1_1_0[qp_i][3] * Node_Value[2][kc[3]] + Geo_Basis_Val_1_1_0[qp_i][4] * Node_Value[2][kc[4]] + Geo_Basis_Val_1_1_0[qp_i][5] * Node_Value[2][kc[5]];
        Map_PHI_Hess_qp_i.m[2][1][0] = Map_PHI_Hess_qp_i.m[2][0][1]; // symmetry!
        Map_PHI_Hess_qp_i.m[2][1][1] = Geo_Basis_Val_0_2_0[qp_i][0] * Node_Value[2][kc[0]] + Geo_Basis_Val_0_2_0[qp_i][1] * Node_Value[2][kc[1]] + Geo_Basis_Val_0_2_0[qp_i][2] * Node_Value[2][kc[2]] + Geo_Basis_Val_0_2_0[qp_i][3] * Node_Value[2][kc[3]] + Geo_Basis_Val_0_2_0[qp_i][4] * Node_Value[2][kc[4]] + Geo_Basis_Val_0_2_0[qp_i][5] * Node_Value[2][kc[5]];
        }
    // compute gradient of metric tensor
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        //
        MAT_2x2x2& Map_PHI_Grad_Metric_qp_i = Map_PHI_Grad_Metric[qp_i];
        // compute d/ds_r d/ds_i vX DOT d/ds_j vX + d/ds_r d/ds_j vX DOT d/ds_i vX
        double HG_DP_A;
        HG_DP_A = Map_PHI_Hess[qp_i].m[0][0][0] * Map_PHI_Grad[qp_i].m[0][0] + Map_PHI_Hess[qp_i].m[1][0][0] * Map_PHI_Grad[qp_i].m[1][0] + Map_PHI_Hess[qp_i].m[2][0][0] * Map_PHI_Grad[qp_i].m[2][0];
        Map_PHI_Grad_Metric_qp_i.m[0][0][0] = 2.0 * HG_DP_A;
        HG_DP_A = Map_PHI_Hess[qp_i].m[0][0][1] * Map_PHI_Grad[qp_i].m[0][1] + Map_PHI_Hess[qp_i].m[1][0][1] * Map_PHI_Grad[qp_i].m[1][1] + Map_PHI_Hess[qp_i].m[2][0][1] * Map_PHI_Grad[qp_i].m[2][1];
        Map_PHI_Grad_Metric_qp_i.m[0][1][1] = 2.0 * HG_DP_A;
        HG_DP_A = Map_PHI_Hess[qp_i].m[0][1][0] * Map_PHI_Grad[qp_i].m[0][0] + Map_PHI_Hess[qp_i].m[1][1][0] * Map_PHI_Grad[qp_i].m[1][0] + Map_PHI_Hess[qp_i].m[2][1][0] * Map_PHI_Grad[qp_i].m[2][0];
        Map_PHI_Grad_Metric_qp_i.m[1][0][0] = 2.0 * HG_DP_A;
        HG_DP_A = Map_PHI_Hess[qp_i].m[0][1][1] * Map_PHI_Grad[qp_i].m[0][1] + Map_PHI_Hess[qp_i].m[1][1][1] * Map_PHI_Grad[qp_i].m[1][1] + Map_PHI_Hess[qp_i].m[2][1][1] * Map_PHI_Grad[qp_i].m[2][1];
        Map_PHI_Grad_Metric_qp_i.m[1][1][1] = 2.0 * HG_DP_A;
        HG_DP_A = Map_PHI_Hess[qp_i].m[0][0][0] * Map_PHI_Grad[qp_i].m[0][1] + Map_PHI_Hess[qp_i].m[1][0][0] * Map_PHI_Grad[qp_i].m[1][1] + Map_PHI_Hess[qp_i].m[2][0][0] * Map_PHI_Grad[qp_i].m[2][1];
        Map_PHI_Grad_Metric_qp_i.m[0][0][1] = HG_DP_A + 0.5*Map_PHI_Grad_Metric_qp_i.m[1][0][0]; // reuse!
        Map_PHI_Grad_Metric_qp_i.m[0][1][0] = Map_PHI_Grad_Metric_qp_i.m[0][0][1]; // symmetric
        HG_DP_A = Map_PHI_Hess[qp_i].m[0][1][1] * Map_PHI_Grad[qp_i].m[0][0] + Map_PHI_Hess[qp_i].m[1][1][1] * Map_PHI_Grad[qp_i].m[1][0] + Map_PHI_Hess[qp_i].m[2][1][1] * Map_PHI_Grad[qp_i].m[2][0];
        Map_PHI_Grad_Metric_qp_i.m[1][1][0] = HG_DP_A + 0.5*Map_PHI_Grad_Metric_qp_i.m[0][1][1]; // reuse!
        Map_PHI_Grad_Metric_qp_i.m[1][0][1] = Map_PHI_Grad_Metric_qp_i.m[1][1][0]; // symmetric
        }
    // compute Christoffel symbols of the 2nd kind \Gamma^k_{i,j}
    // loop through quad points
    for (int qp_i = 0; (qp_i < Num_QP); qp_i++)
        {
        //
        MAT_2x2x2& Map_PHI_Christoffel_2nd_Kind_qp_i = Map_PHI_Christoffel_2nd_Kind[qp_i];
        // compute each component of \Gamma^k_{i,j} individually
        for (int kk = 0; (kk < 2); kk++)
        for (int ii = 0; (ii < 2); ii++)
        for (int jj = 0; (jj < 2); jj++)
            {
            const double Deriv_Metric_Term_0 = Map_PHI_Grad_Metric[qp_i].m[ii][0][jj] + Map_PHI_Grad_Metric[qp_i].m[jj][ii][0] - Map_PHI_Grad_Metric[qp_i].m[0][ii][jj];
            const double Deriv_Metric_Term_1 = Map_PHI_Grad_Metric[qp_i].m[ii][1][jj] + Map_PHI_Grad_Metric[qp_i].m[jj][ii][1] - Map_PHI_Grad_Metric[qp_i].m[1][ii][jj];
            Map_PHI_Christoffel_2nd_Kind_qp_i.m[kk][ii][jj] = (1.0/2.0) *
                 ( Map_PHI_Inv_Metric[qp_i].m[kk][0] * Deriv_Metric_Term_0
                +  Map_PHI_Inv_Metric[qp_i].m[kk][1] * Deriv_Metric_Term_1 );
            }
        }
/*------------   END: Auto Generate ------------*/
}
/***************************************************************************************/
/*------------   END: Auto Generate ------------*/

#undef MGC

#undef MAP_type
#undef GLOBAL_TD
#undef SUB_TD
#undef DOI_TD
#undef GD
#undef NQ
#undef NB
#undef ORIENT

/***/
