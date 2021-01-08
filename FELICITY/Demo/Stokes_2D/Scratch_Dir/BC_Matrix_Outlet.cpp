/***************************************************************************************/
/* Tabulate the tensor for the local element contribution */
void SpecificFEM::Tabulate_Tensor(const CLASS_geom_Outlet_embedded_in_Omega_restricted_to_Outlet& Mesh)
{
    const unsigned int NQ = 2;

    // Compute element tensor using quadrature

    // Loop quadrature points for integral
    for (unsigned int i = 0; i < ROW_NB; i++)
        {
        double  A0_value = 0.0; // initialize
        double  A1_value = 0.0; // initialize
        for (unsigned int qp = 0; qp < NQ; qp++)
            {
            const double  integrand_0 = BC_Out_restricted_to_Outlet->Func_f_Value[0][qp].a*(*Vector_P2_phi_restricted_to_Outlet->Func_f_Value)[i][qp].a;
            const double  integrand_1 = BC_Out_restricted_to_Outlet->Func_f_Value[1][qp].a*(*Vector_P2_phi_restricted_to_Outlet->Func_f_Value)[i][qp].a;
            A0_value += integrand_0 * Mesh.Map_Det_Jac_w_Weight[qp].a;
            A1_value += integrand_1 * Mesh.Map_Det_Jac_w_Weight[qp].a;
            }
        FE_Tensor_0[i] = A0_value;
        FE_Tensor_1[i] = A1_value;
        }
}
/***************************************************************************************/
