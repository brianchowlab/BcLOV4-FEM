/***************************************************************************************/
/* Tabulate the tensor for the local element contribution */
void SpecificFEM::Tabulate_Tensor(const CLASS_geom_Gamma_embedded_in_Global_restricted_to_Gamma& Mesh)
{
    const unsigned int NQ = 3;

    // Compute single entry using quadrature

    // Loop quadrature points for integral
    double  A0_value = 0.0; // initialize
    double  A1_value = 0.0; // initialize
    for (unsigned int qp = 0; qp < NQ; qp++)
        {
        const double  integrand_0 = old_soln_restricted_to_Gamma->Func_f_Value[0][qp].a;
        const double  integrand_1 = geom_Gamma_embedded_in_Global_restricted_to_Gamma->Map_PHI[qp].v[1];
        A0_value += integrand_0 * Mesh.Map_Det_Jac_w_Weight[qp].a;
        A1_value += integrand_1 * Mesh.Map_Det_Jac_w_Weight[qp].a;
        }
    FE_Tensor_0[0] = A0_value;
    FE_Tensor_1[0] = A1_value;
}
/***************************************************************************************/
