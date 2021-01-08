/***************************************************************************************/
/* Tabulate the tensor for the local element contribution */
void SpecificFEM::Tabulate_Tensor(const CLASS_geom_Sigma_embedded_in_Sigma_restricted_to_Sigma& Mesh)
{
    const unsigned int NQ = 4;

    // Compute single entry using quadrature

    // Loop quadrature points for integral
    double  A0_value = 0.0; // initialize
    double  A1_value = 0.0; // initialize
    for (unsigned int qp = 0; qp < NQ; qp++)
        {
        const double  integrand_0 = Sca_restricted_to_Sigma->Func_f_Hess[0][qp].m[0][0];
        const double  integrand_1 = Vec_restricted_to_Sigma->Func_f_d2_ds2[1][qp].a;
        A0_value += integrand_0 * Mesh.Map_Det_Jac_w_Weight[qp].a;
        A1_value += integrand_1 * Mesh.Map_Det_Jac_w_Weight[qp].a;
        }
    FE_Tensor_0[0] = A0_value;
    FE_Tensor_1[0] = A1_value;
}
/***************************************************************************************/
