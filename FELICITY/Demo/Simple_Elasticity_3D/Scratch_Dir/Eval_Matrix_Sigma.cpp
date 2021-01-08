/***************************************************************************************/
/* Tabulate the tensor for the local element contribution */
void SpecificFEM::Tabulate_Tensor(const CLASS_geom_Sigma_embedded_in_Omega_restricted_to_Sigma& Mesh)
{
    const unsigned int NQ = 2;

    // Compute single entry using quadrature

    // Loop quadrature points for integral
    double  A0_value = 0.0; // initialize
    for (unsigned int qp = 0; qp < NQ; qp++)
        {
        const double  integrand_0 = Displace_restricted_to_Sigma->Func_f_Value[0][qp].a*geom_Sigma_embedded_in_Omega_restricted_to_Sigma->Map_Tangent_Vector[0].v[0]+Displace_restricted_to_Sigma->Func_f_Value[1][qp].a*geom_Sigma_embedded_in_Omega_restricted_to_Sigma->Map_Tangent_Vector[0].v[1]+Displace_restricted_to_Sigma->Func_f_Value[2][qp].a*geom_Sigma_embedded_in_Omega_restricted_to_Sigma->Map_Tangent_Vector[0].v[2];
        A0_value += integrand_0 * Mesh.Map_Det_Jac_w_Weight[qp].a;
        }
    FE_Tensor_0[0] = A0_value;
}
/***************************************************************************************/
