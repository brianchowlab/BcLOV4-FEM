/***************************************************************************************/
/* Tabulate the tensor for the local element contribution */
void SpecificFEM::Tabulate_Tensor(const CLASS_geom_Sigma_embedded_in_Global_restricted_to_Sigma& Mesh)
{
    const unsigned int NQ = 3;

    // Compute element tensor using quadrature

    // Loop quadrature points for integral
    for (unsigned int j = 0; j < COL_NB; j++)
        {
        for (unsigned int i = 0; i < ROW_NB; i++)
            {
            double  A0_value = 0.0; // initialize
            for (unsigned int qp = 0; qp < NQ; qp++)
                {
                const double  integrand_0 = E_Space_phi_restricted_to_Sigma->Func_f_Grad[j][qp].v[0]*M_Space_phi_restricted_to_Sigma->Func_f_Grad[i][qp].v[0]+E_Space_phi_restricted_to_Sigma->Func_f_Grad[j][qp].v[1]*M_Space_phi_restricted_to_Sigma->Func_f_Grad[i][qp].v[1]+E_Space_phi_restricted_to_Sigma->Func_f_Grad[j][qp].v[2]*M_Space_phi_restricted_to_Sigma->Func_f_Grad[i][qp].v[2];
                A0_value += integrand_0 * Mesh.Map_Det_Jac_w_Weight[qp].a;
                }
            FE_Tensor_0[j*ROW_NB + i] = A0_value;
            }
        }
}
/***************************************************************************************/
