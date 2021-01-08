/***************************************************************************************/
/* Tabulate the tensor for the local element contribution */
void SpecificFEM::Tabulate_Tensor(const CLASS_geom_Omega_embedded_in_Omega_restricted_to_Omega& Mesh)
{
    const unsigned int NQ = 4;

    // Compute element tensor using quadrature

    // Loop quadrature points for integral
    for (unsigned int j = 0; j < COL_NB; j++)
        {
        for (unsigned int i = 0; i < ROW_NB; i++)
            {
            double  A0_value = 0.0; // initialize
            double  A1_value = 0.0; // initialize
            double  A3_value = 0.0; // initialize
            for (unsigned int qp = 0; qp < NQ; qp++)
                {
                const double  t2 = Vector_P2_phi_restricted_to_Omega->Func_f_Grad[j][qp].v[0]*Vector_P2_phi_restricted_to_Omega->Func_f_Grad[i][qp].v[0];
                const double  t3 = Vector_P2_phi_restricted_to_Omega->Func_f_Grad[j][qp].v[1]*Vector_P2_phi_restricted_to_Omega->Func_f_Grad[i][qp].v[1];
                const double  integrand_0 = t2*2.0+t3;
                const double  integrand_1 = Vector_P2_phi_restricted_to_Omega->Func_f_Grad[j][qp].v[1]*Vector_P2_phi_restricted_to_Omega->Func_f_Grad[i][qp].v[0];
                const double  integrand_3 = t2+t3*2.0;
                A0_value += integrand_0 * Mesh.Map_Det_Jac_w_Weight[qp].a;
                A1_value += integrand_1 * Mesh.Map_Det_Jac_w_Weight[qp].a;
                A3_value += integrand_3 * Mesh.Map_Det_Jac_w_Weight[qp].a;
                }
            FE_Tensor_0[j*ROW_NB + i] = A0_value;
            FE_Tensor_1[j*ROW_NB + i] = A1_value;
            FE_Tensor_3[j*ROW_NB + i] = A3_value;
            }
        }
}
/***************************************************************************************/
