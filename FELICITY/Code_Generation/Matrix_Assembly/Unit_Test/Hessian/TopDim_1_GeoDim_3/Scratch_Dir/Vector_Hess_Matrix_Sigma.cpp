/***************************************************************************************/
/* Tabulate the tensor for the local element contribution */
void SpecificFEM::Tabulate_Tensor(const CLASS_geom_Sigma_embedded_in_Sigma_restricted_to_Sigma& Mesh)
{
    const unsigned int NQ = 4;

    // Compute element tensor using quadrature

    // Loop quadrature points for integral
    for (unsigned int j = 0; j < COL_NB; j++)
        {
        for (unsigned int i = j; i < ROW_NB; i++)
            {
            double  A0_value = 0.0; // initialize
            for (unsigned int qp = 0; qp < NQ; qp++)
                {
                const double  integrand_0 = Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[j][qp].m[0][0]*Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[i][qp].m[0][0]+Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[j][qp].m[0][1]*Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[i][qp].m[0][1]+Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[j][qp].m[0][2]*Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[i][qp].m[0][2]+Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[j][qp].m[1][0]*Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[i][qp].m[1][0]+Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[j][qp].m[1][1]*Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[i][qp].m[1][1]+Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[j][qp].m[1][2]*Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[i][qp].m[1][2]+Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[j][qp].m[2][0]*Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[i][qp].m[2][0]+Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[j][qp].m[2][1]*Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[i][qp].m[2][1]+Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[j][qp].m[2][2]*Vector_P2_phi_restricted_to_Sigma->Func_f_Hess[i][qp].m[2][2];
                A0_value += integrand_0 * Mesh.Map_Det_Jac_w_Weight[qp].a;
                }
            FE_Tensor_0[j*ROW_NB + i] = A0_value;
            }
        }

    // Copy the lower triangular entries to the upper triangular part (by symmetry)
    for (unsigned int j = 0; j < COL_NB; j++)
        {
        for (unsigned int i = j+1; i < ROW_NB; i++)
            {
            FE_Tensor_0[i*ROW_NB + j] = FE_Tensor_0[j*ROW_NB + i];
            }
        }
}
/***************************************************************************************/
