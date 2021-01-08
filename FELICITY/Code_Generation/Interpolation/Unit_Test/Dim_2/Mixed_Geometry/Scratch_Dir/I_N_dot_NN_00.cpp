
    // Compute interpolation

    // only one point for interpolation
    for (unsigned int qp = 0; qp < 1; qp++)
          INTERP = NN_restricted_to_Sigma->Func_f_Value[0][qp].a*geom_Gamma_embedded_in_Gamma_restricted_to_Sigma->Map_Normal_Vector[0].v[0]+NN_restricted_to_Sigma->Func_f_Value[1][qp].a*geom_Gamma_embedded_in_Gamma_restricted_to_Sigma->Map_Normal_Vector[0].v[1]+NN_restricted_to_Sigma->Func_f_Value[2][qp].a*geom_Gamma_embedded_in_Gamma_restricted_to_Sigma->Map_Normal_Vector[0].v[2];
