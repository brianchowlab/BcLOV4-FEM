
    // Compute interpolation

    // only one point for interpolation
    for (unsigned int qp = 0; qp < 1; qp++)
          INTERP = geom_Gamma_embedded_in_Gamma_restricted_to_Sigma->Map_Normal_Vector[0].v[0]*geom_Sigma_embedded_in_Gamma_restricted_to_Sigma->Map_Tangent_Vector[0].v[0]+geom_Gamma_embedded_in_Gamma_restricted_to_Sigma->Map_Normal_Vector[0].v[1]*geom_Sigma_embedded_in_Gamma_restricted_to_Sigma->Map_Tangent_Vector[0].v[1]+geom_Gamma_embedded_in_Gamma_restricted_to_Sigma->Map_Normal_Vector[0].v[2]*geom_Sigma_embedded_in_Gamma_restricted_to_Sigma->Map_Tangent_Vector[0].v[2];
