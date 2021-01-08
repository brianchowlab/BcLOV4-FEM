
    // Compute interpolation

    // only one point for interpolation
    for (unsigned int qp = 0; qp < 1; qp++)
          INTERP = geom_Omega_embedded_in_Omega_restricted_to_Omega->Map_PHI[qp].v[0]*p_restricted_to_Omega->Func_f_Grad[0][qp].v[0]+geom_Omega_embedded_in_Omega_restricted_to_Omega->Map_PHI[qp].v[1]*p_restricted_to_Omega->Func_f_Grad[0][qp].v[1];
