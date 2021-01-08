
    // Compute interpolation

    // only one point for interpolation
    for (unsigned int qp = 0; qp < 1; qp++)
          INTERP = geom_Gamma_embedded_in_Omega_restricted_to_Gamma->Map_Normal_Vector[0].v[0]*v_restricted_to_Gamma->Func_vv_Value[0][qp].v[0]+geom_Gamma_embedded_in_Omega_restricted_to_Gamma->Map_Normal_Vector[0].v[1]*v_restricted_to_Gamma->Func_vv_Value[0][qp].v[1];
