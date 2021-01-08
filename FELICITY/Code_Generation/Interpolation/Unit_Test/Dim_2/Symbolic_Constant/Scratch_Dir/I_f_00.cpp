
    // Compute interpolation

    // only one point for interpolation
    for (unsigned int qp = 0; qp < 1; qp++)
          INTERP = c0_restricted_to_Omega->Constant_C_Value[0].a*f_restricted_to_Omega->Func_f_Grad[0][qp].v[0]+c0_restricted_to_Omega->Constant_C_Value[1].a*f_restricted_to_Omega->Func_f_Grad[0][qp].v[1];
