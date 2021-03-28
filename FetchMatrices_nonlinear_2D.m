function [B_phi,B_psi]=FetchMatrices_2D(params)
    FEM = mex_ReactionDiffusion_Cell_nl_2D_assemble(params.FEM_nl,params.p,params.c,[],...
        params.embed,params.C_DoF,params.M_DoF,params.u_M,params.v_M);
    FEM_Mats = FEMatrixAccessor('Problem',FEM);
    clear FEM

    %Non-linear components
    B_phi = FEM_Mats.Get_Matrix('B_phi');
    B_psi = FEM_Mats.Get_Matrix('B_psi');
end





