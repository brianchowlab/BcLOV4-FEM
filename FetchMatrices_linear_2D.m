function [K_phi,K_psi,K_phi_phi,K_phi_psi,K_psi_phi,A_phi,A_psi]=FetchMatrices_2D(params)
    FEM = mex_ReactionDiffusion_Cell_l_2D_assemble(params.FEM_l,params.p,params.c,[],params.embed,params.C_DoF,params.M_DoF);
    FEM_Mats = FEMatrixAccessor('Problem',FEM);
    clear FEM

    %Bilinear C & C over Omega
    K_phi = FEM_Mats.Get_Matrix('K_phi');

    %Bilinear M & M
    K_psi = FEM_Mats.Get_Matrix('K_psi');
    
    %Bilinear C & C over dOmega
    K_phi_phi = FEM_Mats.Get_Matrix('K_phi_phi');

    %Bilinear C & M
    K_phi_psi = FEM_Mats.Get_Matrix('K_phi_psi');

    %Bilinear M & C
    K_psi_phi = FEM_Mats.Get_Matrix('K_psi_phi');

    %Mass matrix (A)
    A_phi = FEM_Mats.Get_Matrix('A_phi');
    A_psi = FEM_Mats.Get_Matrix('A_psi');
end





