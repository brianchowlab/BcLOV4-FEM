function MATS = MatAssem_ReactionDiffusion_Cell()
    %Define domain (3-D volume)
    Omega = Domain('tetrahedron');
    
    %Define surface
    dOmega = Domain('triangle') < Omega;

    %Define finite element spaces (in bulk and surface)
    C_h = Element(Omega, lagrange_deg1_dim3,1); % piecewise linear
    M_h = Element(dOmega,lagrange_deg1_dim2,1); % piecewise linear
    C_h_t = Element(Omega, lagrange_deg1_dim3,2); % piecewise linear
    M_h_t = Element(dOmega,lagrange_deg1_dim2,2); % piecewise linear

    %Define functions on FE spaces
    phi_h = Test(C_h);
    u_C_h = Trial(C_h);

    psi_h = Test(M_h);
    u_M_h = Trial(M_h);

    %Define FEM matrices
    %Stiffness matrices (K)
    
    %Bilinear C & C over Omega
    K_phi = Bilinear(C_h,C_h);
    K_phi = K_phi + Integral(Omega, phi_h.val' * u_C_h.val );
    
    %Bilinear C & C over dOmega
    K_phi_phi = Bilinear(C_h,C_h);
    K_phi_phi = K_phi_phi + Integral(dOmega, phi_h.val' * u_C_h.val );
    
    %Bilinear M & M over dOmega
    K_psi = Bilinear(M_h,M_h);
    K_psi = K_psi + Integral(dOmega, psi_h.val' * u_M_h.val );
    
    %Bilinear C & M over dOmega
    K_psi_phi = Bilinear(C_h,M_h);
    K_psi_phi = K_psi_phi + Integral(dOmega, phi_h.val' * u_M_h.val );
    
    %Bilinear M & C over dOmega
    K_phi_psi = Bilinear(M_h,C_h);
    K_phi_psi = K_phi_psi + Integral(dOmega, psi_h.val' * u_C_h.val );
    
    %Mass matrix (A)
    A_phi = Bilinear(C_h,C_h);
    A_phi = A_phi + Integral(Omega, phi_h.grad' * u_C_h.grad );
    
    A_psi = Bilinear(M_h,M_h);
    A_psi = A_psi + Integral(dOmega, psi_h.grad' * u_M_h.grad );

    % set the minimum order of accuracy for the quad rule
    Quadrature_Order = 10;
    % define geometry representation - Domain, (default to piecewise linear)
    G1 = GeoElement(Omega);
    % define a set of matrices
    MATS = Matrices(Quadrature_Order,G1);

    % collect all of the matrices together
    MATS = MATS.Append_Matrix(K_phi);
    MATS = MATS.Append_Matrix(K_psi);
    MATS = MATS.Append_Matrix(K_phi_phi);
    MATS = MATS.Append_Matrix(K_psi_phi);
    MATS = MATS.Append_Matrix(K_phi_psi);
    MATS = MATS.Append_Matrix(A_phi);
    MATS = MATS.Append_Matrix(A_psi);
end