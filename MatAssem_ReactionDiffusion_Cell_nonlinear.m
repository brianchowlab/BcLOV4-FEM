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
    
    %Nonlinear terms
    u_M_h_coef = Coef(M_h);
    v_M_h_coef = Coef(M_h);
    
    B_phi = Bilinear(C_h,C_h);
    B_phi = B_phi + Integral(dOmega,(u_M_h_coef.val+v_M_h_coef.val)* phi_h.val * u_C_h.val);
        
    B_psi = Bilinear(M_h,C_h);
    B_psi = B_psi  + Integral(dOmega,(u_M_h_coef.val+v_M_h_coef.val)* psi_h.val * u_C_h.val);

    % set the minimum order of accuracy for the quad rule
    Quadrature_Order = 10;
    % define geometry representation - Domain, (default to piecewise linear)
    G1 = GeoElement(Omega);
    % define a set of matrices
    MATS = Matrices(Quadrature_Order,G1);

    % collect all of the matrices together
    MATS = MATS.Append_Matrix(B_phi);
    MATS = MATS.Append_Matrix(B_psi);
end