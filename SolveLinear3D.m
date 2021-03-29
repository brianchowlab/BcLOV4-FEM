function [Soln] = SolveLinear3D(param,solver_params)
    
    u_h = solver_params.u_h;
    CN = solver_params.CN;
    MN = solver_params.MN;
    idx = solver_params.idx;
    
    %Fetch initial FEM matrices
    solver_params.FEM_l = mex_ReactionDiffusion_Cell_l_assemble(solver_params.FEM,solver_params.p,solver_params.c,...
        [],solver_params.embed,solver_params.C_DoF,solver_params.M_DoF);
    
    [K_phi,K_psi,K_phi_phi,K_phi_psi,K_psi_phi,A_phi,A_psi]=FetchMatrices_linear(solver_params);
    
    k_on_p_store = param.k_on_p;

    %Assign some variables for easier access
    il_phi = K_phi.*sparse(param.il_c);
    il_psi = K_psi.*sparse(param.il_m);
    u_C_1 =  K_phi/param.dt+param.D*A_phi+sparse(param.k_off_p)*K_phi+param.k_on_l*param.S*K_phi_phi;
    u_C_2 = -param.k_on_p*il_phi;
    u_C_3 = -param.k_off_l*K_psi_phi;
    u_C_4 = sparse(CN,MN);
    v_C_1 = -param.k_off_p*K_phi;
    v_C_2 =  K_phi/param.dt+param.D*A_phi+param.k_on_p*il_phi+param.k_on_d*param.S*K_phi_phi;
    v_C_3 = sparse(CN,MN);
    v_C_4 = -param.k_off_d*K_psi_phi;
    u_M_1 = -param.k_on_l*param.S*K_phi_psi;
    u_M_2 = sparse(MN,CN);
    u_M_3 = K_psi/param.dt+param.D_m*A_psi+param.k_off_p*K_psi+param.k_off_l*K_psi;
    u_M_4 = -param.k_on_p*il_psi;
    v_M_1 = sparse(MN,CN);
    v_M_2 = -param.k_on_d*param.S*K_phi_psi;
    v_M_3 = -param.k_off_p*K_psi;
    v_M_4 = K_psi/param.dt+param.D_m*A_psi+param.k_on_p*il_psi+param.k_off_d*K_psi;

    disp('Decomposing matrices')
    if param.debug;start = tic;end
    %Step 1
    Step_1_LHS_lit = [u_C_1,u_C_2,u_C_3,u_C_4;v_C_1,v_C_2,v_C_3,v_C_4;u_M_1,u_M_2,u_M_3,u_M_4;v_M_1,v_M_2,v_M_3,v_M_4];
    Step_1_LHS_decomp_lit = decomposition(Step_1_LHS_lit);

    %Assign some variables for easier access
    u_C_2 = 0*il_phi;
    v_C_2 = K_phi/param.dt+param.D*A_phi+param.k_on_d*param.S*K_phi_phi;
    u_M_4 = 0*il_psi;
    v_M_4 = K_psi/param.dt+param.D_m*A_psi+param.k_off_d*K_psi;

    %Step 1
    Step_1_LHS_dark = [u_C_1,u_C_2,u_C_3,u_C_4;v_C_1,v_C_2,v_C_3,v_C_4;u_M_1,u_M_2,u_M_3,u_M_4;v_M_1,v_M_2,v_M_3,v_M_4];
    Step_1_LHS_decomp_dark = decomposition(Step_1_LHS_dark);
    if param.debug;t_i = toc(start);disp(['Initial decomposition took ',num2str(t_i),'s.']);end
    
    %Time discretisation via backward (implicit) Euler method.
    Soln = zeros(2*MN + 2*CN,1+param.num_steps/param.store_interval);
    Soln(:,1) = u_h;
    k_on_vec = [];
    for ii = 1:param.num_steps  
        disp(['Step ', num2str(ii), ' of ', num2str(param.num_steps)]);
        if param.debug;start = tic;end
        if mod(ii - 1,param.period/param.dt) < param.ex_duration/param.dt
            k_on_vec = [k_on_vec,1];
            param.k_on_p = k_on_p_store;
            Step_1_LHS_decomp = Step_1_LHS_decomp_lit;
        else
            k_on_vec = [k_on_vec,0];
            param.k_on_p = 0;
            Step_1_LHS_decomp = Step_1_LHS_decomp_dark;
        end

        solver_params.u_M = u_h(2*CN+1:2*CN+MN);
        solver_params.v_M = u_h(2*CN+MN+1:end);

        
        RHS = blkdiag(K_phi,K_phi,K_psi,K_psi)*u_h/sparse(param.dt);
        u_h = Step_1_LHS_decomp\RHS;
        
        if param.debug;t_i = toc(start);disp(['Sub-step 1 took ',num2str(t_i),'s.']);end

      
        if mod(ii,param.store_interval) == 0
            Soln(:,ii/param.store_interval+1) = u_h;
        end
    end
    param.k_on_p = k_on_p_store;
    plot(k_on_vec)
end