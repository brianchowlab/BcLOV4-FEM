function [Soln] = SolveNonLinear2D(param,solver_params)
    
    u_h = solver_params.u_h;
    CN = solver_params.CN;
    MN = solver_params.MN;
    idx = solver_params.idx;
    
    %Fetch initial FEM matrices
    solver_params.FEM_l = mex_ReactionDiffusion_Cell_l_2D_assemble(solver_params.FEM,solver_params.p,solver_params.c,...
        [],solver_params.embed,solver_params.C_DoF,solver_params.M_DoF);
    solver_params.FEM_nl = mex_ReactionDiffusion_Cell_nl_2D_assemble(solver_params.FEM,solver_params.p,solver_params.c,[],...
            solver_params.embed,solver_params.C_DoF,solver_params.M_DoF,solver_params.u_M,solver_params.v_M);
        
    [K_phi,K_psi,K_phi_phi,K_phi_psi,K_psi_phi,A_phi,A_psi]=FetchMatrices_linear_2D(solver_params);
    [B_phi,B_psi]=FetchMatrices_nonlinear_2D(solver_params);
    
    k_on_p_store = param.k_on_p;

    %Assign some variables for easier access
    il_phi = K_phi.*sparse(param.il_c);
    il_psi = K_psi.*sparse(param.il_m);
    u_C_1 = @(theta)K_phi/(theta*param.dt)+param.D*A_phi+sparse(param.k_off_p)*K_phi+param.k_on_l*param.S*K_phi_phi;
    u_C_2 = -param.k_on_p*il_phi;
    u_C_3 = -param.k_off_l*K_psi_phi;
    u_C_4 = sparse(CN,MN);
    v_C_1 = -param.k_off_p*K_phi;
    v_C_2 = @(theta) K_phi/(theta*param.dt)+param.D*A_phi+param.k_on_p*il_phi+param.k_on_d*param.S*K_phi_phi;
    v_C_3 = sparse(CN,MN);
    v_C_4 = -param.k_off_d*K_psi_phi;
    u_M_1 = -param.k_on_l*param.S*K_phi_psi;
    u_M_2 = sparse(MN,CN);
    u_M_3 = @(theta) K_psi/(theta*param.dt)+param.D_m*A_psi+param.k_off_p*K_psi+param.k_off_l*K_psi;
    u_M_4 = -param.k_on_p*il_psi;
    v_M_1 = sparse(MN,CN);
    v_M_2 = -param.k_on_d*param.S*K_phi_psi;
    v_M_3 = -param.k_off_p*K_psi;
    v_M_4 = @(theta) K_psi/(theta*param.dt)+param.D_m*A_psi+param.k_on_p*il_psi+param.k_off_d*K_psi;

    u_C_nl_1 = param.k_on_l*B_phi;
    v_C_n1_2 = param.k_on_d*B_phi;
    u_M_nl_1 = -param.k_on_l*B_psi;
    v_M_nl_2 = -param.k_on_d*B_psi;

    disp('Decomposing matrices for linear steps')
    if param.debug;start = tic;end
    %Step 1
    Step_1_LHS_lit = [u_C_1(param.theta),u_C_2,u_C_3,u_C_4;v_C_1,v_C_2(param.theta),v_C_3,v_C_4;u_M_1,u_M_2,u_M_3(param.theta),u_M_4;v_M_1,v_M_2,v_M_3,v_M_4(param.theta)];
    Step_1_LHS_decomp_lit = decomposition(Step_1_LHS_lit);

    %Assign some variables for easier access
    u_C_2 = 0*il_phi;
    v_C_2 = @(theta) K_phi/(theta*param.dt)+param.D*A_phi+param.k_on_d*param.S*K_phi_phi;
    u_M_4 = 0*il_psi;
    v_M_4 = @(theta) K_psi/(theta*param.dt)+param.D_m*A_psi+param.k_off_d*K_psi;

    %Step 1
    Step_1_LHS_dark = [u_C_1(param.theta),u_C_2,u_C_3,u_C_4;v_C_1,v_C_2(param.theta),v_C_3,v_C_4;u_M_1,u_M_2,u_M_3(param.theta),u_M_4;v_M_1,v_M_2,v_M_3,v_M_4(param.theta)];
    Step_1_LHS_decomp_dark = decomposition(Step_1_LHS_dark);
    if param.debug;t_i = toc(start);disp(['Initial decomposition took ',num2str(t_i),'s.']);end
    %Time discretisation via the fractional-step theta method
    Soln = zeros(2*MN + 2*CN,1+param.num_steps/param.store_interval);
    Soln(:,1) = u_h;
    k_on_vec = [];
    for ii = 1:param.num_steps  
        %disp('-------------------------------------------------------');
        disp(['Step ', num2str(ii), ' of ', num2str(param.num_steps)]);
        if param.debug;start = tic;end
        if mod(ii - 1,param.period/param.dt) <param.ex_duration/param.dt
            k_on_vec = [k_on_vec,1];
            param.k_on_p = k_on_p_store;
            Step_1_LHS_decomp = Step_1_LHS_decomp_lit;
            u_C_2 = -param.k_on_p*il_phi;
            v_C_2 = @(theta) K_phi/(theta*param.dt)+param.D*A_phi+param.k_on_p*il_phi+param.k_on_d*param.S*K_phi_phi;
            u_M_4 = -param.k_on_p*il_psi;
            v_M_4 = @(theta) K_psi/(theta*param.dt)+param.D_m*A_psi+param.k_on_p*il_psi+param.k_off_d*K_psi;
        else
            k_on_vec = [k_on_vec,0];
            param.k_on_p = 0;
            Step_1_LHS_decomp = Step_1_LHS_decomp_dark;
            u_C_2 = 0*il_phi;
            v_C_2 = @(theta) K_phi/(theta*param.dt)+param.D*A_phi+param.k_on_d*param.S*K_phi_phi;
            u_M_4 = 0*il_psi;
            v_M_4 = @(theta) K_psi/(theta*param.dt)+param.D_m*A_psi+param.k_off_d*K_psi;
        end

        main = tic;
        solver_params.u_M = u_h(2*CN+1:2*CN+MN);
        solver_params.v_M = u_h(2*CN+MN+1:end);

        [B_phi,B_psi]=FetchMatrices_nonlinear_2D(solver_params);

        u_C_nl_1 = param.k_on_l*B_phi;
        v_C_n1_2 = param.k_on_d*B_phi;
        u_M_nl_1 = -param.k_on_l*B_psi;
        v_M_nl_2 = -param.k_on_d*B_psi;

        RHS = [u_C_nl_1*u_h(1:CN);v_C_n1_2*u_h(CN+1:2*CN);u_M_nl_1*u_h(1:CN);v_M_nl_2*u_h(CN+1:2*CN)] + blkdiag(K_phi,K_phi,K_psi,K_psi)*u_h/sparse(param.theta*param.dt);

        %RHS = blkdiag(K_phi,K_phi,K_psi,K_psi)*u_h/(param.theta*param.dt);
        u_h = Step_1_LHS_decomp\RHS;
        if param.debug;t_i = toc(start);disp(['Sub-step 1 took ',num2str(t_i),'s.']);end

        %Step 2
        %Initial guess
        if param.debug;start = tic;end
        x0 = u_h;
        err = 100;
        counter = 0;
        while err > param.tol
            counter = counter + 1;
            %Jacobian
            solver_params.u_M = u_h(2*CN+1:2*CN+MN);
            solver_params.v_M = u_h(2*CN+MN+1:end);
            [B_phi,B_psi]=FetchMatrices_nonlinear_2D(solver_params);
            u_C_nl_1 = param.k_on_l*B_phi;
            v_C_n1_2 = param.k_on_d*B_phi;
            u_M_nl_1 = -param.k_on_l*B_psi;
            v_M_nl_2 = -param.k_on_d*B_psi;
            u_C = u_h(1:CN);
            v_C = u_h(CN+1:2*CN);


            J1 = @(u_h) [K_phi-(1-2*param.theta)*param.dt*param.k_on_l*B_phi,sparse(CN,CN),-param.k_on_l*(1-2*param.theta*param.dt)*K_phi(:,idx).*u_h(1:CN),-param.k_on_l*(1-2*param.theta*param.dt)*K_phi(:,idx).*u_h(1:CN)];
            J2 = @(u_h) [sparse(CN,CN),K_phi-param.k_on_d*(1-2*param.theta)*param.dt*B_phi,-param.k_on_d*(1-2*param.theta*param.dt)*K_phi(:,idx).*u_h(CN+1:2*CN),-param.k_on_d*(1-2*param.theta*param.dt)*K_phi(:,idx).*u_h(CN+1:2*CN)];
            J3 = @(u_h) [param.k_on_l*(1-2*param.theta)*param.dt*B_psi,sparse(MN,CN),K_psi+param.k_on_l*(1-2*param.theta*param.dt)*K_psi.*u_C(idx),param.k_on_l*(1-2*param.theta*param.dt)*K_psi.*u_C(idx)];
            J4 = @(u_h) [sparse(MN,CN),param.k_on_d*(1-2*param.theta)*param.dt*B_psi,param.k_on_d*(1-2*param.theta*param.dt)*K_psi.*v_C(idx),K_psi+param.k_on_d*(1-2*param.theta*param.dt)*K_psi.*v_C(idx)];

            J = [J1(x0);J2(x0);J3(x0);J4(x0)];

            %Residual

            linear_part = @(u_h) (1-2*param.theta)*param.dt*-1*[u_C_1(Inf),u_C_2,u_C_3,u_C_4;v_C_1,v_C_2(Inf),v_C_3,v_C_4;u_M_1,u_M_2,u_M_3(Inf),u_M_4;v_M_1,v_M_2,v_M_3,v_M_4(Inf)]*u_h;
            nonlinear_part = @(u_h) (1-2*param.theta)*param.dt*[u_C_nl_1*u_h(1:CN);v_C_n1_2*u_h(CN+1:2*CN);u_M_nl_1*u_h(1:CN);v_M_nl_2*u_h(CN+1:2*CN)];
            R = @(x0) blkdiag(K_phi,K_phi,K_psi,K_psi) * (x0 - u_h) - nonlinear_part(x0) - linear_part(u_h);
            RHS = -R(x0);

            if counter > 10
                disp('Convergence issue, using ilu preconditioner. This will cause a slowdown. This likely indicates an excessively large time-step.')
                [L,U] = ilu(J,struct('type','ilutp','droptol',1e-6));
                dx = gmres(J,RHS,5,param.tol/(10^(counter-10)),10,L,U,x0);
            else
                pre = spdiags(spdiags(J,0),0,2*CN+2*MN,2*CN+2*MN);
                dx = gmres(J,RHS,5,param.tol,10,pre);
            end

            %dx = gmres(J,RHS,5,1e-6,10,L,U);
            %pre = spdiags(spdiags(J,0),0,2*CN+2*MN,2*CN+2*MN);
            %dx = gmres(J,RHS,5,param.tol,10,pre);
            %dec = decomposition(J);
            %dx = dec \ -R(x0);
            err = norm(dx);
            x0 = dx + x0;
        end

        u_h = x0;
        if param.debug;t_i = toc(start);disp(['Sub-step 2 took ',num2str(t_i),'s.']);end
        
        %Step 3
        if param.debug;start = tic;end
        solver_params.u_M = u_h(2*CN+1:2*CN+MN);
        solver_params.v_M = u_h(2*CN+MN+1:end);
        [B_phi,B_psi]=FetchMatrices_nonlinear_2D(solver_params);
        u_C_nl_1 = param.k_on_l*B_phi;
        v_C_n1_2 = param.k_on_d*B_phi;
        u_M_nl_1 = -param.k_on_l*B_psi;
        v_M_nl_2 = -param.k_on_d*B_psi;
        RHS = [u_C_nl_1*u_h(1:CN);v_C_n1_2*u_h(CN+1:2*CN);u_M_nl_1*u_h(1:CN);v_M_nl_2*u_h(CN+1:2*CN)] + blkdiag(K_phi,K_phi,K_psi,K_psi)*u_h/(param.theta*param.dt);
        u_h = Step_1_LHS_decomp\RHS;

        %if mod(ii - 1,period/dt) <= ex_duration/dt 
        %    u_h = dMAT_on\RHS;
        %    %Soln(:,ii+1) = gmres(MAT_on,RHS,10,1e-6,100,L_on,U_on,u_h);
        %else
        %    u_h = dMAT_off\RHS;
        %    %Soln(:,ii+1) = gmres(MAT_off,RHS,10,1e-6,100,L_off,U_off,u_h);
        %end
        if mod(ii,param.store_interval) == 0
            Soln(:,ii/param.store_interval+1) = u_h;
        end
        if param.debug;t_i = toc(start);disp(['Sub-step 3 took ',num2str(t_i),'s.']);end
        toc (main)
    end
    param.k_on_p = k_on_p_store;
    plot(k_on_vec)
end