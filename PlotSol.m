function [] = PlotSol(Soln,photo_on_scale,CN,MN,props,param)
    idx_m = unique(props.faces(:));
    idx_excited = find(photo_on_scale > 0.25);
    
    u_C = Soln(1:CN,:);
    v_C = Soln((CN+1):2*CN,:);
    u_M = Soln(2*CN+1:2*CN+MN,:);
    v_M = Soln(2*CN+MN+1:end,:);
    sol_C = u_C + v_C;
    sol_M = u_M + v_M;
    sol_all = sol_C;
    sol_all(idx_m,:) = sol_all(idx_m,:) + sol_M;
    sol_M_embed = zeros(size(sol_C));
    sol_M_embed(idx_m,:) = sol_M_embed(idx_m,:) + sol_M;

    
    figure
    h1 = trisurf(props.faces,props.nodes(:,1),props.nodes(:,2),props.nodes(:,3),sol_M_embed(:,end));
    colorbar
    colormap jet
    view([-1,1,0.5])


    tlist = (0:param.store_interval:param.num_steps)*param.dt;
    figure
    subplot(3,2,1)
    plot(tlist,mean(u_C))
    title('Cytoplasm Lit State')
    subplot(3,2,2)
    plot(tlist,mean(v_C))
    title('Cytoplasm Dark State')
    subplot(3,2,3)
    plot(tlist,mean(u_M))
    title('Membrane Lit State')
    subplot(3,2,4)
    plot(tlist,mean(v_M))
    title('Membrane Dark State')
    subplot(3,2,5)
    plot(tlist,mean(u_C)+mean(v_C))
    title('Mean Cytoplasm')
    subplot(3,2,6)
    plot(tlist,mean(u_M)+mean(v_M))
    title('Mean Membrane')

    figure
    idx_m_excited = intersect(idx_excited,idx_m);
    idx_m_excited_t = zeros(size(idx_m_excited));
    for i = 1:size(idx_m_excited)
        idx_m_excited_t(i) = find(idx_m_excited(i) == idx_m);
    end

    
    subplot(3,2,1)
    plot(tlist,mean(u_C(idx_excited,:)))
    title('Excitation ROI Cytoplasm Lit State')
    subplot(3,2,2)
    plot(tlist,mean(v_C(idx_excited,:)))
    title('Excitation ROI Cytoplasm Dark State')
    subplot(3,2,3)
    plot(tlist,mean(u_M(idx_m_excited_t,:),1))
    title('Excitation ROI Membrane Lit State')
    subplot(3,2,4)
    plot(tlist,mean(v_M(idx_m_excited_t,:),1))
    title('Excitation ROI Membrane Dark State')
    subplot(3,2,5)
    plot(tlist,mean(u_C(idx_excited,:))+mean(v_C(idx_excited,:)))
    title('Excitation ROI Mean Cytoplasm')
    subplot(3,2,6)
    plot(tlist,mean(u_M(idx_m_excited_t,:),1)+mean(v_M(idx_m_excited_t,:),1))
    title('Excitation ROI Mean Membrane')
end