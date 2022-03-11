function [] = PlotSol(Soln,photo_on_scale,CN,MN,props,param)
    idx_m = unique(props.pm_faces(:));
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
    plot(tlist,mean(u_C)/602.2)
    title('Cytoplasm Lit State')
    xlabel('Time (s)')
    ylabel('Concentration (uM)')
    subplot(3,2,2)
    plot(tlist,mean(v_C)/602.2)
    title('Cytoplasm Dark State')
    xlabel('Time (s)')
    ylabel('Concentration (uM)')
    subplot(3,2,3)
    plot(tlist,mean(u_M))
    title('Membrane Lit State')
    ylabel('Density (molecules/um^2)')
    xlabel('Time (s)')
    subplot(3,2,4)
    plot(tlist,mean(v_M))
    title('Membrane Dark State')
    ylabel('Density (molecules/um^2)')
    xlabel('Time (s)')
    subplot(3,2,5)
    plot(tlist,mean(u_C)/602.2+mean(v_C)/602.2)
    xlabel('Time (s)')
    title('Mean Cytoplasm')
    ylabel('Concentration (uM)')
    subplot(3,2,6)
    plot(tlist,mean(u_M)+mean(v_M))
    title('Mean Membrane')
    xlabel('Time (s)')
    ylabel('Density (molecules/um^2)')

    figure
    idx_m_excited = intersect(idx_excited,idx_m);
    idx_m_excited_t = zeros(size(idx_m_excited));
    for i = 1:size(idx_m_excited)
        idx_m_excited_t(i) = find(idx_m_excited(i) == idx_m);
    end

    
    subplot(3,2,1)
    plot(tlist,mean(u_C(idx_excited,:)/602.2))
    title('Excitation ROI Cytoplasm Lit State')
    xlabel('Time (s)')
    ylabel('Concentration (uM)')
    subplot(3,2,2)
    plot(tlist,mean(v_C(idx_excited,:)/602.2))
    title('Excitation ROI Cytoplasm Dark State')
    xlabel('Time (s)')
    ylabel('Concentration (uM)')
    subplot(3,2,3)
    plot(tlist,mean(u_M(idx_m_excited_t,:),1))
    title('Excitation ROI Membrane Lit State')
    xlabel('Time (s)')
    ylabel('Density (molecules/um^2)')
    subplot(3,2,4)
    plot(tlist,mean(v_M(idx_m_excited_t,:),1))
    title('Excitation ROI Membrane Dark State')
    xlabel('Time (s)')
    ylabel('Density (molecules/um^2)')
    subplot(3,2,5)
    plot(tlist,mean(u_C(idx_excited,:))/602.2+mean(v_C(idx_excited,:))/602.2)
    title('Excitation ROI Mean Cytoplasm')
    xlabel('Time (s)')
    ylabel('Concentration (uM)')
    subplot(3,2,6)
    plot(tlist,mean(u_M(idx_m_excited_t,:),1)+mean(v_M(idx_m_excited_t,:),1))
    title('Excitation ROI Mean Membrane')
    xlabel('Time (s)')
    ylabel('Density (molecules/um^2)')
end