function [photo_on_scale,idx_excited] = PlotSol(Soln)
    h1 = trisurf(Faces,Mesh.Points(:,1),Mesh.Points(:,2),Mesh.Points(:,3),sol_M_embed(:,end));
    colorbar
    colormap jet
    view([-1,1,0.5])


    tlist = (0:size(sol_all,2)-1)*param.dt*param.store_interval;
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
    plot(tlist,mean(u_M(idx_m_excited_t,:))+mean(v_M(idx_m_excited_t,:)))
end