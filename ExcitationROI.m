function [photo_on_scale,idx_excited] = ExcitationROI(mesh_c,mask_il,poly,param)
    z = -param.h:param.scale_len:param.h;
    v = poly.il.Vertices;
    c = mean(v);

    [il_3D,coords_il] = ExcitationVolumeDMD(mask_il,z,param);
    coords_il.X = coords_il.X + c(1);
    coords_il.Y = coords_il.Y + c(2);

    x0 = min(coords_il.X,[],'all');x1 = max(coords_il.X,[],'all');
    y0 = min(coords_il.Y,[],'all');y1 = max(coords_il.Y,[],'all');
    z0 = min(coords_il.Z,[],'all');z1 = max(coords_il.Z,[],'all');
    v = mesh_c.Nodes';
    idx_in = find(v(:,1) >= x0 & v(:,1) <= x1 & v(:,2) >= y0 & v(:,2) <= y1 & v(:,3) >= z0 & v(:,3) <= z1);
    photo_on_scale = zeros(size(v,1),1);
    s = [coords_il.X(:),coords_il.Y(:),coords_il.Z(:)];

    yu = unique(coords_il.Y);
    xu = unique(coords_il.X);
    for i=1:size(idx_in,1)
        z_idx = knnsearch(z',v(idx_in(i),3));
        y_idx = knnsearch(yu,v(idx_in(i),2));
        x_idx = knnsearch(xu,v(idx_in(i),1));
        photo_on_scale(idx_in(i)) = il_3D(y_idx,x_idx,z_idx);
    end
    idx_excited = find(photo_on_scale > 0.25);
    %Plot

    if param.plot
        cut=fliplr(logspace(log10(max(max(max(il_3D)))/4),log10(max(max(max(il_3D)))),5));


        hold on
        a = gca;

        h = waitbar(0,'Please wait...');

        for k=1:numel(cut)
            isonormals(coords_il.X,coords_il.Y,coords_il.Z,il_3D,patch(isosurface(coords_il.X,coords_il.Y,coords_il.Z,il_3D,cut(k)),'EdgeColor','none','FaceAlpha',1/k,'FaceColor',[1 (k-1)/max(1,numel(cut)-0.99) 0],'Parent',a));
            waitbar(k / numel(cut))
        end
        close(h);
        view(35,45);
        axis('tight');
        lighting('gouraud');
        grid('on');
        camlight;
        set(gca,'XColor', 'none','YColor','none','ZColor','None')
        set(gca, 'xtick', [])
        set(gca, 'ytick', [])
        set(gca, 'ztick', [])
        set(gcf,'Color','w')
    end
    clear x_il y_il z_il v
end