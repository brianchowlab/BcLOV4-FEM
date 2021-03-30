function [photo_on_scale,idx_excited] = ExcitationROI(mesh_c,mask_il,poly,param)
    z = -param.h:param.scale_len:ceil(param.h/param.scale_len)*param.scale_len;
    v = poly.il.Vertices;
    c = [(max(v(:,1))+min(v(:,1)))/2,(max(v(:,2))+min(v(:,2)))/2];
    
    if param.excitation_type == 1
        [il_3D,coords_il] = ExcitationVolumeDMD(mask_il,z,param);
    else
        %[il_3D,coords_il] = ExcitationVolumeLaser(mask_il,z,param);
        photo_on_scale = ones(size(mesh_c.Nodes,2),1);
        idx_excited = 1:size(mesh_c.Nodes,2);
        return
    end

    
    temp = coords_il.X;
    coords_il.X = coords_il.Y + c(1);
    coords_il.Y = temp + c(2);

    v = mesh_c.Nodes';
    photo_on_scale = interp3(coords_il.Y,coords_il.X,coords_il.Z,il_3D,v(:,2),v(:,1),v(:,3));
    photo_on_scale(isnan(photo_on_scale)) = 0;

    idx_excited = find(photo_on_scale > 0.25);
    if isempty(idx_excited)
       disp('Warning: no excited nodes in ROI') 
    end
    %Plot

    if param.plot
        cut=fliplr(logspace(log10(max(max(max(il_3D)))/4),log10(max(max(max(il_3D)))),5));


        hold on
        a = gca;

        h = waitbar(0,'Please wait...');

        for k=1:numel(cut)
            isonormals(coords_il.Y,coords_il.X,coords_il.Z,il_3D,patch(isosurface(coords_il.X,coords_il.Y,coords_il.Z,il_3D,cut(k)),'EdgeColor','none','FaceAlpha',1/k,'FaceColor',[1 (k-1)/max(1,numel(cut)-0.99) 0],'Parent',a));
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