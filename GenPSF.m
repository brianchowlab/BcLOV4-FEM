function [PSF_3D] = GenPSF(psf_params,param)
    cd MicroscPSF
    PSF_3D = MicroscPSF(psf_params);
    cd ..
    s = sum(PSF_3D,'all');
    PSF_3D = PSF_3D / s;

    if param.plot
        [nx,ny,nz] = size(PSF_3D);

        cut=exp(-1:-1:-5)/s;
        [X,Y,Z]=meshgrid(-(nx/2-1):nx/2,-(ny/2-1):ny/2,-(nz/2-1):nz/2);
        X = X*param.scale_len;
        Y=Y*param.scale_len;
        Z=Z*param.scale_len;

        figure;
        a = gca;

        h = waitbar(0,'Please wait...');

        for k=1:numel(cut)
            isonormals(X,Y,Z,PSF_3D,patch(isosurface(X,Y,Z,PSF_3D,cut(k)),'EdgeColor','none','FaceAlpha',1/k,'FaceColor',[1 (k-1)/max(1,numel(cut)-0.99) 0],'Parent',a));
            waitbar(k / numel(cut))
        end
        close(h);
        view(35,45);
        axis('tight');
        lighting('gouraud');
        grid('on');
        camlight;
        camlight;
        set(gca,'XColor', 'none','YColor','none','ZColor','None')
        set(gca, 'xtick', [])
        set(gca, 'ytick', [])
        set(gca, 'ztick', [])
        set(gcf,'Color','w')
        %figure
        %surf(X(:,:,1),Y(:,:,1),PSF_3D(:,:,ceil(psf_params.size(3)/2)))    
    end
end