function [] = StepwiseInterpolate(folder_name,pde_C,sol_M,props,param,conversion,I)
    mkdir(['./',folder_name]);
    mkdir(['./',folder_name,'_no_PSF']);
    filename = param.PSF;
    param.downsample = 1;
    %param.scale_z = 0.2;
    %param.PSF_axial_ratio = param.scale_z/param.scale_len;
    tstack = Tiff(filename,'r');

    [i,j] = size(tstack.read());
    K = length(imfinfo(filename));
    data = zeros(i,j,K);
    data(:,:,1)  = tstack.read();
    for n = 2:K
        tstack.nextDirectory()
        data(:,:,n) = tstack.read();
    end

    PSF_3D = data(1:end,1:end,1:end);
    s = sum(PSF_3D,'all');
    PSF_3D = PSF_3D / s;

    [nx,ny,nz] = size(PSF_3D);
    cut=exp(-1:-1:-3)/s;


    [X,Y,Z]=meshgrid(-floor(nx/2):nx/2,-floor(ny/2):ny/2,-floor(nz/2):nz/2);
    X = X*param.scale_len;
    Y=Y*param.scale_len;
    Z=Z*param.axial_resolution;%0.2 for WF
    %z_interp = min(Z(:)):param.scale_len:max(Z(:));
    %[Xi,Yi,Zi]=meshgrid(-floor(nx/2):nx/2,-floor(ny/2):ny/2,z_interp);
    %PSF_3D = interp3(X,Y,Z,PSF_3D,Xi,Yi,Zi,'nearest');

    figure;
    a = gca;
    h = waitbar(0,'Please wait...');
    for k=1:numel(cut)
        isonormals(X,Y,Z,PSF_3D,patch(isosurface(X,Y,Z,PSF_3D,cut(k)),'EdgeColor','none','FaceAlpha',1/k,'FaceColor',[1 (k-1)/max(1,numel(cut)-0.99) 0],'Parent',a));
        waitbar(k / numel(cut))
    end
    close(h);
    view(35,45);
    axis('equal');
    lighting('gouraud');
    grid('on');
    camlight;
    set(gca,'XColor', 'none','YColor','none','ZColor','None')
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
    set(gca, 'ztick', [])
    axis equal


    %parfor (i=1:param.interpolation_interval:size(sol_M,2),2)
    iter_vec = [1,2,3,4,5,6,26,51,101,201];
    parfor (j=1:10,3)
        i = iter_vec(j);
        pdem_C = createpde(1);
        gm_C_f = geometryFromMesh(pdem_C,props.nodes',props.elements');
        if i == size(sol_M,2)
            pde_C_i=createPDEResults(pdem_C,pde_C.NodalSolution(:,[i,1]),1:2,'time-dependent');
        else
            pde_C_i=createPDEResults(pdem_C,pde_C.NodalSolution(:,i:i+1),1:2,'time-dependent');
        end
        
        c_intrp = InterpolateCytoplasm(pde_C_i,1,I,param);

        m_intrp = InterpolateMembrane(I,1,props.pm_TR,sol_M(:,i),param);
        m_intrp(isnan(m_intrp)) = 0;


        c_intrp = c_intrp/conversion * param.conc_ratio + param.offset;
        m_intrp = m_intrp * 1.2;
        c_intrp(m_intrp ~= 0) = 0.5*c_intrp(m_intrp ~= 0) + m_intrp(m_intrp~=0);
        %size(c_intrp)
        c_intrp_blurred = ConvolvePSF(c_intrp,PSF_3D);
        %c_intrp_blurred = c_intrp;
        folder_name_no_PSF = [folder_name,'_no_PSF'];
        imwrite(uint16(squeeze(c_intrp_blurred(:,:,ceil(size(c_intrp,3)/2)))),['./',folder_name,'/',num2str(i),'.tif']);
        imwrite(uint16(squeeze(c_intrp(:,:,ceil(size(c_intrp,3)/2)))),['./',folder_name_no_PSF,'/',num2str(i),'.tif']);
    end
end