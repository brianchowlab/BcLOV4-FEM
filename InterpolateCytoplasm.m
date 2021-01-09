function [c_intrp] =InterpolateCytoplasm(pde_C,desired_times,I,param)
    axial_resolution = param.scale_len;
    samp_res = param.scale_len;
    scale=1;

    z = (-7:axial_resolution*scale:7 + 1e-4)+1e-4;

    % Cytoplasm
    x = (0:scale:(size(I,1)-1))*samp_res + samp_res/2;
    y = (0:scale:(size(I,2)-1))*samp_res + samp_res/2;

    [X,Y,Z] = meshgrid(x,y,z);
    v = [X(:),Y(:),Z(:)];
    
    %Interpolate onto the desired grid
    c_intrp = single(interpolateSolution(pde_C,v',1:size(desired_times,2)));
    c_intrp = reshape(c_intrp,size(X,1),size(X,2),size(X,3),size(desired_times,2));
    c_intrp(isnan(c_intrp)) = 0;

    %This gives the concentration of BcLOV at each grid point. We are
    %intersted in the number of BcLOV molecules in each voxel.
    %volume_grid = (samp_res*scale).^2*axial_resolution;
    %idx_c = c_intrp(:) > 0;
    %idx_c = idx_c(1:size(idx_nan,1)/size(c_intrp,4));
    %idx_m = find(~idx_nan(1:size(idx_nan,1)/size(c_intrp,4)));

    
    %vol_int = CytoplasmVolume(v,volume_grid,idx_c,idx_nan,param);
end
