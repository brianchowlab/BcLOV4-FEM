function [c_intrp] =InterpolateCytoplasm(pde_C,desired_idx,I,param)
    axial_resolution = param.scale_len;
    samp_res = param.scale_len;
    scale=param.downsample;

    z = -param.h:param.scale_len*scale:ceil(param.h/param.scale_len/scale)*param.scale_len*scale;

    % Cytoplasm
    x = (0:scale:(size(I,2)-1))*samp_res + scale*samp_res/2;
    y = (0:scale:(size(I,1)-1))*samp_res + scale*samp_res/2;

    [X,Y,Z] = meshgrid(x,y,z);
    v = [Y(:),X(:),Z(:)];
    
    %Interpolate onto the desired grid
    
    if length(desired_idx) < 2
        c_intrp = single(interpolateSolution(pde_C,v',desired_idx));
        c_intrp = reshape(c_intrp,size(X,1),size(X,2),size(X,3),size(desired_idx,2));
        c_intrp(isnan(c_intrp)) = 0;
    else
        c_intrp = single(interpolateSolution(pde_C,v',desired_idx));
        c_intrp = reshape(c_intrp,size(X,1),size(X,2),size(X,3),size(desired_idx,2));
        c_intrp(isnan(c_intrp)) = 0;
    end

    %This gives the concentration of BcLOV at each grid point. We are
    %intersted in the number of BcLOV molecules in each voxel.
    %volume_grid = (samp_res*scale).^2*axial_resolution;
    %idx_c = c_intrp(:) > 0;
    %idx_c = idx_c(1:size(idx_nan,1)/size(c_intrp,4));
    %idx_m = find(~idx_nan(1:size(idx_nan,1)/size(c_intrp,4)));

    
    %vol_int = CytoplasmVolume(v,volume_grid,idx_c,idx_nan,param);
end
