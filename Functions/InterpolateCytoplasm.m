function [c_intrp] =InterpolateCytoplasm(pde_C,desired_idx,I,param)
    axial_resolution = param.axial_resolution;
    samp_res = param.scale_len_x;
    scale=param.downsample;

    z = -param.h:axial_resolution*scale:ceil(param.h/axial_resolution/scale)*axial_resolution*scale;
    
    % Cytoplasm
    x = (0:scale:(size(I,1)-1))*samp_res + scale*samp_res/2;
    y = (0:scale:(size(I,2)-1))*samp_res + scale*samp_res/2;

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
end
