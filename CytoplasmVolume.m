function [vol_int] = CytoplasmVolume(I,interp_props,param)    
    axial_resolution = param.scale_len;
    samp_res = param.scale_len;
    scale=2;

    z = (-7:axial_resolution*scale:7 + 1e-4)+1e-4;

    % Cytoplasm
    x = (0:scale:(size(I,1)-1))*samp_res + samp_res/2;
    y = (0:scale:(size(I,2)-1))*samp_res + samp_res/2;

    [X,Y,Z] = meshgrid(x,y,z);
    v = [X(:),Y(:),Z(:)];
    
    volume_grid = (samp_res*scale).^2*axial_resolution;
    vol_int = zeros(size(v,1),1);
     
    vol_int(interp_props.idx_c) = volume_grid;
    for i = 1:size(interp_props.idx_m)
        if mod(i,1000) == 0
            i
        end
        p = v(idx_m(i),:);
        x = linspace(-samp_res/2 + p(1),samp_res/2 + p(1),2);
        y = linspace(-samp_res/2 + p(2),samp_res/2 + p(2),2);
        z = linspace(-samp_res/2 + p(3),samp_res/2 + p(3),2);
        [X,Y,Z] = meshgrid(x,y,z);
        v_s = [X(:),Y(:),Z(:)];
        fraction_in = nansum(pointLocation(DT,v_s)>0)/8;
        vol_int(idx_m(i)) = fraction_in;
    end
end
