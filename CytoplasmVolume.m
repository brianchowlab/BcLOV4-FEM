function [vol_int] = CytoplasmVolume(v,volume_grid,idx_c,idx_nan,param)    
    samp_res = param.scale_len;
    vol_int = zeros(size(v,1),1);
     
    vol_int(idx_c) = volume_grid;
    for i = 1:size(idx_m)
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
