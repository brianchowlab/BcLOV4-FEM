function [m_intrp] = InterpolateMembrane(I,desired_idx,TR,sol_M,param)
    scale=param.downsample;
    samp_res = param.scale_len;
    axial_resolution = param.axial_resolution;

    y = (0:scale:size(I,1)-1)*samp_res + scale*samp_res/2;
    x = (0:scale:size(I,2)-1)*samp_res + scale*samp_res/2;
    %z = -param.h:param.scale_len*scale:ceil(param.h/param.scale_len/scale)*param.scale_len*scale;
    z = -param.h:axial_resolution*scale:ceil(param.h/axial_resolution/scale)*axial_resolution*scale;

    %size(I)
    m_intrp = NaN(size(y,2),size(x,2),size(z,2),size(desired_idx,2),'single');
    % m_intrp = sparse(size(x,2),size(y,2)*size(z,2)*size(desired_idx,2)); 
    %m_intrp = ndSparse(m_intrp,[size(x,2),size(y,2),size(z,2),size(desired_idx,2)]);
    h = waitbar(0,'Interpolating membrane...');
    for k=1:size(m_intrp,4)
        waitbar(k/size(m_intrp,4));
        F = scatteredInterpolant(TR.Points(:,1),TR.Points(:,2),TR.Points(:,3),sol_M(:,desired_idx(k)));
        planes.n = [zeros(size(z,2),2),ones(size(z,2),1)];
        planes.r = [zeros(size(z,2),2),z'];
        polygons = mesh_xsections( TR.Points,TR.ConnectivityList, planes, 1e-3, 0 );

        P_m_i = {};
        N_m_i = {};
        for i=1:size(z,2)
            if length(polygons{i}) == 0
                continue
            end
            X = polygons{i}{1}(:,1);
            Y = polygons{i}{1}(:,2);
            X = [X;X(1)];
            Y = [Y;Y(1)];

            pt = interparc(size(polygons{i}{1}(:,1),1)*5,X,Y,'linear');

            %pt_sep = pdist(pt(1:2,:),'euclidean');
            pt = [pt,z(i)*ones(size(pt,1),1)];
            P_m_i{end+1} = pt;
            %plot(polygons{1}{1}(:,1),polygons{1}{1}(:,2),'o')
            %plot(pt(:,1),pt(:,2),'*')
            N = F(pt);
            N_m_i{end+1} = N;
            X = pt(:,1);
            Y = pt(:,2);


            %Convert interpolated values to pixels
            X = round((X-scale*samp_res/2)/samp_res/scale)*samp_res*scale + scale*samp_res/2;
            Y = round((Y-scale*samp_res/2)/samp_res/scale)*samp_res*scale + scale*samp_res/2;
            [~,idx_x] = ismembertol(X,x,1e-4);
            [~,idx_y] = ismembertol(Y,y,1e-4);
            

            
            [C,ia,ic] = unique([idx_x,idx_y],'rows');

            C = [C,i*ones(size(C,1),1),k*ones(size(C,1),1)];
            %Average function values when multiple points are in a pixel
            pixel_val = zeros(size(C,1),1);
            for j = 1:size(C,1)
                idx = find(ic == j);
                pixel_val(j) = mean(N(idx));
            end 
            %if k == 1 & i == floor(size(z,2)/2)
            %   plot(C(:,2),pixel_val,'o') 
            %end     
            
            C=sub2ind(size(m_intrp),C(:,2),C(:,1),C(:,3),C(:,4));
            m_intrp(C) = pixel_val;

        end
    end
    close(h);
end
