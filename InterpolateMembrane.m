function [] = InterpolateMembrane(I,c_intrp,TR,sol_M,param)
    scale=1;
    samp_res = param.scale_len;
    axial_resolution = samp_res;

    x = (0:scale:size(I,1)-1)*samp_res + samp_res/2;
    y = (0:scale:size(I,2)-1)*samp_res + samp_res/2;
    z = -param.h:axial_resolution:param.h;

    m_intrp = NaN(size(c_intrp),'single');
    for k=1:size(c_intrp,4)
        F = scatteredInterpolant(TR.Points(:,1),TR.Points(:,2),TR.Points(:,3),sol_M(:,k));
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
            X = round((X-samp_res/2)/samp_res)*samp_res + samp_res/2;
            Y = round((Y-samp_res/2)/samp_res)*samp_res + samp_res/2;
            [~,idx_x] = ismember(X,x);
            [~,idx_y] = ismember(Y,y);
            [C,ia,ic] = unique([idx_x,idx_y],'rows');
            C = [C,i*ones(size(C,1),1),k*ones(size(C,1),1)];

            %Average function values when multiple points are in a pixel
            pixel_val = zeros(size(C,1),1);
            for j = 1:size(C,1)
                idx = find(ic == j);
                pixel_val(j) = mean(N(idx));
            end 
            C=sub2ind(size(m_intrp),C(:,2),C(:,1),C(:,3),C(:,4));
            m_intrp(C) = pixel_val;
        end
        close all
    end
end
