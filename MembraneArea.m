function [] = MembraneArea(x,y,z,TR,param)
    samp_res = param.scale_len;
    axial_resolution = samp_res;
    
    x_box = find(x > (min(TR.Points(:,1))-samp_res) & x < (max(TR.Points(:,1))+samp_res));
    y_box = find(y > (min((TR.Points(:,2))-samp_res) & y < (max((TR.Points(:,2))+samp_res))));
    voxel_mem_area = sparse(size(x,2),size(y,2)*size(z,2)); 
    voxel_mem_area = ndSparse(voxel_mem_area,[size(x,2),size(y,2),size(z,2)]);
    pixel_map = containers.Map;
    for ii = 1:size(x_box,2)
        i = x_box(ii);
        disp(ii/size(x_box,2) * 100)
        pts_vox = [min(y)-samp_res/2,min(z)-axial_resolution/2;...
                    max(y)+samp_res/2,min(z)-axial_resolution/2;...
                    max(y)+samp_res/2,max(z)+axial_resolution/2;...
                    min(y)-samp_res/2,max(z)+axial_resolution/2];
        a = [(x(i)+samp_res/2)*ones(4,1),pts_vox];
        b = [(x(i)-samp_res/2)*ones(4,1),pts_vox];
        c = [a;b];
        DT = delaunayTriangulation(c);
        TR_vox = triangulation(freeBoundary(DT),DT.Points);
        [~,intSurface] = SurfaceIntersection(TR,TR_vox,'debug',false);
        idx_vert = (TR.Points(:,1) > (x(i)-samp_res/2) &...
                    TR.Points(:,1) < (x(i)+samp_res/2));
        intSurface.vertices = [intSurface.vertices;TR.Points(idx_vert,:)];
        y_locs = intSurface.vertices(:,2);
        idx_y = sort(unique(knnsearch(y',y_locs)))';
        idx_y = min(idx_y):max(idx_y);
        if isempty(idx_y)
          continue
        end
        for kk = 1:size(idx_y,2)
            k = idx_y(kk);
            pts_vox = [y(k)-samp_res/2,min(z)-axial_resolution/2;...
                y(k)+samp_res/2,min(z)-axial_resolution/2;...
                y(k)+samp_res/2,max(z)+axial_resolution/2;...
                y(k)-samp_res/2,max(z)+axial_resolution/2];
            a = [(x(i)+samp_res/2)*ones(4,1),pts_vox];
            b = [(x(i)-samp_res/2)*ones(4,1),pts_vox];
            c = [a;b];
            DT = delaunayTriangulation(c);
            TR_vox = triangulation(freeBoundary(DT),DT.Points);
            intMatrix = SurfaceIntersection(TR,TR_vox,'debug',false);
            [i1,i2] = find(intMatrix);
            unique_faces = unique(i1);
            idx_pts_TR = TR.ConnectivityList(unique_faces,:);
            idx_z = [];
            for h=1:size(idx_pts_TR,1)
                pts_TR = TR.Points(idx_pts_TR(h,:),:);
                idx = sort(unique(knnsearch(z',pts_TR(:,3))))';
                idx = (min(idx)-1):(max(idx)+1);
                idx(idx<1 | idx>size(z,2)) = [];
                idx_z = [idx_z,idx];
            end
            idx_z = unique(idx_z);
            %idx_z = min(idx_z):max(idx_z);
            if isempty(idx_z)
                continue
            end
            %if isempty(find(idx_z==36))
            %    continue
            %end
            for jj=1:size(idx_z,2)%jj =1:size(idx_z,2)
                j = idx_z(jj);
                %j=36;
                pts_vox = [y(k)-samp_res/2,z(j)-axial_resolution/2;...
                    y(k)+samp_res/2,z(j)-axial_resolution/2;...
                    y(k)+samp_res/2,z(j)+axial_resolution/2;...
                    y(k)-samp_res/2,z(j)+axial_resolution/2];
                a = [(x(i)+samp_res/2)*ones(4,1),pts_vox];
                b = [(x(i)-samp_res/2)*ones(4,1),pts_vox];
                c = [a;b];
                DT = delaunayTriangulation(c);
                TR_vox = triangulation(freeBoundary(DT),DT.Points);
                [intMatrix, intSurface] = SurfaceIntersection(TR,TR_vox,'debug',false);
                if ~nnz(intMatrix)
                    continue
                end
                %If there is a vertex inside the voxel, we also need to add it
                %to the set
                idx_vert = (TR.Points(:,1) > (x(i)-samp_res/2) &...
                    TR.Points(:,1) < (x(i)+samp_res/2) &...
                    TR.Points(:,2) > (y(k)-samp_res/2) &...
                    TR.Points(:,2) < (y(k)+samp_res/2) &...
                    TR.Points(:,3) > (z(j)-axial_resolution/2) &...
                    TR.Points(:,3) < (z(j)+axial_resolution/2));
                intSurface.vertices = [intSurface.vertices;TR.Points(idx_vert,:)];
                [i1,i2] = find(intMatrix);
                unique_faces = unique(i1);
                idx_pts_TR = TR.ConnectivityList(unique_faces,:);
                pts_TR = TR.Points(idx_pts_TR,:);
                to_store = zeros(size(unique_faces,1),5);
                %t = triangulation(idx_pts_TR,TR.Points);
                %Find which points lie on which triangle of the plasma membrane
                area = 0;
                for l = 1:size(unique_faces,1)
                    m = unique_faces(l);
                    %Find coplanar points
                    idx_pts_TR_l = TR.ConnectivityList(m,:);
                    plane_pnt = TR.Points(idx_pts_TR_l(1),:);
                    fnorm = faceNormal(TR,m);
                    idx_plane = (intSurface.vertices - plane_pnt) * fnorm' < 1e-3;

                    %Throw out some edge cases where the voxel brushes the
                    %triangle (less than 3 points of contact). These are not a
                    %relavent contribution to area anyway.
                    if sum(idx_plane) < 3
                       continue 
                    end
                    %xy = intSurface.vertices(idx_plane,:);
                    %T = delaunayTriangulation(xy);
                    %B = freeBoundary(T);
                    %polyarea(xy(B(:,1),1),xy(B(:,1),2))
                    %Make the normal vector to the plane the new zaxis
                    yaxis = cross([1,0,0],fnorm);
                    yaxis = yaxis/norm(yaxis);
                    xaxis = cross(fnorm,yaxis);
                    xaxis = xaxis/norm(xaxis);
                    x_proj = (intSurface.vertices(idx_plane,:) - plane_pnt) * xaxis';
                    y_proj = (intSurface.vertices(idx_plane,:) - plane_pnt) * yaxis';
                    %Make sure all points are in the triangle;
                    x_t = (TR.Points(idx_pts_TR_l,:) - plane_pnt) * xaxis';
                    y_t = (TR.Points(idx_pts_TR_l,:) - plane_pnt) * yaxis';
                    poly_t = polyshape(x_t,y_t);
                    idx_in = isinterior(poly_t,x_proj,y_proj);
                    if sum(idx_in) < 3
                       to_store(l,:) = -1;
                       continue 
                    end
                    x_proj = x_proj(idx_in);
                    y_proj = y_proj(idx_in);
                    %Check for collinearity
                    xy = [x_proj,y_proj];
                    if rank(xy(2:end,:) - xy(1,:)) == 1
                        to_store(l,:) = -1;
                        continue
                    end

                    [~,area_part] = convhull(x_proj,y_proj);
                    area = area+area_part;
                    to_store(l,1) = m;
                    to_store(l,2) = area_part;
                    to_store(l,3:5) = mean(intSurface.vertices(idx_plane,:));
                end
                %area
                if isempty(to_store)
                    continue
                end
                pixel_map(num2str([i,k,j])) = to_store;
                voxel_mem_area(i,k,j) = area;
                %trimesh(TR);hold on;trimesh(TR_vox);hold on;plot3(pts_TR(:,1),pts_TR(:,2),pts_TR(:,3),'o');
                %hold on;plot3(intSurface.vertices(:,1),intSurface.vertices(:,2),intSurface.vertices(:,3),'o')
                %poly_vox = polyshape(pts_vox);
                %[polyout,shapeID,vertexID] = intersect(poly_1,poly_vox);
                %v1 = polyout.Vertices(shapeID == 0,:);
                %[polyout,shapeID,vertexID] = intersect(poly_2,poly_vox);
                %v2 = polyout.Vertices(shapeID == 0,:);
                %plot(poly_1);hold on;plot(poly_vox) 
            end
        end
    end
    save('1-2_c_area.mat','voxel_mem_area','pixel_map');
end
