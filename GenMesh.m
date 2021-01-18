function [mesh_c,poly,shp_n] = GenMesh(contours,param)
    %Generates a 3D tetrahedral mesh of the cell by first computing surface
    %triangulations and alphaShapes for the cytoplasm and nucleus and using
    %that to seed the MATLAB mesh generation function from the PDE toolbox.

    %-------------Pre-compute polyshapes/centroids-------------%   
    %Polyshapes get rid of co-linear points present in the contours.
    %Illuminated region
    %pg_il = polyshape([contours.il]*param.scale_len);

    %Nucleus/cytoplasm
    pg = polyshape([contours.cytoplasm;contours.nucleus]*param.scale_len);
    pg_c = polyshape([contours.cytoplasm]*param.scale_len);
    pg_n = polyshape([contours.nucleus]*param.scale_len);
    pg_il = polyshape([contours.il]*param.scale_len);
    
    poly.il = pg_il;
    poly.cn = pg;
    poly.c = pg_c;
    poly.n = pg_n;
    
    %Extract centroid and vertices for cytoplasm/nucleus
    cent_n = nanmean(pg_n.Vertices);
    verts_n = [pg_n.Vertices,zeros(size(pg_n.Vertices,1),1)];
    verts_n = verts_n(1:param.downsample:end,:);

    cent_c = cent_n;%nanmedian(pg_c.Vertices);%
    verts_c = [pg_c.Vertices,zeros(size(pg_c.Vertices,1),1)];
    verts_c = verts_c(1:param.downsample:end,:);

    %-------------Compute surface triangulations for cytoplasm/nucleus-------------%
    %Starting from the 2D segemntation of the cytoplasm, compute a 3D
    %surface triangulation ny hemi-ellipsoid extrusion, using param.extrude
    %number of points for the hemi-ellipsoid extrusion.
    new_verts_c = zeros(size(verts_c,1)*(param.extrude-1),3);
    tri_c = [];
    points_per_z_c = zeros(size(verts_c,1),2,param.extrude-1);
    for i=1:size(verts_c,1)
        %Radii of ellipsoid
        r_a = norm(verts_c(i,1:2) - cent_c);
        r_b = param.h;

        %Compute points of ellipsoid
        if verts_c(i,1) < cent_c(1) & verts_c(i,2) < cent_c(2)
            t=fliplr(linspace(-pi/2,pi/2,param.extrude));
        elseif verts_c(i,1) > cent_c(1) & verts_c(i,2) < cent_c(2)
            t=fliplr(linspace(-pi/2,pi/2,param.extrude));
        else
            t=linspace(pi/2,3*pi/2,param.extrude);
        end

        %Location of ellipsoid points in 2D (XY Plane)
        xu=r_a*cos(t);
        yu=r_b*sin(t);

        %Project the 2D ellipsoid points into 3D
        x =  xu;
        y = zeros(size(xu));
        z = yu;
        zc = z;

        %Rotate the points 
        rot_angle = -atan((verts_c(i,1)-cent_c(1))/(verts_c(i,2)-cent_c(2))) - pi/2;
        xr = x*cos(rot_angle) - y*sin(rot_angle) + cent_c(1);
        yr = x*sin(rot_angle) + y*cos(rot_angle) + cent_c(2);
        points_per_z_c(i,:,1:end-1) = [xr(2:end-1);yr(2:end-1)];

        %Add triangles for top and bottom
        if i ~= size(verts_c,1)
            tri_c = [tri_c;size(verts_c,1)*(param.extrude-1) + 1,(param.extrude-1) * (i-1) + 1,(param.extrude-1) * i + 1];
        else
            tri_c = [tri_c;size(verts_c,1)*(param.extrude-1) + 1,(param.extrude-1) * (i-1) + 1,1];
        end

        %Compute connectivity
        if i ~= size(verts_c,1) & param.extrude>2
           %trimesh(tri_c,new_verts_c(:,1),new_verts_c(:,2),new_verts_c(:,3),'FaceAlpha',0.1)
           for j=1:(param.extrude-2)
               %First node of previous is (param.extrude-1) * (i-1) + 1
               %Last node of previous is (param.extrude-1) * i
               %First node of current is (param.extrude-1) * i + 1
               %Top node is size(verts_c,1)*(param.extrude-1) + 1
               %Bottom node is size(verts_c,1)*(param.extrude-1) + 2
               tri_c = [tri_c;(param.extrude-1) * (i-1) + j,(param.extrude-1) * (i-1) + j+1,(param.extrude-1) * i + j];
               tri_c = [tri_c;(param.extrude-1) * i + j,(param.extrude-1) * i + j+1,(param.extrude-1) * (i-1) + j+1];
           end
        elseif param.extrude>2
            for j=1:(param.extrude-2)
               %First node of previous is (param.extrude-1) * (i-1) + 1
               %Last node of previous is (param.extrude-1) * i
               %First node of current is (param.extrude-1) * i + 1
               %Top node is size(verts_c,1)*(param.extrude-1) + 1
               %Bottom node is size(verts_c,1)*(param.extrude-1) + 2
               tri_c = [tri_c;(param.extrude-1) * (i-1) + j,(param.extrude-1) * (i-1) + j+1, j];
               tri_c = [tri_c;j,j+1,(param.extrude-1) * (i-1) + j+1];
           end
        end
        if i ~= size(verts_c,1)
            tri_c = [tri_c;size(verts_c,1)*(param.extrude-1) + 2,(param.extrude-1) * i,(param.extrude-1) * (i+1)];
        else
            tri_c = [tri_c;size(verts_c,1)*(param.extrude-1) + 2,(param.extrude-1) * i,param.extrude-1];
        end

        %Order vertices
        start = (param.extrude - 1) * (i-1) + 1;
        finish = (param.extrude - 1) * i;

        new_verts_c(start:finish,:) = [xr(2:(param.extrude/2))',yr(2:(param.extrude/2))',z(2:(param.extrude/2))';...
        verts_c(i,:);xr((param.extrude/2+1):end-1)',yr((param.extrude/2+1):end-1)',z((param.extrude/2+1):end-1)'];
    end

    %Add top and bottom vertices
    new_verts_c = [new_verts_c;cent_c(1),cent_c(2),param.h;cent_c(1),cent_c(2),-param.h];
    points_per_z_c(:,:,end) = verts_c(:,1:2);
    zc = [zc(2:end-1),0];

    %Compute 3D nucleus
    new_verts_n = zeros(size(verts_n,1)*(param.extrude-1),3);
    tri_n = [];
    points_per_z_n = zeros(size(verts_n,1),2,param.extrude-1);
    for i=1:size(verts_n,1)
        %Radii of ellipsoid
        r_a = norm(verts_n(i,1:2) - cent_n);
        r_b = param.hn;

        %Compute points of ellipsoid
        if verts_n(i,1) < cent_n(1) & verts_n(i,2) < cent_n(2)
            t=fliplr(linspace(-pi/2,pi/2,param.extrude));
        elseif verts_n(i,1) > cent_n(1) & verts_n(i,2) < cent_n(2)
            t=fliplr(linspace(-pi/2,pi/2,param.extrude));
        else
            t=linspace(pi/2,3*pi/2,param.extrude);
        end

        %Ellipsoid nodes in 2D (XY plane)
        xu=r_a*cos(t);
        yu=r_b*sin(t);

        %Project 2D points into 3D
        x =  xu;
        y = zeros(size(xu));
        z = yu;
        zn = z;

        %Rotate points
        rot_angle = -atan((verts_n(i,1)-cent_n(1))/(verts_n(i,2)-cent_n(2))) - pi/2;
        xr = x*cos(rot_angle) - y*sin(rot_angle) + cent_n(1);
        yr = x*sin(rot_angle) + y*cos(rot_angle) + cent_n(2);
        points_per_z_n(i,:,1:end-1) = [xr(2:end-1);yr(2:end-1)];

        %Add triangles for top and bottom
        if i ~= size(verts_n,1)
            tri_n = [tri_n;size(verts_n,1)*(param.extrude-1) + 1,(param.extrude-1) * (i-1) + 1,(param.extrude-1) * i + 1];
        else
            tri_n = [tri_n;size(verts_n,1)*(param.extrude-1) + 1,(param.extrude-1) * (i-1) + 1,1];
        end

        %Compute connectivity
        if i ~= size(verts_n,1) & param.extrude>2
           %trimesh(tri_c,new_verts_c(:,1),new_verts_c(:,2),new_verts_c(:,3))
           for j=1:(param.extrude-2)
               %First node of previous is (param.extrude-1) * (i-1) + 1
               %Last node of previous is (param.extrude-1) * i
               %First node of current is (param.extrude-1) * i + 1
               %Top node is size(verts_c,1)*(param.extrude-1) + 1
               %Bottom node is size(verts_c,1)*(param.extrude-1) + 2
               tri_n = [tri_n;(param.extrude-1) * (i-1) + j,(param.extrude-1) * (i-1) + j+1,(param.extrude-1) * i + j];
               tri_n = [tri_n;(param.extrude-1) * i + j,(param.extrude-1) * i + j+1,(param.extrude-1) * (i-1) + j+1];
           end
        elseif param.extrude>2
            for j=1:(param.extrude-2)
               %First node of previous is (param.extrude-1) * (i-1) + 1
               %Last node of previous is (param.extrude-1) * i
               %First node of current is (param.extrude-1) * i + 1
               %Top node is size(verts_c,1)*(param.extrude-1) + 1
               %Bottom node is size(verts_c,1)*(param.extrude-1) + 2
               tri_n = [tri_n;(param.extrude-1) * (i-1) + j,(param.extrude-1) * (i-1) + j+1, j];
               tri_n = [tri_n;j,j+1,(param.extrude-1) * (i-1) + j+1];
           end
        end
        if i ~= size(verts_n,1)
            tri_n = [tri_n;size(verts_n,1)*(param.extrude-1) + 2,(param.extrude-1) * i,(param.extrude-1) * (i+1)];
        else
            tri_n = [tri_n;size(verts_n,1)*(param.extrude-1) + 2,(param.extrude-1) * i,param.extrude-1];
        end

        %Order vertices
        start = (param.extrude - 1) * (i-1) + 1;
        finish = (param.extrude - 1) * i;

        new_verts_n(start:finish,:) = [xr(2:(param.extrude/2))',yr(2:(param.extrude/2))',z(2:(param.extrude/2))';...
        verts_n(i,:);xr((param.extrude/2+1):end-1)',yr((param.extrude/2+1):end-1)',z((param.extrude/2+1):end-1)'];

    end

    %Add top and bottom nodes
    new_verts_n = [new_verts_n;cent_n(1),cent_n(2),param.hn;cent_n(1),cent_n(2),-param.hn];
    points_per_z_n(:,:,end) = verts_n(:,1:2);
    zn = [zn(2:end-1),0];

    %-------------Construct alphaShapes of cytoplasm and nucleus-------------%
    %First compute a low resolution (non-meshable) point cloud for the cytoplasm
    X = [];
    Y = [];
    Z = [];
    for i=1:size(points_per_z_c,3)
        x = min(new_verts_c(:,1)):0.5:max(new_verts_c(:,1));
        y = min(new_verts_c(:,2)):0.5:max(new_verts_c(:,2));
        [x,y] = meshgrid(x,y);
        x = x(:);
        y = y(:);
        %shp = polyshape(points_per_z_n(:,1,i),points_per_z_n(:,2,i));
        %in = inpolygon(x,y,shp.Vertices(:,1),shp.Vertices(:,2));
        %x = [x(~in);points_per_z_n(:,1,i)];
        %y = [y(~in);points_per_z_n(:,2,i)];
        shp = polyshape(points_per_z_c(:,1,i),points_per_z_c(:,2,i));
        in = inpolygon(x,y,shp.Vertices(:,1),shp.Vertices(:,2));
        x = [x(in);points_per_z_c(:,1,i)];
        y = [y(in);points_per_z_c(:,2,i)];
        X = [X;x];
        Y = [Y;y];
        Z = [Z;ones(size(x)) * zc(i)];
    end
    X = [X;cent_c(1);cent_c(1)];
    Y = [Y;cent_c(2);cent_c(2)];
    Z = [Z;param.h;-param.h];
    shp_c = alphaShape(X,Y,Z,2);
    %plot(shp_c,'FaceAlpha',0.5)
    %hold on
    %plot(pg)

    %Next, compute a low resolution (non-meshable) point cloud for the nucleus
    X = [];
    Y = [];
    Z = [];
    for i=1:size(points_per_z_c,3)
        x = min(new_verts_n(:,1)):0.5:max(new_verts_n(:,1));
        y = min(new_verts_n(:,2)):0.5:max(new_verts_n(:,2));
        [x,y] = meshgrid(x,y);
        x = x(:);
        y = y(:);
        shp = polyshape(points_per_z_n(:,1,i),points_per_z_n(:,2,i));
        in = inpolygon(x,y,shp.Vertices(:,1),shp.Vertices(:,2));
        x = [x(in);points_per_z_n(:,1,i)];
        y = [y(in);points_per_z_n(:,2,i)];
        X = [X;x];
        Y = [Y;y];
        Z = [Z;ones(size(x)) * zn(i)];
    end
    X = [X;cent_n(1);cent_n(1)];
    Y = [Y;cent_n(2);cent_n(2)];
    Z = [Z;param.hn;-param.hn];
    shp_n = alphaShape(X,Y,Z,1);
    %plot(shp_n,'FaceAlpha',0.5)
    %hold on
    %plot(pg)
    %

    %-------------Construct alphaShapes of cytoplasm with nuclear void-------------%
    %High resolution mesh
    res=0.2;
    x = min(new_verts_c(:,1)):0.5:max(new_verts_c(:,1));
    y = min(new_verts_c(:,2)):0.5:max(new_verts_c(:,2));
    z = -param.h:res:param.h;
    [x,y,z] = meshgrid(x,y,z);
    x = x(:);
    y = y(:);
    z = z(:);
    in = inShape(shp_c,x,y,z);
    x = x(in);
    y = y(in);
    z = z(in);
    [~,nodes] = boundaryFacets(shp_c);
    x = [x;nodes(:,1)];
    y = [y;nodes(:,2)];
    z = [z;nodes(:,3)];
    in = inShape(shp_n,x,y,z);
    x = x(~in);
    y = y(~in);
    z = z(~in);
    [~,nodes] = boundaryFacets(shp_n);
    x = [x;nodes(:,1)];
    y = [y;nodes(:,2)];
    z = [z;nodes(:,3)];

    shp_cn = alphaShape(x,y,z,param.alpha_radius);
    numRegions(shp_cn)
    plot(shp_cn,'FaceAlpha',0.5)
    figure
    %-------------Mesh cytoplasm/nuclear void using built in mesh generator-------------%
    [elements,nodes] = boundaryFacets(shp_cn);
    nodes = nodes';
    elements = elements';
    model = createpde();
    gm_c = geometryFromMesh(model,nodes,elements);
    %pdegplot(model,'CellLabels','on','FaceAlpha',0.5)
    mesh_c = generateMesh(model,'GeometricOrder','linear','Hmin',param.min_element_size,'Hmax',param.max_element_size);%mesh_c = generateMesh(model,'Hmin',0.5,'Hmax',30,'GeometricOrder','linear');
    if param.plot
        pdemesh(mesh_c,'FaceAlpha',0.1)
        h = gca;
        h = h.Children;
        h(2).Visible = 'Off';
        h(3).Visible = 'Off';
        h(4).Visible = 'Off';
        h(5).Visible = 'Off';
    end
end