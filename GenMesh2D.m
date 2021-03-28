function [mesh_c,poly,tr_n] = GenMesh2D(contours,param)
    %Generates a 2D triangular mesh using the MATLAB mesh generation function from the PDE toolbox.

    %-------------Pre-compute polyshapes/centroids-------------%   
    %Polyshapes get rid of co-linear points present in the contours.

    %Nucleus/cytoplasm and illuminated region

    pg_c = polyshape([contours.cytoplasm]*param.scale_len);
    pg_n = polyshape([contours.nucleus]*param.scale_len);
    pg = subtract(pg_c,pg_n);
    pg_il = polyshape([contours.il]*param.scale_len);
    
    poly.il = pg_il;
    poly.cn = pg;
    poly.c = pg_c;
    poly.n = pg_n;
    
   
    %-------------Mesh cytoplasm/nuclear void using built in mesh generator-------------%
    tr = triangulation(pg);
    tr_n = triangulation(pg_n);
    nodes = tr.Points';
    elements = tr.ConnectivityList';
    model = createpde();
    gm_c = geometryFromMesh(model,nodes,elements);
    mesh_c = generateMesh(model,'GeometricOrder','linear','Hmin',param.min_element_size,'Hmax',param.max_element_size);%mesh_c = generateMesh(model,'Hmin',0.5,'Hmax',30,'GeometricOrder','linear');

    if param.plot
        pdemesh(mesh_c,'FaceAlpha',0.1)
        h = gca;
        h = h.Children;
        %h(2).Visible = 'Off';
        %h(3).Visible = 'Off';
        %h(4).Visible = 'Off';
        %h(5).Visible = 'Off';
    end
end