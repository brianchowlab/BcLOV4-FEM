function x_verts = crosssections(model,m,plane,plot_flag)
    if ~exist('plot_flag','var')
        plot_flag = 0;
    end
    el = m.Elements;
    verts = m.Nodes;
    x_verts = {};
    for i=1:size(el,2)
        tet = verts(:,el(:,i));
        [int_num, pint] = plane_normal_tetrahedron_intersect(plane.r, plane.n,tet);
        if ( 0 < int_num )
            %hold on
            %plot3 ( pint(1,1:int_num), pint(2,1:int_num), pint(3,1:int_num), 'k.', 'MarkerSize', 20 );
            %color = [1,0,0];
            %patch ( pint(1,1:int_num), pint(2,1:int_num), pint(3,1:int_num), color );
            x_verts{end+1} = pint(:,1:int_num);
        end
    end
    if plot_flag
       pdeplot3D(model,'ColorMapData',ones(size(m.Nodes,2),1),'FaceAlpha',0.15)
       colorbar(gca ,'off' )
       for i =  1:size(x_verts,2)
          v = x_verts{i};
          hold on
          color = [1,0,0];
          patch (v(1,1:end), v(2,1:end), v(3,1:end),color);
       end
       camtarget(plane.r)
       view(plane.n)
    end
end
