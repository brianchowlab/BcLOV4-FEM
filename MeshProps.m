function [props] = MeshProps(Mesh,shp_n)
    %Fetch faces. Note that this is a better way of doing it than using
    %built-in MATLAB subroutines, because it ensures that face winding is
    %correct (i.e. face normals point outwards).
    props.faces = Mesh.freeBoundary();
    props.unique_faces = unique(props.faces(:));
    props.nodes = Mesh.Points;
    props.surface_TR = triangulation(props.faces,props.nodes(props.unique_faces,:));
    props.surface_nodes = props.surface_TR.Points;
    props.surface_vert_normal = vertexNormal(props.surface_TR);
    props.nucleus_surface_nodes_idx = inShape(shp_n,props.surface_nodes + 0.1*props.surface_vert_normal);
    props.pm_surface_nodes = props.surface_nodes;
    props.pm_surface_nodes(props.nucleus_surface_nodes_idx,:) = [];
    props.nucleus_surface_nodes_idx = find(props.nucleus_surface_nodes_idx);
    props.pm_faces = props.faces;
    for i = 1:size(props.nucleus_surface_nodes_idx,1)
        props.pm_faces(any(props.pm_faces == props.nucleus_surface_nodes_idx(i),2),:) = [];
    end
    props.pm_unique_faces = unique(props.pm_faces(:));
    props.pm_unique_faces_renum = props.pm_unique_faces;
    for i = 1:size(props.pm_unique_faces,1)
       props.pm_unique_faces_renum(props.pm_unique_faces_renum == props.pm_unique_faces(i)) = i;
    end
end