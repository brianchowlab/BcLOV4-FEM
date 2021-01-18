function [props] = MeshProps(Mesh,shp_n)
    %Fetch faces. Note that this is a better way of doing it than using
    %built-in MATLAB subroutines, because it ensures that face winding is
    %correct (i.e. face normals point outwards).
    props.faces = Mesh.freeBoundary();
    
    %Indices of nodes located on plasma membrane or nucleus
    props.surface_node_idx = unique(props.faces(:));
    
    %All nodes
    props.nodes = Mesh.Points;
    
    %All elements
    props.elements = Mesh.ConnectivityList;
    
    %Surface triangulation
    props.surface_TR = triangulation(props.faces,props.nodes);
    
    %Surface nodes
    props.surface_nodes = props.nodes(props.surface_node_idx,:);
    
    %Vertex normals of surfaces
    props.surface_vert_normal = vertexNormal(props.surface_TR);
    props.surface_vert_normal = props.surface_vert_normal(props.surface_node_idx,:);
    
    %Indices of nuclear surface nodes - relative to surface nodes
    props.nucleus_surface_nodes_idx = inShape(shp_n,props.surface_nodes + 0.1*props.surface_vert_normal);%1257
    props.nucleus_surface_nodes_idx = find(props.nucleus_surface_nodes_idx);
    
    %Plasma membrane surface nodes
    props.pm_surface_nodes = props.surface_nodes;
    props.pm_surface_nodes(props.nucleus_surface_nodes_idx,:) = [];
    
    %Plasma membrane faces
    props.pm_faces = props.faces;
    temp = props.surface_node_idx(props.nucleus_surface_nodes_idx);%relative to original node numbering
    for i = 1:size(props.nucleus_surface_nodes_idx,1)
        props.pm_faces(any(props.pm_faces == temp(i),2),:) = [];
    end
    
    %Indices of plasma membrane surface nodes
    props.pm_surface_node_idx = unique(props.pm_faces(:));
    
    
    props.pm_faces_renum = props.pm_faces;
    for i = 1:size(props.pm_surface_node_idx,1)
       props.pm_faces_renum(props.pm_faces_renum == props.pm_surface_node_idx(i)) = i;
    end
    props.pm_TR = triangulation(props.pm_faces_renum,props.pm_surface_nodes);
end