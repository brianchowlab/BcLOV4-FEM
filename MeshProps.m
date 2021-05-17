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
    
    %Centroids of elements
    tr = triangulation(Mesh.ConnectivityList,Mesh.Points);
    props.centroids = barycentricToCartesian(tr,[1:size(tr.ConnectivityList,1)]',repmat([1/4,1/4,1/4,1/4],[size(tr.ConnectivityList,1),1]));
    
    %Alphashape of centroids
    props.centroidShape = alphaShape(props.centroids,10,'HoleThreshold',1e6);
    
    %Indices of nuclear surface nodes - relative to surface nodes
    %Use alphaShape of nucleus to determine what surface nodes are on
    %nucleus. This can lead to occasional errors if the tolerance is chosen
    %poorly.
    if size(shp_n.Points,1) > 0
        %The face winding of the nuclear membrane can be confusing. The
        %vertex normals either need to be + or - depending on how they were
        %constructed
        props.nucleus_surface_nodes_idx = inShape(shp_n,props.surface_nodes + .5*props.surface_vert_normal);
        %props.nucleus_surface_nodes_idx = inShape(props.centroidShape,props.surface_nodes);
        %props.nucleus_surface_nodes_idx = [];
        props.nucleus_surface_nodes_idx = find(props.nucleus_surface_nodes_idx);
    else
        props.nucleus_surface_nodes_idx = [];
    end
        
    
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