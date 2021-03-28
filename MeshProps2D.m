function [props] = MeshProps2D(Mesh)

    tr_m = triangulation(Mesh.ConnectivityList,Mesh.Points);
    %Fetch edges.
    props.edges = freeBoundary(tr_m);
    
    %Mesh.Plot;
    
    %Indices of nodes located on plasma membrane or nucleus
    props.surface_node_idx = unique(props.edges(:));
    
    %All nodes
    props.nodes = Mesh.Points;
    
    %All elements
    props.elements = Mesh.ConnectivityList;
    
    %Surface nodes
    props.surface_nodes = props.nodes(props.surface_node_idx,:);
  
    %Centroids of elements
    props.centroids = barycentricToCartesian(tr_m,[1:size(tr_m.ConnectivityList,1)]',repmat([1/3,1/3,1/3],[size(tr_m.ConnectivityList,1),1]));
    
    %Alphashape of centroids
    props.centroidShape = alphaShape(props.centroids,'HoleThreshold',1e6);
    
    %Indices of nuclear surface nodes - relative to surface nodes
    %Use alphaShape of nucleus to determine what surface nodes are on
    %nucleus. This can lead to occasional errors if the tolerance is chosen
    %poorly.
    props.nucleus_surface_nodes_idx = inShape(props.centroidShape,props.surface_nodes(:,1),props.surface_nodes(:,2));
    props.nucleus_surface_nodes_idx = find(props.nucleus_surface_nodes_idx);
    
    %Plasma membrane surface nodes
    props.pm_surface_nodes = props.surface_nodes;
    props.pm_surface_nodes(props.nucleus_surface_nodes_idx,:) = [];
    
    %hold on 
    %plot(props.pm_surface_nodes(:,1),props.pm_surface_nodes(:,2),'o')
    
    %Indices of plasma membrane surface nodes
    props.pm_surface_node_idx = find(ismember(props.nodes,props.pm_surface_nodes,'rows'));
    
    %Plasma membrane edges
    props.pm_edges_idx = ~(ismember(props.edges(:,1),props.nucleus_surface_nodes_idx) | ismember(props.edges(:,2),props.nucleus_surface_nodes_idx));
    props.pm_edges = props.edges(props.pm_edges_idx,:);
    
    %Renumber plasma membrane edges according to plasma membrane nodes
    props.pm_edges_renum = props.pm_edges;
    for i = 1:size(props.pm_surface_node_idx,1)
       props.pm_edges_renum(props.pm_edges_renum == props.pm_surface_node_idx(i)) = i;
    end
end