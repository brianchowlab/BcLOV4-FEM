mesh_c = pdem_C.Mesh;
e = boundaryFacets(mesh_c)';
t = triangulation(e,mesh_c.Nodes');
stlwrite(t,'optics_geom.stl','text')