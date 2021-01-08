%% Set parameters
fid = fopen('params.txt');
params = textscan(fid,'%s %f','delimiter',' ','MultipleDelimsAsOne',1,'CommentStyle','%');
fclose(fid);
param = cell2struct(num2cell(params{2}),string(params{1})');

%Unit conversion
unit_scaling_k_on_l_and_d = 1e15/6.02214076e23;%Convert from M-1 s-1 to um^3 s-1 molecules-1
param.k_on_d = param.k_on_d*unit_scaling_k_on_l_and_d;
param.k_on_l = param.k_on_l*unit_scaling_k_on_l_and_d;

%Excitation Parameters
excitation_frequency = 299792458 * 1e9 / param.excitation_wavelength;%1/s
absorption_cross = param.extinction_coeff * 2303 / 6.0221409e23;%cm-2
photon_excitation_flux_density = param.power_density / excitation_frequency / 6.626e-34;%photons cm-2 s-1
param.k_on_p = param.quantum_yield_signaling_state * absorption_cross * photon_excitation_flux_density;%rate constant of signaling state formation, s-1
param.period = param.excitation_duration * 100/param.duty_cycle;

%% Load image
[contours,mask_il,I] = LoadImages();

%% Make Mesh
[mesh_c,poly,shp_n] = GenMesh(contours,param);

% Compute key values (indices, volumes, etc.)
props.elements = mesh_c.Elements;
props.nodes = mesh_c.Nodes;
props.vol_els = zeros(size(props.elements,2),1);
for i = 1:size(props.elements,2)
    props.vol_els(i) = volume(mesh_c,i);
end
props.vol = zeros(size(props.nodes,2),1);

%Find volumes
for i = 1:size(mesh_c.Nodes,2)
    [loc,attached_el] = ind2sub(size(props.elements),find(props.elements == i));
    props.vol(i) = volume(mesh_c,attached_el)/4;
end

z = -param.h:param.scale_len:param.h;
v = poly.il.Vertices;
c = mean(v);

%% Determine illuminated region
[il_3D,coords_il] = ExcitationVolumeDMD(mask_il,z,param);
coords_il.X = coords_il.X + c(1);
coords_il.Y = coords_il.Y + c(2);

x0 = min(coords_il.X,[],'all');x1 = max(coords_il.X,[],'all');
y0 = min(coords_il.Y,[],'all');y1 = max(coords_il.Y,[],'all');
z0 = min(coords_il.Z,[],'all');z1 = max(coords_il.Z,[],'all');
v = mesh_c.Nodes';
idx_in = find(v(:,1) >= x0 & v(:,1) <= x1 & v(:,2) >= y0 & v(:,2) <= y1 & v(:,3) >= z0 & v(:,3) <= z1);
photo_on_scale = zeros(size(v,1),1);
s = [coords_il.X(:),coords_il.Y(:),coords_il.Z(:)];

yu = unique(coords_il.Y);
xu = unique(coords_il.X);
for i=1:size(idx_in,1)
    z_idx = knnsearch(z',v(idx_in(i),3));
    y_idx = knnsearch(yu,v(idx_in(i),2));
    x_idx = knnsearch(xu,v(idx_in(i),1));
    photo_on_scale(idx_in(i)) = il_3D(y_idx,x_idx,z_idx);
end
idx_excited = find(photo_on_scale > 0.25);
%Plot

cut=fliplr(logspace(log10(max(max(max(il_3D)))/4),log10(max(max(max(il_3D)))),5));


figure;
a = gca;

h = waitbar(0,'Please wait...');

for k=1:numel(cut)
    isonormals(coords_il.X,coords_il.Y,coords_il.Z,il_3D,patch(isosurface(coords_il.X,coords_il.Y,coords_il.Z,il_3D,cut(k)),'EdgeColor','none','FaceAlpha',1/k,'FaceColor',[1 (k-1)/max(1,numel(cut)-0.99) 0],'Parent',a));
    waitbar(k / numel(cut))
end
close(h);
view(35,45);
axis('tight');
lighting('gouraud');
grid('on');
camlight;
set(gca,'XColor', 'none','YColor','none','ZColor','None')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(gca, 'ztick', [])
set(gcf,'Color','w')
clear x_il y_il z_il v
%% Build FEM matrices with FELICITY and solve.
cd FELICITY;FELICITY_paths;cd ..

%Port 3D cytoplasm mesh into FELICITY
Mesh = MeshTetrahedron(props.elements',props.nodes','Omega');


%Fetch faces. Note that this is a better way of doing it than using
%built-in MATLAB subroutines, because it ensures that face winding is
%correct (i.e. face normals point outwards).
Faces = Mesh.freeBoundary();
Faces_unique = unique(Faces(:));
TR = triangulation(Faces,Mesh.Points(Faces_unique,:));
P = TR.Points;
V = vertexNormal(TR);
 
Pt = P + 0.1*V;
V_n = inShape(shp_n,Pt);
P(V_n,:) = [];
idx_n = find(V_n);
for i = 1:size(idx_n,1)
    Faces(any(Faces == idx_n(i),2),:) = [];
end
Faces_unique = unique(Faces(:));
Faces_renum = Faces;
for i = 1:size(Faces_unique,1)
    Faces_renum(Faces_renum == Faces_unique(i)) = i;
end

[~,idx] = ismember(P,Mesh.Points,'rows');
Mesh = Mesh.Append_Subdomain('2D','dOmega',Faces);

%Calculate which nodes are illuminated
param.il_c = photo_on_scale;%inShape(shp_l,Mesh.Points);
%Mapping from membrane nodes to cytoplasmic nodes
param.il_m = photo_on_scale(idx);

%Characteristic length
vol_t = volume(mesh_c);
a = P(Faces_renum(:, 2), :) - P(Faces_renum(:, 1), :);
b = P(Faces_renum(:, 3), :) - P(Faces_renum(:, 1), :);
c = cross(a, b, 2);
area_all = 1/2 * sqrt(sum(c.^2, 2));
area = 1/2 * sum(sqrt(sum(c.^2, 2)));
l_char = vol_t/area;

%Calculate area per node on membrane
TR = triangulation(Faces_renum,Mesh.Points(Faces_unique,:));
vert2face = vertexAttachments(TR);
nodal_areas = cellfun(@(x) sum(area_all(x))/3,vert2face);

% Define FE spaces
C_RefElem_3D = ReferenceFiniteElement(lagrange_deg1_dim3());
C_Space = FiniteElementSpace('C_h',C_RefElem_3D,Mesh,'Omega');
C_Space = C_Space.Set_DoFmap(Mesh,uint32(Mesh.ConnectivityList));

M_RefElem_2D = ReferenceFiniteElement(lagrange_deg1_dim2());
M_Space = FiniteElementSpace('M_h',M_RefElem_2D,Mesh,'dOmega');
M_Space = M_Space.Set_DoFmap(Mesh,uint32(Faces_renum));

Domain_Names = {'Omega'; 'dOmega'};
Omega_Embed = Mesh.Generate_Subdomain_Embedding_Data(Domain_Names);

%FEM = mex_ReactionDiffusion_Cell_assemble([],nodes,uint32(elements),[],...
%    Omega_Embed,C_Space.DoFmap,M_Space.DoFmap);
CN = C_Space.num_dof;
MN = M_Space.num_dof;

%Parameters for matrix construction
solver_params.FEM = [];
solver_params.FEM_l = [];
solver_params.FEM_nl = [];
solver_params.p = Mesh.Points;
solver_params.c = uint32(Mesh.ConnectivityList);
solver_params.embed = Omega_Embed;
solver_params.C_DoF = C_Space.DoFmap;
solver_params.M_DoF = M_Space.DoFmap;
solver_params.CN = CN;
solver_params.MN = MN;
solver_params.idx = idx;

u_h = zeros(2*CN + 2*MN,1);

%Dark, cytosolic state initial condition (assuming no protein in lit
%state initially).
u_h((CN+1):2*CN) = param.conc*ones(CN,1);

u_h(2*CN+MN+1:end,:) = param.conc/(param.k_off_d/param.k_on_d+param.conc) * param.S * ones(MN,1);%k_on_d * 3000/k_off_d*ones(MN,1);
%u_h(2*CN+1:2*CN+MN,:) = conc/(k_off_d/k_on_d+conc) * S * ones(MN,1);%k_on_d * 3000/k_off_d*ones(MN,1);

solver_params.u_h = u_h;
solver_params.u_M = u_h(2*CN+1:2*CN+MN);
solver_params.v_M = u_h(2*CN+MN+1:end);

Soln = SolveNonLinear(param,solver_params);

idx_m = unique(Faces(:));
u_C = Soln(1:CN,:);
v_C = Soln((CN+1):2*CN,:);
u_M = Soln(2*CN+1:2*CN+MN,:);
v_M = Soln(2*CN+MN+1:end,:);
clear Soln A_phi A_psi B_phi B_psi FEM_l FEM_nl J K_phi K_phi_phi K_phi_psi K_psi K_psi_phi
clear Step_1_LHS_dark Step_1_LHS_lit Step_1_LHS_decomp Step_1_LHS_decomp_dark Step_1_LHS_decomp_lit
sol_C = u_C + v_C;
sol_M = u_M + v_M;
sol_all = sol_C;
sol_all(idx_m,:) = sol_all(idx_m,:) + sol_M;
sol_M_embed = zeros(size(sol_C));
desired_times = 1:1/(param.dt*param.store_interval):size(sol_all,2);

sol_M_embed(idx_m,:) = sol_M_embed(idx_m,:) + sol_M;
h1 = trisurf(Faces,Mesh.Points(:,1),Mesh.Points(:,2),Mesh.Points(:,3),sol_M_embed(:,end));
colorbar
colormap jet
view([-1,1,0.5])
sol_M = sol_M(:,desired_times);

tlist = (0:size(sol_all,2)-1)*param.dt*param.store_interval;
figure
subplot(3,2,1)
plot(tlist,mean(u_C))
title('Cytoplasm Lit State')
subplot(3,2,2)
plot(tlist,mean(v_C))
title('Cytoplasm Dark State')
subplot(3,2,3)
plot(tlist,mean(u_M))
title('Membrane Lit State')
subplot(3,2,4)
plot(tlist,mean(v_M))
title('Membrane Dark State')
subplot(3,2,5)
plot(tlist,mean(u_C)+mean(v_C))
title('Mean Cytoplasm')
subplot(3,2,6)
plot(tlist,mean(u_M)+mean(v_M))
title('Mean Membrane')

figure
idx_m_excited = intersect(idx_excited,idx_m);
idx_m_excited_t = zeros(size(idx_m_excited));
for i = 1:size(idx_m_excited)
    idx_m_excited_t(i) = find(idx_m_excited(i) == idx_m);
end
plot(tlist,mean(u_M(idx_m_excited_t,:))+mean(v_M(idx_m_excited_t,:)))



pdem_C = createpde(1);
gm_C_f = geometryFromMesh(pdem_C,Mesh.Points',Mesh.ConnectivityList');
pde_C=createPDEResults(pdem_C,sol_C(:,desired_times),tlist(desired_times),'time-dependent');
save('data.mat','v_M','u_M','v_C','u_C','idx_m','tlist','desired_times','pdem_C','gm_C_f');
clear v_M u_M v_C u_C
clear sol_C sol_all sol_M_embed
clear v_C_1 v_C_2 v_C_3 v_C_4 v_C_nl_2 v_M_1 v_M_2 v_M_3 v_M_4 v_M_nl_2 v_C_nl_2
clear u_M_1 u_M_2 u_M_3 u_M_4 u_M_nl_1 u_C_1 u_C_2 u_C_3 u_C_4
clear shp_n shp_cn shp_c s pre il_phi il_psi il_M il_3D
clear Z params mask
%% Interpolate 
[c_intrp] = InterpolateCytoplasm(pde_C,desired_times,I,param);

% Compute volume of cytoplasm at each grid point
vol_int = zeros(size(v,1),1);
idx_c = c_intrp(:) > 0;
idx_c = idx_c(1:size(idx_nan,1)/size(c_intrp,4));
vol_int(idx_c) = volume_grid;
idx_m = find(~idx_nan(1:size(idx_nan,1)/size(c_intrp,4)));
for i = 1:size(idx_m)
    if mod(i,1000) == 0
        i
    end
    p = v(idx_m(i),:);
    x = linspace(-samp_res/2 + p(1),samp_res/2 + p(1),2);
    y = linspace(-samp_res/2 + p(2),samp_res/2 + p(2),2);
    z = linspace(-samp_res/2 + p(3),samp_res/2 + p(3),2);
    [X,Y,Z] = meshgrid(x,y,z);
    v_s = [X(:),Y(:),Z(:)];
    fraction_in = nansum(pointLocation(DT,v_s)>0)/8;
    vol_int(idx_m(i)) = fraction_in;
end




%% OLD

m_intrp = sparse(size(x,2),size(y,2)*size(z,2)*size(desired_times,2)); 
m_intrp = ndSparse(m_intrp,[size(x,2),size(y,2),size(z,2),size(desired_times,2)]);
m_mat = pixel_map.values;
m_mat = cellfun(@(x) x(find(x(:,2) > 0),:),m_mat,'UniformOutput',false);
keys_keep = cellfun(@(x) ~isempty(x),m_mat);
m_mat = m_mat(keys_keep);
m_keys = cell2mat(cellfun(@str2num,pixel_map.keys,'UniformOutput',false)');
m_keys = m_keys(keys_keep,:);
%cellfun(@(x) x(:,2),m_mat,'UniformOutput',false)
P = pdem_C.Mesh.Nodes(:,idx_m)';
for i=1:size(desired_times,2)
    i
    F = scatteredInterpolant(P(:,1),P(:,2),P(:,3),sol_M(:,desired_times(i)));
    m_val = cell2mat(cellfun(@(x) x(:,2)' * F(x(:,3:end)),m_mat,'UniformOutput',false));
    idx = sub2ind(size(m_intrp),m_keys(:,1),m_keys(:,2),m_keys(:,3),i*ones(size(m_keys,1),1));
    %Mem value is product of value at centroid times area
    m_intrp(idx) = m_val;
end
c_intrp = c_intrp + permute(m_intrp,[2,1,3,4]);

%%
%Plasma membrane

x = (0:scale:size(I,1)-1)*samp_res + samp_res/2;
y = (0:scale:size(I,2)-1)*samp_res + samp_res/2;
z = -7:axial_resolution:7;

m_intrp = NaN(size(c_intrp),'single');
for k=1:size(c_intrp,4)
    F = scatteredInterpolant(P(:,1),P(:,2),P(:,3),sol_M(:,k));
    planes.n = [zeros(size(z,2),2),ones(size(z,2),1)];
    planes.r = [zeros(size(z,2),2),z'];
    polygons = mesh_xsections( P,TR.ConnectivityList, planes, 1e-3, 0 );

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

%Don't include cytoplasmic contribution to membrane fluorescence?
idx_nan = isnan(m_intrp(:));
m_intrp(c_intrp(:) > 0 & idx_nan) = 0;

%% Compute area of membrane 

planes.n = [zeros(size(y,2)+1,1),ones(size(y,2)+1,1),zeros(size(y,2)+1,1)];
planes.r = [zeros(size(y,2)+1,1),[y(1)-samp_res/2,y + samp_res/2]',zeros(size(y,2)+1,1)];
polygons_y = mesh_xsections( P,TR.ConnectivityList, planes, 1e-3, 0 );


planes.n = [ones(size(x,2)+1,1),zeros(size(x,2)+1,1),zeros(size(x,2)+1,1)];
planes.r = [[x(1)-samp_res/2,x + samp_res/2]',zeros(size(x,2)+1,1),zeros(size(x,2)+1,1)];
polygons_x = mesh_xsections( P,TR.ConnectivityList, planes, 1e-3, 0 );

planes.n = [zeros(size(z,2)+1,2),ones(size(z,2)+1,1)];
planes.r = [zeros(size(z,2)+1,2),[z(1)-axial_resolution/2,z + axial_resolution/2]'];
polygons_z = mesh_xsections( P,TR.ConnectivityList, planes, 1e-3, 0 );

close all
%%
[X,Y,Z] = meshgrid(x,y,z);
v = [X(:),Y(:),Z(:)];
v = repmat(v,1,1,8);
count = 0;
for i = -1:2:1
    for j = -1:2:1
        for k = -1:2:1
            count = count + 1;
            v(:,:,count) = v(:,:,8) + [i*samp_res/2,j*samp_res/2,k*axial_resolution/2];
        end
    end
end

vox = squeeze(v(1,:,:))';
DT = delaunayTriangulation(vox);
TR_vox = triangulation(freeBoundary(DT),DT.Points);
%ConnectivityList is the same for each voxel
vox_connect = TR_vox.ConnectivityList;
%%
a = v(1:4000,:,:);
a=permute(a,[3,1,2]);
a = reshape(a,[32000,3]);
a(1:8,1)
plot3(a(1:8,1),a(1:8,2),a(1:8,3),'-')
all_faces = repmat(vox_connect,4000,1);
all_faces = all_faces + 8*(reshape(repmat(1:4000',12,1),[],1)-1);
%%
b = 0:2383;
pts = mod(b,16)<=7;
b=b(pts);
surface2.vertices = a(b+1,:);
surface2.faces = all_faces(1:12*(sum(pts)/8),:);
[intMatrix,intSurface] = SurfaceIntersection(TR,surface2,'debug',false);
%%
FV.faces = TR.ConnectivityList;
FV.vertices = TR.Points;
%d = point2trimesh(FV,'QueryPoints',[X(:),Y(:),Z(:)]);
%%
for i =1:1000
    tic
    surf3.vertices = squeeze(v(1,:,:))';
    surf3.faces = vox_connect;
    [intMatrix,intSurface] = SurfaceIntersection(TR,surf3,'debug',false);
    toc
end

%%
FV.vertices = TR.Points;
FV.faces = TR.ConnectivityList;
IN = inpolyhedron(FV,[x-samp_res/2,x(end)+samp_res/2],[y-samp_res/2,y(end)+samp_res/2]+1e-6,[z-axial_resolution/2,z(end)+axial_resolution/2],'TOL',0,'FLIPNORMALS',0);
voxel_sums = movsum(movsum(movsum(IN,2,1),2,2),2,3);
voxel_sums = voxel_sums(2:end,2:end,2:end);
voxel_sums = voxel_sums > 0 & voxel_sums < 8;
d = v(voxel_sums(:),:);
plot3(d(:,1),d(:,2),d(:,3),'o')
a = delaunayTriangulation(d);
b = triangulation(freeBoundary(a),d);
figure
trimesh(b)
%%
idx = find(voxel_sums);
[ai,bi,ci] = ind2sub(size(voxel_sums),idx);
%touching = v(voxel_sums(:),:);
%i
%i+1
%i+size(y,2)+2
%i+size(y,2)+3
%i+(size(y,2)+1)*(size(x,2)+1)+1
%i+(size(y,2)+1)*(size(x,2)+1)+2
%i+(size(y,2)+1)*(size(x,2)+1)+size(y,2) + 2
%i+(size(y,2)+1)*(size(x,2)+1)+size(y,2) + 3

voxel_mem_area = sparse(size(x,2),size(y,2)*size(z,2)); 
voxel_mem_area = ndSparse(voxel_mem_area,[size(x,2),size(y,2),size(z,2)]);
pixel_map = containers.Map;
for ii = 1:size(ai,1)
    if mod(ii,1000) == 0
        disp(ii/size(ai,1)*100)
    end
    k = ai(ii);
    i = bi(ii);
    j = ci(ii);
    %if j ~= 18
    %    continue
    %end
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
        %Some voxels just brush the surface of the mesh, but do not truly
        %intersect.
        %[i,k,j]
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
            %area
        if isempty(to_store)
            continue
        end
    end
    pixel_map(num2str([i,k,j])) = to_store;
    voxel_mem_area(i,k,j) = area;
end
%%
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

%%
test_m = zeros(422,298);
keys = pixel_map.keys;
for i=1:size(keys,2)
    vals = pixel_map(keys{i});
    f = vals(:,1);
    areas = vals(:,2);
    pts = vals(:,3:5);
    idx = str2num(keys{i});
    test_m(idx(1),idx(2)) = sum(areas);
end

%% Approximate PSF


clear params
solver_params.size = [128,128,64];%[32,32,21]
solver_params.NA = NA;
solver_params.lambda = 610e-9;
solver_params.M = 63;
solver_params.ti0 = 100e-6;
solver_params.resLateral = samp_res * 1e-6;%scale_len * 1e-6;
solver_params.resAxial = axial_resolution *1e-6;%scale_len * 1e-6;
solver_params.pZ = 0;%ceil(params.size(3)/2) * axial_resolution*1e-6;
solver_params.oversampling = 2;
PSF_3D = MicroscPSF(solver_params);
s = sum(PSF_3D,'all')
PSF_3D = PSF_3D / s;

[nx,ny,nz] = size(PSF_3D);

cut=exp(-1:-1:-5)/s;
[X,Y,Z]=meshgrid(-(nx/2-1):nx/2,-(ny/2-1):ny/2,-(nz/2-1):nz/2);
X = X*samp_res;
Y=Y*samp_res;
Z=Z*axial_resolution;

figure;
a = gca;

h = waitbar(0,'Please wait...');

for k=1:numel(cut)
    isonormals(X,Y,Z,PSF_3D,patch(isosurface(X,Y,Z,PSF_3D,cut(k)),'EdgeColor','none','FaceAlpha',1/k,'FaceColor',[1 (k-1)/max(1,numel(cut)-0.99) 0],'Parent',a));
    waitbar(k / numel(cut))
end
close(h);
view(35,45);
axis('tight');
lighting('gouraud');
grid('on');
camlight;
figure
surf(X(:,:,1),Y(:,:,1),PSF_3D(:,:,ceil(solver_params.size(3)/2)))

%% Convolve with PSF
%c_intrp_blurred = zeros(size(c_intrp));
c_intrp_blurred = zeros(size(c_intrp),'single');
%c_intrp(isnan(c_intrp(:))) = 0;
%m_intrp(isnan(m_intrp(:))) = 0;
for i=1:4:size(desired_times,2)%[1,size(desired_times,2)]%
    i
    %c_intrp_blurred(:,:,:,i) = convn(c_intrp(:,:,:,i),PSF_3D,'same');
    c_intrp_blurred(:,:,:,i) = convn(c_intrp(:,:,:,i),single(PSF_3D),'same');
end
%% Line Profile
slice = ceil(size(z,2)/2);
t=21;
[B,L] = bwboundaries(m_intrp(:,:,slice,1),'noholes');
ind = sub2ind(size(m_intrp),B{1}(:,1),B{1}(:,2),slice * ones(size(B{1}(:,1))),t * ones(size(B{1}(:,1))));
prof = m_intrp(ind);
plot(prof)
%%

m_intrp_s = uint8(m_intrp_blurred / max(m_intrp_blurred(:)) * 255);
m_intrp_s = m_intrp_s - min(m_intrp_s(m_intrp_s>0));
%c_intrp_s(isnan(c_intrp_s)) = 0;
%m_intrp_s(m_intrp_s < 0) = 0;
m_intrp_s = m_intrp_s * (255 / max(m_intrp_s(:)));
%%
writerObj = VideoWriter('20s_05x_dark_Kd.avi','Indexed AVI');
writerObj.FrameRate = 30;
writerObj.Colormap = jet;
open(writerObj);
for t=1:size(m_intrp_blurred,4)
  t
  frame = im2frame(imadjust(m_intrp_s(:,:,36,t),[0,1],[0,1],2),jet);
  writeVideo(writerObj, frame);
end
close(writerObj);

%%
%caxis([min(m_intrp(m_intrp>0))-100,max(m_intrp(m_intrp>0))])
slice = ceil(size(z,2)/2);
t=101;
imagesc(c_intrp(:,:,:,t)')
hold on
visboundaries(imresize(mask_il,1/scale))
colormap(jet)
mask = zeros(400,600,'logical'); %<==You can do logical zeros here
position = [300,200];
degree = 45;
L = 50;
x1 = position(1)+ (L*cosd(degree)) ; 
y1 = position(2)+ (L*sind(degree)) ;
set(gca, 'DataAspectRatioMode', 'auto')  %<== Need this to prevent unpredictable behavior, BUT you loose the aspect ratio of the original image!
set(gca, 'unit', 'normalize')
set(gca, 'position', [0 0 1 1]);
%%
%caxis([min(min(uintrp(:,:,slice))),max(max(uintrp(:,:,slice,t)))])
%caxis([0,30])
figure
imagesc(c_intrp_blurred(:,:,slice,t)')
hold on
visboundaries(mask_il)
mask = zeros(400,600,'logical'); %<==You can do logical zeros here
position = [300,200];
degree = 45;
L = 50;
x1 = position(1)+ (L*cosd(degree)) ; 
y1 = position(2)+ (L*sind(degree)) ;
set(gca, 'DataAspectRatioMode', 'auto')  %<== Need this to prevent unpredictable behavior, BUT you loose the aspect ratio of the original image!
set(gca, 'unit', 'normalize')
set(gca, 'position', [0 0 1 1]);
%caxis([0,150*40])
%caxis([min(min(uintrp_blurred(:,:,slice))),max(max(uintrp_blurred(:,:,slice,t)))])
colormap(gray)
figure
imagesc(I)
colormap(gray)
mask = zeros(400,600,'logical'); %<==You can do logical zeros here
position = [300,200];
degree = 45;
L = 50;
x1 = position(1)+ (L*cosd(degree)) ; 
y1 = position(2)+ (L*sind(degree)) ;
set(gca, 'DataAspectRatioMode', 'auto')  %<== Need this to prevent unpredictable behavior, BUT you loose the aspect ratio of the original image!
set(gca, 'unit', 'normalize')
set(gca, 'position', [0 0 1 1]);

%% Background removal
mem_fl = zeros(size(c_intrp,4),1);
se = strel('disk',2);
c = se.Neighborhood;
for i = 1:size(c_intrp,4)
    img = imtophat(m_intrp(:,:,slice,i)',se);

    if mod(i,10) == 0
        figure
        imagesc(img)
    end
    mem_fl(i) = sum(sum(img));
end
plot(tlist(desired_times),(mem_fl-min(mem_fl))/max(mem_fl-min(mem_fl)))
hold on
plot(tlist,(mean(sol_M)-min(mean(sol_M)))/max(mean(sol_M)-min(mean(sol_M))))
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Fluorescence (a.u.)')
%% Initialize video
myVideo = VideoWriter('myVideoFile.avi'); %open video file
myVideo.FrameRate = 30;  %can adjust this, 5 - 10 works well for me
open(myVideo)
% Plot in a loop and grab frames
for i=1:1:size(sol_all,2)
    h1 = trisurf(Faces,Mesh.Points(:,1),Mesh.Points(:,2),Mesh.Points(:,3),sol_all(:,i));
    colorbar
    caxis([0,65])
    colormap jet
    %view(2)
    view([-1,1,0.5])
    pause(0.01) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    clf
end
close(myVideo)

%%
myVideo = VideoWriter('myVideoFile.avi'); %open video file
myVideo.FrameRate = 5;  %can adjust this, 5 - 10 works well for me
open(myVideo)
% Plot in a loop and grab frames
for i=1:1:size(sol_all,2)
    %pdeplot(model,'XYData',uintrp(:,i),'XYStyle','flat','ColorMap','parula')
    imagesc(c_intrp_blurred(:,:,ceil(size(c_intrp_blurred,3)/2),i))
    colormap gray
    mask = zeros(400,600,'logical'); %<==You can do logical zeros here
    position = [300,200];
    degree = 45;
    L = 50;
    x1 = position(1)+ (L*cosd(degree)) ; 
    y1 = position(2)+ (L*sind(degree)) ;
    set(gca, 'DataAspectRatioMode', 'auto')  %<== Need this to prevent unpredictable behavior, BUT you loose the aspect ratio of the original image!
    set(gca, 'unit', 'normalize')
    set(gca, 'position', [0 0 1 1]);
    caxis([0,0.9*max(max(c_intrp_blurred(:,:,ceil(size(c_intrp_blurred,3)/2),end)))])
    pause(0.01) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    clf
end
close(myVideo)
