%% Set parameters
%filename = 'params_cell_1_10s_1dc';
filename = 'params_optics';

fid = fopen([filename,'.txt']);
ims_and_rois = textscan(fid,'%s %s','delimiter',' ','MultipleDelimsAsOne',1,'CommentStyle','%');
fclose(fid);
param = cell2struct(num2cell(ims_and_rois{2}),ims_and_rois{1}');

fn = fieldnames(param);
for k=1:numel(fn)
    if(isnan(str2double((param.(fn{k})))))
        param.(fn{k}) = param.(fn{k}){1};
    else
        param.(fn{k}) = str2double(param.(fn{k}){1});
    end
end


% param.k_off_p = 1/18.5;%0.2632;
% param.k_off_d = 0.0225;
% param.k_on_d = 1130;
% param.offset=104;
% param.dt = 1e-1;
% param.num_steps = 10;
% param.store_interval = 10;
% param.tol = 1e-3;
% param.scale_len = 0.1;
% param.excitation_type = 1;
% param.conc_ratio = 505.52;
% pram.offset = 125.35;

% param.k_off_p = 1/18.5;%0.2632;
% param.k_off_d = 0.0225;
% param.k_on_d = 1130;
% param.offset=104;
param.dt = 2.5e-2;
param.num_steps = 2400;
param.store_interval = 4;
param.tol = 1e-5;

%Unit conversion
unit_scaling_k_on_l_and_d = 1e15/6.02214076e23;%Convert from M-1 s-1 to um^3 s-1 molecules-1
param.k_on_d = param.k_on_d*unit_scaling_k_on_l_and_d;
param.k_on_l = param.k_on_l*unit_scaling_k_on_l_and_d;

%Excitation Parameters
excitation_frequency = 299792458 * 1e9 / param.excitation_wavelength;%1/s
absorption_cross = param.extinction_coeff * 2303 / 6.0221409e23;%cm-2
photon_excitation_flux_density = param.power_density / excitation_frequency / 6.626e-34;%photons cm-2 s-1
param.k_on_p = param.quantum_yield_signaling_state * absorption_cross * photon_excitation_flux_density;%rate constant of signaling state formation, s-1
param.ex_duration = param.period * param.duty_cycle/100;
%Conversion from uM to molecules/um^3
conversion = 1e-6 * 6.022e23 * 1e-15;

param.min_element_size = 0.25;
param.max_element_size = 0.75;


%% Load image
[contours,mask_il,I] = LoadImages(param);

%% Make Mesh
[mesh_c,poly,shp_n] = GenMesh_Optics(contours,param);

% Determine illuminated region
z = 0:0.1:5;
v = poly.il.Vertices;
c = [(max(v(:,1))+min(v(:,1)))/2,(max(v(:,2))+min(v(:,2)))/2];
[il_3D,coords_il] = ExcitationVolumeTIRF(mask_il,z,param);

temp = coords_il.X;
coords_il.X = coords_il.Y + c(1);
coords_il.Y = temp + c(2);

v = mesh_c.Nodes';
photo_on_scale = interp3(coords_il.Y,coords_il.X,coords_il.Z,il_3D,v(:,2),v(:,1),v(:,3));
photo_on_scale(isnan(photo_on_scale)) = 0;

idx_excited = find(photo_on_scale > 0.25);
if isempty(idx_excited)
   disp('Warning: no excited nodes in ROI') 
end
%Plot
%%
if param.plot
    cut=fliplr(logspace(log10(max(max(max(il_3D)))/4),log10(max(max(max(il_3D)))),5));


    hold on
    a = gca;

    h = waitbar(0,'Please wait...');

    for k=1:numel(cut)
        isonormals(coords_il.Y,coords_il.X,coords_il.Z,il_3D,patch(isosurface(coords_il.X,coords_il.Y,coords_il.Z,il_3D,cut(k)),'EdgeColor','none','FaceAlpha',1/k,'FaceColor',[1 (k-1)/max(1,numel(cut)-0.99) 0],'Parent',a));
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
end
clear x_il y_il z_il v temp coords_il il_3D
%%
figure
pdemesh(mesh_c,'FaceAlpha',0.1,'EdgeColor','none','FaceColor',[0.5,0.5,0.5])
h = gca;
h = h.Children;
h(2).Visible = 'Off';
h(3).Visible = 'Off';
view([1,0,0])
cut = 0.1;
a = gca;
isonormals(coords_il.Y,coords_il.X,coords_il.Z,il_3D,patch(isosurface(coords_il.X,coords_il.Y,coords_il.Z,il_3D,cut),'EdgeColor','none','FaceAlpha',1,'FaceColor','b','Parent',a));
ylim([27.5,33.5])
%#648FFF - tirf
%#FE6100 - 1P
%#DC267F - 2p
%%


plot(mesh_c.Nodes(1,idx_excited)',mesh_c.Nodes(2,idx_excited)','ro')
hold on
%plot(mesh_c.Nodes(1,idx_excited)',mesh_c.Nodes(2,idx_excited)','ko')

plot(poly.il)
%% Build FEM matrices with FELICITY and solve.
%cd FELICITY;FELICITY_paths;cd ..;

addpath ./FELICITY
FELICITY_paths

%Port 3D cytoplasm mesh into FELICITY
Mesh = MeshTetrahedron(mesh_c.Elements',mesh_c.Nodes','Omega');

props = MeshProps(Mesh,shp_n);
plot(props.pm_surface_nodes(abs(props.pm_surface_nodes(:,3))<1,1),props.pm_surface_nodes(abs(props.pm_surface_nodes(:,3))<1,2),'o')
[~,idx] = ismember(props.pm_surface_nodes,props.nodes,'rows');
m = find(abs(Mesh.Points(:,3) - param.h)<0.25);
m = intersect(m,idx);
m = Mesh.Points(m,:);
plot(m(:,1),m(:,2),'o')
Mesh = Mesh.Append_Subdomain('2D','dOmega',props.pm_faces);

%Calculate which nodes are illuminated
param.il_c = photo_on_scale;%inShape(shp_l,Mesh.Points);
%Mapping from membrane nodes to cytoplasmic nodes
param.il_m = photo_on_scale(idx);

% Define FE spaces
C_RefElem_3D = ReferenceFiniteElement(lagrange_deg1_dim3());
C_Space = FiniteElementSpace('C_h',C_RefElem_3D,Mesh,'Omega');
C_Space = C_Space.Set_DoFmap(Mesh,uint32(Mesh.ConnectivityList));

M_RefElem_2D = ReferenceFiniteElement(lagrange_deg1_dim2());
M_Space = FiniteElementSpace('M_h',M_RefElem_2D,Mesh,'dOmega');
M_Space = M_Space.Set_DoFmap(Mesh,uint32(props.pm_faces_renum));

Domain_Names = {'Omega'; 'dOmega'};
Omega_Embed = Mesh.Generate_Subdomain_Embedding_Data(Domain_Names);

CN = C_Space.num_dof;
MN = M_Space.num_dof;

%Parameters for matrix construction
solver_params.FEM = [];
solver_params.FEM_l = [];
solver_params.FEM_nl = [];
solver_params.p = props.nodes;
solver_params.c = uint32(props.elements);
solver_params.embed = Omega_Embed;
solver_params.C_DoF = C_Space.DoFmap;
solver_params.M_DoF = M_Space.DoFmap;
solver_params.CN = CN;
solver_params.MN = MN;
solver_params.idx = idx;

mask_cyto_only = poly2mask(contours.cytoplasm(:,1),contours.cytoplasm(:,2),size(I,1),size(I,2))-poly2mask(contours.nucleus(:,1),contours.nucleus(:,2),size(I,1),size(I,2));
u_h = zeros(2*CN + 2*MN,1);

%Dark, cytosolic state initial condition (assuming no protein in lit
%state initially).
%Convert from fluorescence to concentration
%param.conc = (mean(I(logical(mask_cyto_only(:))))-param.offset)/param.conc_ratio;
%param.conc = 1.1552;
%param.conc = 0.5385;

%For loading from file
% data_file = split(filename,'.txt');
% data_file = data_file{1};
% data_file = split(data_file,'params_c');
% data_file = data_file{2};
% 
% data_file = ['Cyto_C',data_file,'.csv'];
% 
% data = readmatrix(data_file);
% data = data(1,2);
% param.conc = (data-param.offset)/param.conc_ratio;

param.conc = 1;

%Convert from uM to molecules/um^3

param.conc = param.conc * conversion;
u_h((CN+1):2*CN) = param.conc*ones(CN,1);

u_h(2*CN+MN+1:end,:) = param.conc/(param.k_off_d/param.k_on_d+param.conc) * param.S * ones(MN,1);
solver_params.u_h = u_h;
solver_params.u_M = u_h(2*CN+1:2*CN+MN);
solver_params.v_M = u_h(2*CN+MN+1:end);

Soln = SolveNonLinear(param,solver_params);

u_C = Soln(1:CN,:);
v_C = Soln((CN+1):2*CN,:);
u_M = Soln(2*CN+1:2*CN+MN,:);
v_M = Soln(2*CN+MN+1:end,:);

sol_C = u_C + v_C;
sol_M = u_M + v_M;

tlist = (0:param.num_steps)*param.dt;
desired_times = (0:param.store_interval:param.num_steps)+1;
PlotSol(Soln,photo_on_scale,CN,MN,props,param)

pdem_C = createpde(1);
gm_C_f = geometryFromMesh(pdem_C,props.nodes',props.elements');
pde_C=createPDEResults(pdem_C,sol_C(:,1:param.interpolation_interval:end),tlist(desired_times(1:param.interpolation_interval:end)),'time-dependent');
save([filename,'-TIRF-hr-1um-nl-3D.mat'],'Soln','sol_M','I','contours','param','props','tlist','desired_times','pdem_C','gm_C_f');
%clear Soln u_C u_M v_C v_M sol_C
%% Interpolate  just at z = 0
z = 0;
scale =1;
samp_res = param.scale_len;
y = (0:scale:size(I,1)-1)*samp_res + scale*samp_res/2;
x = (0:scale:size(I,2)-1)*samp_res + scale*samp_res/2;

[X,Y] = meshgrid(x,y);
bottom_nodes= find(props.pm_surface_nodes(:,3) == 0);
plot(props.pm_surface_nodes(bottom_nodes,1),props.pm_surface_nodes(bottom_nodes,2),'o')
bottom_sol = sol_M(bottom_nodes,:);
m_intrp = NaN(size(X,1),size(X,2),size(bottom_sol,2));
for i = 1:size(bottom_sol,2)
    i
    m_intrp(:,:,i) = griddata(props.pm_surface_nodes(bottom_nodes,1),props.pm_surface_nodes(bottom_nodes,2),bottom_sol(:,i),X,Y);
end
%Max at 301,305
%%
%load('2P_m_intrp.mat')
%Max tirf, 33.31
%Max 2P, 41.2413
%Max 1P, 82.9356
%1-100 are nan, 503-601 are nan
sub_m = m_intrp - m_intrp(:,:,1);
profile = squeeze(sub_m(:,305,:));
x = -30:0.1:30;
clipped_x = x(101:502);
clipped = profile(101:502,:);

%%
figure
plot(x,profile(:,[2,3,6,11,21]),'LineWidth',2)
legend('0.1s','.2s','.5s','1s','2s')
set(gca,'linew',2)
set(gca,'tickdir','out')
ylim([0,1.05*max(profile(:))])
xlim([-20,20])
%xticks(['-30','-15','0','15','30'])
box off
figure
plot(x,profile(:,[31,101,201,401,601]),'LineWidth',2)
legend('2s','10s','20s','40s','60s')
set(gca,'linew',2)
set(gca,'tickdir','out')
box off
ylim([0,1.05*max(profile(:))])
xlim([-20,20])


%%

wid = [];
wid_m = [];
cen = [];
amp=[];
max_height = [];
resamp = min(clipped_x):0.01:max(clipped_x);
for i=1:size(profile,2)
    f = clipped(:,i);
    %halfmax = (min(f) +max(f))/2;
    fr = interp1(clipped_x,f,resamp);
    halfmax = max(f)/2;
    fifthmax = 82.9356/2;
    idx1 = find(fr>= halfmax,1,'first');
    idx2 = find(fr>= halfmax,1,'last');
    fwhm  = resamp(idx2) - resamp(idx1);
    idx1_m = find(fr>= fifthmax,1,'first');
    idx2_m = find(fr>= fifthmax,1,'last');
    fwhm_m  = resamp(idx2_m) - resamp(idx1_m);
    if size(fwhm_m,2)==0
       fwhm_m = NaN; 
    end
    %f(175:230) = max(f);
    f_g = fit(clipped_x',f,'gauss1','Start',[300,0,5],'Lower',[max(f)*1.2,-1,0],'Upper',[1000,1,100]);
    %wid = [wid,f_g.c1];
    wid = [wid,fwhm / 2.355];
    wid_m = [wid_m,fwhm_m];
    cen = [cen,f_g.b1];
    amp = [amp,f_g.a1];
    max_height = [max_height,max(f)];
    if mod(i,60) == 0
        plot(f_g,clipped_x,f)
    end
end
wid(1) = 0;%0.75/sqrt(2*log(2));
figure
t = 0:0.1:60;
plot(t,max_height,'k','LineWidth',2)
xlabel('Time (s)')
ylabel('Amplitude (molecules/um^2)')
set(gca,'FontSize',30)
box on
set(gca,'linew',2)
set(gca,'tickdir','out')
box off
xlim([0,60])
set(gcf,'Position',[1200,1000,800,720])
figure
plot(t,wid,'k','LineWidth',2);


 
xlabel('Time (s)')
ylabel('Standard Deviation ({\mu}m)')
set(gca,'FontSize',30)
box on
set(gca,'linew',2)
set(gca,'tickdir','out')
box off
xlim([0,60])
set(gcf,'Position',[1200,1000,800,720])

figure
plot(t,wid_m,'k','LineWidth',2);


 
xlabel('Time (s)')
ylabel('Width (um)')
set(gca,'FontSize',30)
box on
set(gca,'linew',2)
set(gca,'tickdir','out')
box off
xlim([0,60])
set(gcf,'Position',[1200,1000,800,720])
%%
mem_sub = sol_M - sol_M(:,1);
bottom_nodes= find(props.pm_surface_nodes(:,3) == 0);
top_nodes = find(props.pm_surface_nodes(:,3) > 7);
ratio = sum(mem_sub(bottom_nodes,:)) ./ (sum(mem_sub(bottom_nodes,:))+sum(mem_sub(top_nodes,:)));
ratio(1) = 1;
plot(ratio)
%% Interpolate  just at z = 0
z = 0;
scale =1;
samp_res = param.scale_len;
y = (0:scale:size(I,1)-1)*samp_res + scale*samp_res/2;
x = (0:scale:size(I,2)-1)*samp_res + scale*samp_res/2;
desired_idx = 1:size(Soln,2);
TR = props.pm_TR;
m_intrp = NaN(size(y,2),size(x,2),size(z,2),size(desired_idx,2),'single');
h = waitbar(0,'Interpolating membrane...');
for k=1:size(m_intrp,4)
    waitbar(k/size(m_intrp,4));
    F = scatteredInterpolant(TR.Points(:,1),TR.Points(:,2),TR.Points(:,3),sol_M(:,desired_idx(k)));
    planes.n = [zeros(size(z,2),2),ones(size(z,2),1)];
    planes.r = [zeros(size(z,2),2),z'];
    polygons = mesh_xsections( TR.Points,TR.ConnectivityList, planes, 1e-3, 0 );

    P_m_i = {};
    N_m_i = {};
    for i=1:size(z,2)
        if length(polygons{i}) == 0
            continue
        end
        [~,idx] = max(cell2mat(cellfun(@(c) length(c),polygons{i},'UniformOutput',false)));
        X = polygons{i}{idx}(:,1);
        Y = polygons{i}{idx}(:,2);
        X = [X;X(1)];
        Y = [Y;Y(1)];

        pt = interparc(size(polygons{i}{idx}(:,1),1)*5,X,Y,'linear');

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
        X = round((X-scale*samp_res/2)/samp_res/scale)*samp_res*scale + scale*samp_res/2;
        Y = round((Y-scale*samp_res/2)/samp_res/scale)*samp_res*scale + scale*samp_res/2;
        [~,idx_x] = ismembertol(X,x,1e-4);
        [~,idx_y] = ismembertol(Y,y,1e-4);



        [C,ia,ic] = unique([idx_x,idx_y],'rows');

        C = [C,i*ones(size(C,1),1),k*ones(size(C,1),1)];
        %Average function values when multiple points are in a pixel
        pixel_val = zeros(size(C,1),1);
        for j = 1:size(C,1)
            idx = find(ic == j);
            pixel_val(j) = mean(N(idx));
        end 
        %if k == 1 & i == floor(size(z,2)/2)
        %   plot(C(:,2),pixel_val,'o') 
        %end     

        C=sub2ind(size(m_intrp),C(:,2),C(:,1),C(:,3),C(:,4));
        m_intrp(C) = pixel_val;

    end
end
close(h);
%% Interpolate 
%[voxel_mem_area,pixel_map] = MembraneArea(I,props.surface_TR,param);
param.axial_resolution = 0.175;
c_intrp = InterpolateCytoplasm(pde_C,1:length(desired_times(1:param.interpolation_interval:end)),I,param);

m_intrp = InterpolateMembrane(I,1:param.interpolation_interval:length(desired_times),props.pm_TR,sol_M,param);
%m_intrp = ndSparse(m_intrp,[size(m_intrp)]);
%[c_intrp] = InterpolateCytoplasm(sol_C,1:8:length(desired_times),I,props,param);
m_intrp(isnan(m_intrp)) = 0;
%m_intrp = m_intrp * 1.8448;%1.3869
%c_intrp = c_intrp / conversion * 452.7271;

c_intrp = c_intrp/conversion * param.conc_ratio + param.offset;
%m_intrp = m_intrp * 1.8448 * param.conc_ratio/452.7271;
%m_intrp = m_intrp * 1.85;
%c_intrp(m_intrp ~= 0) = m_intrp(m_intrp~=0);
m_intrp = m_intrp * 1.2;
c_intrp(m_intrp ~= 0) = 0.5*c_intrp(m_intrp ~= 0) + m_intrp(m_intrp~=0);

for i = 1:size(c_intrp,4)
	imwrite(uint16(squeeze(c_intrp(:,:,ceil(size(c_intrp,3)/2),i))),['./',folder_name_no_PSF,'/',num2str(i),'.tif']);
end

%% Unwrap membrane
mem_slice = squeeze(m_intrp(:,:,ceil(size(m_intrp,3)/2),:));
sub_slice = mem_slice - mem_slice(:,:,1);
c_slice = c_intrp_blurred(:,:,ceil(size(m_intrp,3)/2),:);
c_slice = c_slice - c_slice(:,:,1);
[y,x] = find(sub_slice(:,:,100)>0);
contour = polyshape(y,x,'Simplify',false);
contour =bwtraceboundary(sub_slice(:,:,100)>0,[y(1),x(1)],'W');
idx = sub2ind([size(mem_slice,1),size(mem_slice,2)],contour(:,1),contour(:,2));
st = sub_slice(:,:,2);
[~,idx_max] = max(st(idx));
dist = [0;cumsum(sqrt(sum((contour(1:end-1,:) - contour(2:end,:)).^2,2)))] * param.scale_len;
dist = dist - dist(idx_max);

sigma = 5;
sz = 30;    % length of gaussFilter vector
x = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize


profile = [];
for i = 1:size(sub_slice,3)
    temp = sub_slice(:,:,i);
    %profile = [profile,temp(idx)];
    profile = [profile,conv(temp(idx),gaussFilter,'same')];
end

%Circularly permute
offset = idx_max - floor(size(profile,1)/2);
idx_all = [(offset+1):size(profile,1),1:offset];
dist = dist(idx_all);
dist(end-offset+1:end) = dist(end-offset+1:end) - dist(end-offset+1) + dist(end-offset);
profile = profile(idx_all,:);

low = -7.5;
high = 7.5;
[~,idx_min] = min(abs(dist- low));
[~,idx_max] = min(abs(dist- high));
%idx_min = 1;
%idx_max = size(profile,1);
clipped = profile(idx_min:idx_max,:);
clipped_x = dist(idx_min:idx_max);
ti = 0:0.1:60;
to_show = clipped ./ max(clipped(:));
h=heatmap(to_show)
set(h.NodeChildren(3), 'XTickLabelRotation', 0);
h.GridVisible = 'off';
h.Colormap = jet;
%h.ColorLimits = [0,prctile(f_post(:),99.9)];
h.ColorLimits = [0,1];
h.XLabel = 'Time (s)';
h.YLabel = 'Position ({\mu}m)';
cdl = h.XDisplayLabels; 
%xd = num2cell(t_post);
xd = repmat(NaN,size(cdl,1), size(cdl,2));
xd(1:150:end) = string(ti(1:150:end))';
h.XDisplayLabels = xd;
cdl = h.YDisplayLabels; 
%xd = num2cell(t_post);
xd = repmat(NaN,size(cdl,1), size(cdl,2));
xd(1:20:end) = string(clipped_x(1:20:end))';
h.YDisplayLabels = xd;
h = gca;
h.FontSize = 18;

%%
low = -7.5;
high = 7.5;
[~,idx_min] = min(abs(dist- low));
[~,idx_max] = min(abs(dist- high));
%idx_min = 1;
%idx_max = size(profile,1);
clipped = profile(idx_min:idx_max,:);
clipped_x = dist(idx_min:idx_max);


wid = [];
cen = [];
amp=[];
resamp = min(clipped_x):0.01:max(clipped_x);
for i=1:size(clipped,2)
    %f = clipped(:,i) - (clipped(1,i)+clipped(end,i))/2;
    f = clipped(:,i) - clipped(end,i);
    fr = interp1(clipped_x,f,resamp);
    halfmax = max(f)/2;
    idx1 = find(fr>= halfmax,1,'first');
    idx2 = find(fr>= halfmax,1,'last');
    %fwhm  = resamp(idx2) - resamp(idx1);
    %fwhm = resamp(idx2) *2;
    %f(1:103) = flipud(f(105:207));
    %f(105:207) = flipud(f(1:103));
    f_g = fit(clipped_x,f,'gauss1','Start',[300,0,5],'Lower',[0,-1,0],'Upper',[1000,1,100]);
    wid = [wid,f_g.c1];
    %wid = [wid,fwhm / 2.355];
    cen = [cen,f_g.b1];
    amp = [amp,f_g.a1];
    if mod(i,60) == 0
        plot(f_g,clipped_x,f)
    end
end
wid(1) = 0;%0.75/sqrt(2*log(2));
figure
t = 0:0.1:60;
plot(t,wid,'k','LineWidth',2);


 
xlabel('Time (s)')
ylabel('Standard Deviation ({\mu}m)')
set(gca,'FontSize',30)
box on
set(gca,'linew',2)
set(gca,'tickdir','out')
box off
xlim([0,15])
%xticks([0,5,10,15]);
%yticks([0,1,2,3,4,5,6])
%clear sol_M m_intrp
%%%%%%%%%%%%%FOR TOMORROW USE CONCENTRATION CALIBRATION TO GET TO
%%%%%%%%%%%%%FLUORESCENCE
%m_intrp = m_intrp - min(m_intrp(:));
%m_intrp(isnan(m_intrp)) = 0;
%interp_props.idx_nan = isnan(m_intrp(:));
%m_intrp(c_intrp(:) > 0 & nterp_props.idx_nan) = 0;
%interp_props.idx_c = c_intrp(:) > 0;
%interp_props.idx_c = interp_props.idx_c(1:size(interp_props.idx_nan,1)/size(c_intrp,4));
%interp_props.idx_m = find(~interp_props.idx_nan(1:size(interp_props.idx_nan,1)/size(c_intrp,4)));

%vol_int = CytoplasmVolume(I,interp_props,param);

%WriteVideo(uint8(m_intrp/max(m_intrp(:))*255),72,'mem-ex.avi',2,jet)
%WriteVideo(uint8(255*m_intrp/max(m_intrp(:))),72,'mem-ex.avi',2,jet)


%% Approximate PSF
filename = param.PSF;
%param.scale_z = 0.16;
%param.PSF_axial_ratio = param.scale_z/param.scale_len;
tstack = Tiff(filename,'r');

[i,j] = size(tstack.read());
K = length(imfinfo(filename));
data = zeros(i,j,K);
data(:,:,1)  = tstack.read();
for n = 2:K
    tstack.nextDirectory()
    data(:,:,n) = tstack.read();
end

PSF_3D = data;
s = sum(PSF_3D,'all');
PSF_3D = PSF_3D / s;

[nx,ny,nz] = size(PSF_3D);
cut=exp(-1:-1:-3)/s;


[X,Y,Z]=meshgrid(-floor(nx/2):nx/2,-floor(ny/2):ny/2,-floor(nz/2):nz/2);
X = X*param.scale_len;
Y=Y*param.scale_len;
Z=Z*param.axial_resolution;%0.2 for WF
%z_interp = min(Z(:)):param.scale_len:max(Z(:));
%[Xi,Yi,Zi]=meshgrid(-floor(nx/2):nx/2,-floor(ny/2):ny/2,z_interp);
%PSF_3D = interp3(X,Y,Z,PSF_3D,Xi,Yi,Zi,'nearest');

figure;
a = gca;
h = waitbar(0,'Please wait...');
for k=1:numel(cut)
    isonormals(X,Y,Z,PSF_3D,patch(isosurface(X,Y,Z,PSF_3D,cut(k)),'EdgeColor','none','FaceAlpha',1/k,'FaceColor',[1 (k-1)/max(1,numel(cut)-0.99) 0],'Parent',a));
    waitbar(k / numel(cut))
end
close(h);
view(35,45);
axis('equal');
lighting('gouraud');
grid('on');
camlight;
set(gca,'XColor', 'none','YColor','none','ZColor','None')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(gca, 'ztick', [])
axis equal



% psf_params.size = [64,64,64];%[32,32,21]
% psf_params.NA = param.NA;
% psf_params.lambda = 610e-9;
% psf_params.M = param.mag;
% psf_params.ti0 = 100e-6;
% psf_params.resLateral = param.scale_len * 1e-6;%scale_len * 1e-6;
% psf_params.resAxial = param.axial_resolution * 1e-6;%scale_len * 1e-6;
% psf_params.pZ = 0;%ceil(params.size(3)/2) * axial_resolution*1e-6;
% psf_params.oversampling = 2;
% [PSF_3D] = GenPSF(psf_params,param);
%% Convolve with PSF
[c_intrp_blurred] = ConvolvePSF(c_intrp,single(PSF_3D));
%WriteVideo(uint8(c_intrp_blurred/max(c_intrp_blurred(:))*255),ceil(size(c_intrp,3)/2),[filename,'.avi'],10,gray)