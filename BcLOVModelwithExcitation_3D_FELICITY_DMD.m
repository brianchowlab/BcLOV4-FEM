%% Set parameters
filename = 'params_cell_1_1s_10dc';

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
%% Load image
[contours,mask_il,I] = LoadImages(param);

%% Make Mesh
[mesh_c,poly,shp_n] = GenMesh(contours,param);

% Determine illuminated region
[photo_on_scale,idx_excited] = ExcitationROI(mesh_c,mask_il,poly,param);

%% Build FEM matrices with FELICITY and solve.
%cd FELICITY;FELICITY_paths;cd ..;

addpath ./FELICITY
FELICITY_paths

%Port 3D cytoplasm mesh into FELICITY
Mesh = MeshTetrahedron(mesh_c.Elements',mesh_c.Nodes','Omega');

props = MeshProps(Mesh,shp_n);
%plot(props.pm_surface_nodes(abs(props.pm_surface_nodes(:,3))<1,1),props.pm_surface_nodes(abs(props.pm_surface_nodes(:,3))<1,2),'o')
[~,idx] = ismember(props.pm_surface_nodes,props.nodes,'rows');
%m = find(abs(Mesh.Points(:,3))<1);
%m = intersect(m,idx);
%m = Mesh.Points(m,:);
%plot(m(:,1),m(:,2),'o')
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
param.conc = (mean(I(logical(mask_cyto_only(:))))-param.offset)/param.conc_ratio;
%param.conc = 1.1552;
%param.conc = 0.5385;

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
%save([filename,'-nl-3D.mat'],'Soln','sol_M','I','contours','param','props','tlist','desired_times','pdem_C','gm_C_f');
%clear Soln u_C u_M v_C v_M sol_C
%% Interpolate 
%[voxel_mem_area,pixel_map] = MembraneArea(I,props.surface_TR,param);
c_intrp = InterpolateCytoplasm(pde_C,1:length(desired_times(1:param.interpolation_interval:end)),I,param);

m_intrp = InterpolateMembrane(I,1:param.interpolation_interval:length(desired_times),props.pm_TR,sol_M,param);
%m_intrp = ndSparse(m_intrp,[size(m_intrp)]);
%[c_intrp] = InterpolateCytoplasm(sol_C,1:8:length(desired_times),I,props,param);
m_intrp(isnan(m_intrp)) = 0;
%m_intrp = m_intrp * 1.8448;%1.3869
%c_intrp = c_intrp / conversion * 452.7271;

c_intrp = c_intrp/conversion * param.conc_ratio + param.offset;
%m_intrp = m_intrp * 1.8448 * param.conc_ratio/452.7271;
m_intrp = m_intrp * 1.85;
c_intrp(m_intrp ~= 0) = m_intrp(m_intrp~=0);
clear sol_M m_intrp
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
param.scale_z = 0.2;
param.PSF_axial_ratio = param.scale_z/param.scale_len;
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
Z=Z*param.scale_z;%0.2 for WF
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
% psf_params.resAxial = param.scale_len * 1e-6;%scale_len * 1e-6;
% psf_params.pZ = 0;%ceil(params.size(3)/2) * axial_resolution*1e-6;
% psf_params.oversampling = 2;
% [PSF_3D] = GenPSF(psf_params,param);
%% Convolve with PSF
[c_intrp_blurred] = ConvolvePSF(c_intrp,single(PSF_3D));
WriteVideo(uint8(c_intrp_blurred/max(c_intrp_blurred(:))*255),ceil(size(c_intrp,3)/2),[filename,'.avi'],10,gray)