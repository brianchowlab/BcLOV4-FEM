%% Set parameters
fid = fopen('params.txt');
params = textscan(fid,'%s %f','delimiter',' ','HeaderLines',4,'MultipleDelimsAsOne',1,'CommentStyle','%');
fclose(fid);
fid = fopen('params.txt');
ims_and_rois = textscan(fid,'%s %s','delimiter',' ','MultipleDelimsAsOne',1,'CommentStyle','%');
fclose(fid);
param = cell2struct(num2cell(params{2}),string(params{1})');
param.im_file = ims_and_rois{2}{1};
param.il_roi_file = ims_and_rois{2}{2};
param.nucleus_roi_file = ims_and_rois{2}{3};
param.cyto_roi_file = ims_and_rois{2}{4};


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
cd FELICITY;FELICITY_paths;cd ..;

%Port 3D cytoplasm mesh into FELICITY
Mesh = MeshTetrahedron(mesh_c.Elements',mesh_c.Nodes','Omega');

props = MeshProps(Mesh,shp_n);

[~,idx] = ismember(props.surface_nodes,props.nodes,'rows');
Mesh = Mesh.Append_Subdomain('2D','dOmega',props.faces);

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
pde_C=createPDEResults(pdem_C,sol_C,tlist(desired_times),'time-dependent');
save('data.mat','Soln','props','tlist','desired_times','pdem_C','gm_C_f');
%clear v_M u_M v_C u_C
%clear sol_C sol_all sol_M_embed
%clear v_C_1 v_C_2 v_C_3 v_C_4 v_C_nl_2 v_M_1 v_M_2 v_M_3 v_M_4 v_M_nl_2 v_C_nl_2
%clear u_M_1 u_M_2 u_M_3 u_M_4 u_M_nl_1 u_C_1 u_C_2 u_C_3 u_C_4
%clear shp_n shp_cn shp_c s pre il_phi il_psi il_M il_3D
%clear Z params mask
%clear Soln A_phi A_psi B_phi B_psi FEM_l FEM_nl J K_phi K_phi_phi K_phi_psi K_psi K_psi_phi
%clear Step_1_LHS_dark Step_1_LHS_lit Step_1_LHS_decomp Step_1_LHS_decomp_dark Step_1_LHS_decomp_lit
%% Interpolate 
%[voxel_mem_area,pixel_map] = MembraneArea(I,props.surface_TR,param);
[m_intrp]= InterpolateMembrane(I,1:length(desired_times),props.surface_TR,sol_M,param);
[c_intrp] = InterpolateCytoplasm(pde_C,1:length(desired_times),I,param);
m_intrp(isnan(m_intrp)) = 0;
m_intrp = m_intrp * 1.8448 + 436.14;
c_intrp = c_intrp / conversion * 452.7271;
c_intrp = m_intrp + c_intrp;
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

psf_params.size = [128,128,64];%[32,32,21]
psf_params.NA = param.NA;
psf_params.lambda = 610e-9;
psf_params.M = param.mag;
psf_params.ti0 = 100e-6;
psf_params.resLateral = param.scale_len * 1e-6;%scale_len * 1e-6;
psf_params.resAxial = param.scale_len * 1e-6;%scale_len * 1e-6;
psf_params.pZ = 0;%ceil(params.size(3)/2) * axial_resolution*1e-6;
psf_params.oversampling = 2;
[PSF_3D] = GenPSF(psf_params,param);
%% Convolve with PSF
[c_intrp_blurred] = ConvolvePSF(c_intrp,PSF_3D);
WriteVideo(uint8(c_intrp_blurred/max(c_intrp_blurred(:))*255),ceil(size(c_intrp,3)/2),'out.avi',10,jet)