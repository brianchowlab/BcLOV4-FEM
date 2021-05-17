function [] = BcLOVModelwithExcitation_3D_Linear_Loop(filename)
%% Set parameters
fid = fopen(filename);
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

%%%%%%%%%%%%%%%%%%Manually Set
param.k_off_p = 1/18.5;%0.2632;
param.k_off_d = 0.0225;
param.k_on_d = 1130;
param.offset=104;
param.dt = 5e-2;
param.num_steps = 4000;
param.store_interval = 20;
param.tol = 1e-3;

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
[mesh_c,poly,shp_n] = GenMesh2D(contours,param);

% Determine illuminated region
idx_excited = inpolygon(mesh_c.Nodes(1,:),mesh_c.Nodes(2,:),poly.il.Vertices(:,1),poly.il.Vertices(:,2))';
photo_on_scale = double(idx_excited);

%% Build FEM matrices with FELICITY and solve.
%cd FELICITY;FELICITY_paths;cd ..;

addpath ./FELICITY
FELICITY_paths

%Port 3D cytoplasm mesh into FELICITY
Mesh = MeshTriangle(mesh_c.Elements',mesh_c.Nodes','Omega');

props = MeshProps2D(Mesh);
%plot(props.pm_surface_nodes(abs(props.pm_surface_nodes(:,3))<1,1),props.pm_surface_nodes(abs(props.pm_surface_nodes(:,3))<1,2),'o')
Mesh = Mesh.Append_Subdomain('1D','dOmega',props.pm_edges);

%Calculate which nodes are illuminated
param.il_c = photo_on_scale;%inShape(shp_l,Mesh.Points);
%Mapping from membrane nodes to cytoplasmic nodes
param.il_m = photo_on_scale(props.pm_surface_node_idx);

% Define FE spaces
C_RefElem_2D = ReferenceFiniteElement(lagrange_deg1_dim2());
C_Space = FiniteElementSpace('C_h',C_RefElem_2D,Mesh,'Omega');
C_Space = C_Space.Set_DoFmap(Mesh,uint32(Mesh.ConnectivityList));

M_RefElem_1D = ReferenceFiniteElement(lagrange_deg1_dim1());
M_Space = FiniteElementSpace('M_h',M_RefElem_1D,Mesh,'dOmega');
M_Space = M_Space.Set_DoFmap(Mesh,uint32(props.pm_edges_renum));

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
solver_params.idx = props.pm_surface_node_idx;

mask_cyto_only = poly2mask(contours.cytoplasm(:,1),contours.cytoplasm(:,2),size(I,1),size(I,2))-poly2mask(contours.nucleus(:,1),contours.nucleus(:,2),size(I,1),size(I,2));
u_h = zeros(2*CN + 2*MN,1);

%Dark, cytosolic state initial condition (assuming no protein in lit
%state initially).
%Convert from fluorescence to concentration

data_file = split(filename,'.txt');
data_file = data_file{1};
data_file = split(data_file,'params_c');
data_file = data_file{2};

data_file = ['Cyto_C',data_file,'.csv'];

data = readmatrix(data_file);
data = data(1,2);
%param.conc = (mean(I(logical(mask_cyto_only(:))))-param.offset)/param.conc_ratio;
param.conc = (data-param.offset)/param.conc_ratio;

%param.conc = 1.1552;
%param.conc = 0.5385;

%Convert from uM to molecules/um^3

param.conc = param.conc * conversion;
u_h((CN+1):2*CN) = param.conc*ones(CN,1);

u_h(2*CN+MN+1:end,:) = param.conc/(param.k_off_d/param.k_on_d+param.conc) * param.S * ones(MN,1);
solver_params.u_h = u_h;
solver_params.u_M = u_h(2*CN+1:2*CN+MN);
solver_params.v_M = u_h(2*CN+MN+1:end);

Soln = SolveLinear2D(param,solver_params);

u_C = Soln(1:CN,:);
v_C = Soln((CN+1):2*CN,:);
u_M = Soln(2*CN+1:2*CN+MN,:);
v_M = Soln(2*CN+MN+1:end,:);

sol_C = u_C + v_C;
sol_M = u_M + v_M;

tlist = (0:param.num_steps)*param.dt;
desired_times = (0:param.store_interval:param.num_steps)+1;
PlotSol2D(Soln,photo_on_scale,CN,MN,props,param)

pdem_C = createpde(1);
gm_C_f = geometryFromMesh(pdem_C,props.nodes',props.elements');
pde_C=createPDEResults(pdem_C,sol_C(:,1:param.interpolation_interval:end),tlist(desired_times(1:param.interpolation_interval:end)),'time-dependent');
save([filename,'-nl-2D.mat'],'Soln','sol_M','I','contours','param','props','tlist','desired_times','pdem_C','gm_C_f');
end
