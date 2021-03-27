filename = param.z_stack;
tstack = Tiff(filename,'r');

[i,j] = size(tstack.read());
K = length(imfinfo(filename));
data = zeros(i,j,K);
data(:,:,1)  = tstack.read();
for n = 2:K
    tstack.nextDirectory()
    data(:,:,n) = tstack.read();
end
data = double(data)/max(data(:));
I = imread(param.im_file);

corr_v = [];
for i = 1:size(data,3)
    corr_v = [corr_v,corr2(I,data(:,:,i))];
end

[~,idx_mid] = max(corr_v);

nucleus = ReadImageJROI(param.nucleus_roi_file);
nucleus_contour = nucleus.mnCoordinates;
nucleus = poly2mask(nucleus_contour(:,1),nucleus_contour(:,2),size(data,1),size(data,2));


cytoplasm = ReadImageJROI(param.cyto_roi_file);
cytoplasm_contour = cytoplasm.mnCoordinates;
cytoplasm = poly2mask(cytoplasm_contour(:,1),cytoplasm_contour(:,2),size(data,1),size(data,2));
%%
slice = idx_mid;
I = data(:,:,idx_mid);

% I = data(:,:,slice);
% I = imsharpen(I);
% [L,Centers] = imsegkmeans(uint8(I*255),4);
% avg = [];
% for i=1:4
%     avg = [avg,mean(I(L==i))];
% end
% [~,idx] = min(avg);
% 
% L_back = L==idx;
% L_front = ~L_back;
% 
% 
% L(L_back) = 0;
% L(L_front) = 1;
% bw  = logical(L);
% 
% se = strel('disk',5);
% bw = imclose(bw,se);
% cytoplasm = bwareaopen(bw,100);
% cytoplasm_and_nucleus = imfill(cytoplasm,'holes');
% nucleus = cytoplasm_and_nucleus & ~cytoplasm;

%imshow(I)
%hold on
%visboundaries(cytoplasm_and_nucleus,'Color','r');
%visboundaries(nucleus,'Color','y');
%set(gca,'units','pixels') % set the axes units to pixels
%x = get(gca,'position') % get the position of the axes
%set(gcf,'units','pixels') % set the figure units to pixels
%y = get(gcf,'position') % get the figure position
%set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
%set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels

cytoplasm_and_nucleus = cytoplasm | nucleus;

cytoplasm_store = cytoplasm_and_nucleus;
nucleus_store = nucleus;
background = I .* double(~cytoplasm_and_nucleus);
foreground = I .* double(cytoplasm_and_nucleus);
avgs = squeeze(mean(mean(data,1),2));
avgs = avgs - avgs(1);
idx = find(avgs > 0.1*(mean(avgs)));
idx = min(idx):max(idx);
slice=54;
idx=1:134;
%%
clf
T = [1 0 0 0;0 1 0 0;0 0 1.5 0];
for i=idx(10:20:end-10)
    h1 = slice3(uint8(1.5*255*data),T,3,i);
    %  h2 = slice3(squeeze(D),T,2,64);
    %  h3 = slice3(squeeze(D),T,3,14);
    set(h1,'EdgeColor','white','LineStyle','-');
end
colormap gray(88);
view(15,10); axis equal; axis vis3d;
light;
set(gca,'visible','off')
%%

cytoplasm_and_nucleus = cytoplasm_store;
nucleus = nucleus_store;
cyt_3d = zeros(size(data,1),size(data,2),size(idx,2));
nuc_3d = zeros(size(data,1),size(data,2),size(idx,2));
count = 0;

for_v = foreground(foreground>0);
for i = fliplr(idx(1):slice)
    count = count+1;
    I = data(:,:,i);
    I(I > quantile(for_v,0.02)) = 1;
    I(I <= quantile(for_v,0.02)) = 0;
    se = strel('disk',2);
    %I = imclose(I,se);
    bw = activecontour(I,cytoplasm_and_nucleus,500,'Chan-Vese','ContractionBias',.1);
    cytoplasm_and_nucleus = bwareaopen(imfill(imclose(bw,se),'holes'),100);
    cyt_3D(:,:,slice-idx(1)-count+2) = cytoplasm_and_nucleus;
    figure
    imagesc(data(:,:,i))
    hold on;
    visboundaries(cytoplasm_and_nucleus,'Color','r');
    [r, c] = find(nucleus == 1);
    r = floor(mean(r));
    c = floor(mean(c));
    
    I = data(:,:,i);
    I(I > quantile(for_v,0.02)) = 1;
    I(I <= quantile(for_v,0.02)) = 0;
    %imagesc(I)
    mask = zeros(size(nucleus));
    %y = 1:size(nucleus,1);
    %x = 1:size(nucleus,2);
    %[X,Y] = meshgrid(x,y);
    %nucleus = (X-c).^2 + (Y-r).^2 < 4;
    %visboundaries(nucleus,'Color','k');
    bw = activecontour(I,nucleus,100,'Chan-Vese','SmoothFactor',0,'ContractionBias',0);
    se = strel('disk',3);
    nucleus = imclose(bw,se);
    nucleus_filled = imfill(nucleus,'holes');
    if sum(nucleus_filled(:))>sum(nucleus(:))
        nucleus_only = logical(nucleus_filled - nucleus);
    else
        nucleus_only = nucleus_filled;
    end


    %nucleus = nucleus_filled - nucleus;
    %nucleus = bw;
    if mean(I(nucleus_only)) > 0.5
        nucleus_only = zeros(size(nucleus_only));
    end
    %cyt_3D(:,:,count) = bwareaopen(bw,100);
    nuc_3D(:,:,slice-idx(1)-count+2) = nucleus_only;
    nucleus = nucleus_only;
    visboundaries(nucleus_only,'Color','g');
end

cytoplasm_and_nucleus = cytoplasm_store;
nucleus = nucleus_store;
for i = (slice+1):idx(end)
    count = count+1;
    I = data(:,:,i);
    I(I > quantile(for_v,0.02)) = 1;
    I(I <= quantile(for_v,0.02)) = 0;
    %se = strel('disk',2);
    %I = imclose(I,se);
    bw = activecontour(I,cytoplasm_and_nucleus,500,'Chan-Vese','ContractionBias',.1);
    cytoplasm_and_nucleus = bwareaopen(imfill(imclose(bw,se),'holes'),100);
    cyt_3D(:,:,count) = cytoplasm_and_nucleus;
    figure
    imagesc(data(:,:,i))
    hold on;
    visboundaries(cytoplasm_and_nucleus,'Color','r');
    [r, c] = find(nucleus == 1);
    r = floor(mean(r));
    c = floor(mean(c));
    
    I = data(:,:,i);
    I(I > quantile(for_v,0.02)) = 1;
    I(I <= quantile(for_v,0.02)) = 0;
    %imagesc(I)
    mask = zeros(size(nucleus));
    %y = 1:size(nucleus,1);
    %x = 1:size(nucleus,2);
    %[X,Y] = meshgrid(x,y);
    %nucleus = (X-c).^2 + (Y-r).^2 < 4;
    %visboundaries(nucleus,'Color','k');
    bw = activecontour(I,nucleus,100,'Chan-Vese','SmoothFactor',0,'ContractionBias',0.1);
    se = strel('disk',3);
    nucleus = imclose(bw,se);
    nucleus_filled = imfill(nucleus,'holes');
    if sum(nucleus_filled(:))>sum(nucleus(:))
        nucleus_only = logical(nucleus_filled - nucleus);
    else
        nucleus_only = nucleus_filled;
    end
    %nucleus = nucleus_filled - nucleus;
    %nucleus = bw;
    if mean(I(nucleus_only)) > 0.5
        nucleus_only = zeros(size(nucleus_only));
    end
    %cyt_3D(:,:,count) = bwareaopen(bw,100);
    nuc_3D(:,:,count) = nucleus_only;
    nucleus = nucleus_only;
    visboundaries(nucleus_only,'Color','g');
end
%%
nuc_3D_interp = zeros(size(nuc_3D));
for i = 1:size(nuc_3D,1)
    temp = squeeze(nuc_3D(i,:,:));
    nuc_3D_interp(i,:,:) = bwareaopen(imclose(temp,se),5);
end
nuc_3D_interp(cyt_3D == 0) = 0;
nuc_3D_interp_final = bwareaopen(nuc_3D_interp,250);

cell_segment = logical(cyt_3D - nuc_3D_interp);
f = split(filename,'.tif');
f = f{1};
for i = 1:size(nuc_3D,3)
    imwrite(nuc_3D_interp_final(:,:,i),[f,'-nuc.tif'],'WriteMode','append','Compression','none');
    imwrite(cyt_3D(:,:,i),[f,'-cyto.tif'],'WriteMode','append','Compression','none');
    imwrite(cell_segment(:,:,i),[f,'-segment.tif'],'WriteMode','append','Compression','none');
end

%%
temp = cyt_3D- nuc_3D_interp_final;
[x,y,z] = ind2sub(size(cyt_3D),find(temp(:) == 1));
shp_cn = alphaShape(x,y,z,10);
[elements,nodes] = boundaryFacets(shp_cn);
nodes(:,1) = 0.2167 * nodes(:,1);
nodes(:,2) = 0.2167 * nodes(:,2);
nodes(:,3) = 0.150 * nodes(:,3);


model = createpde();
gm_c = geometryFromMesh(model,nodes',elements');
%pdegplot(model,'CellLabels','on','FaceAlpha',0.5)
%mesh_c = generateMesh(model,'GeometricOrder','linear','Hmin',param.min_element_size,'Hmax',param.max_element_size);%mesh_c = generateMesh(model,'Hmin',0.5,'Hmax',30,'GeometricOrder','linear');
mesh_c = generateMesh(model,'GeometricOrder','linear','Hmin',0.5,'Hmax',2);%mesh_c = generateMesh(model,'Hmin',0.5,'Hmax',30,'GeometricOrder','linear');

figure
pdemesh(mesh_c,'FaceAlpha',0.1)
save([f,'-mesh.mat'],'mesh_c')

h = gca;
h = h.Children;
h(2).Visible = 'Off';
h(3).Visible = 'Off';
h(4).Visible = 'Off';
h(5).Visible = 'Off';


%%
[x,y,z] = ind2sub(size(cyt_3D),find(cyt_3D(:) == 1));
shp_c = alphaShape(x,y,z,2);
%plot(sh)
[elements,nodes] = boundaryFacets(shp_c);
nodes = nodes';
elements = elements';
nodes(:,1) = 0.2157 * nodes(:,1);
nodes(:,2) = 0.2157 * nodes(:,2);
nodes(:,2) = 0.150 * nodes(:,3);


shp_cn = alphaShape(x,y,z,2);

model = createpde();
gm_c = geometryFromMesh(model,nodes,elements);
%pdegplot(model,'CellLabels','on','FaceAlpha',0.5)
mesh_c = generateMesh(model,'GeometricOrder','linear','Hmin',param.min_element_size,'Hmax',param.max_element_size);%mesh_c = generateMesh(model,'Hmin',0.5,'Hmax',30,'GeometricOrder','linear');

pdemesh(mesh_c,'FaceAlpha',0.1)
save([f,'-mesh.mat'],'mesh_c')
h = gca;
h = h.Children;
h(2).Visible = 'Off';
h(3).Visible = 'Off';
h(4).Visible = 'Off';
h(5).Visible = 'Off';
