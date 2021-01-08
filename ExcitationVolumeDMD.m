function [il_3D,coords] = ExcitationVolumeDMD(mask_il,z,params)

    mask = mask_il';
    res = params.scale_len;
    alpha = asind(params.NA/params.immersion_n);
    
    within_aperture = @(x,y,z,xo,yo,zo) sqrt((x-xo).^2 + (y-yo).^2) ./ abs(z-zo) < tand(alpha);
    lens_t = @(x,y,z,xo,yo,zo) within_aperture(x,y,z,xo,yo,zo)./((x-xo).^2+(y-yo).^2+(z-zo).^2);
    
    [y,x] = find(mask==1) ;
    x0 = min(x); y0 = min(y) ;
    x1 = max(x); y1 = max(y) ;
    on = mask(y0:y1,x0:x1);
    %on = imresize(on,upsample);
    
    X = tand(alpha) * min(z) - size(on,2)/2 * res:res:tand(alpha) * max(z) + size(on,2)/2 * res;
    Y = tand(alpha) * min(z) - size(on,1)/2 * res:res:tand(alpha) * max(z) + size(on,1)/2 * res;

    [X,Y,Z] = meshgrid(X,Y,z);
    coords.X = X;
    coords.Y = Y;
    coords.Z = Z;
    
    %Embed
    on_uncropped = zeros(size(X,1),size(X,2));
    center_field = size(X) / 2;
    center_field = center_field(1:2);
    y_embed_s = floor(center_field(1) - size(on,1)/2+1);
    y_embed_e = floor(center_field(1) + size(on,1)/2);
    x_embed_s = floor(center_field(2) - size(on,2)/2+1);
    x_embed_e = floor(center_field(2) + size(on,2)/2);
    on_uncropped(y_embed_s:y_embed_e,x_embed_s:x_embed_e) = on;
    
    irradiance_integ = lens_t(X,Y,Z,0,0,0);
    norm_i = mean(sum(sum(irradiance_integ(:,:,1:floor(size(z,2)/2)),1),2));
    irradiance_integ = irradiance_integ / norm_i;
    C = zeros(size(irradiance_integ));
    for i = 1:size(z,2)
        C(:,:,i) = conv2(irradiance_integ(:,:,i),on,'same') / size(on,1) / size(on,2);
    end
    irradiance_integ = C * sum(on_uncropped,'all');
    irradiance_integ(:,:,ceil(size(z,2)/2)) = on_uncropped;
    il_3D = irradiance_integ;%permute(irradiance_integ,[2,1,3]);
    
    
%     %Plot
%     cut=fliplr(logspace(log10(max(max(max(irradiance_integ)))/4),log10(max(max(max(irradiance_integ)))),5));
% 
% 
%     figure;
%     a = gca;
% 
%     h = waitbar(0,'Please wait...');
% 
%     for k=1:numel(cut)
%         isonormals(X,Y,Z,irradiance_integ,patch(isosurface(X,Y,Z,irradiance_integ,cut(k)),'EdgeColor','none','FaceAlpha',1/k,'FaceColor',[1 (k-1)/max(1,numel(cut)-0.99) 0],'Parent',a));
%         waitbar(k / numel(cut))
%     end
%     close(h);
%     view(35,45);
%     axis('tight');
%     lighting('gouraud');
%     grid('on');
%     camlight;
end