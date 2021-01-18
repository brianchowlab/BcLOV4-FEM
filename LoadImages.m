function [contours,mask_il,I] = LoadImages(param)
    I = imread(param.im_file);
    %I = imresize(I,[floor(size(I,1)/2),floor(size(I,2)/2)]);
    I = fliplr(flipud(I));

    %Select nucleus
    if strcmp(param.nucleus_roi_file,'NA')
        imshow(I)
        r = drawrectangle;
        mask = createMask(r);
        nucleus = activecontour(I,mask,500);
        close
        nucleus_contour = bwboundaries(nucleus);
        nucleus_contour = nucleus_contour{1};
    else
        nucleus = ReadImageJROI(param.nucleus_roi_file);
        nucleus_contour = fliplr(size(I)) - nucleus.mnCoordinates;
        nucleus = poly2mask(nucleus_contour(:,1),nucleus_contour(:,2),size(I,1),size(I,2));
    end

    if strcmp(param.cyto_roi_file,'NA')
        mask = ones(size(I));
        bw = activecontour(I,mask,500);
        se = strel('disk',5);
        bw = imclose(bw,se);
        cytoplasm = bwareaopen(bw,100);
        cytoplasm_contour = bwboundaries(cytoplasm);
        cytoplasm_contour = cytoplasm_contour{1};
    else
        cytoplasm = ReadImageJROI(param.cyto_roi_file);
        cytoplasm_contour = fliplr(size(I)) - cytoplasm.mnCoordinates;
        cytoplasm = poly2mask(cytoplasm_contour(:,1),cytoplasm_contour(:,2),size(I,1),size(I,2));
    end
    
    if strcmp(param.il_roi_file,'NA')
        mask_il = ones(size(I));
        il_contour = bwboundaries(mask_il);
        il_contour = il_contour{1};
    else
        roi = ReadImageJROI(param.il_roi_file);
        il_contour = fliplr(size(I)) - roi.mnCoordinates;
        mask_il = poly2mask(il_contour(:,1),il_contour(:,2),size(I,1),size(I,2));
    end
    
    if param.plot
        imshow(I)
        hold on;
        visboundaries(cytoplasm,'Color','r'); 
        visboundaries(nucleus,'Color','y');
        visboundaries(mask_il,'Color','b')
    end

    contours.cytoplasm = cytoplasm_contour;
    contours.nucleus = nucleus_contour;
    contours.il = il_contour;
end