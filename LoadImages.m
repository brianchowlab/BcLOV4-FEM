function [contours,mask_il,I] = LoadImages()
    im_file = input('Location of Image: ');%'C:\Users\Ivan Kuznetsov\Desktop\PhD\BcRGS5\1-2_c.png'
    I = imread(im_file);
    I = imresize(I,[floor(size(I,1)/2),floor(size(I,2)/2)]);
    I = flipud(I);
    imshow(I)
    %Select nucleus
    r = drawrectangle;
    mask = createMask(r);
    nucleus = activecontour(I,mask,500);
    close

    mask = ones(size(I));
    bw = activecontour(I,mask,500);
    se = strel('disk',5);
    bw = imclose(bw,se);
    cytoplasm = bwareaopen(bw,100);
    imshow(I)
    hold on;
    visboundaries(cytoplasm,'Color','r'); 
    visboundaries(nucleus,'Color','g');
    cytoplasm_only = cytoplasm-nucleus;
    domain = cytoplasm_only;

    roi = ReadImageJROI('C:\Users\Ivan Kuznetsov\Desktop\PhD\BcRGS5\1-2.roi');
    il_contour = fliplr(size(I)) - roi.mnCoordinates/2;
    mask_il = poly2mask(il_contour(:,1),il_contour(:,2),size(I,1),size(I,2));
    hold on
    visboundaries(mask_il)
    %r = drawrectangle;
    %mask_il = createMask(r);

    %Contour describing cytoplasmic border
    cytoplasm_contour = bwboundaries(cytoplasm);
    cytoplasm_contour = cytoplasm_contour{1};

    %Contour describing nuclear border
    nucleus_contour = bwboundaries(nucleus);
    nucleus_contour = nucleus_contour{1};
    num_cyto = size(cytoplasm_contour,1);
    num_nuc = size(nucleus_contour,1);

    contours.cytoplasm = cytoplasm_contour;
    contours.nucleus = nucleus_contour;
    contours.il = il_contour;
end