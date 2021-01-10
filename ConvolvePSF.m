function [c_intrp_blurred] = ConvolvePSF(data,PSF)
    c_intrp_blurred = zeros(size(data),'single');
    for i=1:size(data,4)%[1,size(desired_times,2)]%
        i/size(data,4)
        c_intrp_blurred(:,:,:,i) = convolution3D_FFTdomain(data(:,:,:,i),single(PSF));
    end
end