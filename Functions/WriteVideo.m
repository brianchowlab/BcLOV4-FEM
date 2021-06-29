function [] = WriteVideo(data,slice,filename,frame_rate,cmap)
    writerObj = VideoWriter(filename,'Indexed AVI');
    writerObj.FrameRate = frame_rate;
    writerObj.Colormap = cmap;
    open(writerObj);
    for t=1:size(data,4)
      %frame = im2frame(imadjust(data(:,:,slice,t)/max(data(:)),[0,1],[0,1],2),cmap);
      frame = im2frame(data(:,:,slice,t),cmap);
      writeVideo(writerObj, frame);
    end
    close(writerObj);
end