clear all, close all, clc

% filepath = '/home/barracuda/a/lrajendr/turbulent_channel_Mb=0.85_Reb=6900/';
%filepath = '/home/barracuda/a/lrajendr/turbmat-20150108/JHU_Buoyancy_Turb/';
filepath = '/home/barracuda/a/lrajendr/SchlierenRayVis/StereoBOS_pinhole/Processing/BOS/Camera_2/';

files = dir([filepath '*.bin']);
fprintf('No. of files : %d\n', numel(files));

N = numel(files);

%filename = 'Movie';
%writerObj = VideoWriter(strcat(filepath,filename),'Grayscale AVI');
%writerObj.FrameRate = 1;
%open(writerObj)

height = 1024;
width = 1024;
data1 = zeros(height,width);
%N = 1;
for i = 1:N
    i
    if(mod(i,10)==0)
      fprintf('file number : %d\n',i)
    end
    
    fileID = fopen([filepath files(i).name],'r');   
    data1 = fread(fileID,'uint32');
    data1 = reshape(data1,height,width);
    data1 = permute(data1,[2 1]);
    %data1 = data1/(2^32-1)*(2^16-1);
    data1 = uint8(data1);
%     if(i==1)
%         imshow(uint8(data2))
%         title({'Select Region of Interest to be Processed';'First Click Top Left Then Bottom Right of Region of Interest'})
%         coord   =   ginput(2);
%         y1      =   coord(1,2);
%         y2      =   coord(2,2);
%         x1      =   coord(1,1);
%         x2      =   coord(2,1);
%         x1 = round(x1); x2 = round(x2); y1 = round(y1); y2 = round(y2);
%     end
%     
%       data2 = data2(y1:y2,x1:x2);
      
    %writeVideo(writerObj,uint8(data2));
%     save([filepath files(i).name(1:end-4) '.mat'],'data2');
    %     caxis([1*10^6 1*10^7])
    
    %imshow(data1);
    imwrite(data1,[filepath files(i).name(1:end-4) '.tif'],'tif');
    fclose(fileID);
    
    %imwrite(imagesc(data1
end

%close(writerObj);

%implay([filepath 'Movie.avi'])