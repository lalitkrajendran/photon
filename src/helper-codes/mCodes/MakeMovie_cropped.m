% clear; clc; close all;

% value of sc is defined in CombinedAnalysis.m to ensure that a uniform
% scaling factor is used throughout

% sc = 1;
filename = 'Movie';
writerObj = VideoWriter(strcat(filepath,filename),'Grayscale AVI');
writerObj.FrameRate = 1;

load([filepath files(1).name]);
% load(strcat(path,'10.mat'));
temp = uint16(img);
imshow(sc * temp)

title({'Select Region of Interest to be Processed';'First Click Top Left Then Bottom Right of Region of Interest'})
coord   =   ginput(2);
y1      =   coord(1,2);
y2      =   coord(2,2);
x1      =   coord(1,1);
x2      =   coord(2,1);
x1 = round(x1); x2 = round(x2); y1 = round(y1); y2 = round(y2);

open(writerObj);
N = numel(files);

% N = 1;
for i = 1:dn_movie:N
    
    if(mod(i,10)==0) 
        i
    end
        
    load([filepath files(i).name]);
    img = uint16(img*sc);
%     img = img*sc;
%     imshow(img)
    img = double(img);
    
%     img = imresize(img,1/2);
    
%     img = img(y1:y2,x1:x2)/2^16;
    img = img(y1:y2,x1:x2)/(2^16);
    img = imresize(img,1/2);
    img = max(img,0);
    img = min(img,1);
    
    writeVideo(writerObj,img);

end

writerObj.VideoCompressionMethod
writerObj.VideoFormat
writerObj.FileFormat

close(writerObj);

implay([filepath 'Movie.avi'])