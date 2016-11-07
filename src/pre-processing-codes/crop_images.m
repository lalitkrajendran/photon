% this program reads in a set of images and crops them to remove the
% blurring at the edges due to real lens effects

clear 
close all
clc


case_name = '150x150-f16-disp8';
% this is the filepath where the images will be read
read_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/reordered-images/'];

% this is the filepath where the cropped images will be saved
save_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/cropped-images/'];

% these are the boundaries of the cropped region
xmin = 256;
xmax = 768;
ymin = 256;
ymax = 768;

% this loads a list of all the images in the directory
files = dir([read_filepath '*.tif']);

% this is the number of images in the directory
N = length(files);
fprintf('number of images: %d\n', N);

% this loops through all the files in the directory and crops them
for i = 1:N
    % read the image
    img_ref = imread([read_filepath files(i).name]);

    % this crops the image
    img_crop = img_ref(xmin:xmax, ymin:ymax);

    % this is the name of the image file without the path or extension
    [pathstr, name, ext] = fileparts(files(i).name);
    
    % this is the name of the file under which the image will be saved
    save_filename = [name ext];
    
    % this saves the cropped image to file
    imwrite(img_crop, [save_filepath save_filename]);
    
    % this displays the progress to the user
    waitbar(i/N);
end

fprintf('finished cropping images\n');

% % this displays the images before and after cropping
% T = zeros(1024,1024,3);
% T(:,:,1) = img_ref;
% T(xmin:xmax,ymin:ymax,2) = img_crop;
% figure(1)
% imshow(T)