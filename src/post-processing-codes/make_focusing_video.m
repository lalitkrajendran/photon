% program to make a video of the effect of perturbing the camera sensor on
% the image

%% load images

% image read filepath
img_read_filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/focusing/grad_x=5.00/processing/reordered-images/';

% populate list of filenames
files = dir([img_read_filepath '*.tif']);

% sort them in the correct order
% files 2 to 11 are in positive perturbation
% files 12 to 21 are in negative perturation
% file 1 is zero perturbation


% plot them one by one
figure
hold on
for i = 21:12
    img = imread([img_read_filepath files(i).name]);
    
    imshow(img)
    axis image
    
    pause(1)
end

for i = 1:11
    img = imread([img_read_filepath files(i).name]);
    
    imshow(img)
    axis image
    pause(1)
end
% save video