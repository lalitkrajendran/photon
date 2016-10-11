% the purpose is sto study how the various dot sizing schemes in prana work

clear all
close all
clc

case_name = '50x50-f11';

%% set read and write paths

% this is the path to the folder containing the raw images
img_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/reordered-images/'];

% this is the folder where the particle sizing results will be
% stored
particle_size_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/particle-size/'];

i = 10;
image_ctr = 0;
num_analysis = 1;
%% load images

% these are the names of the file to be read
files = dir([img_filepath '*.tif']);

% this is the number of files
N = length(files);
fprintf('number of files in the directory: %d\n', N);

% these are the methods that will be analyzed
size_methods={'GEO','IWC','TPG','FPG','CFPG','LSG','CLSG'};
% size_methods={'FPG','CFPG'};

% this loads the image
img = imread([img_filepath files(i).name]);

% this updates the number of images read
image_ctr = image_ctr + 1;

fprintf('image %d of %d\n', image_ctr, num_analysis);

% these are the paths containing the m-files used for particle
% identification and sizing
addpath('/home/barracuda/a/lrajendr/Software/prana/');
addpath('/home/barracuda/a/lrajendr/Projects/camera_simulation/src/post-processing-codes/from-sayantan/');

% consider only the center parg of the image where teh particles are in
% focus
image_margin = [300 700];
img = img(image_margin(1):image_margin(2), image_margin(1):image_margin(2));

%% estimate diameter

% these are the physical diameters (mm) of the dots for which the images were
% generated
% dot_diameters_physical = [50 75 100 200 300 400 500 600 700 800 900 1000];
% dot_diameters_physical = [10 20 30 40 50 100 150 200 250 300 350 400 450 500];
dot_diameters_physical = 50;

% this is the magnification
M = 123500./700000;

% this is the pixel pitch (mm)
pixel_pitch = 17e-3;

% calculate dot diameter in the image (pixels) from theory
dot_diameters_theory = dot_diameters_physical*1e-3 * M / pixel_pitch;
%% calculate diameter from prana

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% identify the particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define settings for particle identification
IDmethod = {'blob','dynamic','combined'};
Data.ID.method         = IDmethod{2};
Data.ID.run            = 1;
Data.ID.v              = 10;
Data.ID.contrast_ratio = 0;

particleIDprops       = Data.ID;

% call function to id the particles
[ID1.p_matrix,ID1.peaks,ID1.num_p]=particle_ID(img,particleIDprops);
% Data.ID.save_dir       = particle_id_filepath;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% size the particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define settings for particle sizing
%     SIZEmethod={'GEO','IWC','TPG','FPG','CFPG','LSG','CLSG'};
Data.Size.run      = 1;
Data.Size.thresh   = 10;
%     Data.Size.method   = 'CLSG';
Data.Size.p_area   = 0;
Data.Size.sigma    = 4;
Data.Size.errors   = 1;%str2double(Data.Size.errors);

dot_diameters_image = zeros(1,length(size_methods));

for j = 1:length(size_methods)
    fprintf(['running method ' size_methods{j} '\n']);
    Data.Size.method = size_methods{j};
%     Data.Size.save_dir = [particle_size_filepath size_methods{j} '/'];
% 
%     if ~exist(Data.Size.save_dir, 'dir')
%         mkdir(Data.Size.save_dir);
%     end


    sizeprops       = Data.Size;
    % call function to size the particles
    [SIZE1.XYDiameter,SIZE1.mapsizeinfo,SIZE1.locxy]=particle_sizing(img,ID1.p_matrix,...
                        ID1.num_p,sizeprops);

    % extract co-ordinate, size and intensity properties from the results             
    X1=SIZE1.XYDiameter(:,1);
    Y1=SIZE1.XYDiameter(:,2);                
    Dp1=SIZE1.XYDiameter(:,3);
    Ip1=SIZE1.XYDiameter(:,4);                 

    % remove NaN values
    Dp1 = Dp1(~isnan(Dp1));
    dot_diameters_raw{j} = Dp1;
    % calculate the average dot diameter (pixels)
    dot_diameters_image(image_ctr,j) = nanmean(Dp1);
end

%% calculate diameter from imfindcircles

% detect dots and measure their radii (pixels)
[centers, radii] = imfindcircles(img, [1 10], 'Sensitivity', 0.95);

% calculate the representative dot radius (pixels)
radius = nanmean(radii);

% calculate the representative dot diameter (pixels)
diameter = 2 * radius;

% calculate the diameter of the dots (pixels)
dot_diameters_circles = diameter;

%% display results

fprintf('theory: %f\n', dot_diameters_theory);
fprintf('circles: %f\n', dot_diameters_circles);
for j = 1:length(size_methods)
    fprintf(['prana - ' size_methods{j} ': %f \n'], dot_diameters_image(1,j));
end
