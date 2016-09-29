% This program loads an image, measures the diameters of the dots in two
% different ways and compares these values to the reference value

%% clear workspace and close all windows

clear all
close all
clc

figure_ctr = 0;

%% load image

% this is the path to the folder containing the raw images
filepath = '~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/50x50/300/';

% this is the name of the file to be read
filename = 'bos_pattern_image_1.tif';

% this loads the image
img = imread([filepath filename]);

%% calculate theoretical diameter of the dot

% set physical dot diameter (mm)
dp_physical = 500e-3;

% set magnification
M = 123500./700000;

% set pixel pitch (mm)
pixel_pitch = 17e-3;

% calculate dot diameter in the image (pixels)
dp_image = dp_physical * M / pixel_pitch;

% set the calculated dot diameter as the reference
dp_ref = dp_image;

%% calculate diameter from imfindcircles

% detect dots and measure their radii (pixels)
[centers, radii] = imfindcircles(img, [1 5], 'Sensitivity', 0.95);

% calculate the representative dot radius (pixels)
radius = nanmean(radii);

% calculate the diameter of the dots (pixels)
dp_circles = 2 * radius;

%% calculate diameter from prana

% add prana to the path
addpath('/home/barracuda/a/lrajendr/Software/prana/');

% add functions from sayantan to path
addpath('/home/barracuda/a/lrajendr/Projects/camera_simulation/src/post-processing-codes/from-sayantan/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ID the particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define settings for particle identification
IDmethod = {'blob','dynamic','combined'};
Data.ID.method         = IDmethod{2};
Data.ID.run            = 1;
Data.ID.v              = 10;
Data.ID.contrast_ratio = 0;
% Data.ID.s_num          = 0;
% Data.ID.s_name         = 'ID_cam2_';
Data.ID.save_dir       = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/50x50/results/particle-id/';

particleIDprops       = Data.ID;

% call function to id the particles
[ID1.p_matrix,ID1.peaks,ID1.num_p]=particle_ID(img,particleIDprops);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Particle Sizing function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define settings for particle sizing
SIZEmethod={'GEO','IWC','TPG','FPG','CFPG','LSG','CLSG'};
Data.Size.run      = 1;
Data.Size.thresh   = 10;
Data.Size.method   = 'CLSG';
Data.Size.p_area   = 0;
Data.Size.sigma    = 4;
Data.Size.errors   = 1;%str2double(Data.Size.errors);
% Data.Size.s_name   = PTV_Data.Size.savebase;
Data.Size.save_dir = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/50x50/results/particle-size/';

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

% calculate the average dot diameter (pixels)
dp_prana = nanmean(Dp1);

%% compute errors

% bias errors
bias_circles = nanmean(2*radii - dp_ref);
bias_prana = nanmean(Dp1 - dp_ref); 

% random errors
random_circles = std(2*radii);
random_prana = std(Dp1);

% rms error
rms_circles = rms(2*radii - dp_ref);
rms_prana = rms(Dp1 - dp_ref);

%% display results

% display the diameters corresponding to the reference, imfindcircles and
% prana

fprintf('----------- dot diameters (pixels) -----------\n');
fprintf('reference: %f\n', dp_ref);
fprintf('imfindcircles: %f\n', dp_circles);
fprintf('prana: %f\n', dp_prana);

% display the errors

% bias error
fprintf('----------- bias error (pixels) -----------\n');
fprintf('imfindcircles: %f\n', bias_circles);
fprintf('prana: %f\n', bias_prana);

% random error
fprintf('----------- random error (pixels) -----------\n');
fprintf('imfindcircles: %f\n', random_circles);
fprintf('prana: %f\n', random_prana);

% rms error
fprintf('----------- rms error (pixels) -----------\n');
fprintf('imfindcircles: %f\n', rms_circles);
fprintf('prana: %f\n', rms_prana);

% plot identified dots from imfindcircles
figure_ctr = figure_ctr + 1;
figure(figure_ctr)
hold on
imshow(img)
viscircles(centers, radii);
title('imfindcircles');

% plot identified dots from prana
figure_ctr = figure_ctr + 1;
figure(figure_ctr)
hold on
imagesc(img)
colormap('gray');
axis image;
% imagesc(255-im1);colormap(flipud('gray'));axis image;
plot(X1,Y1,'bo', 'Markersize', 10);
