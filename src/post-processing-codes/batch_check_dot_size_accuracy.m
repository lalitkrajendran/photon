% This program loads a set of images, and for each image, it measures the 
% diameters of the dots in two different ways and compares these values 
% to the reference value

% Author - Lalit Rajendran
% Date - 27th September, 2016

%% clear workspace and close all windows

clear all
close all
clc

figure_ctr = 0;

case_name = '150x150-f16-disp5';

%% set read and write paths

% this is the path to the folder containing the raw images
img_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/cropped-images/'];

% this is the folder where the particle identification results will be
% stored
particle_id_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/particle-id/'];

% this is the folder where the particle sizing results will be stored
particle_size_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/particle-size/'];

% this is the folder where the figures will be saved
figure_save_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/plots/'];

% this is the folder where the workspace variables will be saved
workspace_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/'];

%% load images

% these are the names of the file to be read
files = dir([img_filepath '*.tif']);

% this the numbe of files
N = length(files);
fprintf('number of files in the directory: %d\n', N);

%% batch process images

% these are the physical diameters (mm) of the dots for which the images were
% generated
% dot_diameters_physical = [50 75 100 200 300 400 500 600 700 800 900 1000];
% dot_diameters_physical = [10 20 30 40 50 60 70 80 90 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900];
dot_diameters_physical = [10 20 30 40 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950];

% this is the magnification
M = 123500./700000;

% this is the pixel pitch (mm)
pixel_pitch = 17e-3;

% these are the paths containing the m-files used for particle
% identification and sizing
addpath('/home/barracuda/a/lrajendr/Software/prana/');
addpath('/home/barracuda/a/lrajendr/Projects/camera_simulation/src/post-processing-codes/from-sayantan/');

% this is the skip value 
% 0 if all files are analyzed, 1 if alternate files are analyzed
skip = 1;

% this is the number of images that will be analyzed
num_analysis = length(1:1+skip:N);
fprintf('number of files to be analyzed: %d\n', num_analysis);

% this intializes the arrays that store the errors
% first column - empty, 2nd column - imfindcircles, 3rd - prana
bias_error = zeros(num_analysis,3);
random_error = zeros(num_analysis,3);
rms_error = zeros(num_analysis,3);

% intialize arrays that store the diameter values
% first column - ref, 2nd column - imfindcircles, 3rd - prana
dp = zeros(num_analysis,3);

% initialize the counter for the arrays
image_ctr = 0;
for i = 1:1+skip:N
    % display progress
    display_calculation_progress(image_ctr, 1:num_analysis);
    
    % this loads the image
    img = imread([img_filepath files(i).name]);

    % this updates the number of images read
    image_ctr = image_ctr + 1;

    
    %% calculate theoretical diameter of the dot

    % extract physical dot diameter (mm)
    dp_physical = dot_diameters_physical(image_ctr)*1e-3;

    % calculate dot diameter in the image (pixels)
    dp_image = dp_physical * M / pixel_pitch;

    % set the calculated dot diameter as the reference
    dp_ref = dp_image;
    
    % store diameter in array
    dp(image_ctr,1) = dp_image;

    %% calculate diameter from imfindcircles

    % detect dots and measure their radii (pixels)
    [centers, radii] = imfindcircles(img, [1 5], 'Sensitivity', 0.95);

    % calculate the representative dot radius (pixels)
    radius = nanmean(radii);

    % calculate the diameter of the dots (pixels)
    dp_circles = 2 * radius;

    % store diameter in array
    dp(image_ctr,2) = dp_circles;

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
    % Data.ID.s_num          = 0;
    % Data.ID.s_name         = 'ID_cam2_';
    Data.ID.save_dir       = particle_id_filepath;

    particleIDprops       = Data.ID;

    % call function to id the particles
    [ID1.p_matrix,ID1.peaks,ID1.num_p]=particle_ID(img,particleIDprops);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % size the particles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define settings for particle sizing
    SIZEmethod={'GEO','IWC','TPG','FPG','CFPG','LSG','CLSG'};
    Data.Size.run      = 1;
    Data.Size.thresh   = 10;
%     Data.Size.method   = 'CLSG';
    Data.Size.method = 'IWC';
    Data.Size.p_area   = 0;
    Data.Size.sigma    = 4;
    Data.Size.errors   = 1;%str2double(Data.Size.errors);
    % Data.Size.s_name   = PTV_Data.Size.savebase;
    Data.Size.save_dir = particle_size_filepath;

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
    
    % store diameter in array
    dp(image_ctr,3) = dp_prana;

    %% compute errors

%     % bias errors
%     bias_error(image_ctr,2) = nanmean(2*radii - dp_ref);
%     bias_error(image_ctr,3) = nanmean(Dp1 - dp_ref); 
% 
%     % random errors
%     random_error(image_ctr,2) = std(2*radii);
%     random_error(image_ctr,3) = std(Dp1);
% 
%     % rms error
%     rms_error(image_ctr,2) = rms(2*radii - dp_ref);
%     rms_error(image_ctr,3) = rms(Dp1 - dp_ref);
% 
%     [bias_error(image_ctr,2), random_error(image_ctr,2), rms_error(image_ctr,2)] ...
%         = compute_errors(dp_circles, dp_ref);
%     
%     [bias_error(image_ctr,3), random_error(image_ctr,3), rms_error(image_ctr,3)] ...
%         = compute_errors(dp_prana, dp_ref);    

end
%% display results

% plot the diameters
figure_ctr = figure_ctr+1;
figure(figure_ctr);
hold on
plot(dot_diameters_physical, dp(:,1), 'k*-');
plot(dot_diameters_physical, dp(:,2), 'r*-');
plot(dot_diameters_physical, dp(:,3), 'b*-');
xlabel('physical (microns)');
ylabel('image (pixels)');
title('Dot diameters');
legend('theory', 'imfindcircles', 'prana');

% save results to file
savefig(gcf, [figure_save_filepath 'diameters.fig']);
print(gcf, [figure_save_filepath 'diameters.eps'], '-depsc');
print(gcf, [figure_save_filepath 'diameters.png'], '-dpng');

% plot the errors
figure_ctr = figure_ctr+1;
figure(figure_ctr);
hold on

% bias error
subplot(3,1,1)
hold on
plot(dot_diameters_physical, bias_error(:,2), 'r*-');
plot(dot_diameters_physical, bias_error(:,3), 'b*-');
xlabel('physical diameter (microns)');
ylabel('error (pixels)');
title('bias error');
legend('imfindcircles', 'prana');

% random error
subplot(3,1,2)
hold on
plot(dot_diameters_physical, random_error(:,2), 'r*-');
plot(dot_diameters_physical, random_error(:,3), 'b*-');
xlabel('physical diameter (microns)');
ylabel('error (pixels)');
title('random error');
legend('imfindcircles', 'prana');

% rms error
subplot(3,1,3)
hold on
plot(dot_diameters_physical, rms_error(:,2), 'r*-');
plot(dot_diameters_physical, rms_error(:,3), 'b*-');
xlabel('physical diameter (microns)');
ylabel('error (pixels)');
title('rms error');
legend('imfindcircles', 'prana');

% save results to file
savefig(gcf, [figure_save_filepath 'errors-diameter.fig']);
print(gcf, [figure_save_filepath 'errors-diameter.eps'], '-depsc');
print(gcf, [figure_save_filepath 'errors-diameter.png'], '-dpng');

%% save worskpace to file

% this is the name of the current script
script_name_full = mfilename('fullpath');
[pathstr, script_name, ext] = fileparts(script_name_full);
save([workspace_filepath script_name '.mat']);

