% This program reads in the displacements estimated by cross-correlating
% the bos images, and compares it to the theoretical displacement. It then
% computes the mean and rms errors

%% clear workspace and close all windows

clear all
close all
clc

figure_ctr = 0;

%% overall settings

% this specifies whethere the figures will be saved or not
printfig = false;

% this is the case name
case_name = '150x150-f16-disp2';

% this is the dot sizing method ('theory', 'imfindcircles', 'prana')
dot_sizing_method = 'prana';

% this is the base name of the files that contain the results to be analyzed
results_basename =  'BOS_pass1';
%% set read  and write paths

% this is the path to the folder where the results are located
vectors_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/vectors/'];

% this is the path to the folder containing the raw images
img_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/reordered-images/'];

% this is the folder where the particle identification results will be
% stored
particle_id_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/particle-id/'];

% this is the folder where the particle sizing results will be
% stored
particle_size_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/particle-size/'];

% this is the folder where the figures will be saved
figure_save_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/plots/'];

% this is the folder where the workspace variables will be saved
workspace_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/'];

% this creates the write directories if they are not already present
if ~exist(particle_id_filepath,'dir')
    mkdir(particle_id_filepath);
end

if ~exist(particle_size_filepath,'dir')
    mkdir(particle_size_filepath);
end

if ~exist(figure_save_filepath,'dir')
    mkdir(figure_save_filepath);
end

%% set theoretical displacements and other experimental parameters

% these are the theoretical displacements
U_ref = 4.87;
V_ref = 0;

% these are the physical diameters (microns) of the dots for which the images were
% generated
% dot_diameters_physical = [10 20 30 40 50 75 100 200 300 400 500 600 700 800 900 1000];
dot_diameters_physical = [10 20 30 40 50 100 150 200 250 300 350 400 450];

% this is the magnification
M = 123500./700000;

% this is the pixel pitch (mm)
pixel_pitch = 17e-3;

%% load displacement data and compute errors

% this populates all the files in the directory
filenames = dir([vectors_filepath results_basename '*.mat']);

% filenames = filenames(1:end-1);

% this is the number of files present
num_files = length(filenames);
fprintf('No. of files: %d\n', num_files);

%% compute errors
% this is the array containing the bias error
err_U_bias = zeros(1, num_files);
err_V_bias = zeros(1, num_files);

% this is the array containing the random error
err_U_random = zeros(1, num_files);
err_V_random = zeros(1, num_files);

% this is the array containing the rms error
err_U_rms = zeros(1, num_files);
err_V_rms = zeros(1, num_files);

for i = 1:length(filenames)
    % load data
    data = load([vectors_filepath filenames(i).name]);
    
    % calculate errors
    [err_U_bias(i), err_U_random(i), err_U_rms(i)] = compute_errors(data.U(:), U_ref);
    [err_V_bias(i), err_V_random(i), err_V_rms(i)] = compute_errors(data.V(:), V_ref);
end

%% calculate the dot diameter in the image

% these are the names of the file to be read
files = dir([img_filepath '*.tif']);

% this is the number of files
N = length(files);
fprintf('number of files in the directory: %d\n', N);

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

% these are the methods that will be analyzed
% size_methods={'GEO','IWC','TPG','FPG','CFPG','LSG','CLSG'};
size_methods={'IWC','TPG','CLSG'};


% this intializes the arrays that store the errors
% first column - empty, 2nd column - imfindcircles, 3rd - prana
dot_diameters_image = zeros(num_analysis,1);
dot_diameters_circles = zeros(num_analysis, 1);

% initialize the counter for the arrays
image_ctr = 0;

for i = 1:1+skip:N
    % display progress
%     display_calculation_progress(image_ctr, 1:num_analysis);
    % this loads the image
    img = imread([img_filepath files(i).name]);

    % consider only the center parg of the image where teh particles are in
    % focus
    image_margin = [200 800];
    img = img(image_margin(1):image_margin(2), image_margin(1):image_margin(2));

    % this updates the number of images read
    image_ctr = image_ctr + 1;
    waitbar(image_ctr/num_analysis);
    
    fprintf('image %d of %d\n', image_ctr, num_analysis);
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % size the particles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define settings for particle sizing
    Data.Size.run      = 1;
    Data.Size.thresh   = 10;
    Data.Size.method   = 'CLSG';
    Data.Size.p_area   = 0;
    Data.Size.sigma    = 4;
    Data.Size.errors   = 1;%str2double(Data.Size.errors);

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
    dot_diameters_image(image_ctr) = nanmean(Dp1);
    
    %% calculate diameter from imfindcircles

    % detect dots and measure their radii (pixels)
    [centers, radii] = imfindcircles(img, [1 10], 'Sensitivity', 0.95);

    % calculate the representative dot radius (pixels)
    radius = nanmean(radii);

    % calculate the representative dot diameter (pixels)
    diameter = 2 * radius;

    % calculate the diameter of the dots (pixels)
    dot_diameters_circles(image_ctr) = diameter;
end
% % calculate diameter from prana
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % identify the particles
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % define settings for particle identification
% IDmethod = {'blob','dynamic','combined'};
% Data.ID.method         = IDmethod{2};
% Data.ID.run            = 1;
% Data.ID.v              = 10;
% Data.ID.contrast_ratio = 0;
% 
% particleIDprops       = Data.ID;
% 
% % call function to id the particles
% [ID1.p_matrix,ID1.peaks,ID1.num_p]=particle_ID(img,particleIDprops);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % size the particles
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % define settings for particle sizing
% Data.Size.run      = 1;
% Data.Size.thresh   = 10;
% Data.Size.method   = 'CLSG';
% Data.Size.p_area   = 0;
% Data.Size.sigma    = 4;
% Data.Size.errors   = 1;%str2double(Data.Size.errors);
% 
% sizeprops       = Data.Size;
% 
% % call function to size the particles
% [SIZE1.XYDiameter,SIZE1.mapsizeinfo,SIZE1.locxy]=particle_sizing(img,ID1.p_matrix,...
%                     ID1.num_p,sizeprops);
% 
% % extract co-ordinate, size and intensity properties from the results             
% X1=SIZE1.XYDiameter(:,1);
% Y1=SIZE1.XYDiameter(:,2);                
% Dp1=SIZE1.XYDiameter(:,3);
% Ip1=SIZE1.XYDiameter(:,4);                 
% 
% % remove NaN values
% Dp1 = Dp1(~isnan(Dp1));
% 
% % calculate the average dot diameter (pixels)
% dot_diameters_image(j) = nanmean(Dp1);
% if strcmp(dot_sizing_method, 'theory')
%     % calculate dot diameter in the image (pixels) from theory
%     dot_diameters_image = dot_diameters_physical*1e-3 * M / pixel_pitch;
% else
%     dot_diameters_image = measure_dot_diameter_from_image(dot_sizing_method, img_filepath, 'alternate', particle_id_filepath, particle_size_filepath);
% end


% save worskpace to file
save([workspace_filepath 'workspace.mat']);

%% plot results

% plot the errors against the physical dot diameter
figure_ctr = figure_ctr+1;
figure(figure_ctr);
set(gcf, 'Position', [200 200 900 500])
hold on

% bias error
subplot(1,3,1)
hold on
plot(dot_diameters_image, err_U_bias, 'r*-');
plot(dot_diameters_image, err_V_bias, 'b*-');
xlim([0 max(dot_diameters_image)*1.2])
xlabel('dot diameter (pixels)');
ylabel('error (pixels)');
title('bias error');
legend('\Delta x', '\Delta y');
grid on

% random error
subplot(1,3,2)
hold on
plot(dot_diameters_image, err_U_random, 'r*-');
plot(dot_diameters_image, err_V_random, 'b*-');
xlim([0 max(dot_diameters_image)*1.2])
xlabel('dot diameter (pixels)');
ylabel('error (pixels)');
title('random error');
legend('\Delta x', '\Delta y');
grid on

% rms error
subplot(1,3,3)
hold on
plot(dot_diameters_image, err_U_rms, 'r*-');
plot(dot_diameters_image, err_V_rms, 'b*-');
xlim([0 max(dot_diameters_image)*1.2])
xlabel('dot diameter (pixels)');
ylabel('error (pixels)');
title('total error');
legend('\Delta x', '\Delta y', 'location', 'Northwest');
grid on

% save results to file
savefig(gcf, [figure_save_filepath 'errors-displacement-image-diameter.fig']);
print(gcf, [figure_save_filepath 'errors-displacement-image-diameter.eps'], '-depsc');
print(gcf, [figure_save_filepath 'errors-displacement-image-diameter.png'], '-dpng');

% figure
% % bias error
% plot(dot_diameters_physical, err_U_bias, 'r*-');
% plot(dot_diameters_physical, err_V_bias, 'b*-');
% 
% % random error
% plot(dot_diameters_physical, err_U_random, 'ro-');
% plot(dot_diameters_physical, err_V_random, 'bo-');
% 
% % total error
% plot(dot_diameters_physical, err_U_rms, 'rs-');
% plot(dot_diameters_physical, err_V_rms, 'bs-');
% 
% grid on
% 
% legend('bias-U', 'bias-V', 'random-U', 'random-V', 'total-U', 'total-V'); 
% xlabel('physical diameter (microns)');
% ylabel('error (pixels)');
% title('effect of dot diameter on error')

% % save results to file
% savefig(gcf, [figure_save_filepath 'errors-displacement-physical-diameter.fig']);
% print(gcf, [figure_save_filepath 'errors-displacement-physical-diameter.eps'], '-depsc');
% print(gcf, [figure_save_filepath 'errors-displacement-physical-diameter.png'], '-dpng');

figure_ctr = figure_ctr+1;
figure(figure_ctr);
set(gcf, 'Position', [200 200 700 500])
hold on

% bias error
plot(dot_diameters_image, err_U_bias, 'r*-');
plot(dot_diameters_image, err_V_bias, 'b*-');

% random error
plot(dot_diameters_image, err_U_random, 'ro-');
plot(dot_diameters_image, err_V_random, 'bo-');

% total error
plot(dot_diameters_image, err_U_rms, 'rs-');
plot(dot_diameters_image, err_V_rms, 'bs-');

grid on

legend('bias-U', 'bias-V', 'random-U', 'random-V', 'total-U', 'total-V'); 
xlabel('image diameter (pixels)');
ylabel('error (pixels)');
title('effect of dot diameter on error')

% % plot the errors against the image dot diameter
% figure_ctr = figure_ctr+1;
% figure(figure_ctr);
% set(gcf, 'Position', [200 200 800 800])
% hold on
% 
% % bias error
% subplot(3,1,1)
% hold on
% plot(dot_diameters_image, err_U_bias, 'r*-');
% plot(dot_diameters_image, err_V_bias, 'b*-');
% xlabel('dot diameter in image (pixels)');
% ylabel('error (pixels)');
% title('bias error');
% legend('U', 'V');
% grid on
% 
% % random error
% subplot(3,1,2)
% hold on
% plot(dot_diameters_image, err_U_random, 'r*-');
% plot(dot_diameters_image, err_V_random, 'b*-');
% xlabel('dot diameter in image (pixels)');
% ylabel('error (pixels)');
% title('random error');
% legend('U', 'V');
% grid on
% 
% % rms error
% subplot(3,1,3)
% hold on
% plot(dot_diameters_image, err_U_rms, 'r*-');
% plot(dot_diameters_image, err_V_rms, 'b*-');
% xlabel('dot diameter in image (pixels)');
% ylabel('error (pixels)');
% title('rms error');
% legend('U', 'V');
% grid on
% 
% % save results to file
% savefig(gcf, [figure_save_filepath 'errors-displacement-image-diameter.fig']);
% print(gcf, [figure_save_filepath 'errors-displacement-image-diameter.eps'], '-depsc');
% print(gcf, [figure_save_filepath 'errors-displacement-image-diameter.png'], '-dpng');
    