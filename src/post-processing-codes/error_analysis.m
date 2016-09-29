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
case_name = '100x100';

% this is the dot sizing method ('theory', 'imfindcircles', 'prana')
dot_sizing_method = 'prana';

% this is the base name of the files that contain the results to be analyzed
results_basename =  'BOS_pass2';
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
U_ref = 9.78;
V_ref = 0;

% these are the physical diameters (microns) of the dots for which the images were
% generated
% dot_diameters_physical = [10 20 30 40 50 75 100 200 300 400 500 600 700 800 900 1000];
dot_diameters_physical = [10 20 30 40 50 100 150 200 250 300 350 400 450 500];

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

if strcmp(dot_sizing_method, 'theory')
    % calculate dot diameter in the image (pixels) from theory
    dot_diameters_image = dot_diameters_physical*1e-3 * M / pixel_pitch;
else
    dot_diameters_image = measure_dot_diameter_from_image(dot_sizing_method, img_filepath, 'alternate', particle_id_filepath, particle_size_filepath);
end

%% plot results

% plot the errors against the physical dot diameter
figure_ctr = figure_ctr+1;
figure(figure_ctr);
set(gcf, 'Position', [200 200 800 800])
hold on

% bias error
subplot(3,1,1)
hold on
plot(dot_diameters_physical, err_U_bias, 'r*-');
plot(dot_diameters_physical, err_V_bias, 'b*-');
xlabel('physical diameter (microns)');
ylabel('error (pixels)');
title('bias error');
legend('U', 'V');
grid on

% random error
subplot(3,1,2)
hold on
plot(dot_diameters_physical, err_U_random, 'r*-');
plot(dot_diameters_physical, err_V_random, 'b*-');
xlabel('physical diameter (microns)');
ylabel('error (pixels)');
title('random error');
legend('U', 'V');
grid on

% rms error
subplot(3,1,3)
hold on
plot(dot_diameters_physical, err_U_rms, 'r*-');
plot(dot_diameters_physical, err_V_rms, 'b*-');
xlabel('physical diameter (microns)');
ylabel('error (pixels)');
title('rms error');
legend('U', 'V');
grid on

% save results to file
savefig(gcf, [figure_save_filepath 'errors-displacement-physical-diameter.fig']);
print(gcf, [figure_save_filepath 'errors-displacement-physical-diameter.eps'], '-depsc');
print(gcf, [figure_save_filepath 'errors-displacement-physical-diameter.png'], '-dpng');


% plot the errors against the image dot diameter
figure_ctr = figure_ctr+1;
figure(figure_ctr);
set(gcf, 'Position', [200 200 800 800])
hold on

% bias error
subplot(3,1,1)
hold on
plot(dot_diameters_image, err_U_bias, 'r*-');
plot(dot_diameters_image, err_V_bias, 'b*-');
xlabel('dot diameter in image (pixels)');
ylabel('error (pixels)');
title('bias error');
legend('U', 'V');
grid on

% random error
subplot(3,1,2)
hold on
plot(dot_diameters_image, err_U_random, 'r*-');
plot(dot_diameters_image, err_V_random, 'b*-');
xlabel('dot diameter in image (pixels)');
ylabel('error (pixels)');
title('random error');
legend('U', 'V');
grid on

% rms error
subplot(3,1,3)
hold on
plot(dot_diameters_image, err_U_rms, 'r*-');
plot(dot_diameters_image, err_V_rms, 'b*-');
xlabel('dot diameter in image (pixels)');
ylabel('error (pixels)');
title('rms error');
legend('U', 'V');
grid on

% save results to file
savefig(gcf, [figure_save_filepath 'errors-displacement-image-diameter.fig']);
print(gcf, [figure_save_filepath 'errors-displacement-image-diameter.eps'], '-depsc');
print(gcf, [figure_save_filepath 'errors-displacement-image-diameter.png'], '-dpng');
    