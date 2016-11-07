% This program reads in the displacements estimated by cross-correlating
% the bos images, and compares it to the theoretical displacement. It then
% computes the mean and rms errors

%% clear workspace and close all windows

clear 
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
results_basename =  'BOS_pass3';

% this ensures all the figures are docked and displayed from a single
% window
set(0,'DefaultFigureWindowStyle','docked')

%% set theoretical displacements and other experimental parameters

% these are the theoretical displacements
% U_ref = 7.76;
U_ref = 2;
% U_ref = 5;
V_ref = 0;

% this is the magnification
M = 123500./700000;

% this is the pixel pitch (mm)
pixel_pitch = 17e-3;


%% set read  and write paths

% this is the path to the folder where the results are located
vectors_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/vectors/'];

% this is the path to the folder containing the raw images
img_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/reordered-images/'];

% this is the folder where the particle identification results will be
% stored
particle_id_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/particle-id/'];

% this is the folder where the particle sizing results will be
% stored
particle_size_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/particle-size/'];

% this is the folder where the figures will be saved
figure_save_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/plots/'];

% this is the folder where the workspace variables will be saved
workspace_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/'];

% this is the folder containing the raw images with the folder name
% denoting the dot size
dot_diameter_img_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/' case_name '/'];

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

%%

% for each f-number, populate the list of dot sizes for which images
% were generated
folder_list = dir([dot_diameter_img_filepath]);

% ignore entries like .DS_store that are not really folders
folder_list = folder_list([folder_list.isdir]);

% also remove entries like . or ..
folder_list = folder_list(~strncmpi('.',{folder_list.name},1));
    
fprintf('%d folders found\n', length(folder_list));

% initialize an array to store all the physical dot sizes for which
% images were generated
dot_diameters_physical_microns = zeros(length(folder_list),1);

% initialize array to store all the image dot sizes
% dot_diameters_image_pixels = zeros(length(folder_list), 1);

% loop through all the files to retrieve the images and note the dot
% size
for j = 1:length(folder_list)
    fprintf('processing %d out of %d dot sizes\n', j, length(folder_list));

    % get dot size from folder name in microns
    dot_diameters_physical_microns(j) = str2double(folder_list(j).name);
end

% sort the array to ensure that it contains elements in the increasing
% order of diameter. this is because the results of the piv analysis are
% stored this way.
dot_diameters_physical_microns = sort(dot_diameters_physical_microns);
% convert dot size to pixels
dot_diameters_image_pixels = dot_diameters_physical_microns * 1e-3 * M / pixel_pitch;    

%% calculate the dot diameter in the image

% % these are the names of the file to be read
% files = dir([img_filepath '*.tif']);
% 
% % this is the number of files
% N = length(files);
% fprintf('number of files in the directory: %d\n', N);
% 
% % these are the paths containing the m-files used for particle
% % identification and sizing
% addpath('/home/barracuda/a/lrajendr/Software/prana/');
% addpath('/home/barracuda/a/lrajendr/Projects/camera_simulation/src/post-processing-codes/from-sayantan/');
% 
% % this is the skip value 
% % 0 if all files are analyzed, 1 if alternate files are analyzed
% skip = 1;
% 
% % this is the number of images that will be analyzed
% num_analysis = length(1:1+skip:N);
% fprintf('number of files to be analyzed: %d\n', num_analysis);
% 
% % these are the methods that will be analyzed
% % size_methods={'GEO','IWC','TPG','FPG','CFPG','LSG','CLSG'};
% size_methods={'IWC','TPG','CLSG'};
% 
% 
% % this intializes the arrays that store the errors
% % first column - empty, 2nd column - imfindcircles, 3rd - prana
% dot_diameters_image = zeros(num_analysis,1);
% dot_diameters_circles = zeros(num_analysis, 1);
% 
% % initialize the counter for the arrays
% image_ctr = 0;
% 
% for i = 1:1+skip:N
%     % display progress
% %     display_calculation_progress(image_ctr, 1:num_analysis);
%     % this loads the image
%     img = imread([img_filepath files(i).name]);
% 
%     % consider only the center parg of the image where teh particles are in
%     % focus
%     image_margin = [200 800];
%     img = img(image_margin(1):image_margin(2), image_margin(1):image_margin(2));
% 
%     % this updates the number of images read
%     image_ctr = image_ctr + 1;
%     waitbar(image_ctr/num_analysis);
%     
%     fprintf('image %d of %d\n', image_ctr, num_analysis);
%     
%     %% calculate diameter from prana
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % identify the particles
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     % define settings for particle identification
%     IDmethod = {'blob','dynamic','combined'};
%     Data.ID.method         = IDmethod{1};
%     Data.ID.run            = 1;
%     Data.ID.v              = 10;
%     Data.ID.contrast_ratio = 0;
% 
%     particleIDprops       = Data.ID;
% 
%     % call function to id the particles
%     [ID1.p_matrix,ID1.peaks,ID1.num_p]=particle_ID(img,particleIDprops);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % size the particles
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     % define settings for particle sizing
%     Data.Size.run      = 1;
%     Data.Size.thresh   = 10;
%     Data.Size.method   = 'TPG';
%     Data.Size.p_area   = 0;
%     Data.Size.sigma    = 4;
%     Data.Size.errors   = 1;%str2double(Data.Size.errors);
% 
%     sizeprops       = Data.Size;
%     % call function to size the particles
%     [SIZE1.XYDiameter,SIZE1.mapsizeinfo,SIZE1.locxy]=particle_sizing(img,ID1.p_matrix,...
%                         ID1.num_p,sizeprops);
% 
%     % extract co-ordinate, size and intensity properties from the results             
%     X1=SIZE1.XYDiameter(:,1);
%     Y1=SIZE1.XYDiameter(:,2);                
%     Dp1=SIZE1.XYDiameter(:,3);
%     Ip1=SIZE1.XYDiameter(:,4);                 
% 
%     % remove NaN values
%     Dp1 = Dp1(~isnan(Dp1));
% 
%     % calculate the average dot diameter (pixels)
%     dot_diameters_image(image_ctr) = nanmean(Dp1);
%     
% %     %% calculate diameter from imfindcircles
% % 
% %     % detect dots and measure their radii (pixels)
% %     [centers, radii] = imfindcircles(img, [1 10], 'Sensitivity', 0.95);
% % 
% %     % calculate the representative dot radius (pixels)
% %     radius = nanmean(radii);
% % 
% %     % calculate the representative dot diameter (pixels)
% %     diameter = 2 * radius;
% % 
% %     % calculate the diameter of the dots (pixels)
% %     dot_diameters_circles(image_ctr) = diameter;
% end
% % % % calculate diameter from prana
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % identify the particles
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % % define settings for particle identification
% % % IDmethod = {'blob','dynamic','combined'};
% % % Data.ID.method         = IDmethod{2};
% % % Data.ID.run            = 1;
% % % Data.ID.v              = 10;
% % % Data.ID.contrast_ratio = 0;
% % % 
% % % particleIDprops       = Data.ID;
% % % 
% % % % call function to id the particles
% % % [ID1.p_matrix,ID1.peaks,ID1.num_p]=particle_ID(img,particleIDprops);
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % size the particles
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % % define settings for particle sizing
% % % Data.Size.run      = 1;
% % % Data.Size.thresh   = 10;
% % % Data.Size.method   = 'CLSG';
% % % Data.Size.p_area   = 0;
% % % Data.Size.sigma    = 4;
% % % Data.Size.errors   = 1;%str2double(Data.Size.errors);
% % % 
% % % sizeprops       = Data.Size;
% % % 
% % % % call function to size the particles
% % % [SIZE1.XYDiameter,SIZE1.mapsizeinfo,SIZE1.locxy]=particle_sizing(img,ID1.p_matrix,...
% % %                     ID1.num_p,sizeprops);
% % % 
% % % % extract co-ordinate, size and intensity properties from the results             
% % % X1=SIZE1.XYDiameter(:,1);
% % % Y1=SIZE1.XYDiameter(:,2);                
% % % Dp1=SIZE1.XYDiameter(:,3);
% % % Ip1=SIZE1.XYDiameter(:,4);                 
% % % 
% % % % remove NaN values
% % % Dp1 = Dp1(~isnan(Dp1));
% % % 
% % % % calculate the average dot diameter (pixels)
% % % dot_diameters_image(j) = nanmean(Dp1);
% % % if strcmp(dot_sizing_method, 'theory')
% % %     % calculate dot diameter in the image (pixels) from theory
% % %     dot_diameters_image = dot_diameters_physical*1e-3 * M / pixel_pitch;
% % % else
% % %     dot_diameters_image = measure_dot_diameter_from_image(dot_sizing_method, img_filepath, 'alternate', particle_id_filepath, particle_size_filepath);
% % % end
% 
% % load diameters from file
% % previous_workspace = load([workspace_filepath 'workspace.mat']);
% % dot_diameters_image = previous_workspace.dot_diameters_image;
% % dot_diameters_physical = previous_workspace.dot_diameters_physical;
% 
% % % save worskpace to file
% % save([workspace_filepath 'workspace.mat']);

%% plot results

% % plot the errors against the physical dot diameter
% figure_ctr = figure_ctr+1;
% figure(figure_ctr);
% set(gcf, 'Position', [200 200 900 500])
% hold on
% 
% % bias error
% subplot(1,3,1)
% hold on
% plot(dot_diameters_image, err_U_bias, 'r*-');
% plot(dot_diameters_image, err_V_bias, 'b*-');
% xlim([0 max(dot_diameters_image)*1.2])
% xlabel('dot diameter (pixels)');
% ylabel('error (pixels)');
% title('bias error');
% legend('\Delta x', '\Delta y');
% grid on
% 
% % random errorP
% subplot(1,3,2)
% hold on
% plot(dot_diameters_image, err_U_random, 'r*-');
% plot(dot_diameters_image, err_V_random, 'b*-');
% xlim([0 max(dot_diameters_image)*1.2])
% xlabel('dot diameter (pixels)');
% ylabel('error (pixels)');
% title('random error');
% legend('\Delta x', '\Delta y');
% grid on
% 
% % rms error
% subplot(1,3,3)
% hold on
% plot(dot_diameters_image, err_U_rms, 'r*-');
% plot(dot_diameters_image, err_V_rms, 'b*-');
% xlim([0 max(dot_diameters_image)*1.2])
% xlabel('dot diameter (pixels)');
% ylabel('error (pixels)');
% title('total error');
% legend('\Delta x', '\Delta y', 'location', 'Northwest');
% grid on
% 
% % save results to file
% savefig(gcf, [figure_save_filepath 'errors-displacement-image-diameter.fig']);
% print(gcf, [figure_save_filepath 'errors-displacement-image-diameter.eps'], '-depsc');
% print(gcf, [figure_save_filepath 'errors-displacement-image-diameter.png'], '-dpng');

figure_ctr = figure_ctr+1;
figure(figure_ctr);
set(gcf, 'Position', [200 200 900 500])

% bias error
subplot(1,3,1)
hold on
plot(dot_diameters_image_pixels, err_U_bias, 'r*-');
plot(dot_diameters_image_pixels, err_V_bias, 'b*-');
grid on
xlim([0 max(dot_diameters_image_pixels)*1.2])
xlabel('dot diameter (pixels)');
ylabel('error (pixels)');
legend('\Delta x', '\Delta y', 'location', 'Northwest');
title('bias error')

% random error
subplot(1,3,2)
hold on
plot(dot_diameters_image_pixels, err_U_random, 'ro-');
plot(dot_diameters_image_pixels, err_V_random, 'bo-');
grid on
xlim([0 max(dot_diameters_image_pixels)*1.2])
xlabel('dot diameter (pixels)');
ylabel('error (pixels)');
legend('\Delta x', '\Delta y', 'location', 'Northwest');
title('random error')

% total error
subplot(1,3,3)
hold on
plot(dot_diameters_image_pixels, err_U_rms, 'rs-');
plot(dot_diameters_image_pixels, err_V_rms, 'bs-');
grid on
xlim([0 max(dot_diameters_image_pixels)*1.2])
xlabel('dot diameter (pixels)');
ylabel('error (pixels)');
legend('\Delta x', '\Delta y', 'location', 'Northwest');
title('total error')

% save results to file
savefig(gcf, [figure_save_filepath 'errors-displacement-physical-diameter.fig']);
print(gcf, [figure_save_filepath 'errors-displacement-physical-diameter.eps'], '-depsc');
print(gcf, [figure_save_filepath 'errors-displacement-physical-diameter.png'], '-dpng');

% PLOT RANDOM ERROR SEPARATELY

figure_ctr = figure_ctr+1;
figure(figure_ctr);
hold on
plot(dot_diameters_image_pixels, err_U_random, 'ro-');
plot(dot_diameters_image_pixels, err_V_random, 'bo-');
grid on
xlabel('dot diameter (pixels)');
ylabel('error (pixels)');
legend('\Delta x', '\Delta y', 'location', 'Northwest');
title('random error')
legend('\Delta x', '\Delta y', 'location', 'Northwest');

% save results to file
plot_filename = 'random_error';
savefig([figure_save_filepath plot_filename '.fig'])
print([figure_save_filepath plot_filename '.png'], '-dpng');
print([figure_save_filepath plot_filename '.eps'], '-deps');


% figure_ctr = figure_ctr+1;
% figure(figure_ctr);
% set(gcf, 'Position', [200 200 700 500])
% hold on
% 
% % bias error
% plot(dot_diameters_image, err_U_bias, 'r*-');
% plot(dot_diameters_image, err_V_bias, 'b*-');
% 
% % random error
% plot(dot_diameters_image, err_U_random, 'ro-');
% plot(dot_diameters_image, err_V_random, 'bo-');
% 
% % total error
% plot(dot_diameters_image, err_U_rms, 'rs-');
% plot(dot_diameters_image, err_V_rms, 'bs-');
% 
% grid on
% 
% legend('bias-U', 'bias-V', 'random-U', 'random-V', 'total-U', 'total-V'); 
% xlabel('image diameter (pixels)');
% ylabel('error (pixels)');
% title('effect of dot diameter on error')

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

%% save worskpace to file

% this is the name of the current script
script_name_full = mfilename('fullpath');
[pathstr, script_name, ext] = fileparts(script_name_full);
save([workspace_filepath script_name '.mat']);

