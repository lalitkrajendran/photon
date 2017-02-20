% This program reads in the displacements estimated by cross-correlating
% the bos images, and compares it to the theoretical displacement. It then
% computes the mean and rms errors

%% clear workspace and close all windows

clear 
close all
clc

figure_ctr = 0;

starting_index = 1;
%% overall settings

% this specifies whethere the figures will be saved or not
printfig = false;

% these are teh seeding densities considered (particles / 32x32 pix window)
seeding_densities = [1, 4, 8, 10, 12, 16];

% this is the dot sizing method ('theory', 'imfindcircles', 'prana')
dot_sizing_method = 'prana';

% this is the base name of the files that contain the results to be analyzed
results_basename =  'bos_pass3_';

% this is the number of .mat results files corresponding to a single dot
% size
num_files_per_seeding = 20;

% this ensures all the figures are docked and displayed from a single
% window
set(0,'DefaultFigureWindowStyle','docked')

%% set theoretical displacements and other experimental parameters

% these are the theoretical displacements
% U_ref = 7.76;
U_ref = 2.0;
% U_ref = 5;
V_ref = 0;

% this is the magnification
M = 0.111; %123500./700000;

% this is the pixel pitch (mm)
pixel_pitch = 17e-3;

% this is the aperture f number
f_number = 16;

% this is teh wavelength used to calculate the diffraction limited diameter
% (mm)
wavelength = 500e-6;

%% do the following steps for all cases
num_cases = length(seeding_densities);
fprintf('number of cases: %d\n', num_cases);

% this is the array containing the bias error
err_U_bias = zeros(1, num_cases);
err_V_bias = zeros(1, num_cases);

% this is the array containing the random error
err_U_random = zeros(1, num_cases);
err_V_random = zeros(1, num_cases);

% this is the array containing the rms error
err_U_rms = zeros(1, num_cases);
err_V_rms = zeros(1, num_cases);

for case_index = 1:num_cases
    
    % this is the case name
    case_name = num2str(seeding_densities(case_index));
    
    %% set read  and write paths

    % this is the path to the folder where the results are located
    vectors_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/processing/' case_name '/results/vectors/'];

    % this is the path to the folder containing the raw images
    img_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/processing/' case_name '/reordered-images/'];

    % this is the folder where the figures will be saved
    figure_save_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/plots/'];

    % this is the folder where the workspace variables will be saved
    workspace_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/processing/'];

    % this creates the write directories if they are not already present
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
    
    err_U_bias_single = zeros(1,num_files);
    err_V_bias_single = zeros(1,num_files);
    err_U_random_single = zeros(1,num_files);
    err_V_random_single = zeros(1,num_files);
    err_U_rms_single = zeros(1,num_files);
    err_V_rms_single = zeros(1,num_files);

    % this loops through all the files in the foder
    for file_index = 1:num_files
        % load data
        data = load([vectors_filepath filenames(file_index).name]);

        % account for sign convention used in the code
        data.U(:) = -data.U(:);
        
        % calculate errors
        [err_U_bias_single(file_index), err_U_random_single(file_index), ...
            err_U_rms_single(file_index)] = compute_errors(data.U(:)', U_ref);
        [err_V_bias_single(file_index), err_V_random_single(file_index), ...
            err_V_rms_single(file_index)] = compute_errors(data.V(:)', V_ref);
        
    end

    check_convergence
    
    err_U_bias(case_index) = mean(err_U_bias_single);
    err_V_bias(case_index) = mean(err_V_bias_single);
    err_U_random(case_index) = mean(err_U_random_single);
    err_V_random(case_index) = mean(err_V_random_single);
    err_U_rms(case_index) = mean(err_U_rms_single);
    err_V_rms(case_index) = mean(err_V_rms_single);
    
end

%% plot results

figure_ctr = figure_ctr+1;
figure(figure_ctr);
% set(gcf, 'Position', [200 200 900 500])
% set(gca, 'fontsize', 20)

% bias error
subplot(1,3,1)
hold on
plot(seeding_densities, err_U_bias, 'r*-');
plot(seeding_densities, err_V_bias, 'b*-');
grid on
xlim([0 max(seeding_densities)*1.2])
xlabel('dots (32x32 pix)');
ylabel('error (pixels)');
legend('\Delta x', '\Delta y', 'location', 'Northwest');
title('bias error')
set(gca, 'fontsize', 14)

% random error
subplot(1,3,2)
hold on
plot(seeding_densities, err_U_random, 'r*-');
plot(seeding_densities, err_V_random, 'b*-');
grid on
xlim([0 max(seeding_densities)*1.2])
xlabel('dots (32x32 pix)');
ylabel('error (pixels)');
title('random error')
set(gca, 'fontsize', 14)

% total error
subplot(1,3,3)
hold on
plot(seeding_densities, err_U_rms, 'r*-');
plot(seeding_densities, err_V_rms, 'b*-');
grid on
xlim([0 max(seeding_densities)*1.2])
xlabel('dots (32x32 pix)');
ylabel('error (pixels)');
title('total error')
set(gca, 'fontsize', 14)

% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'all-errors-seeding');

% % PLOT RANDOM ERROR SEPARATELY
% 
% figure_ctr = figure_ctr+1;
% figure(figure_ctr);
% semilogy(dot_diameters_image_pixels(starting_index:end), err_U_random(starting_index:end), 'ro-');
% hold on
% semilogy(dot_diameters_image_pixels(starting_index:end), err_V_random(starting_index:end), 'bo-');
% grid on
% xlabel('dot diameter (pixels)');
% ylabel('error (pixels)');
% legend('\Delta x', '\Delta y', 'location', 'Northwest');
% title('random error')
% legend('\Delta x', '\Delta y', 'location', 'Northwest');
% set(gca, 'fontsize', 14)
% 
% % save results to file
% plot_filename = 'random_error';
% % savefig([figure_save_filepath plot_filename '.fig'])
% % print([figure_save_filepath plot_filename '.png'], '-dpng');
% % print([figure_save_filepath plot_filename '.eps'], '-deps');
% save_figure_to_file(gcf, figure_save_filepath, plot_filename);
% 
% % figure_ctr = figure_ctr+1;
% % figure(figure_ctr);
% % set(gcf, 'Position', [200 200 700 500])
% % hold on
% % 
% % % bias error
% % plot(dot_diameters_image, err_U_bias, 'r*-');
% % plot(dot_diameters_image, err_V_bias, 'b*-');
% % 
% % % random error
% % plot(dot_diameters_image, err_U_random, 'ro-');
% % plot(dot_diameters_image, err_V_random, 'bo-');
% % 
% % % total error
% % plot(dot_diameters_image, err_U_rms, 'rs-');
% % plot(dot_diameters_image, err_V_rms, 'bs-');
% % 
% % grid on
% % 
% % legend('bias-U', 'bias-V', 'random-U', 'random-V', 'total-U', 'total-V'); 
% % xlabel('image diameter (pixels)');
% % ylabel('error (pixels)');
% % title('effect of dot diameter on error')
% 
% % % plot the errors against the image dot diameter
% % figure_ctr = figure_ctr+1;
% % figure(figure_ctr);
% % set(gcf, 'Position', [200 200 800 800])
% % hold on
% % 
% % % bias error
% % subplot(3,1,1)
% % hold on
% % plot(dot_diameters_image, err_U_bias, 'r*-');
% % plot(dot_diameters_image, err_V_bias, 'b*-');
% % xlabel('dot diameter in image (pixels)');
% % ylabel('error (pixels)');
% % title('bias error');
% % legend('U', 'V');
% % grid on
% % 
% % % random error
% % subplot(3,1,2)
% % hold on
% % plot(dot_diameters_image, err_U_random, 'r*-');
% % plot(dot_diameters_image, err_V_random, 'b*-');
% % xlabel('dot diameter in image (pixels)');
% % ylabel('error (pixels)');
% % title('random error');
% % legend('U', 'V');
% % grid on
% % 
% % % rms error
% % subplot(3,1,3)
% % hold on
% % plot(dot_diameters_image, err_U_rms, 'r*-');
% % plot(dot_diameters_image, err_V_rms, 'b*-');
% % xlabel('dot diameter in image (pixels)');
% % ylabel('error (pixels)');
% % title('rms error');
% % legend('U', 'V');
% % grid on
% % 
% % % save results to file
% % savefig(gcf, [figure_save_filepath 'errors-displacement-image-diameter.fig']);
% % print(gcf, [figure_save_filepath 'errors-displacement-image-diameter.eps'], '-depsc');
% % print(gcf, [figure_save_filepath 'errors-displacement-image-diameter.png'], '-dpng');

%% save worskpace to file

% this is the name of the current script
script_name_full = mfilename('fullpath');
[pathstr, script_name, ext] = fileparts(script_name_full);
save([workspace_filepath script_name '.mat']);

