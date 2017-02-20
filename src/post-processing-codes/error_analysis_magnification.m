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

% this is the grad_x value that represents the displacement
grad_x_array = 5;

% this is the number of grad_x cases
num_cases_grad_x = length(grad_x_array);

% these are teh focal lengths considered (mm)
focal_length_array = 100:25:300;

% this is the number of magnification cases
num_cases_focal_length = length(focal_length_array);

% this is the base name of the files that contain the results to be analyzed
results_basename =  'bos_pass1_';

% this is the number of .mat results files corresponding to a focal_length
num_files_per_focal_length = 5;

% this ensures all the figures are docked and displayed from a single
% window
set(0,'DefaultFigureWindowStyle','docked')

%% read and write paths
% this is the folder where the figures will be saved
figure_save_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/magnification/plots/';

% this is the folder where the workspace variables will be saved
workspace_save_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/magnification/results/';

% this creates the write directories if they are not already present
if ~exist(figure_save_filepath,'dir')
    mkdir(figure_save_filepath);
end

% this creates the write directories if they are not already present
if ~exist(workspace_save_filepath,'dir')
    mkdir(workspace_save_filepath);
end

%% set theoretical displacements and other experimental parameters

% these are the theoretical displacements (pix.)
% U_ref = 7.76;
%U_ref = 2.0;
U_ref_array = zeros(length(grad_x_array), length(focal_length_array(grad_x_array)));

for grad_x_index = 1:length(grad_x_array)
    for focal_length_index = 1:length(focal_length_array)
        U_ref_array(grad_x_index, focal_length_index) = 2.0 * grad_x_array(grad_x_index)/5.0 * focal_length_array(focal_length_index)/200;
    end
end

% U_ref_array = grad_x_array/5.0 * 2.0 .* focal_length_array/200;
% U_ref = 5;
% V_ref_array = zeros(size(grad_x_array));
V_ref_array = zeros(size(U_ref_array));

% this is threshold for checking whether a vector is valid or not (pix.)
valid_vector_threshold = 5.;

% this is the pixel pitch (mm)
pixel_pitch = 17e-3;

% this is the aperture f number
f_number = 16;

%% do the following steps for all cases
num_cases_total = num_cases_grad_x * num_cases_focal_length;
fprintf('number of cases: %d\n', num_cases_total);

% this is the array containing the bias error
err_U_bias = zeros(num_cases_grad_x, num_cases_focal_length);
err_V_bias = zeros(num_cases_grad_x, num_cases_focal_length);

% this is the array containing the random error
err_U_random = zeros(num_cases_grad_x, num_cases_focal_length);
err_V_random = zeros(num_cases_grad_x, num_cases_focal_length);

% this is the array containing the rms error
err_U_rms = zeros(num_cases_grad_x, num_cases_focal_length);
err_V_rms = zeros(num_cases_grad_x, num_cases_focal_length);

% this is the number of valid vectors
num_valid_vectors = zeros(num_cases_grad_x, num_cases_focal_length);

% this is the total number of vectors
num_total_vectors = zeros(num_cases_grad_x, num_cases_focal_length);

% this is the valid vector detection probability
valid_vector_detection_probability = zeros(num_cases_grad_x, num_cases_focal_length);

for grad_x_index = 1:num_cases_grad_x
    for focal_length_index = 1:num_cases_focal_length
        
        grad_x = grad_x_array(grad_x_index);
        focal_length = focal_length_array(focal_length_index);
    
        % display progress to user
        fprintf('grad_x: %0.2f, focal_length: %d\n', grad_x, focal_length);

        %% set read  and write paths

        % this is the path to the folder where the results are located
        vectors_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/magnification/grad_x=' num2str(grad_x, '%0.2f') ...
                            '/' num2str(focal_length) '/processing/results/vectors/'];


        %% load displacement data and compute errors

        % this populates all the files in the directory
        filenames = dir([vectors_filepath results_basename '*.mat']);

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

            % if U is 3D (if additional peaks are saved), then consider
            % only the final results for error analysis
            if size(size(data.U),2) == 3
                U_temp = data.U(:,:,1);
                U_temp = U_temp(:);
                V_temp = data.V(:,:,1);
                V_temp = V_temp(:);
            else
                U_temp = data.U(:);
                V_temp = data.V(:);
            end

            % remove invalid vectors (that are greater or less than ref
            U_temp_valid = U_temp(abs(U_temp - U_ref_array(grad_x_index, focal_length_index)) < valid_vector_threshold);
            V_temp_valid = V_temp(abs(V_temp - V_ref_array(grad_x_index, focal_length_index)) < valid_vector_threshold);

            fprintf('num_valid_vectors: %d, num_total_vectors: %d\n', length(U_temp_valid), length(U_temp));
            
            % increment total number of valid vectors for this case
            num_valid_vectors(grad_x_index, focal_length_index) = num_valid_vectors(grad_x_index, focal_length_index) ...
                                                                        + length(U_temp_valid);
            num_total_vectors(grad_x_index, focal_length_index) = num_total_vectors(grad_x_index, focal_length_index) ...
                                                                        + length(U_temp);                                                                    

            
            % calculate errors
            [err_U_bias_single(file_index), err_U_random_single(file_index), ...
                err_U_rms_single(file_index)] = compute_errors(U_temp_valid(:)', U_ref_array(grad_x_index, focal_length_index));
            [err_V_bias_single(file_index), err_V_random_single(file_index), ...
                err_V_rms_single(file_index)] = compute_errors(V_temp_valid(:)', V_ref_array(grad_x_index, focal_length_index));
            
%             % count the number of valid vectors
%             [num_valid_temp, num_total_temp] = count_valid_vectors(data);
            
        end
        
%         check_convergence
       

        err_U_bias(grad_x_index, focal_length_index) = mean(err_U_bias_single);
        err_V_bias(grad_x_index, focal_length_index) = mean(err_V_bias_single);
        err_U_random(grad_x_index, focal_length_index) = mean(err_U_random_single);
        err_V_random(grad_x_index, focal_length_index) = mean(err_V_random_single);
        err_U_rms(grad_x_index, focal_length_index) = mean(err_U_rms_single);
        err_V_rms(grad_x_index, focal_length_index) = mean(err_V_rms_single);
        
    
        % calculate valid vector detection probability
        valid_vector_detection_probability(grad_x_index, focal_length_index) = num_valid_vectors(grad_x_index, focal_length_index)/ ...
                                                                                    num_total_vectors(grad_x_index, focal_length_index);
    end
    
    
end

%% plot results

[X,Y] = meshgrid(focal_length_array, U_ref_array);

% max bias error
max_bias_error = max(err_U_bias(:));
max_random_error = max(err_U_random(:));
max_rms_error = max(err_U_rms(:));
max_error = max([max_bias_error, max_random_error, max_rms_error]);


% bias error
figure_ctr = figure_ctr+1;
figure(figure_ctr);
hold on
plot(focal_length_array, err_U_bias, 'r*-')
plot(focal_length_array, err_V_bias, 'b*-')
grid on
legend('\Delta x', '\Delta y')
xlabel('Focal Length (mm)')
ylabel('Error (pixels)')
title('Bias Error')
set(gca, 'fontsize', 14)

% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'bias-errors-magnification-U');

% bias error
figure_ctr = figure_ctr+1;
figure(figure_ctr);
hold on
plot(focal_length_array, err_U_random, 'r*-')
plot(focal_length_array, err_V_random, 'b*-')
grid on
legend('\Delta x', '\Delta y')
xlabel('Focal Length (mm)')
ylabel('Error (pixels)')
title('Random Error')
set(gca, 'fontsize', 14)

% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'random-errors-magnification-U');

% rms error
% bias error
figure_ctr = figure_ctr+1;
figure(figure_ctr);
hold on
plot(focal_length_array, err_U_rms, 'r*-')
plot(focal_length_array, err_V_rms, 'b*-')
grid on
legend('\Delta x', '\Delta y')
xlabel('Focal Length (mm)')
ylabel('Error (pixels)')
title('RMS Error')
set(gca, 'fontsize', 14)

% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'rms-errors-magnification-U');
    
%% save worskpace to file

% this is the name of the current script
script_name_full = mfilename('fullpath');
[pathstr, script_name, ext] = fileparts(script_name_full);
save([workspace_save_filepath script_name '.mat']);
