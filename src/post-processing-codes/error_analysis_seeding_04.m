% This program reads in the displacements estimated by cross-correlating
% the bos images, and compares it to the theoretical displacement. It then
% computes the mean and total errors

%% clear workspace and close all windows

clear 
close all
clc

figure_ctr = 0;

starting_index = 1;
%% overall settings

% this is the grad_x value that represents the displacement
grad_x_array = 0.5:0.5:5;

% this is the number of grad_x cases
num_cases_grad_x = length(grad_x_array);

% these are teh seeding densities considered (particles / 32x32 pix window)
seeding_density_array = 1:20;

% this is the number of seeding density cases
num_cases_seeding_density = length(seeding_density_array);

% this is the base name of the files that contain the results to be analyzed
results_basename =  'bos_pass1_';

% this is the number of .mat results files corresponding to a seeding
num_files_per_seeding = 5;

% this ensures all the figures are docked and displayed from a single
% window
set(0,'DefaultFigureWindowStyle','docked')

%% read and write paths
% this is the folder where the figures will be saved
figure_save_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/plots/';

% this is the folder where the workspace variables will be saved
workspace_save_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/results/';

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
% U_ref_array = grad_x_array/5.0 * 2.0;
% % U_ref = 5;
% V_ref_array = zeros(size(grad_x_array));

% this is threshold for checking whether a vector is valid or not (pix.)
valid_vector_detection_threshold = 0.2;

% this is the magnification
M = 0.0867; %123500./700000;

% this is the pixel pitch (mm)
pixel_pitch = 17e-3;

% this is the aperture f number
f_number = 16;

%% do the following steps for all cases
num_cases_total = num_cases_grad_x * num_cases_seeding_density;
fprintf('number of cases: %d\n', num_cases_total);

% array containing reference displacements along x and y
U_ref = zeros(num_cases_grad_x, num_cases_seeding_density);
V_ref = zeros(num_cases_grad_x, num_cases_seeding_density);

% this is the array containing the bias error
err_U_bias = zeros(num_cases_grad_x, num_cases_seeding_density);
err_V_bias = zeros(num_cases_grad_x, num_cases_seeding_density);

% this is the array containing the random error
err_U_random = zeros(num_cases_grad_x, num_cases_seeding_density);
err_V_random = zeros(num_cases_grad_x, num_cases_seeding_density);

% this is the array containing the total error
err_U_total = zeros(num_cases_grad_x, num_cases_seeding_density);
err_V_total = zeros(num_cases_grad_x, num_cases_seeding_density);

% this is the number of valid vectors
num_valid_vectors = zeros(num_cases_grad_x, num_cases_seeding_density);

% this is the total number of vectors
num_total_vectors = zeros(num_cases_grad_x, num_cases_seeding_density);

% this is the valid vector detection probability
valid_vector_detection_probability = zeros(num_cases_grad_x, num_cases_seeding_density);

all_errors = cell([num_cases_grad_x, num_cases_seeding_density]);

sample_struct = struct('num_total_vectors', [], 'num_valid_vectors', [], 'ref_disp_x', [], 'error', [], 'bias_error', [], 'random_error', [], 'total_error', []);

for grad_x_index = 1:num_cases_grad_x
    for seeding_density_index = 1:num_cases_seeding_density
        
        grad_x = grad_x_array(grad_x_index);
        seeding_density = seeding_density_array(seeding_density_index);
    
        all_errors{grad_x_index, seeding_density_index} = sample_struct;
        
        all_errors{grad_x_index, seeding_density_index}.grad_x = grad_x;
        all_errors{grad_x_index, seeding_density_index}.seeding_density = seeding_density;
        
        % display progress to user
        fprintf('grad_x: %0.2f, seeding_density: %d\n', grad_x, seeding_density);

        %% set read  and write paths

        % this is the path to the folder where the results are located
        vectors_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/grad_x=' num2str(grad_x, '%0.2f') ...
                            '/' num2str(seeding_density) '/processing/results/vectors/'];


        % location where the light ray positions are stored
        light_ray_positions_filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/grad_x=' num2str(grad_x, '%0.2f') ...
                            '/' num2str(seeding_density) '/1/light-ray-positions/'];

        %% calculate reference displacements from light ray positions
        
        % open file containing positions for case w/o density gradients
        fid_1 = fopen([light_ray_positions_filepath 'im1/pos_0000.bin']);
        % open file containing positions for case with density gradients
        fid_2 = fopen([light_ray_positions_filepath 'im2/pos_0000.bin']);
        
        % load position data
        pos_1 = fread(fid_1, 'single');
        pos_2 = fread(fid_2, 'single');
        
        % extract x coordinates (microns)
        pos_1_x = pos_1(1:2:end);
        pos_2_x = pos_2(1:2:end);
        
        % extract y coordinates (microns)
        pos_1_y = pos_1(2:2:end);
        pos_2_y = pos_2(2:2:end);
        
        % find average ray deflection along x and y (microns)
        del_x = nanmean(pos_2_x - pos_1_x);
        del_y = nanmean(pos_2_y - pos_1_y);
        
        % convert to pixel units
        del_x_pixels = del_x*1e-3/pixel_pitch;
        del_y_pixels = del_y*1e-3/pixel_pitch;
        
        % multipy del_x with -1 to account for co-ordinate system in ray
        % tracing code
        del_x_pixels = -del_x_pixels;
        
        % copy the values to the variables already in the code
        U_ref(grad_x_index, seeding_density_index) = del_x_pixels;
        V_ref(grad_x_index, seeding_density_index) = del_y_pixels;

        % close files
        fclose(fid_1);
        fclose(fid_2);
        
        all_errors{grad_x_index, seeding_density_index}.U_ref = del_x_pixels;
        all_errors{grad_x_index, seeding_density_index}.V_ref = del_y_pixels;
        
        
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
        err_U_total_single = zeros(1,num_files);
        err_V_total_single = zeros(1,num_files);

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
            U_temp_valid = U_temp(abs(U_temp - U_ref(grad_x_index, seeding_density_index)) < valid_vector_detection_threshold);
            V_temp_valid = V_temp(abs(V_temp - V_ref(grad_x_index, seeding_density_index)) < valid_vector_detection_threshold);

            fprintf('num_valid_vectors: %d, num_total_vectors: %d\n', length(U_temp_valid), length(U_temp));
            
            % increment total number of valid vectors for this case
            num_valid_vectors(grad_x_index, seeding_density_index) = num_valid_vectors(grad_x_index, seeding_density_index) ...
                                                                        + length(U_temp_valid);
            num_total_vectors(grad_x_index, seeding_density_index) = num_total_vectors(grad_x_index, seeding_density_index) ...
                                                                        + length(U_temp);                                                                    
            all_errors{grad_x_index, seeding_density_index}.num_valid_vectors = num_valid_vectors(grad_x_index, seeding_density_index);
            
            all_errors{grad_x_index, seeding_density_index}.errors_U = U_temp_valid - U_ref(grad_x_index, seeding_density_index);
            all_errors{grad_x_index, seeding_density_index}.errors_V = V_temp_valid - V_ref(grad_x_index, seeding_density_index);
            
            % calculate errors
            [err_U_bias_single(file_index), err_U_random_single(file_index), ...
                err_U_total_single(file_index)] = compute_errors(U_temp_valid(:)', U_ref(grad_x_index, seeding_density_index));
            [err_V_bias_single(file_index), err_V_random_single(file_index), ...
                err_V_total_single(file_index)] = compute_errors(V_temp_valid(:)', V_ref(grad_x_index, seeding_density_index));
            
%             % count the number of valid vectors
%             [num_valid_temp, num_total_temp] = count_valid_vectors(data);
            
        end
        
%         check_convergence
       

        err_U_bias(grad_x_index, seeding_density_index) = mean(err_U_bias_single);
        err_V_bias(grad_x_index, seeding_density_index) = mean(err_V_bias_single);
        err_U_random(grad_x_index, seeding_density_index) = mean(err_U_random_single);
        err_V_random(grad_x_index, seeding_density_index) = mean(err_V_random_single);
        err_U_total(grad_x_index, seeding_density_index) = mean(err_U_total_single);
        err_V_total(grad_x_index, seeding_density_index) = mean(err_V_total_single);
        
%         all_errors{grad_x_index, seeding_density_index}.bias_error.U = err_U_bias(grad_x_index, seeding_density_index);
%         all_errors{grad_x_index, seeding_density_index}.bias_error.V = err_V_bias(grad_x_index, seeding_density_index);
%         all_errors{grad_x_index, seeding_density_index}.random_error.U = err_U_random(grad_x_index, seeding_density_index);
%         all_errors{grad_x_index, seeding_density_index}.random_error.V = err_V_random(grad_x_index, seeding_density_index);
%         all_errors{grad_x_index, seeding_density_index}.total_error.U = err_U_total(grad_x_index, seeding_density_index);
%         all_errors{grad_x_index, seeding_density_index}.total_error.V = err_V_total(grad_x_index, seeding_density_index);
        
    
        % calculate valid vector detection probability
        valid_vector_detection_probability(grad_x_index, seeding_density_index) = num_valid_vectors(grad_x_index, seeding_density_index)/ ...
                                                                                    num_total_vectors(grad_x_index, seeding_density_index);
    end  
    
end

%% plot results

% pause(1)
% 
% [X,Y] = meshgrid(seeding_density_array, U_ref(:,end));
% 
% % max bias error
% max_bias_error = max(err_U_bias(:));
% max_random_error = max(err_U_random(:));
% max_total_error = max(err_U_total(:));
% max_error = max([max_bias_error, max_random_error, max_total_error]);
% 
% 
% % bias error
% figure_ctr = figure_ctr+1;
% figure(figure_ctr);
% surf(X, Y, abs(err_U_bias), 'linestyle', 'none');
% colorbar
% % caxis([0 max_error])
% view(2)
% xlabel('Dots (32x32 pix.)');
% ylabel('Reference Displacement (pix.)');
% title('Bias Error - \Delta x')
% set(gca, 'fontsize', 14)
% 
% % save results to file
% save_figure_to_file(gcf, figure_save_filepath, 'bias-errors-seeding-U');
% 
% % random error
% figure_ctr = figure_ctr+1;
% figure(figure_ctr);
% surf(X, Y, err_U_random, 'linestyle', 'none');
% colorbar
% % caxis([0 max_error])
% view(2)
% xlabel('Dots (32x32 pix.)');
% ylabel('Rerence Displacement (pix.)');
% title('Random Error - \Delta x')
% set(gca, 'fontsize', 14)
% 
% % save results to file
% save_figure_to_file(gcf, figure_save_filepath, 'random-errors-seeding-U');
% 
% % total error
% figure_ctr = figure_ctr+1;
% figure(figure_ctr);
% surf(X, Y, err_U_total, 'linestyle', 'none');
% colorbar
% % caxis([0 max_error])
% view(2)
% xlabel('Dots (32x32 pix.)');
% ylabel('Reference Displacement (pix.)');
% title('Total Error - \Delta x')
% set(gca, 'fontsize', 14)
% 
% % save results to file
% save_figure_to_file(gcf, figure_save_filepath, 'total-errors-seeding-U');
% 
% % number of valid vectors
% figure_ctr = figure_ctr+1;
% figure(figure_ctr);
% surf(X, Y, num_valid_vectors, 'linestyle', 'none');
% colorbar
% view(2)
% xlabel('Dots (32x32 pix.)');
% ylabel('Reference Displacement (pix.)');
% title('Number of valid vectors')
% set(gca, 'fontsize', 14)
% 
% % save results to file
% save_figure_to_file(gcf, figure_save_filepath, 'num-valid-vector-detection-probability-U-contour');
% 
% % valid vector detection probability
% figure_ctr = figure_ctr+1;
% figure(figure_ctr);
% surf(X, Y, valid_vector_detection_probability, 'linestyle', 'none');
% colorbar
% view(2)
% xlabel('Dots (32x32 pix.)');
% ylabel('Reference Displacement (pix.)');
% title('Valid vector detection probability')
% set(gca, 'fontsize', 14)
% 
% % save results to file
% save_figure_to_file(gcf, figure_save_filepath, 'valid-vector-detection-probability-U-contour');

% NIFIFO
IW_size = [32, 32];
N_I = repmat(seeding_density_array, size(valid_vector_detection_probability, 1), 1);
% F_I = repmat((1 - U_ref_array/IW_size(1))', 1, size(valid_vector_detection_probability, 2));
F_I = 1 - U_ref/IW_size(1);
F_O = ones(size(N_I));
NIFIFO = N_I .* F_I .* F_O;

figure_ctr = figure_ctr+1;
figure(figure_ctr);
plot(NIFIFO(:), valid_vector_detection_probability(:) * 100, '*')
xlabel('N_IF_IF_O')
ylabel('Valid Vector Detection Probability (%)')

% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'valid-vector-detection-probability-U-line');

% NIFIFO for different displacements
figure_ctr = figure_ctr+1;
figure(figure_ctr);
hold on
colors = ['r', 'g', 'b', 'k', 'm'];

legend_string_displacement = cell([5,1]);
for i = 2:2:num_cases_grad_x
    plot(NIFIFO(i,:), valid_vector_detection_probability(i,:), [colors(i/2) '*'])
    legend_string_displacement{(i)/2} = ['\Delta x_{ref} = ' num2str(U_ref(i,5), '%.2f') ' pix.'];
end
grid on
legend(legend_string_displacement, 'location', 'southeast')
xlabel('N_IF_IF_O')
ylabel('Valid Vector Detection Probability (%)')
title('Effect of displacement')
% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'NIFIFO-U-displacement');

% NIFIFO for different seedings
figure_ctr = figure_ctr+1;
figure(figure_ctr);
hold on
colors = ['r', 'g', 'b', 'k', 'y'];

legend_string_seeding = cell([4,1]);
for i = 5:5:20
    plot(NIFIFO(:,i), valid_vector_detection_probability(:,i), [colors(i/5) '*'])
    legend_string_seeding{i/5} = ['seeding = ' num2str(seeding_density_array(i)) ' dots/32x32 pix.'];
end
grid on
legend(legend_string_seeding, 'location', 'southeast')
xlabel('N_IF_IF_O')
ylabel('Valid Vector Detection Probability (%)')
title('Effect of seeding')
% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'NIFIFO-U-seeding');

%%%%% error histogram

% set number of bins
num_bins = 20;

% figure_ctr = figure_ctr+1;
% figure(figure_ctr);
% subplot(1,2,1)
% histogram(abs(err_U_bias(:)), num_bins, 'normalization', 'pdf')
% xlabel('error (pix.)')
% title('PDF of bias error')
% subplot(1,2,2)
% histogram(abs(err_U_bias(:)), num_bins, 'normalization', 'cdf')
% xlabel('error (pix.)')
% title('CDF of bias error')
% 
% % save results to file
% save_figure_to_file(gcf, figure_save_filepath, 'error-histogram-bias');
% 
% figure_ctr = figure_ctr+1;
% figure(figure_ctr);
% subplot(1,2,1)
% histogram(err_U_random(:), num_bins, 'normalization', 'pdf')
% xlabel('error (pix.)')
% title('PDF of random error')
% subplot(1,2,2)
% histogram(err_U_random(:), num_bins, 'normalization', 'cdf')
% % H = histogram(err_U_random(:), num_bins, 'normalization', 'cdf');
% % plot(H.BinEdges(1:end-1), H.Values, '*')
% xlabel('error (pix.)')
% title('CDF of random error')
% 
% % save results to file
% save_figure_to_file(gcf, figure_save_filepath, 'error-histogram-random');
% 
% figure_ctr = figure_ctr+1;
% figure(figure_ctr);
% subplot(1,2,1)
% histogram(err_U_total(:), num_bins, 'normalization', 'pdf')
% xlabel('error (pix.)')
% title('PDF of total error')
% subplot(1,2,2)
% histogram(err_U_total(:), num_bins, 'normalization', 'cdf')
% xlabel('error (pix.)')
% title('CDF of total error')
% 
% % save results to file
% save_figure_to_file(gcf, figure_save_filepath, 'error-histogram-total');

% effect of seeding
grad_index = 5;
figure_ctr = figure_ctr+1;
figure(figure_ctr);
subplot(2,1,1)
hold on
for i = 5:5:20
[val,edge] = histcounts(all_errors{grad_x_index,i}.errors_U, 'numbins', num_bins, 'normalization', 'pdf');
center = 0.5 * (edge(1:end-1) + edge(2:end));
area(center, val, 'facecolor', colors(i/5), 'facealpha', 0.25, 'linewidth', 2.0, 'edgecolor', colors(i/5))
end
legend(legend_string_seeding)
grid on
xlabel('error (pix.)')
title({'Effect of displacement, PDF', ['\Delta x_{ref}=' num2str(U_ref(grad_index,5), '%.2f') ' pix.']})

subplot(2,1,2)
hold on
for i = 5:5:20
[val,edge] = histcounts(all_errors{grad_x_index,i}.errors_U, 'numbins', num_bins, 'normalization', 'cdf');
center = 0.5 * (edge(1:end-1) + edge(2:end));
area(center, val, 'facecolor', colors(i/5), 'facealpha', 0.25, 'linewidth', 2.0, 'edgecolor', colors(i/5))
end
legend(legend_string_seeding)
grid on
xlabel('error (pix.)')
title({'Effect of displacement, CDF', ['\Delta x_{ref}=' num2str(U_ref(grad_index,5), '%.2f') ' pix.']})

% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'error-histogram-U-seeding');

% effect of displacement
seeding_density_index = 10;
figure_ctr = figure_ctr+1;
figure(figure_ctr);
subplot(2,1,1)
hold on
for i = 2:2:num_cases_grad_x
[val,edge] = histcounts(all_errors{i, seeding_density_index}.errors_U, 'numbins', num_bins, 'normalization', 'pdf');
center = 0.5 * (edge(1:end-1) + edge(2:end));
area(center, val, 'facecolor', colors(i/2), 'facealpha', 0.25, 'linewidth', 2.0, 'edgecolor', colors(i/2))
end
legend(legend_string_displacement)
grid on
xlabel('error (pix.)')
title({'Effect of seeding, PDF', [num2str(seeding_density_array(seeding_density_index)) ' dots/32x32 pix']})

subplot(2,1,2)
hold on
for i = 2:2:num_cases_grad_x
[val,edge] = histcounts(all_errors{i, seeding_density_index}.errors_U, 'numbins', num_bins, 'normalization', 'cdf');
center = 0.5 * (edge(1:end-1) + edge(2:end));
area(center, val, 'facecolor', colors(i/2), 'facealpha', 0.25, 'linewidth', 2.0, 'edgecolor', colors(i/2))
end
legend(legend_string_displacement)
grid on
xlabel('error (pix.)')
title({'Effect of seeding, CDF', [num2str(seeding_density_array(seeding_density_index)) ' dots/32x32 pix']})

% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'error-histogram-U-displacement');


% figure_ctr = figure_ctr+1;
% figure(figure_ctr);
% % set(gcf, 'Position', [200 200 900 500])
% % set(gca, 'fontsize', 20)
% 
% % bias error
% subplot(1,3,1)
% hold on
% plot(seeding_densities, err_U_bias, 'r*-');
% plot(seeding_densities, err_V_bias, 'b*-');
% grid on
% xlim([0 max(seeding_densities)*1.2])
% xlabel('dots (32x32 pix)');
% ylabel('error (pixels)');
% legend('\Delta x', '\Delta y', 'location', 'Northwest');
% title('bias error')
% set(gca, 'fontsize', 14)
% 
% % random error
% subplot(1,3,2)
% hold on
% plot(seeding_densities, err_U_random, 'r*-');
% plot(seeding_densities, err_V_random, 'b*-');
% grid on
% xlim([0 max(seeding_densities)*1.2])
% xlabel('dots (32x32 pix)');
% ylabel('error (pixels)');
% title('random err                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               