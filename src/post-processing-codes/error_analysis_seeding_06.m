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

% this is the delta_x value that represents the displacement
delta_x_array = 0.1:0.2:1.0;
delta_y_array = delta_x_array;

% this is the number of delta_x cases
num_cases_delta_x = length(delta_x_array);

% these are teh seeding densities considered (particles / 32x32 pix window)
seeding_density_array = 4:4:20;

% this is the number of seeding density cases
num_cases_seeding_density = length(seeding_density_array);

% this is the base name of the files that contain the results to be analyzed
results_basename =  'BOS_pass1_';

% this is the number of .mat results files corresponding to a seeding
num_files_per_seeding = 5;

% this ensures all the figures are docked and displayed from a single
% window
set(0,'DefaultFigureWindowStyle','docked')

%% read and write paths

% this is the folder where the figures will be saved
figure_save_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/plots/16x16/';

% this is the folder where the workspace variables will be saved
workspace_save_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/results/16x16/';

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
% U_ref_array = delta_x_array/5.0 * 2.0;
% % U_ref = 5;
% V_ref_array = zeros(size(delta_x_array));

% this is threshold for checking whether a vector is valid or not (pix.)
valid_vector_detection_threshold = 0.2;

% this is the magnification
M = 0.0867; %123500./700000;

% this is the pixel pitch (mm)
pixel_pitch = 17e-3;

% this is the aperture f number
f_number = 16;

%% do the following steps for all cases
num_cases_total = num_cases_delta_x * num_cases_seeding_density;
fprintf('number of cases: %d\n', num_cases_total);

% array containing reference displacements along x and y
U_ref = zeros(num_cases_delta_x, num_cases_seeding_density);
V_ref = zeros(num_cases_delta_x, num_cases_seeding_density);

% this is the array containing the bias error
err_U_bias = zeros(num_cases_delta_x, num_cases_seeding_density);
err_V_bias = zeros(num_cases_delta_x, num_cases_seeding_density);

% this is the array containing the random error
err_U_random = zeros(num_cases_delta_x, num_cases_seeding_density);
err_V_random = zeros(num_cases_delta_x, num_cases_seeding_density);

% this is the array containing the total error
err_U_total = zeros(num_cases_delta_x, num_cases_seeding_density);
err_V_total = zeros(num_cases_delta_x, num_cases_seeding_density);

% this is the number of valid vectors
num_valid_vectors = zeros(num_cases_delta_x, num_cases_seeding_density);

% this is the total number of vectors
num_total_vectors = zeros(num_cases_delta_x, num_cases_seeding_density);

% this is the valid vector detection probability
valid_vector_detection_probability = zeros(num_cases_delta_x, num_cases_seeding_density);

all_errors = cell([num_cases_delta_x, num_cases_seeding_density]);

sample_struct = struct('num_total_vectors', [], 'num_valid_vectors', [], 'ref_disp_x', [], 'error', [], 'bias_error', [], 'random_error', [], 'total_error', []);

for delta_x_index = 1:num_cases_delta_x
    for seeding_density_index = 1:num_cases_seeding_density
        
        delta_x = delta_x_array(delta_x_index);
        delta_y = delta_y_array(delta_x_index);

        seeding_density = seeding_density_array(seeding_density_index);
    
        all_errors{delta_x_index, seeding_density_index} = sample_struct;
        
        all_errors{delta_x_index, seeding_density_index}.delta_x = delta_x;
        all_errors{delta_x_index, seeding_density_index}.seeding_density = seeding_density;
        
        % display progress to user
        fprintf('delta_x: %0.2f, seeding_density: %d\n', delta_x, seeding_density);

        %% set read  and write paths

        % this is the path to the folder where the results are located
        vectors_filepath = ['/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/seeding/delta_x=' num2str(delta_x, '%0.2f') ...
                            '_delta_y=' num2str(delta_y, '%0.2f') '/' num2str(seeding_density) '/processing-16x16/results/vectors/'];


        % location where the light ray positions are stored
        light_ray_positions_filepath = ['/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/seeding/delta_x=' num2str(delta_x, '%0.2f') ...
                            '_delta_y=' num2str(delta_y, '%0.2f') '/' num2str(seeding_density) '/1/light-ray-positions/'];

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
        
        % TEMP FIX
        del_y_pixels = -del_y_pixels;
        
        % copy the values to the variables already in the code
        U_ref(delta_x_index, seeding_density_index) = del_x_pixels;
        V_ref(delta_x_index, seeding_density_index) = del_y_pixels;

        % close files
        fclose(fid_1);
        fclose(fid_2);
        
        all_errors{delta_x_index, seeding_density_index}.U_ref = del_x_pixels;
        all_errors{delta_x_index, seeding_density_index}.V_ref = del_y_pixels;
        
        
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

        % combine all vectors into a single array
        U_temp = [];
        V_temp = [];
        
        % this loops through all the files in the foder
        for file_index = 1:num_files
            % load data
            data = load([vectors_filepath filenames(file_index).name]);

            % account for sign convention used in the code
            data.U(:) = -data.U(:);

            % if U is 3D (if additional peaks are saved), then consider
            % only the final results for error analysis
            if size(size(data.U),2) == 3
                U_temp_2 = data.U(:,:,1);
                U_temp = [U_temp U_temp_2(:)'];
                V_temp_2 = data.V(:,:,1);
                V_temp = [V_temp V_temp_2(:)'];
            else
                U_temp = [U_temp data.U(:)'];
                V_temp = [V_temp data.V(:)'];
            end
            
        end        
        
        % remove invalid vectors (that are greater or less than ref
        U_temp_valid = U_temp(abs(U_temp - U_ref(delta_x_index, seeding_density_index)) < valid_vector_detection_threshold);
        V_temp_valid = V_temp(abs(V_temp - V_ref(delta_x_index, seeding_density_index)) < valid_vector_detection_threshold);

        fprintf('num_valid_vectors: %d, num_total_vectors: %d\n', length(U_temp_valid), length(U_temp));

        % increment total number of valid vectors for this case
        num_valid_vectors(delta_x_index, seeding_density_index) = length(U_temp_valid);
        num_total_vectors(delta_x_index, seeding_density_index) = length(U_temp);
        
        all_errors{delta_x_index, seeding_density_index}.num_valid_vectors = num_valid_vectors(delta_x_index, seeding_density_index);

        all_errors{delta_x_index, seeding_density_index}.errors_U = U_temp_valid - U_ref(delta_x_index, seeding_density_index);
        all_errors{delta_x_index, seeding_density_index}.errors_V = V_temp_valid - V_ref(delta_x_index, seeding_density_index);

        all_errors{delta_x_index, seeding_density_index}.U = U_temp_valid;
        all_errors{delta_x_index, seeding_density_index}.V = V_temp_valid;

        [err_U_bias(delta_x_index, seeding_density_index), err_U_random(delta_x_index, seeding_density_index), ...
            err_U_total(delta_x_index, seeding_density_index)] = compute_errors(U_temp_valid(:)', U_ref(delta_x_index, seeding_density_index));
        [err_V_bias(delta_x_index, seeding_density_index), err_V_random(delta_x_index, seeding_density_index), ...
            err_V_total(delta_x_index, seeding_density_index)] = compute_errors(V_temp_valid(:)', V_ref(delta_x_index, seeding_density_index));

        all_errors{delta_x_index, seeding_density_index}.bias_error.U = err_U_bias(delta_x_index, seeding_density_index);
        all_errors{delta_x_index, seeding_density_index}.bias_error.V = err_V_bias(delta_x_index, seeding_density_index);
        all_errors{delta_x_index, seeding_density_index}.random_error.U = err_U_random(delta_x_index, seeding_density_index);
        all_errors{delta_x_index, seeding_density_index}.random_error.V = err_V_random(delta_x_index, seeding_density_index);
        all_errors{delta_x_index, seeding_density_index}.total_error.U = err_U_total(delta_x_index, seeding_density_index);
        all_errors{delta_x_index, seeding_density_index}.total_error.V = err_V_total(delta_x_index, seeding_density_index);
        
    
        % calculate valid vector detection probability
        valid_vector_detection_probability(delta_x_index, seeding_density_index) = num_valid_vectors(delta_x_index, seeding_density_index)/ ...
                                                                                    num_total_vectors(delta_x_index, seeding_density_index);
    end  
    
end

%% plot results

% select color map for all the plots
cmap = parula;


%%%% plot absolute bias error %%%%
% create one figure to have scaled ellipse axes (to 1)
figure_ctr = figure_ctr+1;
fig1 = figure(figure_ctr);
hold on

% create another figure to have unscaled ellipse axes (raw values)
figure_ctr = figure_ctr+1;
fig2 = figure(figure_ctr);
hold on

% calculate error magnitude
absolute_bias_error = sqrt(err_U_bias.^2 + err_V_bias.^2);
% calculate ratio of errors
ellipse_axes_ratio = err_V_bias./err_U_bias;

% calculate min and max error
min_error_color = min(absolute_bias_error(:));
max_error_color = max(absolute_bias_error(:));

% set min and max for colorbar
figure(figure_ctr-1)
caxis([min_error_color, max_error_color])

figure(figure_ctr)
caxis([min_error_color, max_error_color])

% iteratore over all the cases and plot the errors as filled ellipses
for delta_x_index = 1:num_cases_delta_x
    for seeding_density_index = 1:num_cases_seeding_density
        
        % get current value of density gradient
        delta_x = 10 * delta_x_array(delta_x_index);
        % get current value of seeding density
        seeding_density = seeding_density_array(seeding_density_index);
        
        % get current error
        current_error = absolute_bias_error(delta_x_index, seeding_density_index);
        % calculate the position in the color map for the current error
        color_id = ceil((current_error-min_error_color)/(max_error_color - min_error_color) * size(cmap,1));

        % account for values at the edge
        if color_id == 0
            color_id = 1;
        elseif isnan(color_id)
            continue
        elseif color_id > size(cmap,1)
            color_id = size(cmap,1);
        end
        
        % assign the color for the current error
        color_spec = cmap(color_id,:);
                        
        figure(figure_ctr-1)
        scaling_factor = 0.6;
        % calculate the major and minor axes of the ellipse
        rx = scaling_factor * 1/sqrt(1+ellipse_axes_ratio(delta_x_index, seeding_density_index)^2);
        ry = scaling_factor * 1/sqrt(1+1/ellipse_axes_ratio(delta_x_index, seeding_density_index)^2);
        
        % draw a filled ellipse at the color value
        create_filled_ellipse(seeding_density, delta_x, rx, ry,color_spec)

        figure(figure_ctr)
        scaling_factor = 2;
     
        rx = scaling_factor * err_U_random(delta_x_index, seeding_density_index);
        ry = scaling_factor * err_V_random(delta_x_index, seeding_density_index);

        % draw a filled ellipse at the color value
        create_filled_ellipse(seeding_density, delta_x, rx, ry,color_spec)
        
    end
end

figure(figure_ctr)
colorbar
xlabel('Dots (32x32 pix.)');
ylabel('Reference Displacement (pix.)');
title('Bias Error')
set(gca, 'fontsize', 14)
axis tight
axis equal

figure(figure_ctr-1)
colorbar
xlabel('Dots (32x32 pix.)');
ylabel('Reference Displacement (pix.)');
title('Bias Error')
set(gca, 'fontsize', 14)
axis tight
axis equal
% save figure to file
save_figure_to_file(fig1, figure_save_filepath, 'bias-errors-seeding-scaled');
save_figure_to_file(fig2, figure_save_filepath, 'bias-errors-seeding');

%%%% plot absolute random error error %%%%
% create one figure to have scaled ellipse axes (to 1)
figure_ctr = figure_ctr+1;
fig1 = figure(figure_ctr);
hold on

% create another figure to have unscaled ellipse axes (raw values)
figure_ctr = figure_ctr+1;
fig2 = figure(figure_ctr);
hold on

% calculate error magnitude
absolute_random_error = sqrt(err_U_random.^2 + err_V_random.^2);
% calculate ratio of errors
ellipse_axes_ratio = err_V_random./err_U_random;

% calculate min and max error
min_error_color = min(absolute_random_error(:));
max_error_color = max(absolute_random_error(:));

% set min and max for colorbar
figure(figure_ctr-1)
caxis([min_error_color, max_error_color])

figure(figure_ctr)
caxis([min_error_color, max_error_color])

% iteratore over all the cases and plot the errors as filled ellipses
for delta_x_index = 1:num_cases_delta_x
    for seeding_density_index = 1:num_cases_seeding_density
        
        % get current value of density gradient
        delta_x = 10 * delta_x_array(delta_x_index);
        % get current value of seeding density
        seeding_density = seeding_density_array(seeding_density_index);
        
        % get current error
        current_error = absolute_random_error(delta_x_index, seeding_density_index);
        % calculate the position in the color map for the current error
        color_id = ceil((current_error-min_error_color)/(max_error_color - min_error_color) * size(cmap,1));

        % account for values at the edge
        if color_id == 0
            color_id = 1;
        elseif isnan(color_id)
            continue
        elseif color_id > size(cmap,1)
            color_id = size(cmap,1);
        end
        
        % assign the color for the current error
        color_spec = cmap(color_id,:);
                        
        figure(figure_ctr-1)
        scaling_factor = 0.2;
        % calculate the major and minor axes of the ellipse
        rx = scaling_factor * 1/sqrt(1+ellipse_axes_ratio(delta_x_index, seeding_density_index)^2);
        ry = scaling_factor * 1/sqrt(1+1/ellipse_axes_ratio(delta_x_index, seeding_density_index)^2);
        
        % draw a filled ellipse at the color value
        create_filled_ellipse(seeding_density, delta_x, rx, ry,color_spec)

        figure(figure_ctr)
        scaling_factor = 1;
     
        rx = scaling_factor * err_U_random(delta_x_index, seeding_density_index);
        ry = scaling_factor * err_V_random(delta_x_index, seeding_density_index);

        % draw a filled ellipse at the color value
        create_filled_ellipse(seeding_density, delta_x, rx, ry,color_spec)
        
    end
end

figure(figure_ctr)
colorbar
xlabel('Dots (32x32 pix.)');
ylabel('Reference Displacement (pix.)');
title('Random Error')
set(gca, 'fontsize', 14)
axis tight
axis equal

figure(figure_ctr-1)
colorbar
xlabel('Dots (32x32 pix.)');
ylabel('Reference Displacement (pix.)');
title('Random Error')
set(gca, 'fontsize', 14)
axis tight
axis equal

% save figure to file
save_figure_to_file(fig1, figure_save_filepath, 'random-errors-seeding-scaled');
save_figure_to_file(fig2, figure_save_filepath, 'random-errors-seeding');

%%%% plot absolute total error %%%%

% create one figure to have scaled ellipse axes (to 1)
figure_ctr = figure_ctr+1;
fig1 = figure(figure_ctr);
hold on

% create another figure to have unscaled ellipse axes (raw values)
figure_ctr = figure_ctr+1;
fig2 = figure(figure_ctr);
hold on

% calculate error magnitude
absolute_total_error = sqrt(err_U_total.^2 + err_V_total.^2);
% calculate ratio of errors
ellipse_axes_ratio = err_V_total./err_U_total;

% calculate min and max error
min_error_color = min(absolute_total_error(:));
max_error_color = max(absolute_total_error(:)) * 0.75;

% set min and max for colorbar
figure(figure_ctr-1)
caxis([min_error_color, max_error_color])

figure(figure_ctr)
caxis([min_error_color, max_error_color])


% iteratore over all the cases and plot the errors as filled ellipses
for delta_x_index = 1:num_cases_delta_x
    for seeding_density_index = 1:num_cases_seeding_density
        
        % get current value of density gradient
        delta_x = 10 * delta_x_array(delta_x_index);
        % get current value of seeding density
        seeding_density = seeding_density_array(seeding_density_index);
        
        % get current error
        current_error = absolute_total_error(delta_x_index, seeding_density_index);
        % calculate the position in the color map for the current error
        color_id = ceil((current_error-min_error_color)/(max_error_color - min_error_color) * size(cmap,1));

        % account for values at the edge
        if color_id == 0
            color_id = 1;
        elseif isnan(color_id)
            continue
        elseif color_id > size(cmap,1)
            color_id = size(cmap,1);
        end
        
        % assign the color for the current error
        color_spec = cmap(color_id,:);

        figure(figure_ctr-1)
        scaling_factor = 0.2;
        
        % calculate the major and minor axes of the ellipse
        rx = scaling_factor * 1/sqrt(1+ellipse_axes_ratio(delta_x_index, seeding_density_index)^2);
        ry = scaling_factor * 1/sqrt(1+1/ellipse_axes_ratio(delta_x_index, seeding_density_index)^2);
        
        % draw a filled ellipse at the color value
        create_filled_ellipse(seeding_density, delta_x, rx, ry,color_spec)

        figure(figure_ctr)
        scaling_factor = 1.0;
        
        rx = scaling_factor * err_U_total(delta_x_index, seeding_density_index);
        ry = scaling_factor * err_V_total(delta_x_index, seeding_density_index);

        % draw a filled ellipse at the color value
        create_filled_ellipse(seeding_density, delta_x, rx, ry,color_spec);
        
    end
end

figure(figure_ctr)
colorbar
xlabel('Dots (32x32 pix.)');
ylabel('Reference Displacement (pix.)');
title('Total Error')
set(gca, 'fontsize', 14)
axis tight
axis equal

figure(figure_ctr-1)
colorbar
xlabel('Dots (32x32 pix.)');
ylabel('Reference Displacement (pix.)');
title('Total Error')
set(gca, 'fontsize', 14)
axis tight
axis equal

% save figure to file
save_figure_to_file(fig1, figure_save_filepath, 'total-errors-seeding-scaled');
save_figure_to_file(fig2, figure_save_filepath, 'total-errors-seeding');


%%%%%%%%%%%%%%%% PLOT NUMBER OF VALID VECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y] = meshgrid(seeding_density_array, U_ref(:,end));

% number of valid vectors
figure_ctr = figure_ctr+1;
figure(figure_ctr);
surf(X, Y, num_valid_vectors, 'linestyle', 'none');
colorbar
view(2)
xlabel('Dots (32x32 pix.)');
ylabel('Reference Displacement (pix.)');
title('Number of valid vectors')
set(gca, 'fontsize', 14)

% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'num-valid-vector-detection-probability-U-contour');

% valid vector detection probability
figure_ctr = figure_ctr+1;
figure(figure_ctr);
surf(X, Y, valid_vector_detection_probability, 'linestyle', 'none');
colorbar
view(2)
xlabel('Dots (32x32 pix.)');
ylabel('Reference Displacement (pix.)');
title('Valid vector detection probability')
set(gca, 'fontsize', 14)

% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'valid-vector-detection-probability-U-contour');

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
for i = 2:2:num_cases_delta_x
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

% % NIFIFO for different seedings
figure_ctr = figure_ctr+1;
figure(figure_ctr);
hold on
colors = ['r', 'g', 'b', 'k', 'y', 'm', 'c'];

legend_string_seeding = cell([5,1]);
for i = 1:2:length(seeding_density_array)
    i
    plot(NIFIFO(:,i), valid_vector_detection_probability(:,i), [colors((i+1)/2) '*'])
    legend_string_seeding{(i+1)/2} = ['seeding = ' num2str(seeding_density_array(i)) ' dots/32x32 pix.'];
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

% effect of seeding
delta_x_index = 5;
figure_ctr = figure_ctr+1;
figure(figure_ctr);
subplot(2,1,1)
hold on
for i = 1:2:length(seeding_density_array)
[val,edge] = histcounts(all_errors{delta_x_index,i}.errors_U, 'numbins', num_bins, 'normalization', 'pdf');
center = 0.5 * (edge(1:end-1) + edge(2:end));
area(center, val, 'facecolor', colors((i+1)/2), 'facealpha', 0.25, 'linewidth', 2.0, 'edgecolor', colors((i+1)/2))
end
legend(legend_string_seeding)
grid on
xlabel('error (pix.)')
title({'Effect of seeding, PDF', ['\Delta x_{ref}=' num2str(U_ref(delta_x_index,5), '%.2f') ' pix.']})

subplot(2,1,2)
hold on
for i = 1:2:length(seeding_density_array)
[val,edge] = histcounts(all_errors{delta_x_index,i}.errors_U, 'numbins', num_bins, 'normalization', 'cdf');
center = 0.5 * (edge(1:end-1) + edge(2:end));
area(center, val, 'facecolor', colors((i+1)/2), 'facealpha', 0.25, 'linewidth', 2.0, 'edgecolor', colors((i+1)/2))
end
legend(legend_string_seeding)
grid on
xlabel('error (pix.)')
title({'Effect of seeding, CDF', ['\Delta x_{ref}=' num2str(U_ref(delta_x_index,5), '%.2f') ' pix.']})

% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'error-histogram-U-seeding');

% effect of displacement
seeding_density_index = 8;
figure_ctr = figure_ctr+1;
figure(figure_ctr);
subplot(2,1,1)
hold on
for i = 2:2:num_cases_delta_x
[val,edge] = histcounts(all_errors{i, seeding_density_index}.errors_U, 'numbins', num_bins, 'normalization', 'pdf');
center = 0.5 * (edge(1:end-1) + edge(2:end));
area(center, val, 'facecolor', colors(i/2), 'facealpha', 0.25, 'linewidth', 2.0, 'edgecolor', colors(i/2))
end
legend(legend_string_displacement)
grid on
xlabel('error (pix.)')
title({'Effect of displacement, PDF', [num2str(seeding_density_array(seeding_density_index)) ' dots/32x32 pix']})

subplot(2,1,2)
hold on
for i = 2:2:num_cases_delta_x
[val,edge] = histcounts(all_errors{i, seeding_density_index}.errors_U, 'numbins', num_bins, 'normalization', 'cdf');
center = 0.5 * (edge(1:end-1) + edge(2:end));
area(center, val, 'facecolor', colors(i/2), 'facealpha', 0.25, 'linewidth', 2.0, 'edgecolor', colors(i/2))
end
legend(legend_string_displacement)
grid on
xlabel('error (pix.)')
title({'Effect of displacement, CDF', [num2str(seeding_density_array(seeding_density_index)) ' dots/32x32 pix']})

% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'error-histogram-U-displacement');

%%% histogram of subpixel displacement %%%
delta_x_index = 5;
figure_ctr = figure_ctr+1;
figure(figure_ctr);
subplot(2,1,1)
hold on
for i = 1:2:length(seeding_density_array)
    [val,edge] = histcounts(all_errors{delta_x_index,i}.U - floor(all_errors{delta_x_index,i}.U), 'numbins', num_bins, 'normalization', 'pdf');
    center = 0.5 * (edge(1:end-1) + edge(2:end));
    area(center, val, 'facecolor', colors((i+1)/2), 'facealpha', 0.25, 'linewidth', 2.0, 'edgecolor', colors((i+1)/2))
end
legend(legend_string_seeding, 'location', 'northwest')
grid on
xlabel('displacement (pix.)')
title({'PDF of sub pixel displacement', ['\Delta x_{ref}=' num2str(U_ref(delta_x_index,5), '%.2f') ' pix.']})
xlim([0 1])

subplot(2,1,2)
hold on
for i = 1:2:length(seeding_density_array)
    [val,edge] = histcounts(all_errors{delta_x_index,i}.U - floor(all_errors{delta_x_index,i}.U), 'numbins', num_bins, 'normalization', 'cdf');
    center = 0.5 * (edge(1:end-1) + edge(2:end));
    area(center, val, 'facecolor', colors((i+1)/2), 'facealpha', 0.25, 'linewidth', 2.0, 'edgecolor', colors((i+1)/2))
end
legend(legend_string_seeding, 'location', 'northwest')
grid on
xlabel('displacement (pix.)')
title({'CDF of sub pixel displacement', ['\Delta y_{ref}=' num2str(V_ref(delta_x_index,5), '%.2f') ' pix.']})
xlim([0 1])

% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'histogram-subpixel-displacement');


%%% plot random error vs displacement %%%

% for each seeding density plot random errors
figure_ctr = figure_ctr+1;
figure(figure_ctr);
hold on
colors = ['r', 'g', 'b', 'k', 'y'];

legend_string_seeding = cell([4,1]);

index_range = 1:2:length(seeding_density_array);

lower_envelope = min(absolute_random_error(:,index_range), [], 2)';
upper_envelope = max(absolute_random_error(:,index_range), [], 2)';

mean_line = nanmean(absolute_random_error(:,index_range), 2)';
std_line = nanstd(absolute_random_error(:,index_range), [], 2)';
% 
% lower_envelope = mean_line - std_line;
% upper_envelope = mean_line + std_line;


X = [U_ref(:,3)', fliplr(U_ref(:,3)')];
Y = [lower_envelope, fliplr(upper_envelope)];

for i = 1:2:length(seeding_density_array)
    i
    plot(U_ref(:,i), absolute_random_error(:,i), '*')
%     plot(U_ref(:,i), err_U_random(:,i), '*')
end
plot(U_ref(:,3), mean_line, 'g', 'linewidth', 2.0)
grid on
% fill(X,Y,'b', 'facealpha', 0.2, 'edgecolor', 'b', 'linewidth', 2.0)
% legend(legend_string_seeding, 'location', 'southeast')
xlabel('Reference Displacement (pix.)')
ylabel('Random Error (pix.)')
title('Effect of seeding')
% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'random-error-displacement');

%%%%%%%% plot pdf of displacements to identify peak locking %%%%%%
all_vectors.U = [];
all_vectors.V = [];

for delta_x_index = 1:num_cases_delta_x
    for seeding_density_index = 1:num_cases_seeding_density
        
        delta_x = delta_x_array(delta_x_index);
        delta_y = delta_y_array(delta_x_index);

        seeding_density = seeding_density_array(seeding_density_index);
        
        all_vectors.U = [all_vectors.U all_errors{delta_x_index, seeding_density_index}.U];
        all_vectors.V = [all_vectors.V all_errors{delta_x_index, seeding_density_index}.V];                        
        
    end
end

num_bins = 100;

figure
[val,edge] = histcounts(all_vectors.U, 'numbins', num_bins, 'normalization', 'pdf');
center = 0.5 * (edge(1:end-1) + edge(2:end));
area(center, val, 'facecolor', 'r', 'facealpha', 0.25, 'linewidth', 2.0, 'edgecolor', 'r')
xlabel('displacement (pixels)')
title('pdf of displacements')
save_figure_to_file(gcf, figure_save_filepath, 'pdf_displacements')


%%% line plot of random error vs displacement to see effect of seeding %%%%
figure
hold on
legend_string = cell([length(seeding_density_array), 1]);
ctr = 0;
for seeding_density_index = 1:2:length(seeding_density_array)
    ctr = ctr + 1;
    plot(U_ref(:,seeding_density_index), abs(err_U_random(:,seeding_density_index)), '-*', 'linewidth', 2.0)
%     plot(delta_x_array, abs(err_V_bias(:,seeding_density_index)), '-o')
    legend_string{ctr,1} = ['seeding=' num2str(seeding_density_array(seeding_density_index)) 'dots/32x32 pix.'];
end
grid on
legend(legend_string{1:ctr})
xlabel('displacement (pix.)')
ylabel('error (pix.)')
title('random error')
set(gca, 'fontsize', 14)
xlim([0.0 1.5])
set(gca, 'xtick', 0:0.1:1.5)
set(gca, 'ytick', 0.06:0.01:0.09)
save_figure_to_file(gcf, figure_save_filepath, 'random-error-U-line')

%% save workspace to file
% this is the name of the current script
script_name_full = mfilename('fullpath');
[pathstr, script_name, ext] = fileparts(script_name_full);
save([workspace_save_filepath script_name '.mat']);

