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

% these are the dot sizes
dotsize_array = [200, 247, 294, 341, 388, 435, 482, 529, 576, 623, 670, 717, 764, 811, 858, 905, 952, 1000];

% this is the number of seeding density cases
num_cases_dotsize = length(dotsize_array);

% this is the base name of the files that contain the results to be analyzed
results_basename =  'bos_pass1_';

% this is the number of .mat results files corresponding to a seeding
num_files_per_seeding = 5;

% this ensures all the figures are docked and displayed from a single
% window
set(0,'DefaultFigureWindowStyle','docked')

%% read and write paths
% this is the folder where the figures will be saved
figure_save_filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/results/dotsize/plots/';

% this is the folder where the workspace variables will be saved
workspace_save_filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/results/dotsize/';

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

dotsize_pixels = dotsize_array * M/(pixel_pitch*1e3)

%% do the following steps for all cases
num_cases_total = num_cases_grad_x * num_cases_dotsize;
fprintf('number of cases: %d\n', num_cases_total);

% array containing reference displacements along x and y
U_ref = zeros(num_cases_grad_x, num_cases_dotsize);
V_ref = zeros(num_cases_grad_x, num_cases_dotsize);

% this is the array containing the bias error
err_U_bias = zeros(num_cases_grad_x, num_cases_dotsize);
err_V_bias = zeros(num_cases_grad_x, num_cases_dotsize);

% this is the array containing the random error
err_U_random = zeros(num_cases_grad_x, num_cases_dotsize);
err_V_random = zeros(num_cases_grad_x, num_cases_dotsize);

% this is the array containing the total error
err_U_total = zeros(num_cases_grad_x, num_cases_dotsize);
err_V_total = zeros(num_cases_grad_x, num_cases_dotsize);

% this is the number of valid vectors
num_valid_vectors = zeros(num_cases_grad_x, num_cases_dotsize);

% this is the total number of vectors
num_total_vectors = zeros(num_cases_grad_x, num_cases_dotsize);

% this is the valid vector detection probability
valid_vector_detection_probability = zeros(num_cases_grad_x, num_cases_dotsize);

all_errors = cell([num_cases_grad_x, num_cases_dotsize]);

sample_struct = struct('num_total_vectors', [], 'num_valid_vectors', [], 'ref_disp_x', [], 'error', [], 'bias_error', [], 'random_error', [], 'total_error', []);

for grad_x_index = 1:num_cases_grad_x
    for dotsize_index = 1:num_cases_dotsize
        
        grad_x = grad_x_array(grad_x_index);
        dotsize = dotsize_array(dotsize_index);
    
        all_errors{grad_x_index, dotsize_index} = sample_struct;
        
        all_errors{grad_x_index, dotsize_index}.grad_x = grad_x;
        all_errors{grad_x_index, dotsize_index}.dotsize = dotsize;
        
        % display progress to user
        fprintf('grad_x: %0.2f, dotsize: %d\n', grad_x, dotsize);

        %% set read  and write paths

        % this is the path to the folder where the results are located
        vectors_filepath = ['/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/dotsize/grad_x=' num2str(grad_x, '%0.2f') ...
                            '/' num2str(dotsize) 'um/processing/results/vectors/'];


        % location where the light ray positions are stored
        light_ray_positions_filepath = ['/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/dotsize/grad_x=' num2str(grad_x, '%0.2f') ...
                            '/' num2str(dotsize) 'um/1/light-ray-positions/'];

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
        U_ref(grad_x_index, dotsize_index) = del_x_pixels;
        V_ref(grad_x_index, dotsize_index) = del_y_pixels;

        % close files
        fclose(fid_1);
        fclose(fid_2);
        
        all_errors{grad_x_index, dotsize_index}.U_ref = del_x_pixels;
        all_errors{grad_x_index, dotsize_index}.V_ref = del_y_pixels;
        
        
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
            U_temp_valid = U_temp(abs(U_temp - U_ref(grad_x_index, dotsize_index)) < valid_vector_detection_threshold);
            V_temp_valid = V_temp(abs(V_temp - V_ref(grad_x_index, dotsize_index)) < valid_vector_detection_threshold);

            fprintf('num_valid_vectors: %d, num_total_vectors: %d\n', length(U_temp_valid), length(U_temp));
            
            % increment total number of valid vectors for this case
            num_valid_vectors(grad_x_index, dotsize_index) = num_valid_vectors(grad_x_index, dotsize_index) ...
                                                                        + length(U_temp_valid);
            num_total_vectors(grad_x_index, dotsize_index) = num_total_vectors(grad_x_index, dotsize_index) ...
                                                                        + length(U_temp);                                                                    
            all_errors{grad_x_index, dotsize_index}.num_valid_vectors = num_valid_vectors(grad_x_index, dotsize_index);
            
            all_errors{grad_x_index, dotsize_index}.errors_U = U_temp_valid - U_ref(grad_x_index, dotsize_index);
            all_errors{grad_x_index, dotsize_index}.errors_V = V_temp_valid - V_ref(grad_x_index, dotsize_index);
            
            % calculate errors
            [err_U_bias_single(file_index), err_U_random_single(file_index), ...
                err_U_total_single(file_index)] = compute_errors(U_temp_valid(:)', U_ref(grad_x_index, dotsize_index));
            [err_V_bias_single(file_index), err_V_random_single(file_index), ...
                err_V_total_single(file_index)] = compute_errors(V_temp_valid(:)', V_ref(grad_x_index, dotsize_index));
            
%             % count the number of valid vectors
%             [num_valid_temp, num_total_temp] = count_valid_vectors(data);
            
        end
        
%         check_convergence
       

        err_U_bias(grad_x_index, dotsize_index) = mean(err_U_bias_single);
        err_V_bias(grad_x_index, dotsize_index) = mean(err_V_bias_single);
        err_U_random(grad_x_index, dotsize_index) = mean(err_U_random_single);
        err_V_random(grad_x_index, dotsize_index) = mean(err_V_random_single);
        err_U_total(grad_x_index, dotsize_index) = mean(err_U_total_single);
        err_V_total(grad_x_index, dotsize_index) = mean(err_V_total_single);
        
%         all_errors{grad_x_index, dotsize_index}.bias_error.U = err_U_bias(grad_x_index, dotsize_index);
%         all_errors{grad_x_index, dotsize_index}.bias_error.V = err_V_bias(grad_x_index, dotsize_index);
%         all_errors{grad_x_index, dotsize_index}.random_error.U = err_U_random(grad_x_index, dotsize_index);
%         all_errors{grad_x_index, dotsize_index}.random_error.V = err_V_random(grad_x_index, dotsize_index);
%         all_errors{grad_x_index, dotsize_index}.total_error.U = err_U_total(grad_x_index, dotsize_index);
%         all_errors{grad_x_index, dotsize_index}.total_error.V = err_V_total(grad_x_index, dotsize_index);
        
    
        % calculate valid vector detection probability
        valid_vector_detection_probability(grad_x_index, dotsize_index) = num_valid_vectors(grad_x_index, dotsize_index)/ ...
                                                                                    num_total_vectors(grad_x_index, dotsize_index);
    end  
    
end

%% plot results

% select color map for all the plots
cmap = parula;

%%%% plot absolute bias error %%%%
figure_ctr = figure_ctr+1;
figure(figure_ctr);
hold on

% calculate error magnitude
absolute_bias_error = sqrt(err_U_bias.^2 + err_V_bias.^2);
% calculate ratio of errors
ellipse_axes_ratio = err_V_bias./err_U_bias;

% calculate min and max error
min_error_color = min(absolute_bias_error(:));
max_error_color = max(absolute_bias_error(:)) * 0.5;

% set min and max for colorbar
caxis([min_error_color, max_error_color])

% iteratore over all the cases and plot the errors as filled ellipses
for grad_x_index = 1:num_cases_grad_x
    for dotsize_index = 1:num_cases_dotsize
        
        % get current value of density gradient
        grad_x = grad_x_array(grad_x_index);
        % get current value of seeding density
        dotsize = dotsize_pixels(dotsize_index);
        
        % get current error
        current_error = absolute_bias_error(grad_x_index, dotsize_index);
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
        
        % set an arbitrary scaling factor the size of the ellipse
        scaling_factor = 0.2;
        
        % calculate the major and minor axes of the ellipse
        rx = scaling_factor * 1/sqrt(1+ellipse_axes_ratio(grad_x_index, dotsize_index)^2);
        ry = scaling_factor * 1/sqrt(1+1/ellipse_axes_ratio(grad_x_index, dotsize_index)^2);
        
        % draw a filled ellipse at the color value
        create_filled_ellipse(dotsize, grad_x, rx, ry,color_spec)
        
    end
end

colorbar
xlabel('Dot size (pix.)');
ylabel('Reference Displacement (pix.)');
title('Bias Error')
set(gca, 'fontsize', 14)
axis tight
axis equal

% save figure to file
save_figure_to_file(gcf, figure_save_filepath, 'bias-errors');

%%%% plot absolute random error error %%%%
figure_ctr = figure_ctr+1;
figure(figure_ctr);
hold on

% calculate error magnitude
absolute_random_error = sqrt(err_U_random.^2 + err_V_random.^2);
% calculate ratio of errors
ellipse_axes_ratio = err_V_random./err_U_random;

% calculate min and max error
min_error_color = min(absolute_random_error(:));
max_error_color = max(absolute_random_error(:)) * 0.8;

% set min and max for colorbar
caxis([min_error_color, max_error_color])

% iteratore over all the cases and plot the errors as filled ellipses
for grad_x_index = 1:num_cases_grad_x
    for dotsize_index = 1:num_cases_dotsize
        
        % get current value of density gradient
        grad_x = grad_x_array(grad_x_index);
        % get current value of seeding density
        dotsize = dotsize_pixels(dotsize_index);
        
        % get current error
        current_error = absolute_random_error(grad_x_index, dotsize_index);
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
        
        % set an arbitrary scaling factor the size of the ellipse
        scaling_factor = 0.2;
        
        % calculate the major and minor axes of the ellipse
        rx = scaling_factor * 1/sqrt(1+ellipse_axes_ratio(grad_x_index, dotsize_index)^2);
        ry = scaling_factor * 1/sqrt(1+1/ellipse_axes_ratio(grad_x_index, dotsize_index)^2);
        
        % draw a filled ellipse at the color value
        create_filled_ellipse(dotsize, grad_x, rx, ry,color_spec)
        
    end
end

colorbar
xlabel('Dot size (pix.)');
ylabel('Reference Displacement (pix.)');
title('Random Error')
set(gca, 'fontsize', 14)
axis tight
axis equal

% save figure to file
save_figure_to_file(gcf, figure_save_filepath, 'random-errors');

%%%% plot absolute total error %%%%

figure_ctr = figure_ctr+1;
figure(figure_ctr);
hold on

% calculate error magnitude
absolute_total_error = sqrt(err_U_total.^2 + err_V_total.^2);
% calculate ratio of errors
ellipse_axes_ratio = err_V_total./err_U_total;

% calculate min and max error
min_error_color = min(absolute_total_error(:));
max_error_color = max(absolute_total_error(:)) * 0.9;

% set min and max for colorbar
caxis([min_error_color, max_error_color])

% iteratore over all the cases and plot the errors as filled ellipses
for grad_x_index = 1:num_cases_grad_x
    for dotsize_index = 1:num_cases_dotsize
        
        % get current value of density gradient
        grad_x = grad_x_array(grad_x_index);
        % get current value of seeding density
        dotsize = dotsize_pixels(dotsize_index);
        
        % get current error
        current_error = absolute_total_error(grad_x_index, dotsize_index);
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
        
        % set an arbitrary scaling factor the size of the ellipse
        scaling_factor = 0.2;
        
        % calculate the major and minor axes of the ellipse
        rx = scaling_factor * 1/sqrt(1+ellipse_axes_ratio(grad_x_index, dotsize_index)^2);
        ry = scaling_factor * 1/sqrt(1+1/ellipse_axes_ratio(grad_x_index, dotsize_index)^2);
        
%         fprintf('rx: %f, ry: %f\n', rx, ry);
        
        % draw a filled ellipse at the color value
        create_filled_ellipse(dotsize, grad_x, rx, ry,color_spec);
        
    end
end

colorbar
xlabel('Dot size (pix.)');
ylabel('Reference Displacement (pix.)');
title('Total Error')
set(gca, 'fontsize', 14)
axis tight
axis equal

% save figure to file
save_figure_to_file(gcf, figure_save_filepath, 'total-errors');

[X,Y] = meshgrid(dotsize_pixels, U_ref(:,end));

% number of valid vectors
figure_ctr = figure_ctr+1;
figure(figure_ctr);
surf(X, Y, num_valid_vectors, 'linestyle', 'none');
colorbar
view(2)
xlabel('Dot size (pix.)');
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
xlabel('Dot size (pix.)');
ylabel('Reference Displacement (pix.)');
title('Valid vector detection probability')
set(gca, 'fontsize', 14)

% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'valid-vector-detection-probability-U-contour');


colors = ['r', 'g', 'b', 'k', 'y'];

legend_string_displacement = cell([5,1]);
for i = 2:2:num_cases_grad_x
    legend_string_displacement{(i)/2} = ['\Delta x_{ref} = ' num2str(U_ref(i,5), '%.2f') ' pix.'];
end

legend_string_dotsize = cell([4,1]);
for i = 5:5:20
    legend_string_dotsize{i/5} = ['dotsize = ' num2str(dotsize_pixels(i)) ' pix.'];
end

%%%%% error histogram

% set number of bins
num_bins = 20;

% effect of seeding
grad_index = 5;
figure_ctr = figure_ctr+1;
figure(figure_ctr);
subplot(2,1,1)
hold on
for i = 5:5:length(dotsize_pixels)
[val,edge] = histcounts(all_errors{grad_x_index,i}.errors_U, 'numbins', num_bins, 'normalization', 'pdf');
center = 0.5 * (edge(1:end-1) + edge(2:end));
area(center, val, 'facecolor', colors(i/5), 'facealpha', 0.25, 'linewidth', 2.0, 'edgecolor', colors(i/5))
end
legend(legend_string_dotsize)
grid on
xlabel('error (pix.)')
title({'Effect of dotsize, PDF', ['\Delta x_{ref}=' num2str(U_ref(grad_index,5), '%.2f') ' pix.']})

subplot(2,1,2)
hold on
for i = 5:5:length(dotsize_pixels)
[val,edge] = histcounts(all_errors{grad_x_index,i}.errors_U, 'numbins', num_bins, 'normalization', 'cdf');
center = 0.5 * (edge(1:end-1) + edge(2:end));
area(center, val, 'facecolor', colors(i/5), 'facealpha', 0.25, 'linewidth', 2.0, 'edgecolor', colors(i/5))
end
legend(legend_string_dotsize)
grid on
xlabel('error (pix.)')
title({'Effect of dotsize, CDF', ['\Delta x_{ref}=' num2str(U_ref(grad_index,5), '%.2f') ' pix.']})

% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'error-histogram-U-seeding');

% effect of displacement
dotsize_index = 10;
figure_ctr = figure_ctr+1;
figure(figure_ctr);
subplot(2,1,1)
hold on
for i = 2:2:num_cases_grad_x
[val,edge] = histcounts(all_errors{i, dotsize_index}.errors_U, 'numbins', num_bins, 'normalization', 'pdf');
center = 0.5 * (edge(1:end-1) + edge(2:end));
area(center, val, 'facecolor', colors(i/2), 'facealpha', 0.25, 'linewidth', 2.0, 'edgecolor', colors(i/2))
end
legend(legend_string_displacement)
grid on
xlabel('error (pix.)')
title({'Effect of displacement, PDF', [num2str(dotsize_pixels(dotsize_index)) ' pix']})

subplot(2,1,2)
hold on
for i = 2:2:num_cases_grad_x
[val,edge] = histcounts(all_errors{i, dotsize_index}.errors_U, 'numbins', num_bins, 'normalization', 'cdf');
center = 0.5 * (edge(1:end-1) + edge(2:end));
area(center, val, 'facecolor', colors(i/2), 'facealpha', 0.25, 'linewidth', 2.0, 'edgecolor', colors(i/2))
end
legend(legend_string_displacement)
grid on
xlabel('error (pix.)')
title({'Effect of displacement, CDF', [num2str(dotsize_pixels(dotsize_index)) ' pix']})

% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'error-histogram-U-displacement');

%%% plot random error vs displacement %%%

% for each seeding density plot random errors
figure_ctr = figure_ctr+1;
figure(figure_ctr);
hold on
colors = ['r', 'g', 'b', 'k', 'y'];

legend_string_dotsize = cell([4,1]);

lower_envelope = min(absolute_random_error(:,3:20), [], 2)';
upper_envelope = max(absolute_random_error(:,3:20), [], 2)';

mean_line = nanmean(absolute_random_error(:,3:20), 2)';
std_line = nanstd(absolute_random_error(:,3:20), [], 2)';
% 
% lower_envelope = mean_line - std_line;
% upper_envelope = mean_line + std_line;


X = [U_ref(:,3)', fliplr(U_ref(:,3)')];
Y = [lower_envelope, fliplr(upper_envelope)];

% for i = 5:5:20
%     plot(U_ref(:,i), absolute_random_error(:,i), [colors(i/5) '*-'])
%     legend_string_seeding{i/5} = ['seeding = ' num2str(dotsize_array(i)) ' dots/32x32 pix.'];
for i = 3:20
    plot(U_ref(:,i), absolute_random_error(:,i), '*')
    
end
plot(U_ref(:,3), mean_line, 'g', 'linewidth', 2.0)
grid on
fill(X,Y,'b', 'facealpha', 0.2, 'edgecolor', 'b', 'linewidth', 2.0)
% legend(legend_string_seeding, 'location', 'southeast')
xlabel('Reference Displacement (pix.)')
ylabel('Random Error (pix.)')
title('Effect of dotsize')
% save results to file
save_figure_to_file(gcf, figure_save_filepath, 'random-error-displacement');

