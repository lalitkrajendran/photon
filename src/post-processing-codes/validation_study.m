% purpose is to check how well the light ray deflections agree with theory.

clear 
close all
clc

figure_ctr = 0;

%% overall settings

% this is the grad_x value that represents the displacement
grad_x_array = 0.5:0.5:5;

% this is the number of grad_x cases
num_cases_grad_x = length(grad_x_array);

% this ensures all the figures are docked and displayed from a single
% window
set(0,'DefaultFigureWindowStyle','docked')

%% read and write paths

% this is the folder where the figures will be saved
figure_save_filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/results/validation/plots/';

% this is the folder where the workspace variables will be saved
workspace_save_filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/results/validation/without-error/';

% this is the top read directory where all the images and light ray
% positions are stored
top_read_filepath{1} = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/validation/with-error/';
top_read_filepath{2} = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/validation/without-error/';

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
displacements_theory = grad_x_array * 2/5;

% size of a single pixel (mm)
pixel_pitch = 17e-3;

%% calculate reference displacements from light ray positions

% define array to store light ray deflections
displacements_ray_tracing = zeros(2,length(grad_x_array));

for grad_x_index = 1:length(grad_x_array)
    
    % get value of grad_x
    grad_x = grad_x_array(grad_x_index);
    fprintf('grad_x=%f\n', grad_x);

    for i = 1:2
        % location where the light ray positions are stored
        light_ray_positions_filepath = [top_read_filepath{i} 'grad_x=' num2str(grad_x, '%0.2f') ...
                                    '/1/light-ray-positions/'];

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
        displacements_ray_tracing(i,grad_x_index) = del_x_pixels;

        % close files
        fclose(fid_1);
        fclose(fid_2);
    end

end

% calculate error
% displacements_error(1,:) = displacements_ray_tracing(1,:) - displacements_theory;

%% plot everything
figure
hold on

plot(displacements_theory, 'r')
plot(displacements_ray_tracing(1,:), 'g')
plot(displacements_ray_tracing(2,:), 'b')

legend('theory', 'with error', 'without error')