% clear
% close all
% clc

% dbstop if errors
% add prana to path
addpath('/home/barracuda/a/lrajendr/Software/prana')

% this is the top filepath
% top_filepath = '/home/sa/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/';
top_filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/dof/';

% these are the displacements
grad_x_array=5.0;

% these are the f numbers
f_number_array = [2.8 4 5.6 8 11 16 22 32 64];

U = zeros(size(f_number_array));

all_errors = cell([1, length(f_number_array)]);

sample_struct = struct('num_total_vectors', [], 'num_valid_vectors', [], 'ref_disp_x', [], 'error', [], 'bias_error', [], 'random_error', [], 'total_error', []);
pixel_pitch = 17e-3;
% loop through all the displacements
for f_number_index = 1:length(f_number_array)

    grad_x = grad_x_array(1);
    f_number = f_number_array(f_number_index);
    grad_x_index=1;
    seeding_density_index = f_number_index;

   % display progress to user
    fprintf('grad_x: %0.2f, f_number: %d\n', grad_x, f_number);

    
    data = load([top_filepath 'grad_x=' num2str(grad_x, '%.2f') '/f_number=' num2str(floor(f_number), '%02d') '/processing/results/vectors/bos_pass1_01.mat']);
    
    data_U = data.U(:,:,1);
    U(f_number_index) = -mean(data_U(:));
    
                
        all_errors{grad_x_index, seeding_density_index} = sample_struct;
        
        all_errors{grad_x_index, seeding_density_index}.grad_x = grad_x;
        
        
        
        %% set read  and write paths



        % location where the light ray positions are stored
        light_ray_positions_filepath = [top_filepath 'grad_x=' num2str(grad_x, '%.2f') '/f_number=' num2str(floor(f_number), '%02d') '/1/light-ray-positions/'];

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
        U_ref(f_number_index) = del_x_pixels;
        V_ref(f_number_index) = del_y_pixels;

        % close files
        fclose(fid_1);
        fclose(fid_2);
        
        all_errors{grad_x_index, f_number_index}.U_ref = del_x_pixels;
        all_errors{grad_x_index, f_number_index}.V_ref = del_y_pixels;
        
  
    
end

U

plot(f_number_array(3:end), U(3:end), 'b*')
hold on
plot(f_number_array(3:end), U_ref(3:end), 'r*')
xlabel('f_number')
ylabel('\Delta x')
