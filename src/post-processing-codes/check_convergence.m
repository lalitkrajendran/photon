% check whether the results of the error analysis are converged

image_pairs = 1:num_files_per_seeding;


mean_bias_error.U = zeros(1,num_files_per_seeding);
mean_bias_error.V = zeros(1,num_files_per_seeding);
mean_random_error.U = zeros(1,num_files_per_seeding);
mean_random_error.V = zeros(1,num_files_per_seeding);
mean_rms_error.U = zeros(1,num_files_per_seeding);
mean_rms_error.V = zeros(1,num_files_per_seeding);

for image_pair_index = 1:num_files_per_seeding
    
    mean_bias_error.U = mean(err_U_bias_single(1:image_pair_index));
    mean_bias_error.V = mean(err_V_bias_single(1:image_pair_index));
    
    mean_random_error.U = mean(err_U_random_single(1:image_pair_index));
    mean_random_error.V = mean(err_V_random_single(1:image_pair_index));
    
    mean_rms_error.U = mean(err_U_rms_single(1:image_pair_index));
    mean_rms_error.V = mean(err_V_rms_single(1:image_pair_index));
    
end

figure

subplot(1,3,1)
hold on
plot(image_pairs, mean_bias_error.U, 'r*-')
plot(image_pairs, mean_bias_error.V, 'b*-')
xlabel('image pairs')
ylabel('error (pixels)')
title('bias error')

subplot(1,3,2)
hold on
plot(image_pairs, mean_random_error.U, 'r*-')
plot(image_pairs, mean_random_error.V, 'b*-')
xlabel('image pairs')
ylabel('error (pixels)')
title('random error')

subplot(1,3,3)
hold on
plot(image_pairs, mean_rms_error.U, 'r*-')
plot(image_pairs, mean_rms_error.V, 'b*-')
xlabel('image pairs')
ylabel('error (pixels)')
title('rms error')

