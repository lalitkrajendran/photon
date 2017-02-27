close all
figure
hold on

cmap = parula;

absolute_bias_error = sqrt(err_U_bias.^2 + err_V_bias.^2);
ellipse_axes_ratio = err_V_bias./err_U_bias;

min_error_color = min(absolute_bias_error(:));
max_error_color = max(absolute_bias_error(:)) * 0.5;


caxis([min_error_color, max_error_color])

for grad_x_index = 1:num_cases_grad_x
    for seeding_density_index = 1:num_cases_seeding_density
        
        grad_x = grad_x_array(grad_x_index);
        seeding_density = seeding_density_array(seeding_density_index);
        
        current_error = absolute_bias_error(grad_x_index, seeding_density_index);
        color_id = ceil((current_error-min_error_color)/(max_error_color - min_error_color) * size(cmap,1));
        
%         color_id
        
        if color_id == 0
            color_id = 1;
        elseif isnan(color_id)
            continue
        elseif color_id > size(cmap,1)
            color_id = size(cmap,1);
        end
        
        
        color_spec = cmap(color_id,:);
        
        scaling_factor = 0.2;
        
        rx = scaling_factor * 1/sqrt(1+ellipse_axes_ratio(grad_x_index, seeding_density_index)^2);
        ry = scaling_factor * 1/sqrt(1+1/ellipse_axes_ratio(grad_x_index, seeding_density_index)^2);
        
        create_filled_ellipse(seeding_density, grad_x, rx, ry,color_spec)
        
    end
end

colorbar
xlabel('Dots (32x32 pix.)');
ylabel('Reference Displacement (pix.)');
title('Bias Error')
set(gca, 'fontsize', 14)
