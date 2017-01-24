clear
close all
clc

figure_ctr = 0;

% set top read location
top_read_directory = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/';

% set figure write directory
figure_write_directory = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/plots/';

% set displacement array
displacements = [1.0, 2.0, 3.0];

% open a figure, and set the hold on or off
figure_ctr = figure_ctr+1;
figure(figure_ctr);
% set(gcf, 'Position', [200 200 900 500])

% array of colors for each plot
colors = ['r', 'g', 'b'];
% loop through each case and do the following
for displacement_index = 1:length(displacements)
    
    % load workspace file
    workspace_filepath = [top_read_directory '50x50-f16-grad_x' num2str(displacements(displacement_index), '%.1f') '/results/'];
    
    error_analysis_results = load([workspace_filepath 'error_analysis.mat']);
    
    % extract parameters required for plotting error
    if(displacement_index==2)
        error_analysis_results.starting_index = 9;
    end
    
    
    % plot errors
 
    % dummy plot for legend
%     subplot(1,4,1)
% %     figure('visible', off)
%     set(gca, 'Visible', 'off')
%     hold on
%     plot(error_analysis_results.dot_diameters_image_pixels(error_analysis_results.starting_index:end),...
%         error_analysis_results.err_U_bias(error_analysis_results.starting_index:end), [colors(displacement_index),'*-']);
%     plot(error_analysis_results.dot_diameters_image_pixels(error_analysis_results.starting_index:end),...
%         error_analysis_results.err_V_bias(error_analysis_results.starting_index:end), 'b*-');
%     grid on
%     xlim([0 max(error_analysis_results.dot_diameters_image_pixels)*1.2])
%     xlabel('dot diameter (pixels)');
%     ylabel('error (pixels)');
% %     legend('\Delta x', '\Delta y', 'location', 'Northwest');
%     title('bias error')
%     set(gca, 'fontsize', 14)
    
    % bias error
    subplot(1,3,1)
    hold on
    plot(error_analysis_results.dot_diameters_image_pixels(error_analysis_results.starting_index:end),...
        error_analysis_results.err_U_bias(error_analysis_results.starting_index:end), [colors(displacement_index),'*-']);
%     plot(error_analysis_results.dot_diameters_image_pixels(error_analysis_results.starting_index:end),...
%         error_analysis_results.err_V_bias(error_analysis_results.starting_index:end), 'b*-');
    grid on
    xlim([0 max(error_analysis_results.dot_diameters_image_pixels)*1.2])
    xlabel('dot diameter (pixels)');
    ylabel('error (pixels)');
%     legend('\Delta x', '\Delta y', 'location', 'Northwest');
    title('bias error')
    set(gca, 'fontsize', 14)
    
    % random error
    subplot(1,3,2)
    hold on
    plot(error_analysis_results.dot_diameters_image_pixels(error_analysis_results.starting_index:end),...
        error_analysis_results.err_U_random(error_analysis_results.starting_index:end), [colors(displacement_index),'*-']);
%     plot(error_analysis_results.dot_diameters_image_pixels(error_analysis_results.starting_index:end),...
%         err_V_random(error_analysis_results.starting_index:end), 'bo-');
    grid on
    xlim([0 max(error_analysis_results.dot_diameters_image_pixels)*1.2])
    xlabel('dot diameter (pixels)');
    ylabel('error (pixels)');
%     legend('\Delta x', '\Delta y', 'location', 'Northwest');
    title('random error')
    set(gca, 'fontsize', 14)
    
    % total error
    subplot(1,3,3)
    hold on
    plot(error_analysis_results.dot_diameters_image_pixels(error_analysis_results.starting_index:end),...
        error_analysis_results.err_U_rms(error_analysis_results.starting_index:end), [colors(displacement_index),'*-']);
%     plot(error_analysis_results.dot_diameters_image_pixels(error_analysis_results.starting_index:end),...
%         error_analysis_results.err_V_rms(error_analysis_results.starting_index:end), 'bs-');
    grid on
    xlim([0 max(error_analysis_results.dot_diameters_image_pixels)*1.2])
    xlabel('dot diameter (pixels)');
    ylabel('error (pixels)');
%     legend('\Delta x', '\Delta y', 'location', 'Northwest');
    title('total error')
    set(gca, 'fontsize', 14)

end

% add legend
subplot(1,3,3)
legend(['\Delta x = ', num2str(displacements(1), '%0.1f'), ' pix.'], ['\Delta x = ', num2str(displacements(2), '%0.1f'), ' pix.'],...
['\Delta x = ', num2str(displacements(3), '%0.1f'), ' pix.']); %, 'location', 'northwestoutside')
% 
% subplot(1,3,2)
% legend(['\Delta x = ', num2str(displacements(1), '%0.1f'), ' pix.'], ['\Delta x = ', num2str(displacements(2), '%0.1f'), ' pix.'],...
% ['\Delta x = ', num2str(displacements(3), '%0.1f'), ' pix.'])

% subplot(1,3,3)
% legend(['\Delta x = ', num2str(displacements(1), '%0.1f'), ' pix.'], ['\Delta x = ', num2str(displacements(2), '%0.1f'), ' pix.'],...
% ['\Delta x = ', num2str(displacements(3), '%0.1f'), ' pix.'], 'location', 'northeastoutside')

set(gca, 'fontsize', 14)

% save figure to file
save_figure_to_file(gcf, figure_write_directory, 'errors_dotsize_all');