figure
hold on
ctr = 0;

% for dotsize_index = 1:4:length(dotsize_array)
%     plot(U_ref(:,dotsize_index), err_U_random(:,dotsize_index), '-*', 'linewidth', 2.0)
% end

legend_string = cell([length(grad_x_array), 1]);
ctr = 0;

for grad_x_index = 1:2:length(grad_x_array)
    ctr = ctr+1;
    plot(dotsize_pixels, err_U_total(grad_x_index,:), '-*', 'linewidth', 2.0)
    legend_string{ctr,1} = ['\Delta x=' num2str(U_ref(grad_x_index, 1)) ' pix.'];
end

legend(legend_string{1:ctr})
% figure
% hold on
% legend_string = cell([length(seeding_density_array), 1]);
% ctr = 0;
% for seeding_density_index = 1:2:length(seeding_density_array)
%     ctr = ctr + 1;
%     plot(U_ref(:,seeding_density_index), abs(err_U_random(:,seeding_density_index)), '-*', 'linewidth', 2.0)
% %     plot(delta_x_array, abs(err_V_bias(:,seeding_density_index)), '-o')
%     legend_string{ctr,1} = ['seeding=' num2str(seeding_density_array(seeding_density_index)) 'dots/32x32 pix.'];
% end
% grid on
% legend(legend_string{1:ctr})
% xlabel('displacement (pix.)')
% ylabel('error (pix.)')
% title('random error')
% set(gca, 'fontsize', 14)
% xlim([0.0 1.5])
% set(gca, 'xtick', 0:0.1:1.5)
% set(gca, 'ytick', 0.06:0.01:0.09)
% save_figure_to_file(gcf, figure_save_filepath, 'random-error-U-line')

% figure, imshow(img(r_min:r_max, c_min:c_max), 'InitialMagnification', 'fit')
% 
% r_N = del_r;
% c_N = del_c;
% 
% r_0 = round(rc - r_min);
% c_0 = round(cc - c_min);
% 
% % plot intensity profile of the selected dot
% figure_ctr = figure_ctr + 1;
% figure(figure_ctr)
% 
% % plot the profile along a row
% subplot(1,2,1)
% hold on
% plot(I_dot(r_0,:),'-*')
% plot(gaussian_fit(r_0,:), 'r-')
% title('center row')
% 
% % plot the profile along a column
% subplot(1,2,2)
% hold on
% plot(I_dot(:,c_0),'-*')
% plot(gaussian_fit(:,c_0),'r-')
% title('center column')
% 
% % % plot the full profile
% % subplot(1,3,3)
% % surf(I_dot_average, 'linestyle','none');
% % view(2)
% % title('full')