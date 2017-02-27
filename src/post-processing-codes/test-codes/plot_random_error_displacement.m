% for each seeding density plot random errors
figure_ctr = figure_ctr+1;
figure(figure_ctr);
hold on
colors = ['r', 'g', 'b', 'k', 'y'];

legend_string_seeding = cell([4,1]);

lower_envelope = min(absolute_random_error(:,3:20), [], 2)';
upper_envelope = max(absolute_random_error(:,3:20), [], 2)';

% mean_line = nanmean(absolute_random_error(:,3:20), 2)';
% std_line = nanstd(absolute_random_error(:,3:20), [], 2)';
% 
% lower_envelope = mean_line - std_line;
% upper_envelope = mean_line + std_line;


X = [U_ref(:,3)', fliplr(U_ref(:,3)')];
Y = [lower_envelope, fliplr(upper_envelope)];

% for i = 5:5:20
%     plot(U_ref(:,i), absolute_random_error(:,i), [colors(i/5) '*-'])
%     legend_string_seeding{i/5} = ['seeding = ' num2str(seeding_density_array(i)) ' dots/32x32 pix.'];
for i = 3:20
    plot(U_ref(:,i), absolute_random_error(:,i), 'r*')
    
end

grid on
fill(X,Y,'b', 'facealpha', 0.2, 'edgecolor', 'b', 'linewidth', 2.0)
% legend(legend_string_seeding, 'location', 'southeast')
xlabel('Reference Displacement (pix.)')
ylabel('Random Error (pix.)')
title('Effect of seeding')