figure, imshow(img(r_min:r_max, c_min:c_max), 'InitialMagnification', 'fit')

r_N = del_r;
c_N = del_c;

r_0 = round(rc - r_min);
c_0 = round(cc - c_min);

% plot intensity profile of the selected dot
figure_ctr = figure_ctr + 1;
figure(figure_ctr)

% plot the profile along a row
subplot(1,2,1)
hold on
plot(I_dot(r_0,:),'-*')
plot(gaussian_fit(r_0,:), 'r-')
title('center row')

% plot the profile along a column
subplot(1,2,2)
hold on
plot(I_dot(:,c_0),'-*')
plot(gaussian_fit(:,c_0),'r-')
title('center column')

% % plot the full profile
% subplot(1,3,3)
% surf(I_dot_average, 'linestyle','none');
% view(2)
% title('full')