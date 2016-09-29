% This program loads the results of a stereo PIV evaluation from prana and
% plots contours of the three components of the velocities, as well as a
% histogram. It also saves these plots to a file

% Author : Lalit Rajendran


%% clear workspace and close all windows

clear all
close all
clc

%% specify filepaths for reading velocity data and saving images.

% this is the path to the folder where the results are located
read_filepath = '~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/50x50/results/';
% this is the name of the file which contains the results
read_filename = 'BOS_pass2_05.mat';

% this is the path to the folder where the plots will be saved
write_filepath = '~/Projects/camera_simulation/results/plots/';
% this is the file extension for saving the image (e.g. tiff, png, jpeg)
write_file_format = 'png';

% these are the reference velocities to be compared against
U_ref = 9.78;
V_ref = 0;


% this is the threshold percentage for the colorbar range. the range will 
% be set to the reference value +/- the threshold percentage
threshold_percentage = 10;

% read data from file
results = load([read_filepath read_filename]);

%% plot contours and histograms

% plot U velocity contour
figure(1)
imagesc(results.U)
colorbar
caxis([U_ref - threshold_percentage/100 * abs(U_ref), U_ref + threshold_percentage/100 * abs(U_ref)])
title(['U (U_{ref} = ' num2str(U_ref) ')'])

% save figure to file
print([write_filepath 'U.' write_file_format], ['-d' write_file_format]);

% plot V velocity contour
figure(2)
imagesc(results.V)
colorbar
caxis([V_ref - threshold_percentage/100 * abs(V_ref), V_ref + threshold_percentage/100 * abs(V_ref)])
title(['V (V_{ref} = ' num2str(V_ref) ')'])

% save figure to file
print([write_filepath 'V.' write_file_format], ['-d' write_file_format]);

% % plot W velocity contour
% figure(3)
% imagesc(results.W)
% colorbar
% caxis([W_ref - threshold_percentage/100 * abs(W_ref), W_ref + threshold_percentage/100 * abs(W_ref)])
% title(['W (W_{ref} = ' num2str(W_ref) ')'])
% 
% % save figure to file
% print([write_filepath 'W.' write_file_format], ['-d' write_file_format]);

% plot histograms of all velocities
figure(4)
subplot(2,1,1), histogram(results.U(:))
title(['U (U_{ref} = ' num2str(U_ref) ')'])
subplot(2,1,2), histogram(results.V(:))
title(['V (V_{ref} = ' num2str(V_ref) ')'])
% subplot(3,1,3), histogram(results.W(:))
% title(['W (W_{ref} = ' num2str(W_ref) ')'])

% save figure to file
print([write_filepath 'histogram.' write_file_format], ['-d' write_file_format]);




