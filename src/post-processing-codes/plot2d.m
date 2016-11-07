% This program loads the results of a PIV evaluation from prana and
% plots contours of the two components of the velocities, as well as a
% histogram. It also saves these plots to a file

% Author : Lalit Rajendran


%% clear workspace and close all windows

clear 
close all
clc

case_name = '150x150-f16-disp12';
%% specify filepaths for reading velocity data and saving images.

% this is the path to the folder where the results are located
read_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/vectors/'];
% this is the name of the file which contains the results
read_filename = 'BOS_pass1_03.mat';

% this is the path to the folder where the plots will be saved
write_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/plots/'];
% this is the file extension for saving the image (e.g. tiff, png, jpeg)
write_file_format = 'png';

% these are the reference velocities to be compared against
U_ref = 11.64;
V_ref = 0;


% this is the threshold percentage for the colorbar range. the range will 
% be set to the reference value +/- the threshold percentage
threshold_percentage = 20;

% read data from file
results = load([read_filepath read_filename]);

%% plot contours and histograms

% plot U velocity contour
figure(1)
imagesc(results.U)
colorbar
caxis([U_ref - threshold_percentage/100 * abs(U_ref), U_ref + threshold_percentage/100 * abs(U_ref)])
title(['\Delta x (\Delta x_{ref} = ' num2str(U_ref) ')'])

% save figure to file
print([write_filepath 'U.' write_file_format], ['-d' write_file_format]);

% plot V velocity contour
figure(2)
imagesc(results.V)
colorbar
caxis([V_ref - threshold_percentage/100 * abs(V_ref), V_ref + threshold_percentage/100 * abs(V_ref)])
title(['\Delta y (\Delta y_{ref} = ' num2str(V_ref) ')'])

% save figure to file
print([write_filepath 'V.' write_file_format], ['-d' write_file_format]);

% plot histograms of all velocities
figure(4)
subplot(2,1,1), histogram(results.U(:))
title(['\Delta x (\Delta x_{ref} = ' num2str(U_ref) ')'])
subplot(2,1,2), histogram(results.V(:))
title(['\Delta y (\Delta y_{ref} = ' num2str(V_ref) ')'])

% save figure to file
print([write_filepath 'histogram.' write_file_format], ['-d' write_file_format]);

figure(5)
quiver(results.X, results.Y, results.U, results.V)




