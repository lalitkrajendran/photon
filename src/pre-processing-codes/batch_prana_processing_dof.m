clear
close all
clc
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

% this is the filepath containing the sample parameter file
sample_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/grad_x=0.50/2/processing/';

% this is the sample parameter filename 
sample_filename = 'seeding_gradx05_dens2.mat';

% load a sample parameter file
sample_bos_params = load([sample_filepath sample_filename]);

% loop through all the displacements
for f_number_index = 1:length(f_number_array)

    grad_x = grad_x_array(1);
    f_number = f_number_array(f_number_index);
    
    % display progress to user
    fprintf('grad_x: %0.2f, f_number: %d\n', grad_x, f_number);

    % copy the parameters
    Data = sample_bos_params.Data;

    % change the input and output filepaths
    Data.imdirec = [top_filepath 'grad_x=' num2str(grad_x, '%0.2f') '/f_number=' num2str(floor(f_number), '%02d') '/processing/cropped-images/'];
    Data.imdirec2 = [top_filepath 'grad_x=' num2str(grad_x, '%0.2f') '/f_number=' num2str(floor(f_number), '%02d') '/processing/cropped-images/'];
    Data.outdirec = [top_filepath 'grad_x=' num2str(grad_x, '%0.2f') '/f_number=' num2str(floor(f_number), '%02d') '/processing/results/vectors/'];

    Data.imbase = 'dof_';
    Data.imbase2 = 'dof_';
    
    % calculate number of cropped images in the folder
    tiff_files = dir([Data.imdirec 'dof*.tif']);

    % change the number of files to process for this setting
    Data.imfend = num2str(length(tiff_files) - 1);

%         %%% change settings for saving peaks %%%
%         % specify that validation is not required for the third pass
%         Data.PIV3.val = '0';
%         % specify that additional peaks have to be saved
%         Data.PIV3.savepeakinfo = '1';
%         % specify that Peaks 1 and 2 have to be saved
%         Data.PIV3.corrpeaknum = '2';
%         % specify that the Peak magnitude should be saved
%         Data.PIV3.savepeakmag = '1';
%         % specify that the resulting velocities should be saved
%         Data.PIV3.savepeakvel = '1';

    % this is the name to save the new parameters
    save_filename = ['f_number_gradx' num2str(grad_x*10, '%02d') '_f' num2str(floor(f_number), '%02d') '.mat'];

    % this is the filepath to save teh file
    save_filepath = [top_filepath 'grad_x=' num2str(grad_x, '%0.2f') '/f_number=' num2str(floor(f_number), '%02d') '/processing/'];

    Data.batchname = ['f_number_gradx' num2str(grad_x*10, '%02d') '_f' num2str(floor(f_number), '%02d')];

    % save the parameters to file
    save([save_filepath save_filename], 'Data');

    % process images
    pranaPIVcode(Data);

end
