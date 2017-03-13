clear
close all
clc
% dbstop if errors
% add prana to path
addpath('/home/barracuda/a/lrajendr/Software/prana')

% this is the top filepath
top_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/';

% these are the displacements
grad_x_array=0.5:0.5:5;

% these are the seeding densities
seeding_density_array = 1:20;

% this is the filepath containing the sample parameter file
sample_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/grad_x=0.50/2/processing/';

% this is the sample parameter filename 
sample_filename = 'seeding_gradx05_dens2.mat';

% load a sample parameter file
sample_bos_params = load([sample_filepath sample_filename]);

% loop through all the displacements
for grad_x_index = 1:length(grad_x_array)
    for seeding_density_index = 1:length(seeding_density_array)
        
        grad_x = grad_x_array(grad_x_index);
        seeding_density = seeding_density_array(seeding_density_index);
    
        % display progress to user
        fprintf('grad_x: %0.2f, seeding_density: %d\n', grad_x, seeding_density);
        
        % copy the parameters
        Data = sample_bos_params.Data;
        
        % change the input and output filepaths
        Data.imdirec = [top_filepath 'grad_x=' num2str(grad_x, '%0.2f') '/' num2str(seeding_density) '/processing/cropped-images/'];
        Data.imdirec2 = [top_filepath 'grad_x=' num2str(grad_x, '%0.2f') '/' num2str(seeding_density) '/processing/cropped-images/'];
        Data.outdirec = [top_filepath 'grad_x=' num2str(grad_x, '%0.2f') '/' num2str(seeding_density) '/processing/results/vectors/'];

        % calculate number of cropped images in the folder
        tiff_files = dir([Data.imdirec 'seeding*.tif']);

        % change the number of files to process for this setting
        Data.imfend = num2str(length(tiff_files) - 1);
        
        Data.PIV1.winres='32,32;32,32';
        Data.PIV1.winsize='64,64';
        Data.PIV1.gridres='32,32';
        
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
        save_filename = ['seeding_gradx' num2str(grad_x*10, '%02d') '_dens' num2str(seeding_density) '.mat'];
        
        % this is the filepath to save teh file
        save_filepath = [top_filepath 'grad_x=' num2str(grad_x, '%0.2f') '/' num2str(seeding_density) '/processing/'];
        
        Data.batchname = ['seeding_gradx' num2str(grad_x*10, '%02d') '_dens' num2str(seeding_density)];
        
        % save the parameters to file
        save([save_filepath save_filename], 'Data');
        
        % process images
        pranaPIVcode(Data);
           
    end

end
