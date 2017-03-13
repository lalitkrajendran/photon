clear
close all
clc
% dbstop if errors
% add prana to path
addpath('/home/barracuda/a/lrajendr/Software/prana')

% this is the top filepath
top_filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/seeding/';

% these are the displacements
delta_x_array = 0.1:0.2:1;
delta_y_array = 0.1:0.2:1;

% these are the seeding densities
seeding_density_array = 4:4:20;

% this is the filepath containing the sample parameter file
sample_filepath = '/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/seeding/delta_x=0.10_delta_y=0.10/3/processing/';

% this is the sample parameter filename 
sample_filename = 'seeding_delta_x_01_delta_y_01_dens_3.mat';

% load a sample parameter file
sample_bos_params = load([sample_filepath sample_filename]);

% loop through all the displacements
for delta_x_index = 1:length(delta_x_array)
    for seeding_density_index = 1:length(seeding_density_array)

        delta_x = delta_x_array(delta_x_index);
        delta_y = delta_y_array(delta_x_index);
        
        seeding_density = seeding_density_array(seeding_density_index);
        if(seeding_density == 19)
            continue
        end
        % display progress to user
        fprintf('delta_x: %0.2f, seeding_density: %d\n', delta_x, seeding_density);
        
        % copy the parameters
        Data = sample_bos_params.Data;
        
        % change the input and output filepaths
        Data.imdirec = [top_filepath 'delta_x=' num2str(delta_x, '%0.2f') '_delta_y=' num2str(delta_x, '%0.2f') '/' num2str(seeding_density) '/processing/cropped-images/'];
        Data.imdirec2 = [top_filepath 'delta_x=' num2str(delta_x, '%0.2f') '_delta_y=' num2str(delta_x, '%0.2f') '/' num2str(seeding_density) '/processing/cropped-images/'];
%         Data.outdirec = [top_filepath 'delta_x=' num2str(delta_x, '%0.2f') '_delta_y=' num2str(delta_x, '%0.2f') '/' num2str(seeding_density) '/processing/results/vectors/'];
        Data.outdirec = [top_filepath 'delta_x=' num2str(delta_x, '%0.2f') '_delta_y=' num2str(delta_x, '%0.2f') '/' num2str(seeding_density) '/processing-16x16/results/vectors/'];

        % calculate number of cropped images in the folder
        tiff_files = dir([Data.imdirec 'seeding*.tif']);

        % change the number of files to process for this setting
        Data.imfend = num2str(length(tiff_files) - 1);
        
        Data.PIV1.winres='16,16;16,16';
        Data.PIV1.winsize='32,32';
        Data.PIV1.gridres='16,16';
        
        
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
        save_filename = ['seeding_delta_x_' num2str(delta_x*10, '%02d') '_delta_y_' num2str(delta_y*10, '%02d') '_dens_' num2str(seeding_density) '.mat'];
        
        % this is the filepath to save teh file
        save_filepath = [top_filepath 'delta_x=' num2str(delta_x, '%0.2f') '_delta_y=' num2str(delta_x, '%0.2f') '/' num2str(seeding_density) '/processing-16x16/'];
        % if the folder does not does not exist, then create it
        if(~exist(save_filepath, 'dir'))
            mkdir(save_filepath);
        end

        Data.batchname = ['seeding_delta_x_' num2str(delta_x*10, '%02d') '_delta_y_' num2str(delta_y*10, '%02d') '_dens_' num2str(seeding_density)];
        
        % save the parameters to file
        save([save_filepath save_filename], 'Data');
        
        % process images
        pranaPIVcode(Data);
           
    end

end
