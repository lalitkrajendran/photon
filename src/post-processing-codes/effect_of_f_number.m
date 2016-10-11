% the purpose of this code is to study the effect of the camera f number on
% the diameter of the dot pattern

% clear workspace
clear 
close all
clc

% this is the case_name
case_name = 'f-number-study';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute dot sizes from theory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% these are the physical diameters (mm) of the dots for which the images were
% generated
% dot_diameters_physical = [50 75 100 200 300 400 500 600 700 800 900 1000];
% dot_diameters_physical = [10 20 30 40 50 100 150 200 250 300 350 400 450 500];
dot_diameters_physical = [10 20 30 40 50 60 70 80 90 100];

% this is the magnification
M = 123500./700000;

% this is the pixel pitch (mm)
pixel_pitch = 17e-3;

% calculate dot diameter in the image (pixels) from theory
dot_diameters_theory = dot_diameters_physical*1e-3 * M / pixel_pitch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute dot sizes from image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is the range of f numbers considered
f_numbers = [4 8 11 16 22];

% this is the top folder where all the images are located
img_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/' case_name '/'];

% this structure holds the image diameter and the corresponding physical
% dot diameter for each image in each f-number considered
all_dot_diameters_physical = cell(length(f_numbers), 1);
all_dot_diameters_image = cell(length(f_numbers), 1);

all_dot_diameters = cell(length(f_numbers), 1);


for i = 1:length(f_numbers)
    
    % this is the name of the sub-folder corresponding to the current
    % f-number
    folder_name = ['f' num2str(f_numbers(i)) '/'];
    
    % for each f-number, populate the list of dot sizes for which images
    % were generated
    folder_list = dir([img_filepath folder_name]);
    
    % ignore entries like .DS_store that are not really folders
    folder_list = folder_list([folder_list.isdir]);
    
    % also remove entries like . or ..
    folder_list = folder_list(~strncmpi('.',{folder_list.name},1));
    
    fprintf('f#: %f, %d folders found\n', f_numbers(i), length(folder_list));
    
    % initialize an array to store all the physical dot sizes for which
    % images were generated
    dot_diameters_physical = zeros(length(folder_list),1);
    
    % initialize array to store all the image dot sizes
    dot_diameters_image = zeros(length(folder_list), 1);
    
    % loop through all the files to retrieve the images and note the dot
    % size
    for j = 1:length(folder_list)
        fprintf('processing %d out of %d dot sizes\n', j, length(folder_list));
        
        % get dot size from folder name in microns
        dot_diameters_physical(j) = str2double(folder_list(j).name);
        
        % load first image from this folder
        img = imread([img_filepath folder_name folder_list(j).name '/bos_pattern_image_1.tif']);
        
        % consider only the center part of the image where the particles are in
        % focus
        image_margin = [300 700];
        img = img(image_margin(1):image_margin(2), image_margin(1):image_margin(2));
        
        % calculate diameter from prana

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % identify the particles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % define settings for particle identification
        IDmethod = {'blob','dynamic','combined'};
        Data.ID.method         = IDmethod{2};
        Data.ID.run            = 1;
        Data.ID.v              = 10;
        Data.ID.contrast_ratio = 0;

        particleIDprops       = Data.ID;

        % call function to id the particles
        [ID1.p_matrix,ID1.peaks,ID1.num_p]=particle_ID(img,particleIDprops);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % size the particles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % define settings for particle sizing
        Data.Size.run      = 1;
        Data.Size.thresh   = 10;
        Data.Size.method   = 'CLSG';
        Data.Size.p_area   = 0;
        Data.Size.sigma    = 4;
        Data.Size.errors   = 1;%str2double(Data.Size.errors);

        sizeprops       = Data.Size;
        
        % call function to size the particles
        [SIZE1.XYDiameter,SIZE1.mapsizeinfo,SIZE1.locxy]=particle_sizing(img,ID1.p_matrix,...
                            ID1.num_p,sizeprops);

        % extract co-ordinate, size and intensity properties from the results             
        X1=SIZE1.XYDiameter(:,1);
        Y1=SIZE1.XYDiameter(:,2);                
        Dp1=SIZE1.XYDiameter(:,3);
        Ip1=SIZE1.XYDiameter(:,4);                 

        % remove NaN values
        Dp1 = Dp1(~isnan(Dp1));

        % calculate the average dot diameter (pixels)
        dot_diameters_image(j) = nanmean(Dp1);

    end

    % store diameters in a cell array
    all_dot_diameters{i} = struct('f_number', f_numbers(i), 'physical', ...
        dot_diameters_physical,'image', dot_diameters_image);    

end

% save data to file
save([img_filepath 'diameter-results.mat'], 'all_dot_diameters');

% plot results
figure
hold on
c = ['r', 'b', 'g', 'c', 'm'];

for i = 1:length(all_dot_diameters)
    % plot theoretical vs physical
    % calculate theoretical diameter
    d_theory = all_dot_diameters{i}.physical*1e-3 * M / pixel_pitch;
    p0 = plot(all_dot_diameters{i}.physical, d_theory, 'k');
    
    % plot image vs physical
    p(i) = plot(all_dot_diameters{i}.physical, all_dot_diameters{i}.image, [c(i) '*']); 
end

grid on

legend([p0 p], 'theory', ['f' num2str(f_numbers(1))], ['f' num2str(f_numbers(2))], ['f' num2str(f_numbers(3))],...
    ['f' num2str(f_numbers(4))], ['f' num2str(f_numbers(5))])
xlabel('physical diameter (microns)');
ylabel('image diameter (pixels)');
title('effect of f-number on dot size');

% save workspace variables to file
save([img_filepath 'workspace.mat']);