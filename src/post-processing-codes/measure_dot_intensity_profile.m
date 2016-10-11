
%% clear workspace and close all windows

clear all
close all
clc

figure_ctr = 0;

case_name = '50x50-f16';

%% set read and write paths

% this is the path to the folder containing the raw images
img_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/reordered-images/'];

% this is the folder where the particle identification results will be
% stored
particle_id_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/particle-id/'];

% this is the folder where the particle sizing results will be
% stored
particle_size_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/particle-size/'];

% this is the folder where the figures will be saved
figure_save_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/plots/'];

% this creates the write directories if they are not already present
if ~exist(particle_id_filepath,'dir')
    mkdir(particle_id_filepath);
end

if ~exist(particle_size_filepath,'dir')
    mkdir(particle_size_filepath);
end

if ~exist(figure_save_filepath,'dir')
    mkdir(figure_save_filepath);
end

%% load images

% these are the names of the file to be read
files = dir([img_filepath '*.tif']);

% this is the number of files
N = length(files);
fprintf('number of files in the directory: %d\n', N);

% this sets which image has to be loaded
for i =15:5:35
%     i = 30;
    
    % this loads the image
    img = imread([img_filepath files(i).name]);

    % consider only the center part of the image where the particles are in
    % focus
    image_margin = [256 768];
    img = img(image_margin(1):image_margin(2), image_margin(1):image_margin(2));

    %% locate dots and measure diameter

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % imfindcircles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % find location and radius of dots
    [centers, radii] = imfindcircles(img, [2 6], 'Sensitivity', 0.95);

    % % remove NAN elements
    % centers = centers(~isnan(centers));
    % radii = radii(~isnan(radii(:)));

    % calculate average diameter of the dot
    average_diameter = 2 * nanmean(radii);

    % calculate the total number of dots that were identified
    number_of_dots = sum(~isnan(radii(:)));

    fprintf('Number of dots: %d\n', number_of_dots);
    fprintf('Average diameter: %f\n', average_diameter);

%     % check if dots are identified correctly
%     figure_ctr = figure_ctr + 1;
%     figure(figure_ctr)
%     hold on
%     imshow(img)
%     viscircles(centers, radii)

    %% extract the intensity profile of a dot

    % set the max size of the matrix that will store the intensity profile of 
    % the region surrounding a dot 
    max_size = round(1.5*average_diameter); %10;
    fprintf('max size: %d\n', max_size);
    
    % initialize an array that will store a running total of the intensity
    % profiles of the selected dots
    I_dot_total = zeros(max_size,max_size);

    % intialize dot counter to zero
    dot_ctr = 0;

    for j = 1:number_of_dots
        % select the dot whose intensity profile will be plotted
        dot_index = j;

        % extract the location of the dot in the image
        dot_loc = centers(dot_index,:);

        % convert the dot location to an integer value
        dot_loc = round(dot_loc);

        % the center is stored as x,y but the array values are stored as r,c. so
        % convert the dot location to an r,c format
        temp = dot_loc;
        dot_loc(1) = temp(2);
        dot_loc(2) = temp(1);

    %     % extract the radius of the dot
    %     dot_radius = radii(dot_index);
    % 
    %     % set the factor by which the radius will be incremented to extract a
    %     % region that is slightly bigger than the size of the dot to ensure that
    %     % the intensity goes to zero at the edges
    %     dot_radius_increment = 2;
    % 
    %     % increment the dot radius
    %     dot_radius = dot_radius + dot_radius_increment;
    % 
    %     % if the dot radius is greater than the maximum allowable size, reduce it
    %     if dot_radius > max_size/2
    %         dot_radius = floor(max_size/2);
    %     end
    % 
    %     % convert the dot radius to an integer value
    %     dot_radius = round(dot_radius);
    % 
    %     % calculate the minimum and maximum indices of the region surrounding the 
    %     % dot that will be extracted from the image
    %     r_min = dot_loc(1) - dot_radius;
    %     r_max = dot_loc(1) + dot_radius;
    % 
    %     c_min = dot_loc(2) - dot_radius;
    %     c_max = dot_loc(2) + dot_radius;


        % calculate the minimum and maximum indices of the region surrounding the 
        % dot that will be extracted from the image
        r_min = dot_loc(1) - max_size/2;
        r_max = dot_loc(1) + max_size/2 - 1;

        c_min = dot_loc(2) - max_size/2;
        c_max = dot_loc(2) + max_size/2 - 1;

        % if any of the indices are out of bounds, ignore this dot
        if r_min <= 0 || c_min <= 0 || r_max > size(img,1) || c_max > size(img,2)
            continue;
        end

        % extract the intensity profile of the selected dot
        I_dot = double(img(r_min:r_max, c_min:c_max));

        % calculate a running total of the intensity profile of the dots
        I_dot_total = I_dot_total + I_dot;

        % update the number of dots whose intensity profiles have been
        % extracted
        dot_ctr = dot_ctr + 1;

    end

    fprintf('number of dot intensity profiles: %d\n', dot_ctr);

    % calculate the average intensity profile of the dot
    I_dot_average = I_dot_total/dot_ctr;

    % find the size of the extracted array
    [r_N,c_N] = size(I_dot_average);

    % fit a gaussian to the row-wise intensity profile
    [fit_r, goodness_r] = fit([1:c_N]',I_dot_average(round(r_N/2),:)','gauss1');
    % fit a gaussian to the column-wise intensity profile
    [fit_c, goodness_c] = fit([1:r_N]',I_dot_average(:,round(c_N/2)),'gauss1');

    fprintf('Goodness of fit details-----------------\n');
    fprintf('i = %d\n', i);
    fprintf('rsquare - row: %f, col: %f\n', goodness_r.rsquare, goodness_c.rsquare);
    fprintf('rmse - row: %f, col: %f\n', goodness_r.rmse, goodness_c.rmse);


    % plot intensity profile of the selected dot
    figure_ctr = figure_ctr + 1;
    figure(figure_ctr)

    % plot the profile along a row
    subplot(1,3,1)
    hold on
    plot(I_dot_average(round(r_N/2),:),'-*')
    plot(fit_r)
    title([num2str(i) ', row'])

    % plot the profile along a column
    subplot(1,3,2)
    hold on
    plot(I_dot_average(:,round(c_N/2)),'-*')
    plot(fit_c)
    title([num2str(i) ', column'])

    % plot the full profile
    subplot(1,3,3)
    surf(I_dot_average, 'linestyle','none');
    view(2)
    title('full')
% 
% % find average and std of intensity profile
% 
% 
% 

end










