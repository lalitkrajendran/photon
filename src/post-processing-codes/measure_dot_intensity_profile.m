% This program checks the accuracy of the continuous least squares guassian
% (CLSG) fit to the intensity distribution of a dot in an image containing 
% a dot pattern used for BOS experiments. Both the dot identification and
% the CLSG based sizing operations are performed by calling functions that 
% are a part of prana. After the sizing operation, the gaussian function
% used in the CLSG fit is reconstructed and the error between the intensity
% distribution given by this function and the actual intensity of the dot
% are compared to evaluate the error. The bias error obtained by averaging
% the individual errors for each pixel in a dot is divided by the average
% intensity of the dot to report a relative bias error. 

% Author: Lalit Rajendran (lrajendr@purdue.edu)

%% clear workspace and close all windows

clear all
close all
clc

% this initializes the figure counter to zero
figure_ctr = 0;

% this is the name for the case
case_name = '150x150-f16-disp5';

% these are the phyiscal dot diameters
% dot_diameters_physical = [10 20 30 40 50 60 70 80 90 100 150 200 250 300 350 400 450 500];
dot_diameters_physical = [10 20 30 40 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950];
%% set read and write paths

% this is the path to the folder containing the raw images
img_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/cropped-images/'];

% this is the folder where the particle identification results will be stored
particle_id_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/particle-id/'];

% this is the folder where the particle sizing results will be stored
particle_size_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/particle-size/'];

% this is the folder where the figures will be saved
figure_save_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/plots/'];

% this is the folder where the workspace variables will be saved
workspace_filepath = ['~/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/results/'];

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

%% load images and perform dot identification and sizing

% these are the paths containing the m-files used for particle
% identification and sizing
addpath('/home/barracuda/a/lrajendr/Software/prana/');
addpath('/home/barracuda/a/lrajendr/Projects/camera_simulation/src/post-processing-codes/from-sayantan/');

% these are the names of the file to be read
files = dir([img_filepath '*.tif']);

% this is the number of files
N = length(files);
fprintf('number of files in the directory: %d\n', N);

% this is the skip value 
% 0 if all files are analyzed, 1 if alternate files are analyzed
skip = 1;

% this is the number of images that will be analyzed
num_analysis = length(1:1+skip:N);
fprintf('number of files to be analyzed: %d\n', num_analysis);

% initialize arrays to hold the errors
bias_error_abs_avg = zeros(num_analysis, 1);
bias_error_rel_avg = zeros(num_analysis, 1);

unsigned_error_abs_avg = zeros(num_analysis, 1);
unsigned_error_rel_avg = zeros(num_analysis, 1);

random_error_abs_avg = zeros(num_analysis, 1);
random_error_rel_avg = zeros(num_analysis, 1);

rms_error_abs_avg = zeros(num_analysis, 1);
rms_error_rel_avg = zeros(num_analysis, 1);

% initialize the counter for the arrays
image_ctr = 0;

% this loops through every other image in the folder and performs dot
% identification and sizing
for i = 1:1+skip:N
%     i = 30;
    
    % this loads the image
    img = imread([img_filepath files(i).name]);

    % this updates the number of images read
    image_ctr = image_ctr + 1;
    waitbar(image_ctr/num_analysis);
    %% locate dots and measure diameter

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % imfindcircles
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     % find location and radius of dots
%     [centers, radii] = imfindcircles(img, [2 6], 'Sensitivity', 0.95);
% 
%     % % remove NAN elements
%     % centers = centers(~isnan(centers));
%     % radii = radii(~isnan(radii(:)));
% 
%     % calculate average diameter of the dot
%     average_diameter = 2 * nanmean(radii);
% 
%     % calculate the total number of dots that were identified
%     number_of_dots = sum(~isnan(radii(:)));
% 
%     fprintf('Number of dots: %d\n', number_of_dots);
%     fprintf('Average diameter: %f\n', average_diameter);
% 
% %     % check if dots are identified correctly
% %     figure_ctr = figure_ctr + 1;
% %     figure(figure_ctr)
% %     hold on
% %     imshow(img)
% %     viscircles(centers, radii)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify the particles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define settings for particle identification
    IDmethod = {'blob','dynamic','combined'};
    Data.ID.method         = IDmethod{1};
    Data.ID.run            = 1;
    Data.ID.v              = 10;
    Data.ID.contrast_ratio = 0;

    particleIDprops       = Data.ID;

    % call function to id the particles
    [ID1.p_matrix,ID1.peaks,ID1.num_p]=particle_ID(img,particleIDprops);
    Data.ID.save_dir       = particle_id_filepath;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % size the particles [intial guess]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define settings for particle sizing
    Data.Size.run      = 1;
    Data.Size.thresh   = 10;
    Data.Size.method   = 'IWC';
    Data.Size.p_area   = 0;
    Data.Size.sigma    = 4;
    Data.Size.errors   = 1;%str2double(Data.Size.errors);

    sizeprops       = Data.Size;
    % call function to size the particles
    [SIZE1.XYDiameter,SIZE1.mapsizeinfo,SIZE1.locxy]=particle_sizing(img,ID1.p_matrix,...
                        ID1.num_p,sizeprops);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % size the particles 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define settings for particle sizing
    Data.Size.run      = 1;
    Data.Size.thresh   = 10;
    Data.Size.method   = 'CLSG';
    Data.Size.p_area   = mean(SIZE1.XYDiameter(:,3));
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
%     dot_diameters_image(image_ctr) = nanmean(Dp1);

    % calculate average diameter of the dot
    average_diameter = nanmean(Dp1);

    % calculate the total number of dots that were identified
    number_of_dots = sum(~isnan(Dp1(:)));

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

    % intialize percentage rms error to zero
    bias_error_abs = 0.0;
    bias_error_rel = 0.0;
    unsigned_error_abs = 0.0;
    unsigned_error_rel = 0.0;
    random_error_abs = 0.0;
    random_error_rel = 0.0;
    rms_error_abs = 0.0;
    rms_error_rel = 0.0;
    
    for j = 1:number_of_dots
%         j = 15;
        % select the dot whose intensity profile will be plotted
        dot_index = j;

        % extract the location of the dot in the image
%         dot_loc = centers(dot_index,:);
        dot_loc = [X1(j,1) Y1(j,1)];
        
        % convert the dot location to an integer value
        dot_loc = round(dot_loc);

        % the center is stored as x,y but the array values are stored as r,c. so
        % convert the dot location to an r,c format
        temp = dot_loc;
        dot_loc(1) = temp(2);
        dot_loc(2) = temp(1);
  
        % calculate the minimum and maximum indices of the region surrounding the 
        % dot that will be extracted from the image
        r_min = SIZE1.locxy(j,2);
        r_max = SIZE1.locxy(j,2) + SIZE1.mapsizeinfo(j,2);

        c_min = SIZE1.locxy(j,1);
        c_max = SIZE1.locxy(j,1) + SIZE1.mapsizeinfo(j,1);
        
        % if any of the indices are out of bounds, ignore this dot
        if r_min <= 0 || c_min <= 0 || r_max > size(img,1) || c_max > size(img,2)
            continue;
        end

        % extract the intensity profile of the selected dot
        I_dot = double(img(r_min:r_max, c_min:c_max));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate a 2D guassian with the given parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % top left corner 
        r_00 = r_min;
        c_00 = c_min;

        % center
        rc = SIZE1.XYDiameter(j,2);
        cc = SIZE1.XYDiameter(j,1);

        % extent
        del_r = r_max - r_min + 1;
        del_c = c_max - c_min + 1;

        % max intensity
        I0 = SIZE1.XYDiameter(j,4);
        if(isnan(I0))
            continue
        end
        
        % standard deviation
        sigma1 = SIZE1.XYDiameter(j,3)/4;

        % get equation for gaussian
        gaussian_fit = zeros(del_r, del_c); 
        for l = 1:del_r
            for m = 1:del_c
                gaussian_fit(l,m) = I0*exp(-((l + r_min - 1 - rc)^2 + (m + c_min - 1 - cc)^2)/(2*sigma1^2));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate error in the gaussian fit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % calculate the rms error in the fit (%)
%         error_gaussian_fit = error_gaussian_fit + norm(gaussian_fit - I_dot)/norm(I_dot)*100;
%         [bias_error, random_error, rms_error] = compute_errors(gaussian_fit(:), I_dot(:));
        bias_error = mean(gaussian_fit(:) - I_dot(:));
        unsigned_error = mean(abs(gaussian_fit(:) - I_dot(:)));
        
        % NOTE: random error does not make sense here, because it is just
        % going to be the standard deviation of the intensities of the
        % pixels that make up the image of the dot. And since this is also
        % a measure of the dot diameter, it would never be zero in a real
        % case, and hence does not have any meaningful information stored
        % in it. For the same reason, I am also abandoning the rms error,
        % because it includes the effect of the random error. (RMS = bias^2
        % + random^2)
        
        % calculate total absolute bias error
        bias_error_abs = bias_error_abs + bias_error;
        unsigned_error_abs = unsigned_error_abs + unsigned_error;
%         random_error_abs = random_error_abs + random_error;
%         rms_error_abs = rms_error_abs + rms_error;
        
        % calculate total relative bias error
        bias_error_rel = bias_error_rel + bias_error/mean(I_dot(:)) * 100;
        unsigned_error_rel = unsigned_error_rel + unsigned_error/mean(I_dot(:)) * 100;
%         random_error_rel = random_error_rel + random_error/mean(I_dot(:))*100;
%         rms_error_rel = rms_error_rel + rms_error/mean(I_dot(:))*100;
        
        % update the number of dots whose intensity profiles have been
        % extracted
        dot_ctr = dot_ctr + 1;

    end

    fprintf('number of dot intensity profiles: %d\n', dot_ctr);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate average error in the gaussian fit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bias_error_abs_avg(image_ctr) = bias_error_abs/dot_ctr;
    unsigned_error_abs_avg(image_ctr) = unsigned_error_abs/dot_ctr;
%     random_error_abs_avg(image_ctr) = random_error_abs/dot_ctr;
%     rms_error_abs_avg(image_ctr) = rms_error_abs/dot_ctr;

    bias_error_rel_avg(image_ctr) = bias_error_rel/dot_ctr;
    unsigned_error_rel_avg(image_ctr) = unsigned_error_rel/dot_ctr;

%     random_error_rel_avg(image_ctr) = random_error_rel/dot_ctr;
%     rms_error_rel_avg(image_ctr) = rms_error_rel/dot_ctr;

    fprintf('image: %d\n', image_ctr);
    fprintf('average signed bias error of gaussian fit: %f\n', bias_error_abs_avg(image_ctr));
    fprintf('average unsigned error of gaussian fit: %f\n', unsigned_error_rel_avg(image_ctr));

%     % plot intensity profile of the selected dot
%     figure_ctr = figure_ctr + 1;
%     figure(figure_ctr)
% 
%     % plot the profile along a row
%     subplot(1,3,1)
%     hold on
%     plot(I_dot_average(round(r_N/2),:),'-*')
%     plot(fit_r)
%     title([num2str(i) ', row'])
% 
%     % plot the profile along a column
%     subplot(1,3,2)
%     hold on
%     plot(I_dot_average(:,round(c_N/2)),'-*')
%     plot(fit_c)
%     title([num2str(i) ', column'])
% 
%     % plot the full profile
%     subplot(1,3,3)
%     surf(I_dot_average, 'linestyle','none');
%     view(2)
%     title('full')
% 
% % find average and std of intensity profile
% 
% 
% 

end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the ABSOLUTE errors against the physical dot diameter
figure(1)

% bias error
subplot(1,2,1)
hold on
plot(dot_diameters_physical, bias_error_abs_avg, '*-');
xlabel('dot diameter (microns)');
ylabel('error');
title('bias error (absolute)');
grid on

% % random error
% subplot(1,3,2)
% hold on
% plot(dot_diameters_physical, random_error_abs_avg, '*-');
% xlabel('dot diameter (microns)');
% ylabel('error');
% title('random error (absolute)');
% grid on

% rms error
subplot(1,2,2)
hold on
plot(dot_diameters_physical, unsigned_error_abs_avg, '*-');
xlabel('dot diameter (microns)');
ylabel('error');
title('E[|x-x\hat|] error (absolute)');
grid on


% save figure to file
savefig(gcf, [figure_save_filepath 'abs-error-gaussian-fit.fig']);
print(gcf, [figure_save_filepath 'abs-error-gaussian-fit.eps'], '-depsc');
print(gcf, [figure_save_filepath 'abs-error-gaussian-fit.png'], '-dpng');


% plot the relative errors against the physical dot diameter
figure(2)

% bias error
subplot(1,2,1)
hold on
plot(dot_diameters_physical, bias_error_rel_avg, '*-');
xlabel('dot diameter (microns)');
ylabel('% error');
title('bias error (relative)');
grid on

% % random error
% subplot(1,3,2)
% hold on
% plot(dot_diameters_physical, random_error_rel_avg, '*-');
% xlabel('dot diameter (microns)');
% ylabel('% error');
% title('random error (relative)');
% grid on

% rms error
subplot(1,2,2)
hold on
plot(dot_diameters_physical, unsigned_error_rel_avg, '*-');
xlabel('dot diameter (microns)');
ylabel('% error');
title('E[|x-x\hat|] error (relative)');
grid on

% save figure to file
savefig(gcf, [figure_save_filepath 'rel-error-gaussian-fit.fig']);
print(gcf, [figure_save_filepath 'rel-error-gaussian-fit.eps'], '-depsc');
print(gcf, [figure_save_filepath 'rel-error-gaussian-fit.png'], '-dpng');

%% save worskpace to file

% this is the name of the current script
script_name_full = mfilename('fullpath');
[pathstr, script_name, ext] = fileparts(script_name_full);
save([workspace_filepath script_name '.mat']);








