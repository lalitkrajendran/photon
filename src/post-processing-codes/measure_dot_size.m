function dot_diameters = measure_dot_size(method, filepath, read_order)

% these are the paths containing the m-files used for particle
% identification and sizing
addpath('/home/barracuda/a/lrajendr/Software/prana/');
addpath('/home/barracuda/a/lrajendr/Projects/camera_simulation/src/post-processing-codes/from-sayantan/');

% these are the names of the file to be read
files = dir([filepath '*.tif']);

% this the number of files
N = length(files);
fprintf('number of files in the directory: %d\n', N);

% this is the skip value 
% 0 if all files are analyzed, 1 if alternate files are analyzed
if strcmp(read_order, 'alternate')
    skip = 1;
else
    skip = 0;
end

% this is the number of images that will be analyzed
num_analysis = length(1:1+skip:N);
fprintf('number of files to be analyzed: %d\n', num_analysis);

%% measure dot diameters
 
% intialize the array that stores the diameter values
dot_diameters = zeros(num_analysis);

% initialize the counter for the arrays
image_ctr = 0;

for i = 1:1+skip:N
    % display progress
    display_calculation_progress(image_ctr, 1:num_analysis);
    
    % this loads the image
    img = imread([img_filepath files(i).name]);

    % this updates the number of images read
    image_ctr = image_ctr + 1;
    
    if strcmp(method, 'imfindcircles')
        %% calculate diameter from imfindcircles

        % detect dots and measure their radii (pixels)
        [centers, radii] = imfindcircles(img, [1 5], 'Sensitivity', 0.95);

        % calculate the representative dot radius (pixels)
        radius = nanmean(radii);

        % calculate the diameter of the dots (pixels)
        dot_diameter(image_ctr) = 2 * radius;

    elseif strcmp(method, 'prana')
        %% calculate diameter from prana

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % identify the particles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % define settings for particle identification
        IDmethod = {'blob','dynamic','combined'};
        Data.ID.method         = IDmethod{2};
        Data.ID.run            = 1;
        Data.ID.v              = 10;
        Data.ID.contrast_ratio = 0;
        % Data.ID.s_num          = 0;
        % Data.ID.s_name         = 'ID_cam2_';
        Data.ID.save_dir       = particle_id_filepath;

        particleIDprops       = Data.ID;

        % call function to id the particles
        [ID1.p_matrix,ID1.peaks,ID1.num_p]=particle_ID(img,particleIDprops);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % size the particles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % define settings for particle sizing
        SIZEmethod={'GEO','IWC','TPG','FPG','CFPG','LSG','CLSG'};
        Data.Size.run      = 1;
        Data.Size.thresh   = 10;
        Data.Size.method   = 'CLSG';
        Data.Size.p_area   = 0;
        Data.Size.sigma    = 4;
        Data.Size.errors   = 1;%str2double(Data.Size.errors);
        % Data.Size.s_name   = PTV_Data.Size.savebase;
        Data.Size.save_dir = particle_size_filepath;

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
        dp_prana = nanmean(Dp1);

        % store diameter in array
        dot_diameters(image_ctr) = dp_prana;

    else
        fprintf('invalid method entered, exiting');
        return;
    end

end

end