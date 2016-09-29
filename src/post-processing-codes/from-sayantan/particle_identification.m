% Particle identification
clear;
camnum=[2 4];
%% Read particle images
imdir='/home/shannon/a/bhattac3/Tomographic_Reconstruction/Images/Uniform_flow/sheet_4mm_ang_35/camera_images/particle_images/';
im1=imread(fullfile(imdir,'camera_02','particle_image_frame_0001.tif'));
im2=imread(fullfile(imdir,'camera_04','particle_image_frame_0001.tif'));


%% Run particle id function
% --- ID ---
IDmethod = {'blob','dynamic','combined'};
Data.ID.method         = IDmethod{2};
Data.ID.run            = 1;
Data.ID.v              = 10;
Data.ID.contrast_ratio = 0;
% Data.ID.s_num          = 0;
% Data.ID.s_name         = 'ID_cam2_';
Data.ID.save_dir       = '/home/shannon/a/bhattac3/Tomographic_Reconstruction/Processing/Particle_Identification/Try2/';

particleIDprops       = Data.ID;

[ID1.p_matrix,ID1.peaks,ID1.num_p]=particle_ID(im1,particleIDprops);
[ID2.p_matrix,ID2.peaks,ID2.num_p]=particle_ID(im2,particleIDprops);

%% Run Particle Sizing function

% --- Sizing ---
SIZEmethod={'GEO','IWC','TPG','FPG','CFPG','LSG','CLSG'};
Data.Size.run      = 1;
Data.Size.thresh   = 10;
Data.Size.method   = 'CLSG';
Data.Size.p_area   = 0;
Data.Size.sigma    = 4;
Data.Size.errors   = 1;%str2double(Data.Size.errors);
% Data.Size.s_name   = PTV_Data.Size.savebase;
% Data.Size.save_dir = PTV_Data.Size.save_dir;

sizeprops       = Data.Size;

[SIZE1.XYDiameter,SIZE1.mapsizeinfo,SIZE1.locxy]=particle_sizing(im1,ID1.p_matrix,...
                    ID1.num_p,sizeprops);
[SIZE2.XYDiameter,SIZE2.mapsizeinfo,SIZE2.locxy]=particle_sizing(im2,ID2.p_matrix,...
                    ID2.num_p,sizeprops);     
                
%% Plot identified particle on image                
                
X1=SIZE1.XYDiameter(:,1);
Y1=SIZE1.XYDiameter(:,2);                
Dp1=SIZE1.XYDiameter(:,3);
Ip1=SIZE1.XYDiameter(:,4); 

X2=SIZE2.XYDiameter(:,1);
Y2=SIZE2.XYDiameter(:,2);                
Dp2=SIZE2.XYDiameter(:,3);
Ip2=SIZE2.XYDiameter(:,4); 

figure(1);hold on;
imagesc(255-im1);colormap(flipud('gray'));axis image;
plot(X1,Y1,'bo');

figure(2);hold on;
imagesc(255-im2);colormap(flipud('gray'));axis image;
plot(X2,Y2,'bo');

% save('particle_id_im2.mat','X1','Y1','Dp1','Ip1','im1');
% save('particle_id_im4.mat','X2','Y2','Dp2','Ip2','im2');


                