clear
close all
clc

% this is the path where the data is located
read_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/blur/50x50-f8-dotsize300-blur0.05/';

% this is the filepath where the image is located
img_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/blur/50x50-f8-dotsize300-blur0.05/';

% this is the path where the figures will be stored
figure_write_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/blur/50x50-f8-dotsize300-blur0.05/';

% this is the folder where the workspace variables will be saved
workspace_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/blur/50x50-f8-dotsize300-blur0.05/';

% this is the file extension for saving the image (e.g. tiff, png, jpeg)
% NOTE: In addition, a .fig file will also be saved
write_file_format = 'png';

% load PIV results
piv_results = load([read_filepath 'PIV_pass1_1.mat']);

%% calculate theoretical displacements

% these are the coefficients of the refractive index distrbution
K = [1.0003312515937128, 5.44e-13, 0.00e+00];

% this is the distance of the center of the refractive index medium from the camera (microns)
Z_D = 25e3;

% this is the dimension of a pixel on the camera sensor (microns)
l_p = 17;

% this is the Magnification
M = 0.1761;

% this is the depth of the refractive index medium (microns)
W = 50e3;

% this is the displacement (pixels)
% D = Z_D * M / l_p * 2 * K(3) * (x - b) * W;

% this is the minimum x co-ordinate of the refractive index field
x_min = -5e4;

% this is the maximum x co-ordinate of the refractive index field
x_max = 5e4;

% this is the number of grid points along x
x_num = 200;

% this is the minimum x co-ordinate of the refractive index field that is parabolic
x_parabolic_min = -10e3;

% this is the maximum x co-ordinate of the refractive index field that is parabolic
x_parabolic_max = 10e3;

% this is the number of grid points along x in the parbolic region
x_parabolic_num = round(x_num * 1.0 * (x_parabolic_max - x_parabolic_min) / (x_max - x_min));

% this is the total field of view
del_x = x_parabolic_max - x_parabolic_min;

% this is the f number of the camera 
f_number = 8;

% this is the Gladstone-Dale constant for air (m^3/kg)
K_gd = 0.226*1e-3;

% this is the density of air at STP (kg/m^3)
rho_air = 1.225;

% this is the refractive index of air at STP
n_air = 1 + K_gd*rho_air;

% generate the co-ordinate along the x direction
x_1 = linspace(x_min, x_parabolic_min, round((x_num - x_parabolic_num)/2.)+1);
x_1 = x_1(1:end-1);
x_parabolic = linspace(x_parabolic_min, x_parabolic_max, x_parabolic_num + 1);
x_parabolic = x_parabolic(1:end-1);
x_2 = linspace(x_parabolic_max, x_max, round((x_num - x_parabolic_num)/2.) + 1);
x_2 = x_2(1:end-1);

x = [x_1 x_parabolic x_2];

n_1 = n_air * ones(size(x_1));
n_parabolic = K(1) - K(2) * (x_parabolic - K(3)).^2;
n_2 = n_air * ones(size(x_2));

n = [n_1, n_parabolic, n_2];

% generate the displacement field (pixels)
d = Z_D*M/l_p * W * diff(n)./diff(x);
d = [0 d];

% convert image x to physical x
image_x = 17 * (1:1024);
physical_x = (1:1024) * 17/M;
physical_x_plot = physical_x(1:8:end);

%% plot contours and histograms

% plot U velocity contour
figure()
imagesc(piv_results.U)
colorbar
title('\Delta x (pixels)')

% save figure to file
print([figure_write_filepath 'U.' write_file_format], ['-d' write_file_format]);

% plot V velocity contour
figure()
imagesc(piv_results.V)
colorbar
title('\Delta y (pixels)')

% save figure to file
print([figure_write_filepath 'V.' write_file_format], ['-d' write_file_format]);

% plot histograms of all velocities
figure()
subplot(2,1,1), histogram(piv_results.U(:))
title('\Delta x )')
subplot(2,1,2), histogram(piv_results.V(:))
title('\Delta y )')

% save figure to file
print([figure_write_filepath 'histogram.' write_file_format], ['-d' write_file_format]);

figure()
quiver(piv_results.X, piv_results.Y, piv_results.U, piv_results.V)
% save figure to file
print([figure_write_filepath 'quiver.' write_file_format], ['-d' write_file_format]);


%% plot averaged displacement, differential and integral along a line

U_avg = mean(piv_results.U, 1);

figure()
plot(piv_results.X(64,:), U_avg)
xlabel('X (pixels)')
ylabel('\Delta x (pixels)')
title('\Delta x_{avg} (pixels)')

% save figure to file
print([figure_write_filepath 'U_avg.' write_file_format], ['-d' write_file_format]);
savefig([figure_write_filepath 'U_avg.fig']);

V_avg = mean(piv_results.V, 1);

figure()
plot(piv_results.X(64,:), V_avg)
xlabel('X (pixels)')
ylabel('\Delta y (pixels)')
title('\Delta y_{avg} (pixels)')

% save figure to file
print([figure_write_filepath 'V_avg.' write_file_format], ['-d' write_file_format]);
savefig([figure_write_filepath 'V_avg.fig']);

U_diff = diff(U_avg)/(piv_results.X(1,2) - piv_results.X(1,1));

figure()
plot(piv_results.X(64,:), [0 U_diff])
xlabel('X (pixels)')
ylabel('d/dx(\Delta x)')
title('d/dx(\Delta x_{avg})')

% save figure to file
print([figure_write_filepath 'U_avg_diff.' write_file_format], ['-d' write_file_format]);
savefig([figure_write_filepath 'U_avg_diff.fig']);

U_integ = zeros(size(U_avg));
for i = 1:length(U_integ)
    U_integ(i) = trapz(U_avg(1:i));
end

figure()
plot(piv_results.X(64,:), U_integ)
xlabel('X (pixels)')
ylabel('integral(\Delta x)')
title('integral(\Delta x_{avg})')

% save figure to file
print([figure_write_filepath 'U_avg_integ.' write_file_format], ['-d' write_file_format]);
savefig([figure_write_filepath 'U_avg_integ.fig']);

%% plot correlation plane

% load image
img = imread([img_filepath 'bos_pattern_image_2.tif']);

% load data
corr_data = load([read_filepath 'PIV_pass1_corrplanes_1.mat']);

% this is the number of correlation planes that are stored
N = length(corr_data.Xloc);
fprintf('no. of planes: %d\n', N);

% these are the locations where a plot of the correlation plane should be
% saved
X1 = 128;
Y1 = 128;

X2 = 512;
Y2 = 512;

plane_id_1 = find(corr_data.Xloc(find(corr_data.Yloc==Y1))==X1);
plane_id_2 = find(corr_data.Xloc(find(corr_data.Yloc==Y2))==X2);

% load correlation plane
cplane = corr_data.C_planes(:,:,plane_id_1);

% % find maximum and its location
% [maxval, maxloc] = max(cplane(:));
% 
% % find maximum location as R,C
% r = fix(maxloc/size(cplane,1));
% c = rem(maxloc,size(cplane,1));

figure()
h1 = subplot(1,2,1);
imshow(img(X1-32:X1+31, Y1-32:Y1+31))
colormap(h1, gray)
h2 = subplot(1,2,2);
surf(corr_data.C_planes(:, :, plane_id_1), 'linestyle', 'none');
colormap(h2, jet)
colorbar(h2)
view([-180 0])
axis([1 64 1 64 0 maxval])
title('no blurring');

% save figure to file
print([figure_write_filepath 'corr_plane_ref.' write_file_format], ['-d' write_file_format]);
savefig([figure_write_filepath 'corr_plane_ref.fig']);


% plane_id_2 = 8000;
% load correlation plane
cplane = corr_data.C_planes(:,:,plane_id_2);

% find maximum and its location
[maxval, maxloc] = max(cplane(:));

% find maximum location as R,C
r = fix(maxloc/size(cplane,2));
c = rem(maxloc,size(cplane,2));

figure()
h1 = subplot(1,2,1);
imshow(img(X2-32:X2+31, Y2-32:Y2+31))
colormap(h1, gray)
h2 = subplot(1,2,2);
surf(corr_data.C_planes(:, :, plane_id_2), 'linestyle', 'none');
colormap(h2, jet)
colorbar(h2)
view([-180 0])
axis([1 64 1 64 0 maxval])
% colormap([gray; jet])

title('blurred');

% save figure to file
print([figure_write_filepath 'corr_plane_blur.' write_file_format], ['-d' write_file_format]);
savefig([figure_write_filepath 'corr_plane_blur.fig']);

%% save worskpace to file

% % this is the name of the current script
% script_name_full = mfilename('fullpath');
% [pathstr, script_name, ext] = fileparts(script_name_full);
% save([workspace_filepath script_name '.mat']);
