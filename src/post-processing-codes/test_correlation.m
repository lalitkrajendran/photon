% this program loads a set of correlation planes and plots them as a time
% sequence

clear
close all
clc

figure_ctr = 0;

% this is the path where the correlation plane data is located
read_filepath = '/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/150x150-f16/test-processing/';

% load data
data = load([read_filepath 'PIV_pass1_corrplanes_01.mat']);

% this is the number of correlation planes that are stored
N = length(data.Xloc);
fprintf('no. of planes: %d\n', N);

% % plot mesh video
% figure_ctr = figure_ctr + 1;
% figure(figure_ctr)
% hold on
% for i = 1:100
%     % load correlation plane
%     cplane = data.C_planes(:,:,i);
%     
%     % find maximum and its location
%     [maxval, maxloc] = max(cplane(:));
%     
%     % find maximum location as R,C
%     r = fix(maxloc/size(cplane,1));
%     c = rem(maxloc,size(cplane,1));
%     
%     mesh(data.C_planes(:, :, i));
%     view(3)
%     colorbar
%     
%     str = sprintf('i = %d, r,c = %d, %d', i, r,c);
%     title(str);
%     pause(0.1);
% end


% plot surf video
figure_ctr = figure_ctr + 1;
figure(figure_ctr)
hold on

r0 = 32;
c0 = 32;
for i = 1:100:N
    % load correlation plane
    cplane = data.C_planes(:,:,i);
    
    % find maximum and its location
    [maxval, maxloc] = max(cplane(:));
    
    % find maximum location as R,C
    r = fix(maxloc/size(cplane,1));
    c = rem(maxloc,size(cplane,1));
    
    surf(data.C_planes(:, :, i), 'linestyle', 'none');
    view(3)
    colorbar
    axis([1 64 1 64])
    
    
    str = sprintf('i = %d, del-x,del-y = %d, %d\nx,y = %d,%d', i,r-r0,c-c0,data.Xloc(i),data.Yloc(i));
    title(str);
    pause(0.1);
end

% 
% 
