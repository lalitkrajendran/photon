clear all
close all
clc

% filepath = '/home/barracuda/a/lrajendr/SchlierenRayVis/StereoBOS_pinhole/Processing/singleCameraResults/NormShock_1_5_5deg/';
filepath = '/home/barracuda/a/lrajendr/SchlierenRayVis/StereoBOS_random/Processing/shock_3/';
filename = 'BOS_pass1_1.mat';
load([filepath filename]);

figure(1)
quiver(U,V,10)

N = 1024;
%U2 = U;
U2 = U;
[m,n] = size(U2);

U2 = [zeros(m,1) U2 zeros(m,1)];

rho = ones(m,n+2);
grid_res = 8;
del_x = X(1,2)-X(1,1);

c = 0.55; % principal distance
R = 1.25; % distance of camera from origin
r = -0.5; % location of target from origin
scale = c/(R-r-c); % percentage of screen occupied by the image
del_x = del_x/N*1/scale;


for i = 1:m
    for j = 1:n
      %if(U2(i,j)>=0)
        rho(i,j+1) = rho(i,j) + U2(i,j)*del_x; 
      %end
      %rho(i,j+1) = rho(i,j) + U2(i,j)*grid_res;
  end
end

rho(:,n+2) = rho(:,n+1);

rho_avg = mean(rho,1);
%rho_avg = sum(rho,1)/ctr_non_zero*127;
%rho_avg = rho(63,:);

%rho_avg_scale = 1+(rho_avg-1)/(max(rho_avg)-1)*2;
rho_avg_scale = rho_avg;
x_avg = (0:grid_res:N)*1/N;
% hold on
figure(2)
plot(x_avg,rho_avg_scale,'r','LineWidth',1.5)
xlabel('x (arbitrary units) \rightarrow')
ylabel('Density (\rho) \rightarrow')
title('Density Ratio = 3')
% % Reference Density
% Ny = N;
% scale = 0.1;
% scale = 0.025;
% 
% %filepath = '/home/barracuda/a/lrajendr/SchlierenRayVis/StereoBOS/';
% %filename = 'NormShock_2.nrrd';
% 
% rho_ref = ones(1,Ny);
% rho_ref = 2+erf((-Ny/2:Ny/2-1)*scale);
% 
% x_ref = (1:N)*1/N;
% 
% %plot(x_avg,rho_avg_scale,'b',x_ref,rho_ref,'r')
% plot(x_avg,rho_avg_scale)
% legend('Reconstructed','Original')
% title('Density Variation across a Normal Shock')
% xlabel('x (pixels)')
% ylabel('\rho')
% 
% 
% % y = 2+erf((-5:5)*0.5);
% % rho_ref = [ones(1,119) y 3*ones(1,126)];
% % x_ref = 0:255;
% % plot(x_avg,rho_avg,'b',x_ref,rho_ref,'r')
% % rho = ones(31,33);
% % plot(x_avg,rho_avg,'b',x_ref,rho_ref,'r')
% % rho = ones(31,33);
% % for i = 1:31
% % for j = 1:31
% % rho(i,j+1) = rho(i,j) + U2(i,j)*8;
% % end
% % end
% % rho_avg = mean(rho,1);
% % rho_avg_scale = min(rho_avg)+(rho_avg-min(rho_avg))/(max(rho_avg)-min(rho_avg))*2;
% % plot(x_avg,rho_avg,'b',x_ref,rho_ref,'r')
% % rho(:,33) = rho(:,32);
% % rho_avg = mean(rho,1);
% % rho_avg_scale = min(rho_avg)+(rho_avg-min(rho_avg))/(max(rho_avg)-min(rho_avg))*2;
% % x_ref = 0:255;
% % 
