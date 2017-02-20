% This code creates the density profile for a normal shock using erf()

clear all
close all
clc

Nx = 256; Ny = 256; Nz = 128;
scale = 0.2;  % default is 0.2
r = 6; % density ratio
rho_min = 1; % density ahead of the shock
rho_max = r*rho_min;
delta_rho = rho_max - rho_min;

%filepath = '/home/barracuda/a/lrajendr/SchlierenRayVis/StereoBOS_pinhole/';
filepath = '~/SchlierenRayVis/Data/';
filename = 'NormShock_6_scale_0_2.nrrd';

rho = ones(Nx,Ny,Nz);

for i = 1:Nx
    for k = 1:Nz
        % step by step
        x = (-Ny/2:Ny/2-1)*scale; % scale along x
        y = 0.5*delta_rho*erf(x); % scale along y
        rho(i,:,k) = rho_min + (y-min(y)); % shift up
        
        % condensed expression
        %rho(i,:,k) = rho_min/2*(erf((-Ny/2:Ny/2-1)*scale)*(r-1) + r+1); %
    end
end
K = 0.226; % gladstone dale constant for air (cm^3/g)
n = 1+K/1000*rho;
h = surf(rho(:,:,1));
set(h,'LineStyle','none')
nrrdWriter([filepath filename],single(n),[1 1 1],[0 0 0],'raw');