clear 
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set densities before and after the shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is the upstream density (kg/m^3)
rho_1 = 1.225;

% this is the density ratio (max value for perfect gas is 6)
density_ratio = 6;

% set shock thickness (cm)
shock_thickness = 100e-7;

% set gladstone dale constant (m^3/kg)
gladstone_dale = 0.226e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate density and refractive index profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate the spatial coordinates (mm)
x = linspace(-5*shock_thickness,5*shock_thickness, 1e3);

% generate the density profile
rho = rho_1 + (1 + tanh(4*x/shock_thickness)) * rho_1/2 * (density_ratio -1);

% generate the refractive index profile
n = rho*gladstone_dale + 1;

% calcualte 1st derivative of density
rho_x = diff(rho)./diff(x*1e-2);

% add a zero element at the end to make the dimensions consistent with x
rho_x = [rho_x 0];

% calculate 1st derivative of refractive index
n_x = rho_x/gladstone_dale;

% calculate 2nd derivative of density
rho_xx = diff(rho_x)./diff(x*1e-2);

% add a zero element at the end to make the dimensions consistent with x
rho_xx = [rho_xx 0];

% calculate 2nd derivative of refractive index
n_xx = rho_xx/gladstone_dale;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)

% plot density and refractive index 
subplot(3,1,1)
yyaxis left
plot(x, rho)
ylabel('Density (kg/m^3)')

yyaxis right
plot(x, n)
ylabel('Refractive Index')
xlabel('x (cm)')

% plot 1st derivatives
subplot(3,1,2)
yyaxis left
plot(x, rho_x)
ylabel('$d\rho/dx (kg/m^4)$', 'Interpreter', 'LaTex')

yyaxis right
plot(x, n_x)
ylabel('$dn/dx$', 'interpreter', 'latex')
xlabel('x (cm)')

% plot second derivatives
subplot(3,1,3)
yyaxis left
plot(x, rho_xx)
ylabel('$d^2\rho/dx^2 (kg/m^5)$', 'interpreter', 'latex')

yyaxis right
plot(x, n_xx)
ylabel('$d^2n/dx^2$', 'interpreter', 'latex')
xlabel('x (cm)')

