% Plot variation of density ratio with Mach number for a normal shock

clear all
close all
clc

% Define specific heat ratio
gamm = 1.4;

% Define Mach number range
M = linspace(1,50,100);

% Calculate density ratio
numerator = (gamm+1)*M.^2;
denominator = (gamm-1)*M.^2 + 2;
densityRatio = numerator./denominator;

plot(M,densityRatio)
xlabel('M_1 \rightarrow')
ylabel('\rho_2/\rho_1 \rightarrow')
title('Variation of Density ratio with Mach Number for a Normal Shock')
xlim([min(M) max(M)])

