clear all
close all
clc

table_filename = 'mie_table_oil_air.dat';
scattering_data=dlmread(table_filename,'',5,0);

% This extracts the independent variable giving the angle for the rays
angle_data=scattering_data(:,1);
% This converts the angular data from degrees to radians
angle_data=pi*angle_data/180;

% This extracts the dependent variable giving the unpolarized scattering
% magnitude
s11_data=scattering_data(:,2);
% This extracts the dependent variable giving the quantity of polarization
% (ie this will be zero for unpolarized scattering . . . I think)
pol_data=scattering_data(:,3);

% This computes the dependent variable giving the differences between the
% perpendicular and parallel polarization magnitudes
s12_data=-s11_data.*pol_data;

% This computes the dependent variable giving the scattering magnitude
% perpendicular to the scattering plane
s1_data=(s11_data-s12_data);
% This computes the dependent variable giving the scattering magnitude
% parallel to the scattering plane
s2_data=(s11_data+s12_data);


%% plot results

var = s11_data;
max_val = max(var);


% P = polar(angle_data,log(max_val)*ones(size(angle_data)));
% set(P, 'Visible', 'off')
% hold on
h = polar(angle_data,log(var));
set(h,'LineWidth',1.0);
title('\bf 10\mu m oil particle in air');
% axis tight
% print('plots/mie_oil_air_10.png','-dpng');

