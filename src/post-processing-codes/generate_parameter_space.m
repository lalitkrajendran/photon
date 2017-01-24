% generate a plot of the parameter space for a BOS experiment, when the
% real lens effects like depth of field are taken into account

clear 
close all
clc

% set the pixel pitch (mm)
pixel_pitch = 17e-3;

% set the circle of confusion on the image plane (in pixels)
c_i = 3;

% set the range of f-number
f_numbers = [1.4, 2, 2.8, 4, 5.6, 8, 11, 16, 22];

% set the magnification range
magnification = linspace(0.05,0.5,11);

