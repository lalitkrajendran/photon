% this program loads the results of a PRANA analysis and computes the
% average x displacement

clear
close all
clc

% specify the f number
f_number = 64;

% set the filepath where the results are located
filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dof/f' num2str(f_number) '/'];

% set the name of the file containing the results
filename = dir([filepath '/*3_1.mat']);

% load results into workspace
load([filepath filename.name])

% display f number and average x displacement to the user
fprintf('f number: %d\n', f_number);
fprintf('average displacement: %.2f\n', mean(U(:)));