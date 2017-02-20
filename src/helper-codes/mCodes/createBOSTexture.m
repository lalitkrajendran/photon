% This code generates a BOS Texture taking image size (pixels) and 

clear all, close all, clc

N = 1024; % Generate image has NxN pixels
thresh = 0.99; % higher threshold generates fewer points on the final image

filepath = '/home/barracuda/a/lrajendr/SchlierenRayVis/StereoBOS_random/';
filename = ['BOS_texture_' num2str(N) '_fewer.bin'];


A = 255*rand([N,N]); % Generate Matrix

A2 = A - thresh*max(A(:)); 
A3 = 0.5*(A2 + abs(A2)); 
A4 = uint8(A3./A3)*255; 

figure, imshow(A4)

% Save image to file
fID = fopen([filepath filename],'wb');
fwrite(fID,uint32(A4'),'uint32');
