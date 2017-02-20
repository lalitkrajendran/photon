% This program reads in a set of files in H5 or HDF5 format and saves them
% in NRRD format.

clear all, close all, clc

% rho = h5read('channel_isothermal_n104551_t1.500e+02.h5','/rho');
% x = h5read('channel_isothermal_n104551_t1.500e+02.h5','/x');
% y = h5read('channel_isothermal_n104551_t1.500e+02.h5','/y');
% z = h5read('channel_isothermal_n104551_t1.500e+02.h5','/z');

% specify filepath and filenames
filepath = '/home/barracuda/a/lrajendr/SchlierenRayVis/Data/turbulent_channel_M_5/';
files = dir([filepath '*.h5']);
fprintf('No. of files : %d\n', numel(files));

% calculate no. of files
N = numel(files);
N = 1;

% open a shell batch file
fileID = fopen([filepath 'batchfile.sh'],'w');

% Loop through all files, reading them in, converting them to nrrd and
% adding an execution command to the shell batch file.
for i = 1:N
  if(mod(i,10)==0)
      fprintf('file number : %d\n',i)
  end
  
  % Read variables from file
  rho = h5read([filepath files(i).name], '/rho');
  x = h5read([filepath files(i).name], '/x');
  y = h5read([filepath files(i).name], '/y');
  z = h5read([filepath files(i).name], '/z');
  u = h5read([filepath files(i).name], '/u');
  v = h5read([filepath files(i).name], '/v');
  w = h5read([filepath files(i).name], '/w');
    
  % Original Data has the grid stretched in the y-direction. Re-map points
  % to have uniform grid spacing along y
%   rho = remap_points_on_uniform_grid(x,y,z,rho);
  
  % Permute it to alter permutation in nrrdWriter
  rho = permute(rho, [2 1 3]);
 
  % specify file name for the nrrd file as well as the tiff file where the
  % rendered schlieren image will be saved
  temp_name = files(i).name;
  nrrd_name = [temp_name(1:numel(temp_name)-2) 'nrrd'];
  tiff_name = [temp_name(1:numel(temp_name)-2) 'tiff'];
  % save file as nrrd
  x_diff = diff(x);
  y_diff = diff(y);
  z_diff = diff(z);
  
  del_x = x_diff(1);
  del_y = mean(y_diff);
  del_z = z_diff(1);
  origin_x = min(x);
  origin_y = min(y);
  origin_z = min(z);
  
  nrrdWriter([filepath nrrd_name],single(rho),[del_x del_y del_z],[origin_x,origin_y,origin_z],'raw');
  % add execution command for the current file to the batch file
  fprintf(fileID,['./schlieren ' filepath nrrd_name ' ' tiff_name '\n']);
end

fclose(fileID);


