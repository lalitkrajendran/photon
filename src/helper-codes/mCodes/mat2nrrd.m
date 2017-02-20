% THIS CODE HAS TO BE RUN IN THE SAME COMPUTER WHEN SCHLIERENRAY IS USED,
% so that the raw files can be read consistently

clear all, close all, clc

% rho = h5read('channel_isothermal_n104551_t1.500e+02.h5','/rho');
% x = h5read('channel_isothermal_n104551_t1.500e+02.h5','/x');
% y = h5read('channel_isothermal_n104551_t1.500e+02.h5','/y');
% z = h5read('channel_isothermal_n104551_t1.500e+02.h5','/z');


filepath = '/home/barracuda/a/lrajendr/JHU_Buoyancy_Turb/';
files = dir([filepath '*.mat']);
fprintf('No. of files : %d\n', numel(files));

N = numel(files);
% N = 1;

fileID = fopen([filepath 'batchfile.sh'],'w');
spacing = 4*pi/1024; % For JHU data
for i = 1:N
  
  if(mod(i,10)==0)
      fprintf('file number : %d\n',i)
  end
  
  A = load([filepath files(i).name]);
  rho = A.rho;
  
 
  name = ['jhu_buoy_turb_' num2str(i,'%02d') '.nrrd'];
  nrrdWriter([filepath name],single(rho),[1 1 1],[0,0,0],'raw');
  fprintf(fileID,['./schlieren ' filepath name '\n']);
  fprintf(['./schlieren ' filepath name '\n'])
end

fclose(fileID);


