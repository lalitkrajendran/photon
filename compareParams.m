clear all
close all
clc

ref_var = load('reference_params.mat');
names = fieldnames(ref_var);

py_var = load('start_particle.mat');

for i = 1:length(names)
    field = names(i);
    field = field{:};
    disp(field);
    
%     disp_string = ['ref_var' ref_var.(field) 'py_var' py_var.(field)];
%     disp(disp_string)
    
    disp('ref_var');
    disp(ref_var.(field));
    disp('py_var');
    disp(py_var.(field));
    
%     fprintf('ref_var : %s, py_var : %s',getfield(ref_var,names(i)),getfield(py_var,names(i)));
end
