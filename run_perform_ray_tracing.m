% this script reads in the parameters for rendering the particle
% and calibration field and passes them to perform_ray_tracing_03.m

% clear all
close all
clc

%% particle field

for frame_index = 1:1
    particle_filename = ['particle_' num2str(frame_index,'%02d') '.mat'];
    load(particle_filename)
    % convert cell arrays to struct 
    a = struct(optical_system.design.optical_element{:});
    optical_system.design.optical_element = a;

    b = struct(optical_system.design.optical_element(1).optical_element{:});
    optical_system.design.optical_element(1).optical_element = b;
    % reconvert NAN to Null to stay consistent with the code

    optical_system.design.optical_element(1).optical_element(1).element_number=[ ];
    optical_system.design.optical_element(1).optical_element(1).elements_coplanar=[];
    optical_system.design.optical_element(1).optical_element(1).element_properties.abbe_number=[];
    optical_system.design.optical_element(1).element_geometry = [];
    optical_system.design.optical_element(1).element_properties = [];
    scattering_data = [];
    field_type = 'particle';
    I = perform_ray_tracing_03(piv_simulation_parameters,optical_system,pixel_gain,scattering_data,scattering_type,lightfield_source,field_type);
end

% %% calibration
% 
% for plane_index = 1:1
%     calibration_filename = ['calibration_' num2str(plane_index,'%02d') '.mat'];
%     load(calibration_filename)
%     % convert cell arrays to struct 
%     a = struct(optical_system.design.optical_element{:});
%     optical_system.design.optical_element = a;
% 
%     b = struct(optical_system.design.optical_element(1).optical_element{:});
%     optical_system.design.optical_element(1).optical_element = b;
%     % reconvert NAN to Null to stay consistent with the code
% 
%     optical_system.design.optical_element(1).optical_element(1).element_number=[ ];
%     optical_system.design.optical_element(1).optical_element(1).elements_coplanar=[];
%     optical_system.design.optical_element(1).optical_element(1).element_properties.abbe_number=[];
%     optical_system.design.optical_element(1).element_geometry = [];
%     optical_system.design.optical_element(1).element_properties = [];
%     scattering_data = [];
%     field_type = 'calibration';
%     perform_ray_tracing_03(piv_simulation_parameters,optical_system,pixel_gain,scattering_data,scattering_type,lightfield_source,field_type);
% end
%     