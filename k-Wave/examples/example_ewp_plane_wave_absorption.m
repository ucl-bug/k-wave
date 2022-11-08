% Plane Wave Absorption Example
%
% This example illustrates the characteristics of the Kelvin-Voigt
% absorption model used in the k-Wave simulation functions pstdElastic2D,
% pstdElastic3D. It builds on the Explosive Source In A Layered Medium
% Example.
%
% author: Bradley Treeby
% date: 17th January 2014
% last update: 25th July 2019
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2019 Bradley Treeby

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

clearvars;

% =========================================================================
% SET GRID PARAMETERS
% =========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 32;            % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium    
medium.sound_speed_compression = 1800;  % [m/s]
medium.sound_speed_shear       = 1200;  % [m/s]
medium.density                 = 1000;  % [kg/m^3] 

% set the absorption properties
medium.alpha_coeff_compression = 1;     % [dB/(MHz^2 cm)]
medium.alpha_coeff_shear       = 1;     % [dB/(MHz^2 cm)]

% define binary sensor mask with two sensor positions
sensor.mask = zeros(Nx, Ny);
pos1 = 45;                              % [grid points]
pos2 = 65;                              % [grid points]
sensor.mask(pos1, Ny/2) = 1;
sensor.mask(pos2, Ny/2) = 1;

% calculate the distance between the sensor positions
d_cm = (pos2 - pos1) * dx * 100;        % [cm]

% set sensor to record to particle velocity
sensor.record = {'u'};

% define source mask
source_mask = ones(Nx, Ny);
source_pos = 35;                        % [grid points]

% set the CFL
cfl = 0.05;

% define the properties of the PML to allow plane wave propagation
pml_alpha = 0;
pml_size  = [30, 2];

% set the input arguments
input_args = {'PlotScale', 'auto', 'PMLSize', pml_size,...
    'PMLAlpha', pml_alpha, 'PlotPML', false};

% =========================================================================
% COMPRESSIONAL PLANE WAVE SIMULATION
% =========================================================================
        
% define source
source.u_mask = source_mask;
source.ux = zeros(Nx, Ny);
source.ux(source_pos, :) = 1;
source.ux = smooth(source.ux, true);
source.ux = 1e-6.*reshape(source.ux, [], 1);
        
% set end time
t_end = 3.5e-6;
        
% create the time array
kgrid.makeTime(max(medium.sound_speed_compression, medium.sound_speed_shear), cfl, t_end);

% run the simulation
sensor_data_comp = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});

% calculate the amplitude spectrum at the two sensor positions
[~, as1] = spect(sensor_data_comp.ux(1, :), 1/kgrid.dt);
[f_comp, as2] = spect(sensor_data_comp.ux(2, :), 1/kgrid.dt);

% calculate the attenuation from the amplitude spectrums
attenuation_comp = -20 * log10(as2 ./ as1) ./ d_cm;

% calculate the corresponding theoretical attenuation in dB/cm
attenuation_th_comp = medium.alpha_coeff_compression .* (f_comp .* 1e-6).^2;

% calculate the maximum supported frequency
f_max_comp = medium.sound_speed_compression / (2 * dx);

% find the maximum frequency in the frequency vector
[~, f_max_comp_index] = findClosest(f_comp, f_max_comp);

% =========================================================================
% SHEAR PLANE WAVE SIMULATION
% =========================================================================        
                
% define source
clear source
source.u_mask = source_mask;
source.uy = zeros(Nx, Ny);
source.uy(source_pos, :) = 1;
source.uy = smooth(source.uy, true);
source.uy = 1e-6 .* reshape(source.uy, [], 1);
                
% set end time
t_end = 4e-6;
        
% create the time array
kgrid.makeTime(max(medium.sound_speed_compression, medium.sound_speed_shear), cfl, t_end);

% run the simulation
sensor_data_shear = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});

% calculate the amplitude at the two sensor positions
[~, as1] = spect(sensor_data_shear.uy(1, :), 1/kgrid.dt);
[f_shear, as2] = spect(sensor_data_shear.uy(2, :), 1/kgrid.dt);

% calculate the attenuation from the amplitude spectrums
attenuation_shear = -20 * log10(as2 ./ as1) ./ d_cm;

% calculate the corresponding theoretical attenuation in dB/cm
attenuation_th_shear = medium.alpha_coeff_shear .* (f_shear .* 1e-6).^2;

% calculate the maximum supported frequency
f_max_shear = medium.sound_speed_shear / (2 * dx);

% find the maximum frequency in the frequency vector
[~, f_max_shear_index] = findClosest(f_shear, f_max_shear);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot compressional wave traces
figure;
subplot(4, 1, 1);
t_axis = (0:length(sensor_data_comp.ux) - 1) * kgrid.dt * 1e6;
plot(t_axis, sensor_data_comp.ux.', 'k-');
axis tight;
xlabel('Time [\mus]');
ylabel('Particle Velocity');
title('Compressional Wave');

% plot compressional wave absorption
subplot(4, 1, 2);
plot(f_comp .* 1e-6, attenuation_comp, 'ko', ...
     f_comp .* 1e-6, attenuation_th_comp, 'k-');
set(gca, 'XLim', [0, f_max_comp * 1e-6], ...
         'YLim', [0, attenuation_th_comp(f_max_comp_index) * 1.1]);
box on;
xlabel('Frequency [MHz]');
ylabel('\alpha [dB/cm]');

% plot shear wave traces
subplot(4, 1, 3);
t_axis = (0:length(sensor_data_shear.uy)-1) * kgrid.dt * 1e6;
plot(t_axis, sensor_data_shear.uy.', 'k-');
axis tight;
xlabel('Time [\mus]');
ylabel('Particle Velocity');
title('Shear Wave');

% plot shear wave absorption
subplot(4, 1, 4);
plot(f_shear .* 1e-6, attenuation_shear, 'ko', ...
     f_shear .* 1e-6, attenuation_th_shear, 'k-');
set(gca, 'XLim', [0, f_max_shear * 1e-6], ...
         'YLim', [0, attenuation_th_shear(f_max_shear_index) * 1.1]);
box on;
xlabel('Frequency [MHz]');
ylabel('\alpha [dB/cm]');

scaleFig(1, 1.5);