% Modelling Sensor Frequency Response Example
%
% This example shows how to account for the frequency response of detectors
% when this response has a Gaussian shape (e.g., piezoelectric
% transducers). It builds on the Homogeneous Propagation Medium Example. 
%
% author: Bradley Treeby
% date: 6th September 2010
% last update: 9th June 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2017 Bradley Treeby

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
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;

% define the time array
kgrid.makeTime(medium.sound_speed);

% create initial pressure distribution using makeDisc
disc_magnitude = 5; % [Pa]
disc_x_pos = 50;    % [grid points]
disc_y_pos = 50;    % [grid points]
disc_radius = 8;    % [grid points]
disc_1 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_magnitude = 3; % [Pa]
disc_x_pos = 80;    % [grid points]
disc_y_pos = 60;    % [grid points]
disc_radius = 5;    % [grid points]
disc_2 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

source.p0 = disc_1 + disc_2;

% define a centered circular sensor
sensor_radius = 4e-3;   % [m]
num_sensor_points = 50;
sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% define the frequency response of the sensor elements
center_freq = 3e6;      % [Hz]
bandwidth = 80;         % [%]
sensor.frequency_response = [center_freq, bandwidth];

% re-run the simulation
sensor_data_filtered = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% calculate the frequency spectrum at the first sensor element
[~, sensor_data_as] = spect(sensor_data(1, :), 1/kgrid.dt);
[f, sensor_data_filtered_as] = spect(sensor_data_filtered(1, :), 1/kgrid.dt);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the simulated sensor data at the first sensor element
figure;
[t_sc, scale, prefix] = scaleSI(max(kgrid.t_array(:)));
plot(scale * kgrid.t_array, sensor_data(1, :), 'k-', ...
     scale * kgrid.t_array, sensor_data_filtered(1, :), 'r-');
ylabel('Pressure [au]');
xlabel(['Time [' prefix 's]']);
legend('Original Time Series', 'Filtered Time Series');

% plot the amplitude spectrum
figure;
[f_sc, scale, prefix] = scaleSI(max(f(:)));
plot(scale * f, sensor_data_as, 'k-', ...
     scale * f, sensor_data_filtered_as, 'r-');
ylabel('Amplitude Spectrum');
xlabel(['Frequency [' prefix 'Hz]']);
legend('Original Time Series', 'Filtered Time Series');
set(gca, 'XLim', [0 10]);