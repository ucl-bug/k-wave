% Comparison Of Modelling Functions Example
%
% This example provides a short comparison between the simulation functions
% kspaceFirstOrder2D and kspaceSecondOrder. It builds on the Homogeneous
% Propagation Medium and Using A Binary Sensor Mask examples. 
%
% author: Bradley Treeby
% date: 27th October 2010
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

example_number = 1;
% 1: non-absorbing medium, no absorbing boundary layer
% 2: non-absorbing medium, using PML and ExpandGrid
% 3: absorbing medium, no absorbing boundary layer
% 4: absorbing medium, using PML and ExpandGrid

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
if example_number > 2
    medium.alpha_power = 1.5;   % [dB/(MHz^y cm)]
    medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
end

% create the time array
t_end = 6e-6;
kgrid.makeTime(medium.sound_speed, 0.3, t_end);

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

% define a centered circular sensor pushed right to the edge of the grid
sensor_radius = 6.3e-3;   % [m]
num_sensor_points = 50;
sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);

% convert the cartesian sensor mask to a binary sensor mask
sensor.mask = cart2grid(kgrid, sensor.mask);

if example_number == 1  || example_number == 3      % no absorbing boundary layer

    % run the simulation using the first order code
    sensor_data_first_order = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLAlpha', 0);

    % run the simulation using the second order code
    sensor_data_second_order = kspaceSecondOrder(kgrid, medium, source, sensor, 'ExpandGrid', false);
    
elseif example_number == 2  || example_number == 4  % using PML and ExpandGrid
    
    % run the simulation using the first order code
    sensor_data_first_order = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLInside', false);

    % run the simulation using the second order code
    sensor_data_second_order = kspaceSecondOrder(kgrid, medium, source, sensor, 'ExpandGrid', true);
    
end

% =========================================================================
% VISUALISATION
% =========================================================================

% plot a single time series
figure;
[t_sc, t_scale, t_prefix] = scaleSI(kgrid.t_array(end));
subplot(2, 1, 1);
plot(kgrid.t_array*t_scale, sensor_data_second_order(1, :), 'k-', ...
     kgrid.t_array*t_scale, sensor_data_first_order(1, :), 'bx');
set(gca, 'YLim', [-1, 1.5]);
xlabel(['Time [' t_prefix 's]']);
ylabel('Pressure [au]');
legend('kspaceSecondOrder', 'kspaceFirstOrder2D', 'Location', 'NorthWest');
title('Recorded Signals');

subplot(2, 1, 2);
plot(kgrid.t_array*t_scale, sensor_data_second_order(1, :) - sensor_data_first_order(1, :), 'k-');
xlabel(['Time [' t_prefix 's]']);
ylabel('Pressure [au]');
title('Difference In Recorded Signals');