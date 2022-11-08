% Explosive Source In A Layered Medium Example
%
% This example provides a simple demonstration of using k-Wave for the
% simulation and detection of compressional and shear waves in elastic and
% viscoelastic media within a two-dimensional heterogeneous medium. It
% builds on the Homogenous Propagation Medium and Heterogeneous Propagation
% Medium examples.
%
% author: Bradley Treeby
% date: 11th February 2014
% last update: 29th May 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2017 Bradley Treeby

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

% define the properties of the upper layer of the propagation medium
medium.sound_speed_compression = 1500 * ones(Nx, Ny);   % [m/s]
medium.sound_speed_shear       = zeros(Nx, Ny);         % [m/s]
medium.density                 = 1000 * ones(Nx, Ny);   % [kg/m^3]

% define the properties of the lower layer of the propagation medium
medium.sound_speed_compression(Nx/2:end, :) = 2000;     % [m/s]
medium.sound_speed_shear(Nx/2:end, :)       = 800;      % [m/s]
medium.density(Nx/2:end, :)                 = 1200;     % [kg/m^3]

% define the absorption properties
medium.alpha_coeff_compression = 0.1; % [dB/(MHz^2 cm)]
medium.alpha_coeff_shear       = 0.5; % [dB/(MHz^2 cm)]

% create the time array
cfl   = 0.1;    % Courant-Friedrichs-Lewy number
t_end = 8e-6;   % [s]
kgrid.makeTime(max(medium.sound_speed_compression(:)), cfl, t_end);

% create initial pressure distribution using makeDisc
disc_magnitude = 5; % [Pa]
disc_x_pos = 30;    % [grid points]
disc_y_pos = 64;    % [grid points]
disc_radius = 5;    % [grid points]
source.p0 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

% define a centered circular sensor
sensor.mask = makeCircle(Nx, Ny, Nx/2, Ny/2, 20);

% define a custom display mask showing the position of the interface from
% the fluid side
display_mask = false(Nx, Ny);
display_mask(Nx/2 - 1, :) = 1;

% define input arguments
input_args = {'PlotScale', [-0.75, 0.75, -0.15, 0.15], 'PlotPML', false,...
    'DisplayMask', display_mask, 'DataCast', 'single'};

% run the simulation
sensor_data = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});

% reorder the simulation data
sensor_data_reordered = reorderSensorData(kgrid, sensor, sensor_data);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot source, sensor, and position of the interface
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, ...
    double(source.p0 | sensor.mask | display_mask), [-1 1]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;

% plot the re-ordered sensor data
figure;
imagesc(sensor_data_reordered, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;