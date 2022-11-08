% Controlling The Absorbing Boundary Layer Example
%
% The first-order simulation codes included within the k-Wave Toolbox
% (kspaceFirstOrder1D, kspaceFirstOrder2D, and kspaceFirstOrder3D) use a
% special type of anisotropic absorbing boundary layer known as a perfectly
% matched layer (PML) to absorb acoustic waves when they reach the edge of
% the computational domain. By default, this layer occupies a strip of 20
% grid points (10 in 3D) around the edge of the domain. Without this
% boundary layer, the computation of the spatial derivates via the FFT
% causes waves leaving one side of the domain to reappear at the opposite
% side. The use of the PML thus facilitates infinite domain simulations
% without the need to increase in the size of the computational grid. This
% example demonstrates how to control the parameters of the PML within
% k-Wave via optional input parameters.
%
% author: Bradley Treeby
% date: 30th June 2009
% last update: 4th June 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby

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

% modify this parameter to run the different examples
example_number = 1;
% 1: PML with no absorption
% 2: PML with the absorption value set too high 
% 3: partially effective PML 
% 4: PML set to be outside the computational domain

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

% set the input arguments
switch example_number
    case 1
        input_args = {'PMLAlpha', 0};
    case 2
        input_args = {'PMLAlpha', 1e6};        
    case 3
        input_args = {'PMLSize', 2};
    case 4
        input_args = {'PMLInside', false};
end

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the simulated sensor data
figure;
imagesc(sensor_data, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;