% Sensor Element Directivity In 2D Example
%
% This example shows how to attribute a directional response to a
% single-element sensor, or individual elements of a multi-element sensor
% array. Directionality can be included without a separate function through
% explicit averaging, as shown in the examples Modelling Sensor Directivity
% In 2D and Focussed Detector In 2D, but the functionality described here
% allows greater flexibility. Note that directivity defined in this way is
% currently only supported in 2D. This example builds on the Homogeneous
% Propagation Medium and Using A Binary Sensor Mask examples.
%
% author: Ben Cox and Bradley Treeby
% date: 20th January 2010
% last update: 4th June 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2017 Ben Cox and Bradley Treeby

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
% DEFINE THE GRID AND MEDIUM PROPERTIES
% =========================================================================

% create the computational grid
Nx = 64;        % number of grid points in the x (row) direction
Ny = 64;        % number of grid points in the y (column) direction
dx = 1e-3/Nx;   % grid point spacing in the x direction [m]
dy = dx;     	% grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]

% define the array of temporal points 
t_end = 600e-9;      % [s]
kgrid.makeTime(medium.sound_speed, [], t_end);

% =========================================================================
% DEFINE THE DIRECTIONAL SENSOR ARRAY
% =========================================================================

% define a line of sensor points
sensor.mask = zeros(Nx, Ny);
sensor.mask(24, 2:2:63) = 1;

% define the angle of max directivity for each sensor point:
%    0             = max sensitivity in x direction (up/down)
%    pi/2 or -pi/2 = max sensitivity in y direction (left/right)
dir_angles = (-1:1/15:1).' * pi/2;

% assign to the directivity mask
sensor.directivity_angle = zeros(Nx, Ny);
sensor.directivity_angle(sensor.mask == 1) = dir_angles;

% define the directivity pattern
sensor.directivity_pattern = 'pressure';

% define the directivity size
sensor.directivity_size = 16 * kgrid.dx;

% =========================================================================
% SIMULATION AND VISUALISATION FOR AN INITIAL VALUE PROBLEM
% =========================================================================

% define the initial pressure distribution
source.p0 = zeros(Nx, Ny);
source.p0(39:41, :) = 2;
 
% turn off the PML in the y-direction
input_args = {'PMLAlpha', [2, 0]};

% run the simulation
sensor_data1 = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% extract the maximum pressure after the simulation has reach steady state
max_data = max(sensor_data1, [], 2);

% plot the largest value of the output for each sensor
figure;
plot(dir_angles, max_data ./ max(max_data), 'o');
xlabel('Sensor Directivity Angle (radians)');
ylabel('Maximum Recorded Pressure');
axis([-pi/2, pi/2, 0, 1.05]);

% =========================================================================
% SIMULATION AND VISUALISATION FOR TIME-VARYING SOURCE PROBLEM
% =========================================================================

% define a time varying sinusoidal source (instead of an initial pressure)
source_freq = 12e6;     % [Hz]
source_mag = 0.25;      % [Pa]
source.p = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

% define source mask and force to be binary
source.p_mask = source.p0;
source.p_mask(source.p_mask ~= 0) = 1; 

% remove initial pressure field
source = rmfield(source, 'p0');

% run the simulation
sensor_data2 = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% extract the maximum pressure after the simulation has reach steady state
max_data = max(sensor_data2(:, 150:end), [], 2);

% plot the largest value of the output for each sensor
figure;
plot(dir_angles, max_data ./ max(max_data), 'o');
xlabel('Sensor Directivity Angle (radians)');
ylabel('Maximum Recorded Pressure');
axis([-pi/2, pi/2, 0, 1.05]);

% plot the same data using a polar plot
figure;
polar(dir_angles, max_data ./ max(max_data));