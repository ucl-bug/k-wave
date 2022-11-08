% Focussed Detector Example In 2D
%
% This example shows how k-Wave can be used to model the output of a
% focussed semicircular detector where the directionality arises from
% spatially averaging across the detector surface. It builds on the
% Homogeneous Propagation Medium and Using A Binary Sensor Mask examples.
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
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 180;           % number of grid points in the x (row) direction
Ny = 180;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;    % [m/s]

% define a sensor as part of a circle centred on the grid
sensor_radius = 65; % [grid points]
arc_angle = pi;     % [rad]
sensor.mask = makeCircle(Nx, Ny, Nx/2, Ny/2, sensor_radius, arc_angle);

% define the array of temporal points 
t_end = 11e-6;      % [s]
kgrid.makeTime(medium.sound_speed, [], t_end);

% place a disc-shaped source near the focus of the detector
source.p0 = 2 * makeDisc(Nx, Ny, Nx/2, Ny/2, 4);

% run the first simulation
sensor_data1 = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% place a disc-shaped source away from the focus of the detector
source.p0 = 2 * makeDisc(Nx, Ny, Nx/2, Ny/2 + 20, 4);

% run the second simulation
sensor_data2 = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% average the outputs of the simulations over the points making up the sensor
sensor_output1 = sum(sensor_data1, 1) ./ sum(sensor.mask(:));
sensor_output2 = sum(sensor_data2, 1) ./ sum(sensor.mask(:));

% =========================================================================
% VISUALISATION
% =========================================================================

% get a suitable scaling factor for the time axis
[~, t_scale, t_prefix] = scaleSI(kgrid.t_array(end));

% plot the recorded time signals
figure
plot(kgrid.t_array * t_scale, sensor_output1, 'k', ...
     kgrid.t_array * t_scale, sensor_output2, 'r');
xlabel(['Time [' t_prefix 's]']);
ylabel('Average Pressure Measured Over Detector [au]') ;
legend(['Source on focus, sum(output^2) = ' num2str(round(sum(sensor_output1.^2) * 100) / 100)], ...
      ['Source off focus, sum(output^2) = ' num2str(round(sum(sensor_output2.^2) * 100) / 100)]);