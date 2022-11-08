% Focussed Detector In 3D Example
%
% This example shows how k-Wave can be used to model the output of a
% focussed bowl detector in 3D where the directionality arises from
% spatially averaging across the detector surface. It builds on the
% Focussed Detector In 2D Example. 
%
% author: Ben Cox and Bradley Treeby
% date: 29th October 2010
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
Nx = 64;            % number of grid points in the x direction
Ny = 64;            % number of grid points in the y direction
Nz = 64;            % number of grid points in the z direction
dx = 100e-3/Nx;     % grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]
dz = dx;            % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the propagation medium
medium.sound_speed = 1500;	% [m/s]

% create the time array
kgrid.makeTime(medium.sound_speed);

% create a concave sensor
sphere_offset = 10;                             % [grid points]
diameter = Nx/2 + 1;                            % [grid points]
radius = Nx/2;                                  % [grid points]
bowl_pos = [1 + sphere_offset, Ny/2, Nz/2];     % [grid points]
focus_pos = [Nx/2, Ny/2, Nz/2];                 % [grid points]
sensor.mask = makeBowl([Nx, Ny, Nz], bowl_pos, radius, diameter, focus_pos);

% define a time varying sinusoidal source
source_freq = 0.25e6;                           % [Hz]
source_mag = 1;                                 % [Pa]
source.p = source_mag * sin(2 * pi * source_freq * kgrid.t_array);
source.p = filterTimeSeries(kgrid, medium, source.p);

% place the first point source near the focus of the detector
source1 = zeros(Nx, Ny, Nz);
source1(1 + sphere_offset + radius, Ny/2, Nz/2) = 1;

% place the second point source off axis
source2 = zeros(Nx, Ny, Nz);
source2(1 + sphere_offset + radius, Ny/2 + 6, Nz/2 + 6) = 1;

% run the first simulation
source.p_mask = source1;
input_args = {'PMLSize', 10, 'DataCast', 'single', 'PlotSim', false};
sensor_data1 = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% average the data recorded at each grid point to simulate the measured
% signal from a single element focussed detector
sensor_data1 = sum(sensor_data1, 1);

% run the second simulation
source.p_mask = source2;
sensor_data2 = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% average the data recorded at each grid point to simulate the measured
% signal from a single element focussed detector
sensor_data2 = sum(sensor_data2, 1);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the detector and on-axis and off-axis point sources
voxelPlot(sensor.mask + source1 + source2);
view([130, 40]);

% get a suitable scaling factor for the time axis
[~, t_scale, t_prefix] = scaleSI(kgrid.t_array(end));

% plot the time series corresponding to the on-axis and off-axis sources
figure
plot(kgrid.t_array * t_scale, sensor_data1, '-');
hold on
plot(kgrid.t_array * t_scale, sensor_data2, 'r-');
xlabel(['Time [' t_prefix 's]']);
ylabel('Average Pressure Measured By Focussed Detector [Pa]');
legend('Source on axis', 'Source off axis');