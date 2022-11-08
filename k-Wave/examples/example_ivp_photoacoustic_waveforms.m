% Photoacoustic Waveforms in 1D, 2D and 3D Example
%
% The time-varying pressure signals recorded from a photoacoustic source
% look different depending on the number of dimensions used in the
% simulation. This difference occurs because a point source in 1D
% corresponds to a plane wave in 3D, and a point source in 2D corresponds
% to an infinite line source in 3D. This examples shows the difference
% between the signals recorded in each dimension. It builds on the
% Simulations in One Dimension, Homogeneous Propagation Medium, and
% Simulations in Three Dimensions examples.
%
% author: Bradley Treeby and Ben Cox
% date: 29th January 2011
% last update: 20th June 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011-2017 Bradley Treeby and Ben Cox

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
% SETTINGS
% =========================================================================

% size of the computational grid
Nx = 64;    % number of grid points in the x (row) direction
x = 1e-3;   % size of the domain in the x direction [m]
dx = x/Nx;  % grid point spacing in the x direction [m]

% define the properties of the propagation medium
medium.sound_speed = 1500;      % [m/s]

% size of the initial pressure distribution
source_radius = 2;              % [grid points]

% distance between the centre of the source and the sensor
source_sensor_distance = 10;    % [grid points]

% time array
dt = 2e-9;                      % [s]
t_end = 300e-9;                 % [s]

% computation settings
input_args = {'DataCast', 'single'};

% =========================================================================
% ONE DIMENSIONAL SIMULATION
% =========================================================================

% create the computational grid
kgrid = kWaveGrid(Nx, dx);

% create the time array
kgrid.setTime(round(t_end / dt) + 1, dt);

% create initial pressure distribution
source.p0 = zeros(Nx, 1);
source.p0(Nx/2 - source_radius:Nx/2 + source_radius) = 1;

% define a single sensor point
sensor.mask = zeros(Nx, 1);
sensor.mask(Nx/2 + source_sensor_distance) = 1;

% run the simulation
sensor_data_1D = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% TWO DIMENSIONAL SIMULATION
% =========================================================================

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Nx, dx);

% create the time array
kgrid.setTime(round(t_end / dt) + 1, dt);

% create initial pressure distribution
source.p0 = makeDisc(Nx, Nx, Nx/2, Nx/2, source_radius);

% define a single sensor point
sensor.mask = zeros(Nx, Nx);
sensor.mask(Nx/2 - source_sensor_distance, Nx/2) = 1;

% run the simulation
sensor_data_2D = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% THREE DIMENSIONAL SIMULATION
% =========================================================================

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Nx, dx, Nx, dx);

% create the time array
kgrid.setTime(round(t_end / dt) + 1, dt);

% create initial pressure distribution
source.p0 = makeBall(Nx, Nx, Nx, Nx/2, Nx/2, Nx/2, source_radius);

% define a single sensor point
sensor.mask = zeros(Nx, Nx, Nx);
sensor.mask(Nx/2 - source_sensor_distance, Nx/2, Nx/2) = 1;

% run the simulation
sensor_data_3D = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the time signals recorded in each dimension
figure;
[t_sc, t_scale, t_prefix] = scaleSI(t_end);
plot(kgrid.t_array * t_scale, sensor_data_1D ./ max(abs(sensor_data_1D)), 'b-');
hold on;
plot(kgrid.t_array * t_scale, sensor_data_2D ./ max(abs(sensor_data_2D)), 'r-');
plot(kgrid.t_array * t_scale, sensor_data_3D ./ max(abs(sensor_data_3D)), 'k-');
xlabel(['Time [' t_prefix 's]']);
ylabel('Recorded Pressure [au]');
legend('1D', '2D', '3D');
axis tight;