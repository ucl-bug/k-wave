% Modelling Sensor Directivity In 2D Example
%
% This example demonstrates how the sensitivity of a large single element
% detector varies with the angular position of a point-like source. It
% builds on Monopole Point Source In A Homogeneous Propagation Medium and
% Focussed Detector In 2D examples.
%
% author: Ben Cox and Bradley Treeby
% date: 28th October 2010
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
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 50e-3/Nx;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]

% create the time array
Nt = 350;
dt = 7e-8;      % [s]
kgrid.setTime(Nt, dt)

% define a large area detector
sz = 20;        % [grid points]
sensor.mask = zeros(Nx, Ny);
sensor.mask(Nx/2 + 1, (Ny/2 - sz/2 + 1):(Ny/2 + sz/2 + 1)) = 1;

% define equally spaced point sources lying on a circle centred at the
% centre of the detector face
radius = 30;    % [grid points]
points = 11;
circle = makeCartCircle(radius * dx, points, [0, 0], pi);

% find the binary sensor mask most closely corresponding to the Cartesian
% coordinates from makeCartCircle
circle = cart2grid(kgrid, circle);

% find the indices of the sources in the binary source mask
source_positions = find(circle == 1);

% define a time varying sinusoidal source
source_freq = 0.25e6;   % [Hz]
source_mag = 1;         % [Pa]
source.p = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

% filter the source to remove high frequencies not supported by the grid
source.p = filterTimeSeries(kgrid, medium, source.p);

% pre-allocate array for storing the output time series
single_element_data = zeros(Nt, points);

% run a simulation for each of these sources to see the effect that the
% angle from the detector has on the measured signal
for source_loop = 1:points
    
    % select a point source
    source.p_mask = zeros(Nx, Ny);
    source.p_mask(source_positions(source_loop)) = 1;

    % create a display mask to display the transducer
    display_mask = source.p_mask + sensor.mask;

    % run the simulation
    input_args = {'DisplayMask', display_mask, 'PlotScale', [-0.5, 0.5]};
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

    % average the data recorded for each grid point to simulate the
    % measured signal from a large aperture, single element, detector
    single_element_data(:, source_loop) = sum(sensor_data, 1);

end

% =========================================================================
% VISUALISATION
% =========================================================================

% plot source points and sensor mask
figure
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, circle + sensor.mask, [-1, 1])
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image

% get a suitable scaling factor for the time axis
[~, t_scale, t_prefix] = scaleSI(kgrid.t_array(end));

% plot the time series recorded for each of the sources
figure;
plot(kgrid.t_array * t_scale, single_element_data);
xlabel(['Time [' t_prefix 's]']);
ylabel('Pressure [au]');
title('Time Series For Each Source Direction');

% calculate angle between source and centre of detector face
angles = atan((kgrid.y(source_positions)) ./ kgrid.x(source_positions));

% plot the maximum amplitudes for each of the sources, showing that the
% detector sensitivity falls off at low angles as expected.
figure;
plot(angles, max(single_element_data(200:350, :)), 'o')
xlabel('Angle Between Source and Centre of Detector Face [rad]');
ylabel('Maximum Detected Pressure [au]');