% 3D Time Reversal Reconstruction For A Spherical Sensor Example
%
% This example demonstrates the use of k-Wave for the time-reversal
% reconstruction of a three-dimensional photoacoustic wave-field recorded
% over a spherical sensor. The sensor data is simulated and then
% time-reversed using kspaceFirstOrder3D. It builds on the 3D Time Reversal
% Reconstruction For A Planar Sensor Example. 
%
% author: Bradley Treeby
% date: 15th July 2009
% last update: 25th July 2019
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2019 Bradley Treeby

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
PML_size = 10;              % size of the PML in grid points
Nx = 64 - 2 * PML_size;     % number of grid points in the x direction
Ny = Nx;                    % number of grid points in the y direction
Nz = Nx;                    % number of grid points in the z direction
dx = 0.2e-3;                % grid point spacing in the x direction [m]
dy = dx;                    % grid point spacing in the y direction [m]
dz = dx;                    % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the propagation medium
medium.sound_speed = 1500;	% [m/s]

% create initial pressure distribution using makeBall
ball_magnitude = 10;        % [Pa]
ball_x_pos = 16;            % [grid points]
ball_y_pos = 26;            % [grid points]
ball_z_pos = 22;            % [grid points]
ball_radius = 3;            % [grid points]
p0_binary = ball_magnitude * makeBall(Nx, Ny, Nz, ball_x_pos, ball_y_pos, ball_z_pos, ball_radius);

% smooth the initial pressure distribution and restore the magnitude
p0 = smooth(p0_binary, true);

% assign to the source structure
source.p0 = p0;

% define a Cartesian spherical sensor
sensor_radius = 4e-3;       % [m]
center_pos = [0, 0, 0];     % [m]
num_sensor_points = 100;
sensor_mask = makeCartSphere(sensor_radius, num_sensor_points, center_pos, true);

% assign to the sensor structure
sensor.mask = sensor_mask;

% create the time array
kgrid.makeTime(medium.sound_speed);

% set the input arguements
input_args = {'PMLSize', PML_size, 'PMLInside', false, 'PlotPML', false, ...
    'Smooth', false, 'DataCast', 'single', 'CartInterp', 'nearest'};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% reset the initial pressure
source.p0 = 0;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

% run the time-reversal reconstruction
p0_recon = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% create a binary sensor mask of an equivalent continuous sphere 
sensor_radius_grid_points = round(sensor_radius / kgrid.dx);
binary_sensor_mask = makeSphere(kgrid.Nx, kgrid.Ny, kgrid.Nz, sensor_radius_grid_points);

% assign to the sensor structure
sensor.mask = binary_sensor_mask;

% interpolate data to remove the gaps and assign to time reversal data
sensor.time_reversal_boundary_data = interpCartData(kgrid, sensor_data, sensor_mask, binary_sensor_mask);

% run the time-reversal reconstruction
p0_recon_interp = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the initial pressure and sensor surface in voxel form
voxelPlot(double(p0_binary | cart2grid(kgrid, sensor_mask)));
view([60, 20]);

% plot the initial pressure
figure;
plot_scale = [-10, 10];
subplot(2, 2, 1);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, squeeze(p0(:, :, ball_z_pos)), plot_scale);
title('x-y plane');
axis image;

subplot(2, 2, 2);
imagesc(kgrid.z_vec * 1e3, kgrid.x_vec * 1e3, squeeze(p0(:, ball_y_pos, :)), plot_scale);
title('x-z plane');
axis image;
xlabel('(All axes in mm)');

subplot(2, 2, 3);
imagesc(kgrid.z_vec * 1e3, kgrid.y_vec * 1e3, squeeze(p0(ball_x_pos, :, :)), plot_scale);
title('y-z plane');
axis image;
colormap(getColorMap);

% plot the reconstructed initial pressure
figure;
plot_scale = [-0.5, 0.5];
subplot(2, 2, 1);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, squeeze(p0_recon(:, :, ball_z_pos)), plot_scale);
title('x-y plane');
axis image;

subplot(2, 2, 2);
imagesc(kgrid.z_vec * 1e3, kgrid.x_vec * 1e3, squeeze(p0_recon(:, ball_y_pos, :)), plot_scale);
title('x-z plane');
axis image;
xlabel('(All axes in mm)');

subplot(2, 2, 3);
imagesc(kgrid.z_vec * 1e3, kgrid.y_vec * 1e3, squeeze(p0_recon(ball_x_pos, :, :)), plot_scale);
title('y-z plane');
axis image;
colormap(getColorMap);

% plot the reconstructed initial pressure
figure;
plot_scale = [-10, 10];
subplot(2, 2, 1);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, squeeze(p0_recon_interp(:, :, ball_z_pos)), plot_scale);
title('x-y plane');
axis image;

subplot(2, 2, 2);
imagesc(kgrid.z_vec * 1e3, kgrid.x_vec * 1e3, squeeze(p0_recon_interp(:, ball_y_pos, :)), plot_scale);
title('x-z plane');
axis image;
xlabel('(All axes in mm)');

subplot(2, 2, 3);
imagesc(kgrid.z_vec * 1e3, kgrid.y_vec * 1e3, squeeze(p0_recon_interp(ball_x_pos, :, :)), plot_scale);
title('y-z plane');
axis image;
colormap(getColorMap);