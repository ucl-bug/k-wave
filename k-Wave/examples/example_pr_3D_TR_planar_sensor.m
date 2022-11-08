% 3D Time Reversal Reconstruction For A Planar Sensor Example
%
% This example demonstrates the use of k-Wave for the reconstruction of a
% three-dimensional photoacoustic wave-field recorded  over a planar array
% of sensor elements.  The sensor data is simulated and then time-reversed
% using kspaceFirstOrder3D. It builds on the 3D FFT Reconstruction For A
% Planar Sensor and 2D Time Reversal Reconstruction For A Line Sensor
% examples.
%
% author: Bradley Treeby
% date: 10th July 2009
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

% change scale to 2 to reproduce the higher resolution figures used in the
% help file
scale = 1;

% create the computational grid
PML_size = 10;                  % size of the PML in grid points
Nx = 32 * scale - 2 * PML_size; % number of grid points in the x direction
Ny = 64 * scale - 2 * PML_size; % number of grid points in the y direction
Nz = 64 * scale - 2 * PML_size; % number of grid points in the z direction
dx = 0.2e-3 / scale;            % grid point spacing in the x direction [m]
dy = 0.2e-3 / scale;            % grid point spacing in the y direction [m]
dz = 0.2e-3 / scale;            % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the propagation medium
medium.sound_speed = 1500;      % [m/s]

% create initial pressure distribution using makeBall
ball_magnitude = 10;            % [Pa]
ball_radius = 3 * scale;        % [grid points]
p0 = ball_magnitude * makeBall(Nx, Ny, Nz, Nx/2, Ny/2, Nz/2, ball_radius);

% smooth the initial pressure distribution and restore the magnitude
p0 = smooth(p0, true);

% assign to the source structure
source.p0 = p0;

% define a binary planar sensor
sensor.mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
sensor.mask(1, :, :) = 1;

% create the time array
kgrid.makeTime(medium.sound_speed);

% set the input arguements
input_args = {'PMLSize', PML_size, 'PMLInside', false, ...
    'PlotPML', false, 'Smooth', false, 'DataCast', 'single'};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% reset the initial pressure
source.p0 = 0;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

% run the time-reversal reconstruction
p0_recon = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

% add first order compensation for only recording over a half plane
p0_recon = 2 * p0_recon;

% apply a positivity condition
p0_recon(p0_recon < 0) = 0;

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the initial pressure
figure;
plot_scale = [-10, 10];
subplot(2, 2, 1);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, squeeze(p0(:, :, Nz/2)), plot_scale);
title('x-y plane');
axis image;

subplot(2, 2, 2);
imagesc(kgrid.z_vec * 1e3, kgrid.x_vec * 1e3, squeeze(p0(:, Ny/2, :)), plot_scale);
title('x-z plane');
axis image;
xlabel('(All axes in mm)');

subplot(2, 2, 3);
imagesc(kgrid.z_vec * 1e3, kgrid.y_vec * 1e3, squeeze(p0(Nx/2, :, :)), plot_scale);
title('y-z plane');
axis image;
colormap(getColorMap);

% plot the reconstructed initial pressure
figure;
subplot(2, 2, 1);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, squeeze(p0_recon(:, :, Nz/2)), plot_scale);
title('x-y plane');
axis image;

subplot(2, 2, 2);
imagesc(kgrid.z_vec * 1e3, kgrid.x_vec * 1e3, squeeze(p0_recon(:, Ny/2, :)), plot_scale);
title('x-z plane');
axis image;
xlabel('(All axes in mm)');

subplot(2, 2, 3);
imagesc(kgrid.z_vec * 1e3, kgrid.y_vec * 1e3, squeeze(p0_recon(Nx/2, :, :)), plot_scale);
title('y-z plane');
axis image;
colormap(getColorMap);