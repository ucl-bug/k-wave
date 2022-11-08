% 3D FFT Reconstruction For A Planar Sensor Example
%
% This example demonstrates the use of k-Wave for the reconstruction of a
% three-dimensional photoacoustic wave-field recorded over a planar array
% of sensor elements. The sensor data is simulated using kspaceFirstOrder3D
% and reconstructed using kspacePlaneRecon. It builds on the Simulations In
% Three Dimensions and 2D FFT Reconstruction For A Line Sensor examples.
%
% author: Bradley Treeby
% date: 3rd July 2009
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
medium.sound_speed = 1500;	% [m/s]

% create initial pressure distribution using makeBall
ball_magnitude = 10;        % [Pa]
ball_radius = 3 * scale;    % [grid points]
p0 = ball_magnitude * makeBall(Nx, Ny, Nz, Nx/2, Ny/2, Nz/2, ball_radius);

% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(p0, true);

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

% reshape sensor data to y, z, t
sensor_data_rs = reshape(sensor_data, Ny, Nz, kgrid.Nt);

% reconstruct the initial pressure
p_xyz = kspacePlaneRecon(sensor_data_rs, kgrid.dy, kgrid.dz, kgrid.dt, ...
    medium.sound_speed, 'DataOrder', 'yzt', 'PosCond', true, 'Plot', true);

% =========================================================================
% VISUALISATION
% =========================================================================

% define a k-space grid using the dimensions of p_xyz
[Nx_recon, Ny_recon, Nz_recon] = size(p_xyz);
kgrid_recon = kWaveGrid(Nx_recon, kgrid.dt * medium.sound_speed, Ny_recon, kgrid.dy, Nz_recon, kgrid.dz);

% define a k-space grid with the same z-spacing as p0
[Nx_p0, Ny_p0, Nz_p0] = size(source.p0);
kgrid_interp = kWaveGrid(Nx_p0, kgrid.dx, Ny_p0, kgrid.dy, Nz_p0, kgrid.dz);

% resample the p_xyz to be the same size as p0; for a matrix indexed as 
% [M, N, P], the axis variables passed to interp3 must be given in the 
% order N, M, P
p_xyz_rs = interp3(kgrid_recon.y - min(kgrid_recon.y(:)), ...
                   kgrid_recon.x - min(kgrid_recon.x(:)), ...
                   kgrid_recon.z - min(kgrid_recon.z(:)), ...
                   p_xyz, ...
                   kgrid_interp.y - min(kgrid_interp.y(:)), ...
                   kgrid_interp.x - min(kgrid_interp.x(:)), ...
                   kgrid_interp.z - min(kgrid_interp.z(:)));

% plot the initial pressure and sensor surface in voxel form
voxelPlot(double(p0 | sensor.mask));
set(gca, 'Projection', 'perspective');
view([0, 99]);

% plot the initial pressure
figure;
plot_scale = [-10, 10];
subplot(2, 2, 1);
imagesc(kgrid_interp.y_vec * 1e3, kgrid_interp.x_vec * 1e3, squeeze(source.p0(:, :, Nz/2)), plot_scale);
title('x-y plane');
axis image;

subplot(2, 2, 2);
imagesc(kgrid_interp.z_vec * 1e3, kgrid_interp.x_vec * 1e3, squeeze(source.p0(:, Ny/2, :)), plot_scale);
title('x-z plane');
axis image;
xlabel('(All axes in mm)');

subplot(2, 2, 3);
imagesc(kgrid_interp.z_vec * 1e3, kgrid_interp.y_vec * 1e3, squeeze(source.p0(Nx/2, :, :)), plot_scale);
title('y-z plane');
axis image;
colormap(getColorMap);

% plot the reconstructed initial pressure
figure;
subplot(2, 2, 1);
imagesc(kgrid_interp.y_vec * 1e3, kgrid_interp.x_vec * 1e3, squeeze(p_xyz_rs(:, :, Nz/2)), plot_scale);
title('x-y plane');
axis image;

subplot(2, 2, 2);
imagesc(kgrid_interp.z_vec * 1e3, kgrid_interp.x_vec * 1e3, squeeze(p_xyz_rs(:, Ny/2, :)), plot_scale);
title('x-z plane');
axis image;
xlabel('(All axes in mm)');

subplot(2, 2, 3);
imagesc(kgrid_interp.z_vec * 1e3, kgrid_interp.y_vec * 1e3, squeeze(p_xyz_rs(Nx/2, :, :)), plot_scale);
title('y-z plane');
axis image;
colormap(getColorMap);

% view reconstruction slice by slice
flyThrough(p_xyz_rs);