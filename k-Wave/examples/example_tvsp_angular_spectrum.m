% Holographic Projection Using The Angular Spectrum Method Example
%
% This example demonstrates how to project time-domain data between
% parallel planes using the angular spectrum method. It builds on the
% Simulations In Three Dimensions Example. 
%
% author: Bradley Treeby
% date: 12th February 2019
% last update: 13th February 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019 Bradley Treeby

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
Nx = 44;                    % number of grid points in the x direction
Ny = 44;                    % number of grid points in the y direction
Nz = 44;                    % number of grid points in the z direction
dx = 1e-4;                  % grid point spacing [m]
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% define the properties of the propagation medium
c0 = 1500;                  % [m/s]
medium.sound_speed = c0;    % [m/s]

% define two planar sensor masks using opposing corners of a cuboid
z_pos_1 = 10;               % [grid points]
z_pos_2 = Nz;               % [grid points]
sensor.mask = [1, 1, z_pos_1, Nx, Ny, z_pos_1; ...
               1, 1, z_pos_2, Nx, Ny, z_pos_2].';

% calculate distance between the two planes
proj_dist = (z_pos_2 - z_pos_1) * dx;

% define time axis, modifying dt so that the z/c0 offset between the two
% planes is an integer number of time steps
kgrid.makeTime(c0);
time_offset = proj_dist / c0;
time_offset_steps = round(time_offset / kgrid.dt);
kgrid.setTime(kgrid.Nt, time_offset / time_offset_steps);

% define source mask as focused bowl
radius_of_curvature = 33;       % [grid points]
aperture_diameter = 25;         % [grid points]
source.p_mask = makeBowl([Nx, Ny, Nz], [Nx/2, Ny/2, 1], ...
    radius_of_curvature, aperture_diameter, [Nx/2, Ny/2, Nz]);

% define source signal as a three-cycle tone burst at 2 MHz
source_freq = 2e6;              % [Hz]
source_cycles = 3;
source.p = toneBurst(1/kgrid.dt, source_freq, source_cycles);

% run k-Wave simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
    'PMLInside', false', ...
    'PMLSize', 'auto', ...
    'DataCast', 'single');

% assign the sensor data from the two planes
plane_1_kw = squeeze(sensor_data(1).p);
plane_2_kw = squeeze(sensor_data(2).p);

% run projection from the first plane to the second plane using the angular
% spectrum method 
[~, plane_2_as] = angularSpectrum(plane_1_kw, dx, ...
    kgrid.dt, proj_dist, c0, 'GridExpansion', 50);

% trim the time dimension of the plane 2 data calculated using
% kspaceFirstOrder3D and angularSpectrum to the same length
plane_2_kw = plane_2_kw(:, :, time_offset_steps + 1:end);
plane_2_as = plane_2_as(:, :, 1:size(plane_2_kw, 3));

% compute percentage difference
diff = 100 * abs(plane_2_kw - plane_2_as) ./ max(abs(plane_2_kw(:)));

% =========================================================================
% VISUALISATION
% =========================================================================

% plot maximum intensity projection (MIP) through the data
figure;

subplot(2, 2, 1);
imagesc(1e3 * kgrid.y_vec, 1e3 * kgrid.x_vec, max(plane_2_kw, [], 3));
axis image;
ylabel('x-position [mm]');
xlabel('y-position [mm]');
cb = colorbar;
ylabel(cb, '[Pa]');
title('MIP kspaceFirstOrder3D');

subplot(2, 2, 2);
imagesc(1e3 * kgrid.y_vec, 1e3 * kgrid.x_vec, max(plane_2_as, [], 3));
axis image;
ylabel('x-position [mm]');
xlabel('y-position [mm]');
cb = colorbar;
ylabel(cb, '[Pa]');
title('MIP angularSpectrum');

subplot(2, 2, 3);
imagesc(1e3 * kgrid.y_vec, 1e3 * kgrid.x_vec, max(diff, [], 3));
axis image;
ylabel('x-position [mm]');
xlabel('y-position [mm]');
cb = colorbar;
ylabel(cb, '[%]');
title('MIP Difference');

colormap(parula(1024));

% plot complete time series data
figure;

subplot(3, 1, 1);
imagesc(reshape(plane_2_kw, Nx * Ny, []));
xlabel('Time Index');
ylabel('Linear Grid Index');
cb = colorbar;
ylabel(cb, '[Pa]');
title('Simulated data using kspaceFirstOrder3D');

subplot(3, 1, 2);
imagesc(reshape(plane_2_as, Nx * Ny, []));
xlabel('Time Index');
ylabel('Linear Grid Index');
cb = colorbar;
ylabel(cb, '[Pa]');
title('Projected data using angularSpectrum');

subplot(3, 1, 3);
imagesc(reshape(diff, Nx * Ny, []));
xlabel('Time Index');
ylabel('Linear Grid Index');
cb = colorbar;
ylabel(cb, '[%]');
title('Difference');

colormap(parula(1024));
scaleFig(1, 1.5);

% plot single time trace in the middle of the plane
figure;

subplot(2, 1, 1);
plot(squeeze(plane_2_kw(Nx/2, Ny/2, :)));
hold on;
plot(squeeze(plane_2_as(Nx/2, Ny/2, :)), '--');
xlabel('Time Index');
ylabel('Pressure');
legend('kspaceFirstOrder3D', 'angularSpectrum', 'Location', 'Best');
title('Time signal from centre of projected plane');

subplot(2, 1, 2);
plot(100 * abs(squeeze(plane_2_kw(Nx/2, Ny/2, :)) - squeeze(plane_2_as(Nx/2, Ny/2, :))) ./ max(abs(squeeze(plane_2_kw(Nx/2, Ny/2, :)))));
xlabel('Time Index');
ylabel('Difference [%]');