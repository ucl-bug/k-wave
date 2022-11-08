% Image Reconstruction With Bandlimited Sensors Example
%
% This example demonstrates how the bandwidth of sensor elements can give
% rise to artefacts in time reversal photoacoustic image reconstruction. It
% builds on the previous 2D time reversal examples. 
%
% author: Bradley Treeby
% date: 27th August 2010
% last update: 4th June 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2017 Bradley Treeby

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
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]

% compute the time array
kgrid.makeTime(medium.sound_speed);

% create initial pressure distribution using makeDisc
disc_magnitude = 5; % [Pa]
disc_x_pos = 50;    % [grid points]
disc_y_pos = 50;    % [grid points]
disc_radius = 8;    % [grid points]
disc_1 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_magnitude = 4; % [Pa]
disc_x_pos = 65;    % [grid points]
disc_y_pos = 85;    % [grid points]
disc_radius = 7;    % [grid points]
disc_2 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_magnitude = 3; % [Pa]
disc_x_pos = 80;    % [grid points]
disc_y_pos = 60;    % [grid points]
disc_radius = 5;    % [grid points]
disc_3 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

source.p0 = disc_1 + disc_2 + disc_3;

% define a centered circular sensor
sensor_radius = 44; % [grid points]
sensor.mask = makeCircle(Nx, Ny, Nx/2, Ny/2, sensor_radius);

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% filter the sensor data using a high pass filter
Fs = 1/kgrid.dt;        % [Hz]
cutoff_freq = 1e6;      % [Hz]
sensor_data_high_pass = zeros(size(sensor_data));
for index = 1:sum(sensor.mask(:))
    sensor_data_high_pass(index, :) = applyFilter(sensor_data(index, :), Fs, cutoff_freq, 'HighPass', 'ZeroPhase', true);
end

% filter the sensor data using a low pass filter
Fs = 1/kgrid.dt;        % [Hz]
cutoff_freq = 1e6;      % [Hz]
sensor_data_low_pass = zeros(size(sensor_data));
for index = 1:sum(sensor.mask(:))
    sensor_data_low_pass(index, :) = applyFilter(sensor_data(index, :), Fs, cutoff_freq, 'LowPass', 'ZeroPhase', true);
end

% filter the sensor data using a Gaussian filter
Fs = 1/kgrid.dt;        % [Hz]
center_freq = 3e6;      % [Hz]
bandwidth = 100;        % [%]
sensor_data_gaussian = gaussianFilter(sensor_data, Fs, center_freq, bandwidth);

% =========================================================================
% TIME REVERSAL IMAGE RECONSTRUCTION
% =========================================================================

% reset the initial pressure
source.p0 = 0;

% assign the unfiltered time reversal data
sensor.time_reversal_boundary_data = sensor_data;

% run the time reversal reconstruction
p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% assign the high pass filtered time reversal data
sensor.time_reversal_boundary_data = sensor_data_high_pass;

% re-run the time reversal reconstruction
p0_recon_high_pass = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% assign the low pass filtered time reversal data
sensor.time_reversal_boundary_data = sensor_data_low_pass;

% re-run the time reversal reconstruction
p0_recon_low_pass = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% assign the gaussian filtered time reversal data
sensor.time_reversal_boundary_data = sensor_data_gaussian;

% re-run the time reversal reconstruction
p0_recon_gaussian = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the reconstructed initial pressure with no filtering
figure;
mx = max(abs(p0_recon(:)));
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0_recon + mx*sensor.mask, [-mx, mx]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;

% plot the reconstructed initial pressure with high pass filtering
figure;
mx = max(abs(p0_recon_high_pass(:)));
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0_recon_high_pass + mx*sensor.mask, [-mx, mx]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;

% plot the reconstructed initial pressure with low pass filtering
figure;
mx = max(abs(p0_recon_low_pass(:)));
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0_recon_low_pass + mx*sensor.mask, [-mx, mx]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;

% plot the reconstructed initial pressure with Gaussian filtering
figure;
mx = max(abs(p0_recon_gaussian(:)));
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0_recon_gaussian + mx*sensor.mask, [-mx, mx]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;

% plot a profile through the top disc in each reconstruction
figure;
disc_x_pos = 50;
y = kgrid.y_vec * 1e3;
plot(y, p0_recon(disc_x_pos, :), 'k-',...
     y, p0_recon_high_pass(disc_x_pos, :), 'b-',...
     y, p0_recon_low_pass(disc_x_pos, :), 'g-',...
     y, p0_recon_gaussian(disc_x_pos, :), 'r-');
set(gca, 'XLim', [min(y), max(y)], 'YLim', [-1, 5.1]);
ylabel('Pressure [au]');
xlabel('y-position [mm]');
legend('Not Filtered', 'High Pass Filtered', 'Low Pass Filtered', 'Gaussian Filtered');