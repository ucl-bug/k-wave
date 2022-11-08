% Image Reconstruction With Directional Sensors Example
%
% This example demonstrates how the directionality of sensor elements can
% give rise to artefacts in time reversal photoacoustic image
% reconstruction. It builds on the Sensor Element Directivity in 2D and 2D
% Time Reversal Reconstruction For A Line Sensor examples.
%
% author: Bradley Treeby & Ben Cox
% date: 18th January 2010
% last update: 25th July 2019
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2019 Bradley Treeby and Ben Cox

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
PML_size = 20;              % size of the PML in grid points
Nx = 128 - 2 * PML_size;    % number of grid points in the x direction
Ny = 256 - 2 * PML_size;    % number of grid points in the y direction
dx = 0.1e-3;                % grid point spacing in the x direction [m]
dy = 0.1e-3;                % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;	% [m/s]

% create initial pressure distribution using makeDisc
disc_magnitude = 5;         % [Pa]
disc_x_pos = 60;            % [grid points]
disc_y_pos = 140;           % [grid points]
disc_radius = 5;            % [grid points]
disc_2 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_x_pos = 30;            % [grid points]
disc_y_pos = 110;           % [grid points]
disc_radius = 8;            % [grid points]
disc_1 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

% smooth the initial pressure distribution and restore the magnitude
p0 = smooth(disc_1 + disc_2, true);

% assign to the source structure
source.p0 = p0;

% define a four-sided, square sensor
sensor.mask = zeros(kgrid.Nx, kgrid.Ny);
sensor.mask(1, :) = 1;
sensor.mask(end, :) = 1;
sensor.mask(:, 1) = 1;
sensor.mask(:, end) = 1;

% create the time array
kgrid.makeTime(medium.sound_speed);

% set the input arguements
input_args = {'PMLInside', false, 'PMLSize', PML_size, 'PlotPML', false, 'Smooth', false};

% run the simulation for omnidirectional detector elements
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% define the directionality of the sensor elements
sensor.directivity_angle = zeros(kgrid.Nx, kgrid.Ny);
sensor.directivity_angle(1, :) = 0;    	 % max sensitivity in x direction
sensor.directivity_angle(end, :) = 0;  	 % max sensitivity in x direction
sensor.directivity_angle(:, 1) = pi/2;   % max sensitivity in y direction
sensor.directivity_angle(:, end) = pi/2; % max sensitivity in y direction

% define the directivity size
sensor.directivity_size = 20 * kgrid.dx;

% run the simulation with directional elements
sensor_data_directional = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% reset the initial pressure
source.p0 = 0;

% assign the time reversal data for the omnidirectional case
sensor.time_reversal_boundary_data = sensor_data;

% run the time reversal reconstruction
p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% assign the time reversal data for the directional case
sensor.time_reversal_boundary_data = sensor_data_directional;

% run the time reversal reconstruction with directional elements
p0_recon_directional = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the initial pressure and sensor distribution
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0 + sensor.mask * disc_magnitude, [-disc_magnitude, disc_magnitude]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
colorbar;
scaleFig(1, 0.65);

% plot the reconstructed initial pressure 
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0_recon, [-disc_magnitude, disc_magnitude]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
colorbar;
scaleFig(1, 0.65);

% plot the reconstructed initial pressure with directivity
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p0_recon_directional, [-disc_magnitude, disc_magnitude]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
colorbar;
scaleFig(1, 0.65);

% plot a profile for comparison
figure;
plot(kgrid.y_vec * 1e3, p0(disc_x_pos, :), 'k-',...
     kgrid.y_vec * 1e3, p0_recon(disc_x_pos, :), 'r--',...
     kgrid.y_vec * 1e3, p0_recon_directional(disc_x_pos, :), 'b:');
xlabel('y-position [mm]');
ylabel('Pressure');
legend('True', 'Omnidirectional','Directional');
axis tight;
set(gca, 'YLim', [0, 5.1]);