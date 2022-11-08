% Simulations In Three Dimensions Example
%
% This example provides a simple demonstration of using k-Wave to model
% elastic waves in a three-dimensional heterogeneous propagation medium. It
% builds on the Explosive Source In A Layered Medium and Simulations In
% Three-Dimensions examples.
%
% author: Bradley Treeby
% date: 14th February 2014
% last update: 29th May 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2017 Bradley Treeby

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
PML_size = 10;
Nx = 64;            % number of grid points in the x direction
Ny = 64;            % number of grid points in the y direction
Nz = 64;            % number of grid points in the z direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]
dz = dx;            % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the upper layer of the propagation medium
medium.sound_speed_compression = 1500*ones(Nx, Ny, Nz); % [m/s]
medium.sound_speed_shear       = zeros(Nx, Ny, Nz);     % [m/s]
medium.density                 = 1000*ones(Nx, Ny, Nz); % [kg/m^3]

% define the properties of the lower layer of the propagation medium
medium.sound_speed_compression(Nx/2:end, :, :) = 2000;  % [m/s]
medium.sound_speed_shear(Nx/2:end, :, :)       = 800;   % [m/s]
medium.density(Nx/2:end, :, :)                 = 1200;  % [kg/m^3]

% define the absorption properties
medium.alpha_coeff_compression = 0.1; % [dB/(MHz^2 cm)]
medium.alpha_coeff_shear       = 0.5; % [dB/(MHz^2 cm)]

% create the time array
cfl = 0.1;
t_end = 5e-6;
kgrid.makeTime(max(medium.sound_speed_compression(:)), cfl, t_end);

% define source mask to be a square piston
source_x_pos = 11;      % [grid points]
source_radius = 15;     % [grid points]
source.u_mask = zeros(Nx, Ny, Nz);
source.u_mask(source_x_pos, Ny/2 - source_radius + 1:Ny/2 + source_radius, Nz/2 - source_radius + 1:Nz/2 + source_radius) = 1;

% define source to be a velocity source
source_freq = 2e6;      % [Hz]
source_cycles = 3;
source_mag = 1e-6;      % [m/s]
source.ux = source_mag * toneBurst(1/kgrid.dt, source_freq, source_cycles);

% set source focus
source.ux = focus(kgrid, source.ux, source.u_mask, [0, 0, 0], 1500);

% define sensor mask in x-y plane using cuboid corners, where a rectangular
% mask is defined using the xyz coordinates of two opposing corners in the
% form [x1, y1, z1, x2, y2, z2].'
sensor.mask = [1 + PML_size, 1 + PML_size, Nz/2, Nx - PML_size, Ny - PML_size, Nz/2].';

% record the maximum pressure in the plane
sensor.record = {'p_max'};

% define input arguments
input_args = {'PlotScale', [-2, 2, -0.1, 0.1], 'DataCast', 'single',...
    'PMLSize', PML_size, 'DisplayMask', source.u_mask};

% run the simulation with PML inside
sensor_data = pstdElastic3D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the sensor data
figure;
imagesc(sensor_data.p_max);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;