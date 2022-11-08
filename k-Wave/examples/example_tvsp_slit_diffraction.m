% Diffraction Through A Slit Example Example
%
% This example illustrates the diffraction of a plane acoustic wave through
% a slit. It builds on the Monopole Point Source In A Homogeneous
% Propagation Medium and Simulating Transducer Field Patterns examples.  
%
% author: Bradley Treeby
% date: 1st November 2010
% last update: 31st July 2019
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2019 Bradley Treeby

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

example_number = 1;
% 1: single slit with the source wavelength = slit size
% 2: single slit with the source wavelength = slit size / 4
% 3: double slit with the source wavelength = slit size

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid and define the bulk medium properties
scale = 1;                          % change to 2 to produce higher resolution images
pml_size = 10;                      % size of the perfectly matched layer [grid points]
Nx = scale * 128;                   % number of grid points in the x (row) direction
Ny = Nx;                            % number of grid points in the y (column) direction
dx = 50e-3 / Nx;                    % grid point spacing in the x direction [m]
dy = dx;                            % grid point spacing in the y direction [m]
c0 = 1500;                          % [m/s]
rho0 = 1000;                        % [kg/m^3]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the ratio between the barrier and background sound speed and
% density
barrier_scale = 20;

% define the barrier and the source wavelength depending on the example
switch example_number
    case 1
        
        % create a mask of a barrier with a slit
        slit_thickness = scale * 2;             % [grid points]
        slit_width = scale * 10;                % [grid points]
        slit_x_pos = Nx - Nx/4;                 % [grid points]
        slit_offset = Ny/2 - slit_width/2 - 1;  % [grid points]
        slit_mask = zeros(Nx, Ny);
        slit_mask(slit_x_pos:slit_x_pos + slit_thickness, 1:1 + slit_offset) = 1;
        slit_mask(slit_x_pos:slit_x_pos + slit_thickness, end - slit_offset:end) = 1;
        
        % define the source wavelength to be the same as the slit size
        source_wavelength = slit_width * dx;    % [m]
        
    case 2
        
        % create a mask of a barrier with a slit
        slit_thickness = scale * 2;             % [grid points]
        slit_width = scale * 30;                % [grid points]
        slit_x_pos = Nx - Nx/4;                 % [grid points]
        slit_offset = Ny/2 - slit_width/2 - 1;  % [grid points]
        slit_mask = zeros(Nx, Ny);
        slit_mask(slit_x_pos:slit_x_pos + slit_thickness, 1:1 + slit_offset) = 1;
        slit_mask(slit_x_pos:slit_x_pos + slit_thickness, end - slit_offset:end) = 1;
        
        % define the source wavelength to be a quarter of the slit size
        source_wavelength = 0.25 * slit_width * dx;     % [m]
        
    case 3
        
        % create a mask of a barrier with a double slit
        slit_thickness = scale * 2;             % [grid points]
        slit_width = scale * 10;                % [grid points]
        slit_spacing = scale * 20;              % [grid points]
        slit_x_pos = Nx - Nx/4;                 % [grid points]
        slit_offset = Ny/2 - slit_width - slit_spacing/2 - 1; 
        slit_mask = zeros(Nx, Ny);
        slit_mask(slit_x_pos:slit_x_pos + slit_thickness, 1:1 + slit_offset) = 1;
        slit_mask(slit_x_pos:slit_x_pos + slit_thickness, slit_offset + slit_width + 2:slit_offset + slit_width + slit_spacing + 1) = 1;
        slit_mask(slit_x_pos:slit_x_pos + slit_thickness, end - slit_offset:end) = 1;
        
        % define the source wavelength to be the same as the slit size
        source_wavelength = slit_width * dx;    % [m]
        
end

% assign the slit to the properties of the propagation medium
medium.sound_speed = c0 * ones(Nx, Ny);
medium.density = rho0 * ones(Nx, Ny);
medium.sound_speed(slit_mask == 1) = barrier_scale * c0;
medium.density(slit_mask == 1) = barrier_scale * rho0;

% assign the reference sound speed to the background medium
medium.sound_speed_ref = c0;

% find the time step at the stability limit
c_ref = medium.sound_speed_ref;
c_max = barrier_scale * c0;
k_max = max(kgrid.k(:));
dt_limit = 2 / (c_ref * k_max) * asin(c_ref / c_max);

% create the time array, with the time step just below the stability limit
dt = 0.95 * dt_limit;   % [s]
t_end = 40e-6;          % [s]
kgrid.setTime(round(t_end / dt) + 1, dt);

% create a source mask of a single line
source.u_mask = zeros(Nx, Ny);
source.u_mask(end - pml_size, :) = 1;

% create and filter the time varying sinusoidal source
source_mag = 2 / (c0 * rho0);
source_freq = c0 / source_wavelength;
source.ux = source_mag * sin(2 * pi * source_freq * kgrid.t_array);
source.ux = filterTimeSeries(kgrid, medium, source.ux);

% define the field parameters to record
sensor.mask = ones(Nx, Ny);
sensor.record = {'u_final', 'p_final'};

% set the input options
input_args = {...
    'PMLSize', pml_size, ...
    'PlotPML', false, ...
    'DisplayMask', slit_mask, ...
    'DataCast', 'single', ...
    };

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% VISUALISATION
% =========================================================================

% remove the PML from the view
sensor_data.p_final  = sensor_data.p_final (1 + pml_size:end - pml_size, 1 + pml_size:end - pml_size);
sensor_data.ux_final = sensor_data.ux_final(1 + pml_size:end - pml_size, 1 + pml_size:end - pml_size);
sensor_data.uy_final = sensor_data.uy_final(1 + pml_size:end - pml_size, 1 + pml_size:end - pml_size);
slit_mask = slit_mask(1 + pml_size:end - pml_size, 1 + pml_size:end - pml_size);

% get plot axes
y_ax = kgrid.y_vec(1 + pml_size:end - pml_size) * 1e3;
x_ax = kgrid.x_vec(1 + pml_size:end - pml_size) * 1e3;

% plot the final wave-field
figure;
mx = max(abs(sensor_data.p_final(:)));
sensor_data.p_final(slit_mask == 1) = mx;
imagesc(y_ax, x_ax, sensor_data.p_final, [-mx, mx]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('p');

% plot the final wave-field
figure;
mx = max(abs(sensor_data.ux_final(:)));
sensor_data.ux_final(slit_mask == 1) = mx;
imagesc(y_ax, x_ax, sensor_data.ux_final, [-mx, mx]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('ux');

% plot the final wave-field
figure;
mx = max(abs(sensor_data.uy_final(:)));
sensor_data.uy_final(slit_mask == 1) = mx;
imagesc(y_ax, x_ax, sensor_data.uy_final, [-mx, mx]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('uy');