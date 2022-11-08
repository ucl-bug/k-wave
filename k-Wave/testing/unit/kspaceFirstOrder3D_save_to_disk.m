function test_pass = kspaceFirstOrder3D_save_to_disk(~, ~)
% DESCRIPTION:
%       Unit test to save input file to disk using kspaceFirstOrder3D.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 21st February 2017
%       last update - 21st February 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017 Bradley Treeby

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

%#ok<*INUSD>

% set pass variable
test_pass = true;

% pathname for the input file
pathname = tempdir;

% input filename (this must have the .h5 extension)
filename_1  = 'example_input_1.h5';
filename_2  = 'example_input_2.h5';

% =========================================================================
% PRESSURE SOURCE
% =========================================================================

% create the computational grid
Nx = 256;                   % number of grid points in the x direction
Ny = 128;                   % number of grid points in the y direction
Nz = 64;                    % number of grid points in the z direction
dx = 0.1e-3;                % grid point spacing in the x direction [m]
dy = 0.1e-3;                % grid point spacing in the y direction [m]
dz = 0.1e-3;                % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% set the size of the PML
pml_size = 10;              % [grid points]

% define a scattering ball
ball_radius = 20;           % [grid points]
ball_x      = Nx/2 + 40;    % [grid points]
ball_y      = Ny/2;         % [grid points]
ball_z      = Nz/2;         % [grid points]
ball        = makeBall(Nx, Ny, Nz, ball_x, ball_y, ball_z, ball_radius);

% define the properties of the propagation medium
medium.sound_speed            = 1500*ones(Nx, Ny, Nz);	 % [m/s]
medium.sound_speed(ball == 1) = 1800;                    % [m/s]
medium.density                = 1000*ones(Nx, Ny, Nz);   % [kg/m^3]
medium.density(ball == 1)     = 1200;                    % [kg/m^3]
medium.alpha_coeff            = 0.75;                    % [dB/(MHz^y cm)]
medium.alpha_power            = 1.5;

% create the time array
Nt = 1200;                  % number of time steps
dt = 15e-9;                 % time step [s]
kgrid.t_array = (0:(Nt-1))*dt;

% define a square source element facing in the x-direction
source_y_size = 60;         % [grid points]
source_z_size = 30;         % [grid points]
source.p_mask = zeros(Nx, Ny, Nz);
source.p_mask(1 + pml_size, Ny/2 - source_y_size/2:Ny/2 + source_y_size/2, Nz/2 - source_z_size/2:Nz/2 + source_z_size/2) = 1;

% define a time varying sinusoidal source
source_freq     = 2e6;      % [Hz]
source_strength = 0.5e6;    % [Pa]
source.p        = source_strength*sin(2*pi*source_freq*kgrid.t_array);

% filter the source to remove high frequencies not supported by the grid
source.p = filterTimeSeries(kgrid, medium, source.p);

% define a sensor mask through the central plane
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(:, :, Nz/2) = 1;

% set the input arguments
input_args = {'PMLSize', pml_size};

% save the input data to disk and then exit
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', [pathname filename_1]);

% =========================================================================
% VELOCITY SOURCE
% =========================================================================

% define a square source element facing in the x-direction
clear source;
source.u_mask = zeros(Nx, Ny, Nz);
source.u_mask(1 + pml_size, Ny/2 - source_y_size/2:Ny/2 + source_y_size/2, Nz/2 - source_z_size/2:Nz/2 + source_z_size/2) = 1;

% define a time varying sinusoidal source
source.ux = source_strength*sin(2*pi*source_freq*kgrid.t_array);

% filter the source to remove high frequencies not supported by the grid
source.ux = filterTimeSeries(kgrid, medium, source.ux);

% save the input data to disk and then exit
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', [pathname filename_2]);

% =========================================================================
% DELETE FILES
% =========================================================================

% check input files were created, if so, delete them
if exist([pathname filename_1], 'file')
    delete([pathname filename_1]);
else
    test_pass = false;
end
if exist([pathname filename_2], 'file')
    delete([pathname filename_2]);
else
    test_pass = false;
end