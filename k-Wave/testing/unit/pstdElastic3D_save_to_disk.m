function test_pass = pstdElastic3D_save_to_disk(~, ~)
% DESCRIPTION:
%     Unit test to save input file to disk using pstdElastic3D.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 1st May 2017
%     last update - 1st May 2017
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
filename  = 'example_input_1.h5';

% =========================================================================
% PRESSURE SOURCE
% =========================================================================

% create the computational grid
Nx = 96;                   % number of grid points in the x direction
Ny = 80;                   % number of grid points in the y direction
Nz = 64;                    % number of grid points in the z direction
dx = 0.1e-3;                % grid point spacing in the x direction [m]
dy = 0.1e-3;                % grid point spacing in the y direction [m]
dz = 0.1e-3;                % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% set the size of the PML
pml_size = 10;              % [grid points]

% define the properties of the upper layer of the propagation medium
medium.sound_speed_compression = 1500 * ones(Nx, Ny, Nz);   % [m/s]
medium.sound_speed_shear       = zeros(Nx, Ny, Nz);         % [m/s]
medium.density                 = 1000 * ones(Nx, Ny, Nz);   % [kg/m^3]
% medium.alpha_coeff_compression = 0.75 * ones(Nx, Ny, Nz);   % [dB/(MHz^y cm)]
% medium.alpha_coeff_shear       = zeros(Nx, Ny, Nz);         % [dB/(MHz^y cm)]

% define the properties of the lower layer of the propagation medium
medium.sound_speed_compression(Nx/2:end, :, :) = 2000;  % [m/s]
medium.sound_speed_shear(Nx/2:end, :, :)       = 800;   % [m/s]
medium.density(Nx/2:end, :, :)                 = 1200;  % [kg/m^3]
% medium.alpha_coeff_compression(Nx/2:end, :, :) = 1;
% medium.alpha_coeff_shear(Nx/2:end, :, :)       = 3;

% create the time array
Nt = 1200;                  % number of time steps
dt = 15e-9;                 % time step [s]
kgrid.t_array = (0:(Nt-1))*dt;

% define source mask to be a square piston
source_x_pos = 11;      % [grid points]
source_radius = 15;     % [grid points]
source.u_mask = zeros(Nx, Ny, Nz);
source.u_mask(source_x_pos, Ny/2 - source_radius + 1:Ny/2 + source_radius, Nz/2 - source_radius + 1:Nz/2 + source_radius) = 1;

% define source to be a velocity source
source_freq = 2e6;      % [Hz]
source_cycles = 3;
source_mag = 1e-6;
source.ux = source_mag*toneBurst(1/kgrid.dt, source_freq, source_cycles);

% set source focus
source.ux = focus(kgrid, source.ux, source.u_mask, [0, 0, 0], 1500);

% define a sensor mask through the central plane
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(:, :, Nz/2) = 1;

% set the input arguments
input_args = {'PMLSize', pml_size};

% save the input data to disk and then exit
pstdElastic3D(kgrid, medium, source, sensor, input_args{:}, 'SaveToDisk', [pathname filename]);

% =========================================================================
% DELETE FILES
% =========================================================================

% check input files were created, if so, delete them
if exist([pathname filename], 'file')
    delete([pathname filename]);
else
    test_pass = false;
end