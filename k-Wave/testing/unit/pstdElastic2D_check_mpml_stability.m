function test_pass = pstdElastic2D_check_mpml_stability(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to test the stability of the pml and m-pml.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 12 August 2014
%       last update - 12 August 2014
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox

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

% set pass variable
test_pass = true;

% create the computational grid
PML_SIZE = 20;
Nx = 128 - 2*PML_SIZE;      % number of grid points in the x (row) direction
Ny = 128 - 2*PML_SIZE;      % number of grid points in the y (column) direction
dx = 0.1e-3;                % grid point spacing in the x direction [m]
dy = 0.1e-3;                % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the upper layer of the propagation medium
medium.sound_speed_compression = 1500*ones(Nx, Ny); % [m/s]
medium.sound_speed_shear       = zeros(Nx, Ny);     % [m/s]
medium.density                 = 1000*ones(Nx, Ny); % [kg/m^3]

% define the properties of the lower layer of the propagation medium
medium.sound_speed_compression(Nx/2:end, :) = 2000; % [m/s]
medium.sound_speed_shear(Nx/2:end, :)       = 1000;  % [m/s]
medium.density(Nx/2:end, :)                 = 1200; % [kg/m^3]

% create the time array
cfl   = 0.1;        % Courant-Friedrichs-Lewy number
t_end = 10e-6;      % [s]
kgrid.t_array = makeTime(kgrid, max(medium.sound_speed_compression(:)), cfl, t_end);

% define the source
source.s_mask = makeCircle(Nx, Ny, Nx/2 - 10, Ny/2, 20, [1, 3]);
source.s_mask(50:end, :) = 0;
source.sxx = toneBurst(1/kgrid.dt, 1e6, 3);
source.syy = source.sxx;

% define sensor
sensor.record = {'u_final'};
sensor.mask = [];

% define a custom display mask showing the position of the interface from
% the fluid side
display_mask = false(Nx, Ny);
display_mask(Nx/2 - 1, :) = 1;

% define input arguments
input_args = {'PlotScale', [-0.5, 0.5, -0.01, 0.01],...
    'DisplayMask', display_mask, 'PMLAlpha', 2, 'PMLInside', false,...
    'PMLSize', PML_SIZE, 'PlotSim', plot_simulations};

% run the simulations
sensor_data_pml = pstdElastic2D(kgrid, medium, source, sensor, input_args{:}, 'MultiAxialPMLRatio', 0);
sensor_data_mpml = pstdElastic2D(kgrid, medium, source, sensor, input_args{:}, 'MultiAxialPMLRatio', 0.1);

% check magnitudes
pml_max = max([sensor_data_pml.ux_final(:); sensor_data_pml.uy_final(:)]);
mpml_max = max([sensor_data_mpml.ux_final(:); sensor_data_mpml.uy_final(:)]);

% set reference magnitude (initial source)
ref_max = 1 / ( max(medium.sound_speed_shear(:)) * max(medium.density(:)) );

% check results - the test should fail if the pml DOES work (i.e., it
% doesn't become unstable), or if the m-pml DOESN'T work (i.e., it does
% become unstable)
if pml_max < ref_max || mpml_max > ref_max
    test_pass = false;
end

% plot if required
if plot_comparisons
    figure;
    subplot(2, 2, 1)
    imagesc(sensor_data_pml.ux_final);
    axis image;
    colorbar;
    title('u_x with PML');
    
    subplot(2, 2, 2)
    imagesc(sensor_data_pml.uy_final);
    axis image;
    colorbar;
    title('u_y with PML');
    
    subplot(2, 2, 3)
    imagesc(sensor_data_mpml.ux_final);
    axis image;
    colorbar;
    title('u_x with M-PML');
    
    subplot(2, 2, 4)
    imagesc(sensor_data_mpml.uy_final);
    axis image;
    colorbar;
    title('u_y with M-PML');
end