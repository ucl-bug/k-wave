function test_pass = pstdElastic2D_compare_time_reversal_with_kspaceFirstOrder2D(plot_comparisons, plot_simulations)
% DESCRIPTION:
%       Unit test to compare time reversal image reconstruction using
%       elastic and fluid codes assuming the shear wave speed is zero
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 13th August 2014
%       last update - 25th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2019 Bradley Treeby and Ben Cox

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

% set literals
COMPARISON_THRESH   = 1e-14;

% set pass variable
test_pass = true;

% create the computational grid
PML_size = 20;          % size of the PML in grid points
Nx = 128 - 2*PML_size;  % number of grid points in the x (row) direction
Ny = 256 - 2*PML_size;  % number of grid points in the y (column) direction
dx = 0.1e-3;            % grid point spacing in the x direction [m]
dy = 0.1e-3;            % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium_elastic.sound_speed_compression = 1500;	% [m/s]
medium_elastic.sound_speed_shear = 0;
medium_elastic.density = 1000;

medium_fluid.sound_speed = 1500;
medium_fluid.density = 1000;

% create initial pressure distribution using makeDisc
disc_magnitude = 5; % [au]
disc_x_pos = 60;    % [grid points]
disc_y_pos = 140;  	% [grid points]
disc_radius = 5;    % [grid points]
disc_2 = disc_magnitude*makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_x_pos = 30;    % [grid points]
disc_y_pos = 110; 	% [grid points]
disc_radius = 8;    % [grid points]
disc_1 = disc_magnitude*makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

% smooth the initial pressure distribution and restore the magnitude
p0 = smooth(disc_1 + disc_2, true);

% assign to the source structure
source.p0 = p0;

% define a binary line sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, :) = 1;

% create the time array
CFL = 0.3;
kgrid.t_array = makeTime(kgrid, medium_fluid.sound_speed, CFL);

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PMLSize', PML_size, 'Smooth', false, ...
    'PlotSim', plot_simulations, 'PlotPML', false};

% run the forward simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium_fluid, source, sensor, input_args{:});

% reset the initial pressure
source.p0 = 0;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

% run the time reversal reconstructions
p0_recon_elastic = pstdElastic2D(kgrid, medium_elastic, source, sensor, input_args{:}, 'MultiAxialPMLRatio', 0);
p0_recon_fluid = kspaceFirstOrder2D(kgrid, medium_fluid, source, sensor, input_args{:}, 'UsekSpace', false);

% compute L2 error
L2_error = sqrt( sum( (p0_recon_elastic(:) - p0_recon_fluid(:)).^2 ) / sum( p0_recon_fluid(:).^2 ) );

% set pass flag
if L2_error > COMPARISON_THRESH
    test_pass = false;
end

% plot the reconstructed initial pressures
if plot_comparisons

    figure;
    subplot(3, 1, 1)
    imagesc(p0_recon_elastic);
    axis image;
    colorbar

    subplot(3, 1, 2)
    imagesc(p0_recon_fluid);
    axis image;
    colorbar

    subplot(3, 1, 3)
    imagesc(p0_recon_elastic - p0_recon_fluid);
    axis image;
    colorbar
    
end