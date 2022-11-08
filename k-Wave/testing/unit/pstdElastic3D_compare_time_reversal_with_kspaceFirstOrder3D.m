function test_pass = pstdElastic3D_compare_time_reversal_with_kspaceFirstOrder3D(plot_comparisons, plot_simulations)
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
PML_size = 10;              % size of the PML in grid points
Nx = 32 - 2*PML_size; % number of grid points in the x direction
Ny = 64 - 2*PML_size; % number of grid points in the y direction
Nz = 64 - 2*PML_size; % number of grid points in the z direction
dx = 0.2e-3;          % grid point spacing in the x direction [m]
dy = 0.2e-3;          % grid point spacing in the y direction [m]
dz = 0.2e-3;          % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the propagation medium
medium_elastic.sound_speed_compression = 1500;	% [m/s]
medium_elastic.sound_speed_shear = 0;
medium_elastic.density = 1000;

medium_fluid.sound_speed = 1500;
medium_fluid.density = 1000;

% create initial pressure distribution using makeBall
ball_magnitude = 10;        % [au]
ball_radius = 3;    	% [grid points]
p0 = ball_magnitude*makeBall(Nx, Ny, Nz, Nx/2, Ny/2, Nz/2, ball_radius);

% smooth the initial pressure distribution and restore the magnitude
p0 = smooth(p0, true);

% assign to the source structure
source.p0 = p0;

% define a binary planar sensor
sensor.mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
sensor.mask(1, :, :) = 1;

% create the time array
kgrid.t_array = makeTime(kgrid, medium_fluid.sound_speed);

% set the input arguements
input_args = {'PMLSize', PML_size, 'PMLInside', false, 'PlotPML', false, ...
    'Smooth', false, 'PlotSim', plot_simulations};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium_fluid, source, sensor, input_args{:});

% reset the initial pressure
source.p0 = 0;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;

% run the time-reversal reconstructions
p0_recon_elastic = pstdElastic3D(kgrid, medium_elastic, source, sensor, input_args{:}, 'MultiAxialPMLRatio', 0);
p0_recon_fluid = kspaceFirstOrder3D(kgrid, medium_fluid, source, sensor, input_args{:}, 'UsekSpace', false);

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
    imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, squeeze(p0_recon_fluid(:, :, Nz/2)));
    axis image;
    colorbar;

    subplot(3, 1, 2)
    imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, squeeze(p0_recon_elastic(:, :, Nz/2)));
    axis image;
    colorbar;

    subplot(3, 1, 3)
    imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, squeeze(p0_recon_elastic(:, :, Nz/2)) - squeeze(p0_recon_fluid(:, :, Nz/2)));
    axis image;
    colorbar;
    
end