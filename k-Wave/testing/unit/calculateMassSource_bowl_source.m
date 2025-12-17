function test_pass = calculateMassSource_bowl_source(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to check the calculation of a mass source using
%     calculateMassSource.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 9th January 2019
%     last update - 22nd January 2020
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019- Bradley Treeby

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

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
    plot_simulations = true;
end

% set comparison threshold
comparison_threshold_proj = 2;
comparison_threshold_source_plane = 3;

% set pass variable
test_pass = true;

% =========================================================================
% PROJECT TRUE SOURCE
% =========================================================================

% simulation settings
pml_size = 10;
Nx = 64 - 2*pml_size;
Ny = 64 - 2*pml_size;
Nz = 64 - 2*pml_size;
dx = 1e-4;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% define medium properties
c0 = 1500;
medium.sound_speed = c0;

% define sensor mask
input_plane_pos = 10;
sensor.mask = [1, 1, input_plane_pos, Nx, Ny, input_plane_pos; ...
               1, Ny/2, 1 + input_plane_pos, Nx, Ny/2, Nz].';

% define 
kgrid.makeTime(c0);

% define source mask as focused bowl
radius = 33;
diameter = 25;
source.p_mask = makeBowl([Nx, Ny, Nz], [Nx/2, Ny/2, 1], radius, diameter, [Nx/2, Ny/2, Nz]);

% define input as toneburst
f0 = 2e6;
cycles = 3;
source.p = toneBurst(1/kgrid.dt, f0, cycles);

% run simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
    'PMLInside', false, ...
    'PMLSize', pml_size, ...
    'PlotSim', plot_simulations, ...
    'DataCast', 'single', ...
    'DisplayMask', 'off');

% assign input and output planes
input_plane_ref = squeeze(sensor_data(1).p);
lateral_plane_ref = squeeze(sensor_data(2).p);

% =========================================================================
% CALCULATE EQUIVALENT SOURCE
% =========================================================================

% calculate equivalent source
num_iterations = 10;
grid_expansion = 0;
[source_estimate, output] = calculateMassSource(input_plane_ref, ...
    kgrid.dx, kgrid.dt, c0, input_plane_pos - 1, grid_expansion, ...
    'NumSteps', num_iterations, ...
    'Plot', plot_simulations);

% =========================================================================
% PROJECT EQUIVALENT SOURCE
% =========================================================================

% assign source plane
source.p_mask = zeros(Nx, Ny, Nz);
source.p_mask(:, :, 1) = 1;
source.p = reshape(source_estimate, [], size(source_estimate, 3));

% run simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
    'PMLInside', false', ...
    'PMLSize', pml_size, ...
    'PlotSim', plot_simulations, ...
    'DataCast', 'single', ...
    'DisplayMask', 'off');

% assign input and output planes
input_plane_eq = squeeze(sensor_data(1).p);
lateral_plane_eq = squeeze(sensor_data(2).p);

% =========================================================================
% ERROR
% =========================================================================

% trim sources to account for propagation time between places
Nt_trim = round(input_plane_pos * kgrid.dx / (kgrid.dt * c0));
input_plane_eq      = input_plane_eq(:, :, Nt_trim + 1:end);
lateral_plane_eq    = lateral_plane_eq(:, :, Nt_trim + 1:end);
input_plane_ref     = input_plane_ref(:, :, 1:size(input_plane_eq, 3));
lateral_plane_ref   = lateral_plane_ref(:, :, 1:size(lateral_plane_eq, 3));

% compute error in the projection
err_proj = 100 * max(abs(abs(lateral_plane_ref(:)) - abs(lateral_plane_eq(:)))) / max(abs(lateral_plane_ref(:)));

% check for test pass
if (err_proj > comparison_threshold_proj) || (output.l2_error(end) > comparison_threshold_source_plane)
    test_pass = false;
end

% =========================================================================
% PLOT
% =========================================================================

if plot_comparisons
   
    % plot the measured and modelled time series
    figure;
    subplot(3, 1, 1);
    imagesc(reshape(input_plane_ref, [], size(input_plane_eq, 3)))
    colorbar;
    title('Reference');
    
    subplot(3, 1, 2);
    imagesc(reshape(input_plane_eq, [], size(input_plane_eq, 3)))
    colorbar;
    title('Equivalent Source Projection');
    
    subplot(3, 1, 3);
    imagesc(100 * abs(...
        reshape(input_plane_ref, [], size(input_plane_eq, 3)) - ...
        reshape(input_plane_eq, [], size(input_plane_eq, 3)) ...
        ) ./ max(input_plane_ref(:)));
    colorbar;
    title('Error (%)');
    
    % plot the projected time series
    figure;
    subplot(3, 1, 1);
    imagesc(reshape(lateral_plane_ref, [], size(lateral_plane_eq, 3)))
    colorbar;
    title('Reference');
    
    subplot(3, 1, 2);
    imagesc(reshape(lateral_plane_eq, [], size(lateral_plane_eq, 3)))
    colorbar;
    title('Equivalent Source Projection');
    
    subplot(3, 1, 3);
    imagesc(100 * abs(...
        reshape(lateral_plane_ref, [], size(lateral_plane_eq, 3)) - ...
        reshape(lateral_plane_eq, [], size(lateral_plane_eq, 3)) ...
        ) ./ max(lateral_plane_ref(:)));
    colorbar;
    title('Error (%)');
    
    % plot the projected plane
    figure;
    subplot(3, 1, 1);
    imagesc(max(lateral_plane_ref, [], 3));
    axis image;
    colorbar;
    title('Reference');
    
    subplot(3, 1, 2);
    imagesc(max(lateral_plane_eq, [], 3));
    axis image;
    colorbar;
    title('Equivalent Source Projection');
    
    subplot(3, 1, 3);
    imagesc(100 * max(abs(lateral_plane_ref - lateral_plane_eq), [], 3) ./ max(abs(lateral_plane_ref(:))));
    axis image;
    colorbar;
    title('Error (%)');
    
end