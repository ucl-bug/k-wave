function test_pass = calculateMassSourceCW_disc_source(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to check the calculation of a mass source using
%     calculateMassSourceCW.
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
comparison_threshold_proj = 0.05;
comparison_threshold_source_plane = 1;

% set pass variable
test_pass = true;

% =========================================================================
% DEFINE LITERALS
% =========================================================================

% transducer properties 
f0 = 1.1e6;        % [Hz]
c0 = 1500;         % [m/s]
source_radius = 15;

% position of input plane on grid [grid points]
input_plane_pos = 16;

% grid properties 
Nx = 64;
Ny = 64;
Nz = 128;
dx = 0.5e-3;

% =========================================================================
% PROJECT TRUE SOURCE
% =========================================================================

% create source plane in the shape of a disc
disc = makeDisc(Nx, Ny, Nx/2, Ny/2, source_radius);

% assign input plane
amp_in            = zeros(Nx, Ny, Nz);
amp_in(:, :, 1)   = disc;
phase_in          = 0;

% compute projected field
projected_field_ref = acousticFieldPropagator(amp_in, phase_in, dx, f0, c0);

% assign input and output planes
input_plane = projected_field_ref(:, :, input_plane_pos);
lateral_plane_ref = squeeze(projected_field_ref(:, Ny/2, 1 + input_plane_pos:end));

% =========================================================================
% CALCULATE EQUIVALENT SOURCE
% =========================================================================

% calculate equivalent source
num_iterations = 10;
grid_expansion = 0;
[source_estimate, output] = calculateMassSourceCW(input_plane, dx, f0, c0, input_plane_pos - 1, grid_expansion, ...
    'NumSteps', num_iterations, ...
    'Plot', plot_simulations);

% =========================================================================
% PROJECT EQUIVALENT SOURCE
% =========================================================================

% assign input plane
amp_in            = zeros(Nx, Ny, Nz);
amp_in(:, :, 1)   = abs(source_estimate);
phase_in          = zeros(Nx, Ny, Nz);
phase_in(:, :, 1) = angle(source_estimate);

% compute projected field
projected_field_eq = acousticFieldPropagator(amp_in, phase_in, dx, f0, c0);

% assign output plane
lateral_plane_eq = squeeze(projected_field_eq(:, Ny/2, 1 + input_plane_pos:end));

% =========================================================================
% ERROR
% =========================================================================

% compute error in the projection
err_proj = max(abs(abs(lateral_plane_ref(:)) - abs(lateral_plane_eq(:)))) / max(abs(lateral_plane_ref(:)));

% check for test pass
if (err_proj > comparison_threshold_proj) || (output.l2_error(end) > comparison_threshold_source_plane)
    test_pass = false;
end

% =========================================================================
% PLOT
% =========================================================================

if plot_comparisons
   
    figure;
    
    subplot(3, 2, 1);
    imagesc(abs(lateral_plane_ref));
    axis image;
    colorbar;
    title('Reference');
    
    subplot(3, 2, 2);
    imagesc(angle(lateral_plane_ref));
    axis image;
    colorbar;
    
    subplot(3, 2, 3);
    imagesc(abs(lateral_plane_eq));
    axis image;
    colorbar;
    title('Equivalent Source Projection');
    
    subplot(3, 2, 4);
    imagesc(angle(lateral_plane_eq));
    axis image;
    colorbar;
    
    subplot(3, 2, 5);
    imagesc(abs(abs(lateral_plane_ref) - abs(lateral_plane_eq)) ./ max(abs(lateral_plane_ref(:))));
    axis image;
    colorbar;
    title('Error');
    
    subplot(3, 2, 6);
    imagesc(abs(angle(lateral_plane_ref) - angle(lateral_plane_eq)) ./ (2*pi));
    axis image;
    colorbar;
    
end