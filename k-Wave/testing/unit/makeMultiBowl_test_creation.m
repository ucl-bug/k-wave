function test_pass = makeMultiBowl_test_creation(plot_comparisons, ~)
%MAKEMULTIBOWL_TEST_CREATION Check multiple bowls are created successfully.
%
% DESCRIPTION:
%     makeMultiBowl_test_creation checks whether the created bowls and
%     labels match those created manually.
%
% INPUTS:
%     plot_comparisons - Boolean controlling whether any comparisons (e.g.,
%                        error norms) are plotted 
%
% OUTPUTS:
%     test_pass        - Boolean denoting whether the test passed
%
% ABOUT:
%     author           - Bradley Treeby
%     date             - 12th April 2017
%     last update      - 12th April 2017
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

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
end

% set pass variable
test_pass = true;

% define grid parameters
x_size              = 300e-3;
Nx                  = 256;
dx                  = x_size / Nx;
grid_size           = [Nx, Nx, Nx];

% create a Cartesian sphere with the x, y, z positions of the bowls
sphere_radius       = 125e-3;
num_bowls           = 32;
bowl_pos            = makeCartSphere(sphere_radius, num_bowls, [1, 1, 1] * x_size / 2).';

% convert the Cartesian bowl positions to grid points
bowl_pos            = round(bowl_pos/dx);

% =========================================================================
% SINGLE DIAMETER, RADIUS, FOCUS_POS
% =========================================================================

% define element parameters
radius              = round(x_size / (2 * dx));
diameter            = 21;
focus_pos           = [1, 1, 1] * Nx/2;

% create bowls
[bowl_bin, bowl_lab] = makeMultiBowl(grid_size, bowl_pos, radius, diameter, focus_pos);

% allocate comparison matrices
bowl_bin_check = zeros(Nx, Nx, Nx);
bowl_lab_check = zeros(Nx, Nx, Nx);

% create arcs manually
for index = 1:size(bowl_pos, 1)
    
    % create arc
    bowl = makeBowl(grid_size, bowl_pos(index, :), radius, diameter, focus_pos);
    
    % add to matrices
    bowl_bin_check = bowl_bin_check + bowl;
    bowl_lab_check(bowl == 1) = index;
    
end

% compare
if sum(bowl_bin(:) - bowl_bin_check(:)) || sum(bowl_lab(:) - bowl_lab_check(:))
    test_pass = false;
end

% plot
if plot_comparisons
   
    figure;
    subplot(3, 1, 1);
    imagesc(bowl_bin(:, :, end/2));
    axis image;
    colorbar;
    
    subplot(3, 1, 2);
    imagesc(bowl_bin_check(:, :, end/2));
    axis image;
    colorbar;
    
    subplot(3, 1, 3);
    imagesc(abs(bowl_bin(:, :, end/2) - bowl_bin_check(:, :, end/2)));
    axis image;
    colorbar;
    
end

% =========================================================================
% MIULTIPLE DIAMETER, RADIUS, FOCUS_POS
% =========================================================================

% define element parameters
radius    = (1:num_bowls) + 20;
diameter  = (1:num_bowls)*2 + 1;
focus_pos = [Nx/2 + (1:num_bowls); Nx/2 + (1:num_bowls); Nx/2 + (1:num_bowls)].';

% create bowls
[bowl_bin, bowl_lab] = makeMultiBowl(grid_size, bowl_pos, radius, diameter, focus_pos);

% allocate comparison matrices
bowl_bin_check = zeros(Nx, Nx, Nx);
bowl_lab_check = zeros(Nx, Nx, Nx);

% create arcs manually
for index = 1:size(bowl_pos, 1)
    
    % create arc
    bowl = makeBowl(grid_size, bowl_pos(index, :), radius(index), diameter(index), focus_pos(index, :));
    
    % add to matrices
    bowl_bin_check = bowl_bin_check + bowl;
    bowl_lab_check(bowl == 1) = index;
    
end

% compare
if sum(bowl_bin(:) - bowl_bin_check(:)) || sum(bowl_lab(:) - bowl_lab_check(:))
    test_pass = false;
end

% plot
if plot_comparisons
   
    figure;
    subplot(3, 1, 1);
    imagesc(bowl_bin(:, :, end/2));
    axis image;
    colorbar;
    
    subplot(3, 1, 2);
    imagesc(bowl_bin_check(:, :, end/2));
    axis image;
    colorbar;
    
    subplot(3, 1, 3);
    imagesc(abs(bowl_bin(:, :, end/2) - bowl_bin_check(:, :, end/2)));
    axis image;
    colorbar;
    
end