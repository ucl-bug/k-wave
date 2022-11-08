function test_pass = makeMultiArc_test_creation(plot_comparisons, ~)
%MAKEMULTIARC_TEST_CREATION Check multiple arcs are created successfully.
%
% DESCRIPTION:
%     makeMultiArc_test_creation checks whether the created arcs and labels
%     match those created manually.
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
grid_size           = [Nx, Nx];

% create a Cartesian circle with the x, y positions of the arcs
circle_radius       = 125e-3;
num_arcs            = 16;
arc_pos             = makeCartCircle(circle_radius, num_arcs, [1, 1] * x_size / 2).';
 
% convert the Cartesian arc positions to grid points
arc_pos             = round(arc_pos/dx);

% =========================================================================
% SINGLE DIAMETER, RADIUS, FOCUS_POS
% =========================================================================

% define element parameters
radius    = round(x_size / (2 * dx));
diameter  = 21;
focus_pos = [1, 1] * Nx/2;

% create arcs
[arc_bin, arc_lab] = makeMultiArc(grid_size, arc_pos, radius, diameter, focus_pos);

% allocate comparison matrices
arc_bin_check = zeros(Nx, Nx);
arc_lab_check = zeros(Nx, Nx);

% create arcs manually
for index = 1:size(arc_pos, 1)
    
    % create arc
    arc = makeArc(grid_size, arc_pos(index, :), radius, diameter, focus_pos);
    
    % add to matrices
    arc_bin_check = arc_bin_check + arc;
    arc_lab_check(arc == 1) = index;
    
end

% compare
if sum(arc_bin(:) - arc_bin_check(:)) || sum(arc_lab(:) - arc_lab_check(:))
    test_pass = false;
end

% plot
if plot_comparisons
   
    figure;
    subplot(3, 1, 1);
    imagesc(arc_bin);
    axis image;
    colorbar;
    
    subplot(3, 1, 2);
    imagesc(arc_bin_check);
    axis image;
    colorbar;
    
    subplot(3, 1, 3);
    imagesc(abs(arc_bin - arc_bin_check));
    axis image;
    colorbar;
    
end

% =========================================================================
% MIULTIPLE DIAMETER, RADIUS, FOCUS_POS
% =========================================================================

% define element parameters
radius    = (1:num_arcs) + 20;
diameter  = (1:num_arcs)*2 + 1;
focus_pos = [Nx/2 + (1:num_arcs); Nx/2 + (1:num_arcs)].';

% create arcs
[arc_bin, arc_lab] = makeMultiArc(grid_size, arc_pos, radius, diameter, focus_pos);

% allocate comparison matrices
arc_bin_check = zeros(Nx, Nx);
arc_lab_check = zeros(Nx, Nx);

% create arcs manually
for index = 1:size(arc_pos, 1)
    
    % create arc
    arc = makeArc(grid_size, arc_pos(index, :), radius(index), diameter(index), focus_pos(index, :));
    
    % add to matrices
    arc_bin_check = arc_bin_check + arc;
    arc_lab_check(arc == 1) = index;
    
end

% compare
if sum(arc_bin(:) - arc_bin_check(:)) || sum(arc_lab(:) - arc_lab_check(:))
    test_pass = false;
end

% plot
if plot_comparisons
   
    figure;
    subplot(3, 1, 1);
    imagesc(arc_bin);
    axis image;
    colorbar;
    
    subplot(3, 1, 2);
    imagesc(arc_bin_check);
    axis image;
    colorbar;
    
    subplot(3, 1, 3);
    imagesc(abs(arc_bin - arc_bin_check));
    axis image;
    colorbar;
    
end