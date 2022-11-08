function test_pass = makeArc_test_rotation(plot_comparisons, ~)
%MAKEARC_TEST_ROTATION Check rotated arc.
%
% DESCRIPTION:
%     makeArc_test_diameter creates a series of arcs with radius = inf at
%     different angles and checks whether the arcs are all contained within
%     a disc of the same diameter.
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
%     date             - 3rd February 2017
%     last update      - 9th April 2017
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

% specify the arc parameters
grid_size       = [64, 64];
arc_pos         = [32, 32];
radius          = 20;
diameter        = 21;
plot_arc        = false;

% create list of focus positions to rotate line
focus_pos_array = [1, 1; 8, 1; 16, 1; 24, 1; 32, 1; 40, 1; 48, 1; 56, 1; 64, 1; 64, 8; 64, 16; 64, 24; 64, 32; 64, 40; 64, 48; 64, 56];

% run test for both discs and arcs
for loop_index = 1:2
    
    % set radius to inf second time around
    if loop_index == 2
        radius = inf;
    end
    
    % create empty matrix
    arc_sum = zeros(grid_size);

    % loop through positions
    figure;
    for focus_index = 1:length(focus_pos_array)

        % get current focus position
        focus_pos = focus_pos_array(focus_index, :);

        % create arc
        arc = makeArc(grid_size, arc_pos, radius, diameter, focus_pos, 'Plot', plot_arc);

        % add to total arc
        arc_sum = double(arc_sum | arc);

        if plot_comparisons
            imagesc(arc_sum);
            axis image;
            colormap(flipud(gray));
            drawnow;
            pause(0.1);
        end

    end

    % get bounding box
    [x, y] = ind2sub(size(arc_sum), find(arc_sum == 1));

    % check if it is the same as the diameter in each direction
    if ((max(x) - min(x) + 1) ~= diameter) || ((max(y) - min(y) + 1) ~= diameter)
        test_pass = false;
    end

    % create a disc with the same radius
    disc = makeDisc(grid_size(1), grid_size(2), arc_pos(1), arc_pos(2), ceil(diameter/2));

    % check all points are within disc
    if any( (disc == 0) & (arc_sum == 1))
        test_pass = false;
    end

    if plot_comparisons
        figure;
        imagesc(disc + arc_sum);
        axis image;
        colormap(flipud(gray));
    end
    
end
    