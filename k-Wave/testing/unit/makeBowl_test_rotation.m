function test_pass = makeBowl_test_rotation(plot_comparisons, ~)
%MAKEBOWL_TEST_ROTATION Check rotated bowl.
%
% DESCRIPTION:
%     makeBowl_test_diameter creates a series of bowls with radius = inf at
%     different angles and checks whether the bowls are all contained
%     within a ball of the same diameter.
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

% specify the bowl parameters
grid_size       = [64, 64, 64];
bowl_pos        = [32, 32, 32];
radius          = 20;
diameter        = 21;
plot_bowl       = plot_comparisons;
binary          = false;

% create list of focus positions to rotate disc
focus_pos_array = [1, 1, 32; 8, 1, 32; 16, 1, 32; 24, 1, 32; 32, 1, 32; 40, 1, 32; 48, 1, 32; 56, 1, 32; 64, 1, 32; ...
                   64, 8, 32; 64, 16, 32; 64, 24, 32; 64, 32, 32; 64, 40, 32; 64, 48, 32; 64, 56, 32];

% run test for both discs and bowls
for loop_index = 1:2
    
    % set radius to inf second time around
    if loop_index == 2
        radius = inf;
    end
               
    % create empty matrix
    bowl_sum = zeros(grid_size);

    % loop through positions
    figure;
    for focus_index = 1:length(focus_pos_array)

        % get current focus position
        focus_pos = focus_pos_array(focus_index, :);

        % create bowl
        bowl = makeBowl(grid_size, bowl_pos, radius, diameter, focus_pos, 'Plot', plot_bowl, 'Binary', binary);

        % add to total bowl
        bowl_sum = double(bowl_sum | bowl);

    end

    % get bounding box
    [x, y, z] = ind2sub(size(bowl_sum), find(bowl_sum == 1));

    % check if it is the same as the diameter in each direction
    if ((max(x) - min(x) + 1) ~= diameter) || ((max(y) - min(y) + 1) ~= diameter) || ((max(z) - min(z) + 1) ~= diameter)
        test_pass = false;
    end

    % create a ball with the same radius
    ball = makeBall(grid_size(1), grid_size(2), grid_size(3), bowl_pos(1), bowl_pos(2), bowl_pos(3), ceil(diameter/2));

    % check all points are within disc
    if any( (ball == 0) & (bowl_sum == 1))
        test_pass = false;
    end

    % plot if required
    if plot_comparisons

        % plot 3D sum
        voxelPlot(bowl_sum);

        % plot plane
        figure;
        imagesc(squeeze(ball(end/2, :, :)) + squeeze(bowl_sum(end/2, :, :)));
        axis image;
        colormap(flipud(gray));

    end
    
end