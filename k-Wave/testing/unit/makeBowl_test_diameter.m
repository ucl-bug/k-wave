function test_pass = makeBowl_test_diameter(plot_comparisons, ~)
%MAKEBOWL_TEST_DIAMETER Check bowl diameter.
%
% DESCRIPTION:
%     makeBowl_test_diameter checks whether the bowl returned by makeBowl
%     has the correct diameter.
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
bowl_pos        = [1, 32, 32];
radius          = 20;
diams           = [9, 15, 21, 35];
focus_pos       = grid_size ./ 2;
plot_bowl       = plot_comparisons;
binary          = false;

% run test for both discs and bowls
for loop_index = 1:2
    
    % set radius to inf second time around
    if loop_index == 2
        radius = inf;
    end

    % loop through diameters
    for index = 1:length(diams)

        % set diameter
        diameter = diams(index);

        % create bowl
        bowl = makeBowl(grid_size, bowl_pos, radius, diameter, focus_pos, 'Plot', plot_bowl, 'Binary', binary);

        % take 2D projection
        proj = squeeze(max(bowl, [], 1));

        % get bounding box
        [x, y] = ind2sub(size(proj), find(proj == 1));

        % check if it is the same as the diameter in each direction
        if ((max(x) - min(x) + 1) ~= diameter) || ((max(y) - min(y) + 1) ~= diameter)
            test_pass = false;
        end

        % plot if required
        if plot_comparisons
            figure;
            imagesc(proj(min(x):max(x), min(y):max(y)));
            axis image;
            title(['Diameter input = ' num2str(diameter)]);
        end

    end
    
end