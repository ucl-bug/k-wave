function test_pass = makeArc_test_arc_pos(plot_comparisons, ~)
%MAKEARC_TEST_ARC_POS Check arc position lies on the arc.
%
% DESCRIPTION:
%     makeArc_test_arc_pos checks whether the specified arc position (arc
%     midpoint) actually lies on the created arc.
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
%     date             - 4th February 2017
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
arc_pos_array   = [1, 32; 3, 27; 9, 40; 10, 43];
radius          = 20;
diameter        = 9;
focus_pos       = grid_size ./ 2;
plot_arc        = false;

% run test for both lines and arcs
for loop_index = 1:2
    
    % set radius to inf second time around
    if loop_index == 2
        radius = inf;
    end

    % loop through arc positions
    for index = 1:length(arc_pos_array)

        % set arc position
        arc_pos = arc_pos_array(index, :);

        % create arc
        arc = makeArc(grid_size, arc_pos, radius, diameter, focus_pos, 'Plot', plot_arc);
                
        % check arc_pos is part of the arc
        if arc(arc_pos(1), arc_pos(2)) ~= 1
            test_pass = false;
        end
        
        % plot 
        if plot_comparisons
            arc(arc_pos(1), arc_pos(2)) = 2;
            figure;
            imagesc(arc);
            axis image;
            colormap(flipud(gray));
        end
        
    end
    
end