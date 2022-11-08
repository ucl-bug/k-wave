function [arcs, arcs_labelled] = makeMultiArc(grid_size, arc_pos, radius, diameter, focus_pos, varargin)
%MAKEMULTIARC Create a binary map of multiple arcs within a 2D grid.
%
% DESCRIPTION:
%     makeMultiArc creates a binary map of multiple arcs within a
%     two-dimensional grid using makeArc. The position of the arcs is
%     denoted by 1's in the matrix with 0's elsewhere. A labelled matrix
%     can also be returned, where the position of the first arc is denoted
%     by 1's, the position of the second arc by 2's, and so on.
%
% USAGE:
%     [arcs, arcs_labelled] = makeMultiArc(grid_size, arc_pos, radius, diameter, focus_pos)
%     [arcs, arcs_labelled] = makeMultiArc(grid_size, arc_pos, radius, diameter, focus_pos, ...)
%
% INPUTS:
%     grid_size       - size of the 2D grid given as a two element vector
%                       [Nx, Ny] [grid points]
%     arc_pos         - midpoint of each arc given as a matrix of
%                       dimensions N x 2, with each row specifying the
%                       midpoint of each arc as a two element vector 
%                       [ax, ay] [grid points]
%     radius          - radius of curvature of each arc given as either
%                       a single number (if the arcs have the same radius),
%                       or an N-element vector containing the radius for
%                       each arc [grid points] 
%     diameter        - diameter of each arc (length of straight line
%                       between the end points) given as either a single 
%                       number (if the arcs have the same diameter), or an
%                       N-element vector containing the diameter for each
%                       arc [grid points] 
%     focus_pos       - any point on the beam axis of the arc given as
%                       either a two element vector (if the arcs have the
%                       same focus_pos), or as a matrix of dimensions 
%                       N x 2, with each row specifying the focus_pos for
%                       each arc given as a two element vector
%                       [fx, fy] [grid points]
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%
%     'Plot'          - Boolean controlling whether the arc array is
%                        plotted using imagesc (default = false). 
%
% OUTPUTS:
%     arcs            - 2D binary map of arcs
%     arcs_labelled   - 2D labelled matrix of arcs
%
% ABOUT:
%     author          - Yan To Ling and Bradley Treeby
%     date            - 19th January 2015
%     last update     - 13th January 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2015-2019 Yan To Ling and Bradley Treeby
%
% See also makeArc, makeMultiBowl

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

% =========================================================================
% DEFINE LITERALS
% =========================================================================

% optional input defaults
NUM_REQUIRED_INPUTS = 5;
PLOT_ARC_DEF        = false;

% =========================================================================
% INPUT CHECKING
% =========================================================================

% assign optional inputs with defaults
plot_arc = PLOT_ARC_DEF;

% replace with user defined values if provided
if nargin < NUM_REQUIRED_INPUTS
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}                
            case 'Plot'
                
                % assign input
                plot_arc = varargin{input_index + 1};
                
                % check value
                if ~islogical(plot_arc)
                    error('Optional input ''Plot'' must be Boolean.');
                end
          
            otherwise
                error('Unknown optional input.');
        end
    end
end

% check for plot_arc input
if nargin < 6 || isempty(plot_arc)
    plot_arc = false;
end

% check inputs
if size(arc_pos, 2) ~= 2
    error('arc_pos should contain 2 columns, with [ax, ay] in each row.');
end
if length(radius) ~= 1 && length(radius) ~= size(arc_pos,1)
    error('The number of rows in arc_pos and radius does not match.');
end
if length(diameter) ~= 1 && length(diameter) ~= size(arc_pos,1)
    error('The number of rows in arc_pos and diameter does not match.');
end

% force integer grid size values
grid_size = round(grid_size);
arc_pos   = round(arc_pos);
focus_pos = round(focus_pos);
diameter  = round(diameter);
radius    = round(radius);

% =========================================================================
% CREATE ARCS
% =========================================================================

% create empty matrix
arcs = zeros(grid_size);
if nargout == 2
    arcs_labelled = zeros(grid_size);
end

% loop for calling makeArc
for k = 1:size(arc_pos,1)
    
    % get parameters for current arc
    if size(arc_pos, 1) > 1
        arc_pos_k = arc_pos(k, :);
    else
        arc_pos_k = arc_pos;
    end
    if length(radius) > 1
        radius_k = radius(k);
    else
        radius_k = radius;
    end
    if length(diameter) > 1
        diameter_k = diameter(k);
    else
        diameter_k = diameter;
    end      
    if size(focus_pos, 1) > 1
        focus_pos_k = focus_pos(k, :);
    else
        focus_pos_k = focus_pos;
    end      
    
    % create new arc
    new_arc = makeArc(grid_size, arc_pos_k, radius_k, diameter_k, focus_pos_k);
    
    % add arc to arc matrix
    arcs = arcs + new_arc;    
    
    % add new arc to labelling matrix
    if nargout == 2
        arcs_labelled(new_arc == 1) = k;
    end
    
end

% check if any of the arcs are overlapping
if maxND(arcs) > 1
    
    % display warning
    disp(['WARNING: ' num2str(maxND(arcs) - 1) ' arcs are overlapping']);
    
    % force the output to be binary
    arcs(arcs ~= 0) = 1;
    
end

% create the figure
if plot_arc
    figure;
    imagesc(arcs, [-1 1]);
    colormap(getColorMap);
    axis image;
    xlabel('y-position [grid points]');
    ylabel('x-position [grid points]');
end