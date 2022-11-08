function arc = makeArc(grid_size, arc_pos, radius, diameter, focus_pos, varargin)
%MAKEARC Create a binary map of an arc within a 2D grid.
%
% DESCRIPTION:
%     makeArc creates a binary map of an arc within a two-dimensional grid.
%     The arc position is denoted by 1's in the matrix with 0's elsewhere.
%     The arc is found using the mid-point circle algorithm, and is simply
%     connected so that any grid point on the arc will have at most 2
%     neighbours.
%
%     The midpoint of the arc is set by arc_pos. The orientation of the arc
%     is set by focus_pos, which corresponds to any point on the axis of
%     the arc (note, this must not be equal to arc_pos). It is assumed that
%     the arc angle is equal to or less than pi radians. If the radius is
%     set to inf, a line is generated.
%
% USAGE:
%     arc = makeArc(grid_size, arc_pos, radius, diameter, focus_pos)
%     arc = makeArc(grid_size, arc_pos, radius, diameter, focus_pos, ...)
%
% INPUTS:
%     grid_size     - size of the 2D grid given as a two element vector
%                     [Nx, Ny] [grid points]  
%     arc_pos       - midpoint of the arc given as a two element vector
%                     [ax, ay] [grid points]
%     radius        - radius of curvature of the arc [grid points]
%     diameter      - aperture diameter (length of line connecting arc
%                     endpoints) [grid points]  
%     focus_pos     - any point on the beam axis of the arc given as a two
%                     element vector [fx, fy] [grid points]
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%
%     'Plot'        - Boolean controlling whether the arc is plotted using
%                     imagesc (default = false). 
%
% OUTPUTS:
%     arc           - 2D binary map of an arc
%
% ABOUT:
%     author        - Bradley Treeby and Yan To Ling
%     date          - 25th January 2015
%     last update   - 9th April 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2015-2017 Bradley Treeby and Yan To Ling
%
% See also makeBowl, makeCartArc, makeLine, makeMultiArc

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

% force integer input values
grid_size = round(grid_size);
arc_pos   = round(arc_pos);
focus_pos = round(focus_pos);
diameter  = round(diameter);
radius    = round(radius);

% check the input ranges
if any(grid_size < 1)
    error('The grid size must be positive.');
end
if radius <= 0
    error('The radius must be positive.');
end
if diameter <= 0
    error('The diameter must be positive.');
end
if any(arc_pos < 1) || any(arc_pos > grid_size)
    error('The centre of the arc must be within the grid.');
end
if diameter > 2*radius
    error('The diameter of the arc must be less than twice the radius of curvature.');
end
if ~rem(diameter, 2)
    error('The diameter must be an odd number of grid points.');
end
if all(arc_pos == focus_pos)
    error('The focus_pos must be different to the arc_pos.');
end

% assign variable names to vector components
Nx = grid_size(1);
Ny = grid_size(2);
ax = arc_pos(1);
ay = arc_pos(2);
fx = focus_pos(1);
fy = focus_pos(2);

% =========================================================================
% CREATE ARC
% =========================================================================
   
if ~isinf(radius)

    % find half the arc angle
    half_arc_angle = asin(diameter / 2 / radius);

    % find centre of circle on which the arc lies
    distance_cf = sqrt( (ax - fx).^2 + (ay - fy).^2 );
    cx = round(radius ./ distance_cf .* (fx - ax) + ax);
    cy = round(radius ./ distance_cf .* (fy - ay) + ay);
    c = [cx, cy];

    % create circle
    arc = makeCircle(Nx, Ny, cx, cy, radius);
    
    % form vector from the geometric arc centre to the arc midpoint
    v1 = arc_pos - c;
    
    % calculate length of vector
    l1 = sqrt(sum((arc_pos - c).^2));
    
    % extract all points that form part of the arc
    arc_ind = find(arc == 1);    
    
    % loop through the arc points
    for ind = 1:length(arc_ind)
        
        % extract the indices of the current point
        [x_ind, y_ind] = ind2sub([Nx, Ny], arc_ind(ind));
        p = [x_ind, y_ind];

        % form vector from the geometric arc centre to the current point
        v2 = p - c;

        % calculate length of vector
        l2 = sqrt(sum((p - c).^2));

        % find the angle between the two vectors using the dot product,
        % normalised using the vector lengths
        theta = acos(sum( v1 .* v2 ./ (l1 .* l2) ));
        
        % if the angle is greater than the half angle of the arc, remove
        % it from the arc
        if theta > half_arc_angle
            arc(arc_ind(ind)) = 0;
        end
        
    end     
    
else
    
    % calculate arc direction angle, then rotate by 90 degrees
    ang = atan( (fx - ax) ./ (fy - ay) ) + pi/2;
    
    % draw lines to create arc with infinite radius
    arc = double(makeLine(Nx, Ny, arc_pos, ang, (diameter - 1)/2) | ...
        makeLine(Nx, Ny, arc_pos, ang + pi, (diameter - 1)/2));

end

% =========================================================================
% PLOT
% =========================================================================

% create the figure
if plot_arc
    figure;
    imagesc(arc, [-1 1]);
    colormap(getColorMap);
    axis image;
    xlabel('y-position [grid points]');
    ylabel('x-position [grid points]');
end