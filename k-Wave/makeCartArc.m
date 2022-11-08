function arc = makeCartArc(arc_pos, radius, diameter, focus_pos, num_points, plot_arc)
%MAKECARTARC Create evenly distributed Cartesian points covering an arc.
%
% DESCRIPTION:
%     makeCartArc creates a 2 x num_points array of the Cartesian
%     coordinates of points evenly distributed over an arc. The midpoint of
%     the arc is set by arc_pos. The orientation of the arc is set by
%     focus_pos, which corresponds to any point on the axis of the arc
%     (note, this must not be equal to arc_pos). It is assumed that the arc
%     angle is equal to or less than pi radians. If the radius is set to
%     inf, a line is generated.
%     
%     Note, the first and last points do not lie on the end-points of the
%     underlying canonical arc, but are offset by half the angular spacing
%     between the points.
%
% USAGE:
%     arc = makeCartArc(arc_pos, radius, diameter, focus_pos, num_points, plot_arc)
%
% INPUTS:
%     arc_pos     - Centre of the rear surface (midpoint) of the arc given
%                   as a two element vector [bx, by] [m].
%     radius      - Radius of curvature of the arc [m].
%     diameter    - Diameter of the opening of the arc (length of line
%                   connecting arc endpoints) [m].
%     focus_pos   - Any point on the beam axis of the arc given as a two
%                   element vector [fx, fy] [m].
%     num_points  - Number of points on the arc.
%
% OPTIONAL INPUTS
%     plot_arc    - Boolean controlling whether the Cartesian points are
%                   plotted (default = false).
%
% OUTPUTS:
%     points      - 2 x num_points array of Cartesian coordinates.
%
% ABOUT:
%     author      - Elliott Wise and Bradley Treeby
%     date        - 15th May 2017
%     last update - 26th June 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017-2019 Elliott Wise and Bradley Treeby
%
% See also makeArc, makeCartDisc

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

% check for plot_arc input
if nargin < 6
    plot_arc = false;
end

% check input values
if radius <= 0
    error('The radius must be positive.');
end
if diameter <= 0
    error('The diameter must be positive.');
end
if diameter > 2*radius
    error('The diameter of the arc must be less than twice the radius of curvature.');
end
if all(arc_pos == focus_pos)
    error('The focus_pos must be different to the arc_pos.');
end

% check for infinite radius of curvature, and make a large finite number
if isinf(radius)
    radius = 1e10 * diameter;
end

% compute arc angle from chord
varphi_max = asin(diameter ./ (2 * radius));

% angle between points
dvarphi = 2 * varphi_max ./ num_points;

% compute canonical arc points where arc is centered on the origin, and its
% back is placed at a distance of radius along the positive y-axis
t = linspace(-varphi_max + dvarphi/2, varphi_max - dvarphi/2, num_points);  
p0 = radius * [sin(t); cos(t)];

% linearly transform canonical points to give arc in correct orientation
[R, b] = computeLinearTransform2D(arc_pos, radius, focus_pos);
arc = bsxfun(@plus, R * p0, b);

% plot results
if plot_arc
    
    % select suitable axis scaling factor
    [~, scale, prefix] = scaleSI(max(abs(arc(:)))); 
    
    % create the figure
    figure;
    plot(arc(2, :) .* scale, arc(1, :) .* scale, 'b.');
    set(gca, 'YDir', 'reverse');
    xlabel(['y-position [' prefix 'm]']);
    ylabel(['x-position [' prefix 'm]']);
    axis equal;
    
end

end

function [R, b] = computeLinearTransform2D(arc_pos, radius, focus_pos)
%COMPUTELINEARTRANSFORM Compute a linear transformation.
%
% DESCRIPTION:
%     computeLinearTransform calculates a rotation matrix to tranform the
%     computed arc points to the orientation specified by the arc and focus
%     positions.
%     
% ABOUT:
%     author      - Elliott Wise
%     date        - 15th May 2017
%     last update - 28th January 2018

% vector pointing from arc_pos to focus_pos
beam_vec = focus_pos - arc_pos;

% normalise to give unit beam vector
beam_vec = beam_vec ./ norm(beam_vec);

% canonical normalised beam_vec (canonical arc_pos is [0 1])
beam_vec0 = [0, -1];

% find the angle between canonical and specified beam_vec
theta = atan2(beam_vec(2), beam_vec(1)) - atan2(beam_vec0(2), beam_vec0(1));

% convert to a rotation matrix
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

% compute an offset for the arc, where arc_centre = move from arc_pos
% towards focus by radius 
b = arc_pos.' + radius * beam_vec.';

end