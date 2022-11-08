function circle = makeCartCircle(radius, num_points, center_pos, arc_angle, plot_circle)
%MAKECARTCIRCLE Create a 2D Cartesian circle or arc.
%
% DESCRIPTION:
%     MakeCartCircle creates a 2 x num_points array of the Cartesian
%     coordinates of points evenly distributed over a circle or arc (if
%     arc_angle is given).
%
% USAGE:
%     circle = makeCartCircle(radius, num_points)
%     circle = makeCartCircle(radius, num_points, center_pos)
%     circle = makeCartCircle(radius, num_points, center_pos, arc_angle)
%     circle = makeCartCircle(radius, num_points, center_pos, arc_angle, plot_circle)
%
% INPUTS:
%     radius          - circle radius [m]
%     num_points      - number of points in the circle
%
% OPTIONAL INPUTS:
%     center_pos      - [x, y] position of the circle center [m] 
%                       (default = [0, 0])
%     arc_angle       - arc angle for incomplete circle [radians]
%                       (default = 2*pi)
%     plot_circle     - Boolean controlling whether the Cartesian points
%                       are plotted (default = false)
%
% OUTPUTS:
%     circle          - 2 x num_points array of Cartesian coordinates
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 5th June 2009
%     last update     - 8th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby
%
% See also cart2grid, makeCartSphere, makeCircle

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

% check for plot_circle input
if nargin < 5 || isempty(plot_circle)
    plot_circle = false;
end

% check for arc_angle input
if nargin < 4 || isempty(arc_angle)
    arc_angle = 2*pi;
    full_circle = true;
elseif arc_angle == 2*pi
    full_circle = true;
else
    full_circle = false;
end

% check for center_pos input
if nargin < 3 || isempty(center_pos)
    cx = 0;
    cy = 0;
else
    cx = center_pos(1);
    cy = center_pos(2);
end

% ensure there is only a total of num_points including the endpoints when
% arc_angle is not equal to 2*pi
if ~full_circle
    num_points = num_points - 1;
end

% create angles
angles = (0:(num_points)) * arc_angle / num_points + pi/2;

% discard repeated final point if arc_angle is equal to 2*pi
if full_circle
    angles = angles(1:end - 1);
end

% create cartesian grid
% circle = flipud([radius*cos(angles); radius*sin(-angles)]);   % B.0.3
circle = ([radius*cos(angles); radius*sin(-angles)]);           % B.0.4

% offset if needed
circle(1, :) = circle(1, :) + cx;
circle(2, :) = circle(2, :) + cy;

% plot results
if plot_circle
    
    % select suitable axis scaling factor
    [~, scale, prefix] = scaleSI(max(abs(circle(:)))); 
    
    % create the figure
    figure;
    plot(circle(2, :) .* scale, circle(1, :) .* scale, 'b.');
    set(gca, 'YDir', 'reverse');
    xlabel(['y-position [' prefix 'm]']);
    ylabel(['x-position [' prefix 'm]']);
    axis equal;
    
end