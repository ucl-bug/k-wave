function bowl = makeCartBowl(bowl_pos, radius, diameter, focus_pos, num_points, plot_bowl)
%MAKECARTBOWL Create evenly distributed Cartesian points covering a bowl.
%
% DESCRIPTION:
%     makeCartBowl creates a 3 x num_points array of the Cartesian
%     coordinates of points evenly distributed over the surface of a bowl.
%     The bowl points are calculated using Fermat's spiral. The position of
%     the bowl is set by bowl_pos, which corresponds to the center of the
%     rear bowl surface. The orientation of the bowl is set by focus_pos,
%     which corresponds to any point on the axis of the bowl (note, this
%     must not be equal to bowl_pos). It is assumed that the solid angle of
%     the bowl is equal to or less than 2*pi steradians. If the radius is
%     set to inf, a disc is generated.
%
% USAGE:
%     bowl = makeCartBowl(bowl_pos, radius, diameter, focus_pos, num_points)
%     bowl = makeCartBowl(bowl_pos, radius, diameter, focus_pos, num_points, plot_bowl)
%
% INPUTS:
%     bowl_pos    - Cartesian position of the centre of the rear surface of
%                   the bowl given as a three element vector [bx, by, bz]
%                   [m].
%     radius      - Radius of curvature of the bowl [m].
%     diameter    - Diameter of the opening of the bowl [m].
%     focus_pos   - Any point on the beam axis of the bowl given as a three
%                   element vector [fx, fy, fz] [m].
%     num_points  - Number of points on the bowl.
%
% OPTIONAL INPUTS
%     plot_bowl   - Boolean controlling whether the Cartesian points are
%                   plotted (default = false).
%
% OUTPUTS:
%     points      - 3 x num_points array of Cartesian coordinates.
%
% ABOUT:
%     author      - Elliott Wise and Bradley Treeby
%     date        - 5th December 2016
%     last update - 4th February 2018
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2016-2018 Elliott Wise and Bradley Treeby
%
% See also makeBowl, makeCartArc, makeCartSphere, makeCartSphericalSection

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
%
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
% more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

% define literals (ref: http://www.wolframalpha.com/input/?i=golden+angle)
GOLDEN_ANGLE = 2.39996322972865332223155550663361385312499901105811504;

% check for plot_bowl input
if nargin < 6
    plot_bowl = false;
end

% check input values
if radius <= 0
    error('The radius must be positive.');
end
if diameter <= 0
    error('The diameter must be positive.');
end
if diameter > 2*radius
    error('The diameter of the bowl must be equal or less than twice the radius of curvature.');
end
if all(bowl_pos == focus_pos)
    error('The focus_pos must be different to the bowl_pos.');
end

% check for infinite radius of curvature, and call makeCartDisc instead
if isinf(radius)
    bowl = makeCartDisc(bowl_pos, diameter / 2, focus_pos, num_points, plot_bowl);
    return
end

% compute arc angle from chord (ref:
% https://en.wikipedia.org/wiki/Chord_(geometry))
varphi_max = asin(diameter ./ (2*radius));

% compute spiral parameters
theta = @(t) GOLDEN_ANGLE .* t;
C = 2*pi .* (1 - cos(varphi_max)) ./ (num_points - 1);
varphi = @(t) acos(1 - C .* t ./ (2*pi));

% compute canonical spiral points
t = linspace(0, num_points - 1, num_points);
p0 = [cos(theta(t)) .* sin(varphi(t));...
      sin(theta(t)) .* sin(varphi(t));...
                       cos(varphi(t))];
p0 = radius * p0;

% linearly transform the canonical spiral points to give bowl in correct
% orientation
[R, b] = computeLinearTransform(bowl_pos, focus_pos, radius);
bowl = bsxfun(@plus, R * p0, b);

% plot results
if plot_bowl
    
    % select suitable axis scaling factor
    [~, scale, prefix] = scaleSI(max(bowl(:)));
    
    % create the figure
    figure;
    plot3(bowl(1, :) * scale, bowl(2, :) * scale, bowl(3, :) * scale, '.');
    xlabel(['[' prefix 'm]']);
    ylabel(['[' prefix 'm]']);
    zlabel(['[' prefix 'm]']);
    axis equal;
    grid on;
    box on;
    
end