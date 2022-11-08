function segment = makeCartSphericalSegment(bowl_pos, radius, inner_diameter, outer_diameter, focus_pos, num_points, plot_bowl, num_points_inner)
%MAKECARTSPHERICALSEGMENT Create evenly distributed Cartesian points covering a spherical segment.
%
% DESCRIPTION:
%     makeCartSphericalSegment creates a 3 x num_points array of the
%     Cartesian coordinates of points evenly distributed over the surface
%     of a spherical segment. This function is closely related to
%     makeCartBowl, however, only a segment of the bowl is used (equivalent
%     to an annulus if the bowl was a disc).
% 
%     The points are calculated using Fermat's spiral. The position of the
%     segment is set by bowl_pos, which corresponds to the center of the
%     rear surface of the bowl from which the annulus is derived. The
%     orientation of the annulus is set by focus_pos, which corresponds to
%     any point on the axis of the annulus (note, this must not be equal to
%     bowl_pos).
%
%     It is assumed that the solid angle of the bowl from which the annulus
%     is derived is equal to or less than 2*pi steradians. If the radius is
%     set to inf, an annular disc is generated.
%
% USAGE:
%     segment = makeCartSphericalSegment(bowl_pos, radius, inner_diameter, outer_diameter, focus_pos, num_points)
%     segment = makeCartSphericalSegment(bowl_pos, radius, inner_diameter, outer_diameter, focus_pos, num_points, plot_bowl)
%
% INPUTS:
%     bowl_pos       - Cartesian position of the centre of the rear surface
%                      of the underlying bowl on which the spherical
%                      segment lies given as a three element vector [bx,
%                      by, bz] [m].
%     radius         - Radius of curvature of the underlying bowl [m].
%     inner_diameter - Inner aperture diameter of the spherical segment [m].
%     outer_diameter - Outer aperture diameter of the spherical segment [m].
%     focus_pos      - Any point on the beam axis of the underlying bowl
%                      given as a three element vector [fx, fy, fz] [m].
%     num_points     - Number of points on the spherical segment.
%
% OPTIONAL INPUTS
%     plot_bowl      - Boolean controlling whether the Cartesian points are
%                      plotted (default = false).
%     num_points_inner 
%                    - If constructing an annular array with contiguous
%                      elements (no kerf), the positions of the points will
%                      not exactly match makeCartBowl, as each element has
%                      no knowledge of the number of points on the internal
%                      elements. To force the points to match, specify the
%                      total number of points used on all annular segments
%                      internal to the current one.
%
% OUTPUTS:
%     points      	 - 3 x num_points array of Cartesian coordinates.
%
% ABOUT:
%     author         - Bradley Treeby
%     date           - 30th March 2021
%     last update    - 1st April 2021
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2021 Bradley Treeby
%
% See also makeCartBowl

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
if nargin < 7
    plot_bowl = false;
end

% check input values
if radius <= 0
    error('The radius must be positive.');
end
if inner_diameter < 0
    error('The inner diameter must be positive.');
end
if inner_diameter >= outer_diameter
    error('The inner diameter must be less than the outer diameter.');
end
if outer_diameter <= 0
    error('The outer diameter must be positive.');
end
if outer_diameter > 2*radius
    error('The outer diameter of the bowl must be equal or less than twice the radius of curvature.');
end
if all(bowl_pos == focus_pos)
    error('The focus_pos must be different to the bowl_pos.');
end

% check for infinite radius of curvature
if isinf(radius)
    error('Annular disc (infinite radius) not yet supported.');
end

% compute arc angle from chord (ref:
% https://en.wikipedia.org/wiki/Chord_(geometry))
varphi_min = asin(inner_diameter ./ (2 * radius));
varphi_max = asin(outer_diameter ./ (2 * radius));

% compute spiral parameters over annulus - if the number of points used so
% far in the inside of the annulus is given, use this so the points follow
% exactly the same trajectory, otherwise, just distribution them evenly
theta = @(t) GOLDEN_ANGLE .* t;
if nargin == 8
    C = (1 - cos(varphi_max)) ./ ((num_points + num_points_inner) - 1);
    varphi = @(t) (acos(1 - C .* t));
    t_start = ceil((1 - cos(varphi_min)) / C);
    t = linspace(t_start, num_points_inner + num_points - 1, num_points);
else
    C = (1 - cos(varphi_max)) ./ (num_points - 1);
    varphi = @(t) (acos(1 - C .* t));
    t_start = ceil((1 - cos(varphi_min)) / C);
    t = linspace(t_start, num_points - 1, num_points);
end

% compute canonical spiral points
p0 = [cos(theta(t)) .* sin(varphi(t));...
      sin(theta(t)) .* sin(varphi(t));...
                       cos(varphi(t))];
p0 = radius * p0;

% linearly transform the canonical spiral points to give bowl in correct
% orientation
[R, b] = computeLinearTransform(bowl_pos, focus_pos, radius);
segment = bsxfun(@plus, R * p0, b);

% plot results
if plot_bowl
    
    % select suitable axis scaling factor
    [~, scale, prefix] = scaleSI(max(segment(:)));
    
    % create the figure
    figure;
    plot3(segment(1, :) * scale, segment(2, :) * scale, segment(3, :) * scale, '.');
    xlabel(['[' prefix 'm]']);
    ylabel(['[' prefix 'm]']);
    zlabel(['[' prefix 'm]']);
    axis equal;
    grid on;
    box on;
    
end