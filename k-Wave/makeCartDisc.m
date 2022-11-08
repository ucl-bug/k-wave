function disc = makeCartDisc(disc_pos, radius, focus_pos, num_points, plot_disc, use_spiral)
%MAKECARTDISC Create evenly distributed Cartesian points covering a disc.
%
% DESCRIPTION:
%     makeCartDisc creates an array of the Cartesian coordinates of points
%     evenly distributed over a disc. The disc points are calculated using
%     either Fermat's spiral or a set of concentric circles. The position
%     of the disc is set by disc_pos. If given as a 2 element vector, the
%     Cartesian disc points are returned in 2D. If given as a 3 element
%     vector, the Cartesian disc points are returned in 3D assuming the
%     disc lies in the x-y plane. In 3D, the orientation of the disc can
%     optionally be defined using the focus_pos input.
%
% USAGE:
%     disc = makeCartDisc(disc_pos, radius, focus_pos, num_points)
%     disc = makeCartDisc(disc_pos, radius, focus_pos, num_points, plot_disc)
%     disc = makeCartDisc(disc_pos, radius, focus_pos, num_points, plot_disc, use_spiral)
%
% INPUTS:
%     disc_pos    - Cartesian position of the centre of the disc given as a
%                   two (2D) or three (3D) element vector [m].
%     radius      - Radius of the disc [m].
%     focus_pos   - Any point on the beam axis of the disc given as a three
%                   element vector [fx, fy, fz] [m]. Can be set to [] to
%                   define a disc in the x-y plane.
%     num_points  - Number of points on the disc.
%
% OPTIONAL INPUTS
%     plot_disc   - Boolean controlling whether the Cartesian points are
%                   plotted (default = false).
%     use_spiral  - Boolean controlling whether the Cartesian points are
%                   chosen using a spiral sampling pattern (default =
%                   false). Concentric sampling is used by default.
%
% OUTPUTS:
%     disc        - 2 x num_points or 3 x num_points array of Cartesian
%                   coordinates.
%
% ABOUT:
%     author      - Elliott Wise and Bradley Treeby
%     date        - 30th March 2017
%     last update - 2nd July 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017-2019 Elliott Wise and Bradley Treeby
%
% See also makeDisc, makeCartBowl

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

% define literals
% (ref: http://www.wolframalpha.com/input/?i=golden+angle)
% (packing number refers to the algebraic increase in samples with radius)
GOLDEN_ANGLE = 2.39996322972865332223155550663361385312499901105811504;
PACKING_NUMBER = 7;%2*pi;

if nargin < 6
    use_spiral = false;
end

% check for plot_disc input
if nargin < 5
    plot_disc = false;
end

% check input values
if radius <= 0
    error('The radius must be positive.');
end

% compute Cartesian points using spiral sampling
if use_spiral

    % compute spiral parameters
    theta = @(t) GOLDEN_ANGLE .* t;
    C = pi .* radius.^2 ./ (num_points - 1);
    r = @(t) sqrt(C .* t ./ pi);

    % compute canonical spiral points
    t = linspace(0, num_points - 1, num_points);
    p0 = bsxfun(@times, r(t), [cos(theta(t)); sin(theta(t))]);
    
% otherwise use concentric circles (note that the num_points is increased
% to ensure a full set of concentric rings)
else
    
    % compute radial sampling parameters
    num_radial = ceil(sqrt(num_points ./ pi));
    d_radial = radius ./ (num_radial - 1);
    r = (0:num_radial - 1) .* (radius - d_radial/2) ./ (num_radial - 1);
    
    % recompute the number of points that will be created below
    num_points = 1;
    for k = 2:num_radial
        num_theta = round((k - 1) * PACKING_NUMBER);
        num_points = num_points + num_theta;
    end

    % compute canonical concentric circle points
    p0 = NaN(2, num_points);
    p0(:, 1) = [0; 0];
    i_left = 2;
    for k = 2:num_radial
        num_theta = round((k - 1) * PACKING_NUMBER);
        thetas = (0:num_theta - 1) .* 2*pi ./ num_theta;
        p = r(k) .* [cos(thetas); sin(thetas)];
        i_right = i_left + num_theta - 1;
        p0(:, i_left:i_right) = p;
        i_left = i_left + num_theta;
    end
    
end

% add z-dimension points if in 3D
if length(disc_pos) == 3
    p0 = [p0; zeros(1, num_points)];
end

% if 3D and focus_pos is defined, rotate the canonical points to give the
% specified disc 
if (length(disc_pos) == 3) && (nargin > 3) && (~isempty(focus_pos))
    
    % check the focus position isn't coincindent with the disc position
    if all(disc_pos == focus_pos)
        error('The focus_pos must be different to the disc_pos.');
    end
    
    % compute rotation matrix and apply
    R = computeLinearTransform(disc_pos, focus_pos);
    p0 = R * p0;
    
end

% shift the disc to the appropriate centre
disc = bsxfun(@plus, p0, disc_pos.');

% plot results
if plot_disc
    
    % select suitable axis scaling factor
    [~, scale, prefix] = scaleSI(max(disc(:)));
    
    % create the figure
    if length(disc_pos) == 2
        figure;
        plot(disc(2, :) .* scale, disc(1, :) .* scale, '.');
        set(gca, 'YDir', 'reverse');
        xlabel(['y-position [' prefix 'm]']);
        ylabel(['x-position [' prefix 'm]']);
        axis equal;
    else
        figure;
        plot3(disc(1, :) * scale, disc(2, :) * scale, disc(3, :) * scale, '.');
        xlabel(['[' prefix 'm]']);
        ylabel(['[' prefix 'm]']);
        zlabel(['[' prefix 'm]']);
        axis equal;
        grid on;
        box on;
    end
        
end