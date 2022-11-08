function rect = makeCartRect(rect_pos, Lx, Ly, theta, num_points, plot_rect)
%MAKECARTRECT Create evenly distributed Cartesian points covering a rectangle.
%
% DESCRIPTION:
%     makeCartRect creates an array of the Cartesian coordinates of points
%     evenly distributed over a rectangle. The position of the rectangle is
%     set by rect_pos, and its height and width are specified by Lx and Ly.
%     If rect_pos is given as a 2 element vector, the Cartesian rectangle
%     points are returned in 2D. These points are rotated within the plane
%     if an angle theta is given. If rect_pos is given as a 3 element
%     vector, the Cartesian rectangle points are returned in 3D assuming
%     the rectangle lies in the x-y plane. The orientation of the 3D points
%     can be adjusted by giving theta as a 3 element vector specifying the
%     yaw, pitch, and roll of the rectangle.
%
% USAGE:
%     rect = makeCartRect(rect_pos, Lx, Ly, theta, num_points)
%     rect = makeCartRect(rect_pos, Lx, Ly, theta, num_points, plot_rect)
%
% INPUTS:
%     rect_pos    - Cartesian position of the centre of the rectangle given
%                   as a two (2D) or three (3D) element vector [m].
%     Lx          - Height of the rectangle (along the x-coordinate before 
%                   rotation) [m].
%     Ly          - Width of the rectangle (along the y-coordinate before
%                   rotation) [m].
%     theta       - Either a scalar (2D) or three element vector (3D) 
%                   [tx, ty, tz] specifying the orientation of the
%                   rectangle [deg]. In 3D, the rotations are specified
%                   about x-y'-z'' (intrinsic rotations) or z-y-x
%                   (extrinsic rotations). All rotations are
%                   counter-clockwise. Can be set to [] if no rotation.
%     num_points  - Approximate number of points on the rectangle. This is
%                   rounded upwards to ensure the number of points along
%                   each side is an integer. 
%
% OPTIONAL INPUTS
%     plot_rect   - Boolean controlling whether the Cartesian points are
%                   plotted (default = false).
%
% OUTPUTS:
%     rect        - 2 x num_points* or 3 x num_points* array of Cartesian
%                   coordinates.
%                   *Note: num_points is adjusted, as described above.
%
% ABOUT:
%     author      - Elliott Wise and Bradley Treeby
%     date        - 30th April 2018
%     last update - 16th July 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2019 Elliott Wise and Bradley Treeby
%
% See also makeCartDisc, makeCartBowl

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

% check for theta input
if nargin < 5
    theta = [];
end

% check for plot_rect input
if nargin < 6
    plot_rect = false;
end

% find number of points in along each axis
npts_x = ceil(sqrt(num_points .* Lx ./ Ly));
npts_y = ceil(num_points ./ npts_x);

% recalculate the true number of points
num_points = npts_x .* npts_y;

% distance between points in each dimension
d_x = 2 ./ npts_x;
d_y = 2 ./ npts_y;

% compute canonical rectangle points ([-1, 1] x [-1, 1], z=0 plane)
p_x = linspace(-1 + d_x/2, 1 - d_x/2, npts_x);
p_y = linspace(-1 + d_y/2, 1 - d_y/2, npts_y);
[P_x, P_y] = meshgrid(p_x, p_y);
p0 = [P_x(:), P_y(:)].';

% add z-dimension points if in 3D
if length(rect_pos) == 3
    p0 = [p0; zeros(1, num_points)];
end

% transform the canonical rectangle points to give the specified rectangle
if length(rect_pos) == 2
    
    % scaling transformation
    S = [Lx, 0; 0, Ly] ./ 2;
    
    % rotation
    if isempty(theta)
        R = eye(2);
    else
        R = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
    end
    
else
    
    % scaling transformation
    S = [Lx, 0, 0; 0, Ly, 0; 0, 0, 2] ./ 2;
    
    % rotation
    if isempty(theta)
        
        % no rotation
        R = eye(3);
        
    else
        
        % ref: https://en.wikipedia.org/wiki/Rotation_matrix#General_rotations
        % Intrinsic rotations chain from right to left, extrinsic rotations
        % chain from left to right. Here, we use x-y'-z'' (intrinsic
        % rotations) which is equivalent to z-y-x (extrinsic rotations).
        R = Rz(theta(3)) * Ry(theta(2)) * Rx(theta(1));
        
    end
end

% combine scaling and rotation matrices
A = R * S;

% apply this transformation to the canonical points
p0 = A * p0;

% shift the rectangle to the appropriate centre
rect = bsxfun(@plus, p0, rect_pos.');

% plot results
if plot_rect
    
    % select suitable axis scaling factor
    [~, scale, prefix] = scaleSI(max(rect(:)));
    
    % create the figure
    if length(rect_pos) == 2
        figure;
        plot(rect(2, :) .* scale, rect(1, :) .* scale, '.');
        set(gca, 'YDir', 'reverse');
        xlabel(['y-position [' prefix 'm]']);
        ylabel(['x-position [' prefix 'm]']);
        axis equal;
    else
        figure;
        plot3(rect(2, :) * scale, rect(1, :) * scale, rect(3, :) * scale, '.');
        xlabel(['y-position [' prefix 'm]']);
        ylabel(['x-position [' prefix 'm]']);
        zlabel(['z-position [' prefix 'm]']);
        axis equal;
        grid on;
        box on;
    end
        
end

end

% generate 3D rotation matrices
function R = Rx(theta)
    R = [1, 0, 0; 0, cosd(theta), -sind(theta); 0, sind(theta), cosd(theta)];
end

function R = Ry(theta)
    R = [cosd(theta), 0, sind(theta); 0, 1, 0; -sind(theta), 0, cosd(theta)];
end

function R = Rz(theta)
    R = [cosd(theta), -sind(theta), 0; sind(theta), cosd(theta), 0; 0, 0, 1];
end