function [rotMat, offsetPos] = computeLinearTransform(pos1, pos2, offset)
%COMPUTELINEARTRANSFORM Compute a linear transformation matrix from two points.
%
% DESCRIPTION:
%     computeLinearTransform calculates a linear transformation matrix that
%     maps from the origin to a point in space, where the orientation of
%     the transformation is defined by a second point in space. It can be
%     used to generate a transformation matrix for bowl (or disc)
%     transducers, where pos1 is the Cartesian position of the centre of
%     the rear surface of the bowl or disc, and pos2 is any point on the
%     beam axis.
%
%     An offset can optionally be defined which specifies a distance along
%     the vector from pos1 to pos2. The position of this offset in
%     transformed space is returned as offsetPos. Note, giving a value for
%     offset does not modify the rotation matrix.
%
% USAGE:
%     [rotMat, offsetPos] = computeLinearTransform(pos1, pos2)
%     [rotMat, offsetPos] = computeLinearTransform(pos1, pos2, offset)
%
% INPUTS:
%     pos1        - Cartesian position given as a three element vector [bx,
%                   by, bz] [m].
%     pos2        - Cartesian position given as a three element vector [fx,
%                   fy, fz] [m].
%
% OPTIONAL INPUTS:
%     offset      - Scalar distance along the vector pos1->pos2 [m].
%
% OUTPUTS:
%     rotMat      - Rotation matrix.
%     offsetPos   - Vector position of offset input [m].
%
% ABOUT:
%     author      - Elliott Wise and Bradley Treeby
%     date        - 5th December 2016
%     last update - 20th October 2022
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2016-2022 Elliott Wise and Bradley Treeby
%
% See also getAffineMatrix

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

% compute vector pointing from bowl_pos to focus_pos
beam_vec = pos2 - pos1;

% normalise to give unit beam vector
beam_vec = beam_vec ./ norm(beam_vec); 

% canonical normalised beam_vec (canonical bowl_pos is [0, 0, 1])
beam_vec0 = [0, 0, -1];                

% find the rotation matrix for the bowl
u = cross(beam_vec0, beam_vec);

% normalise the rotation matrix if not zero
if any(u ~= 0)
    u = u.' ./ norm(u);
end

% find the axis-angle transformation between beam_vec and e1
theta = acos(dot(beam_vec0, beam_vec));
    
% convert axis-angle transformation to a rotation matrix
% https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
A = [0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
rotMat = cos(theta) * eye(3) + sin(theta) * A + (1 - cos(theta)) * (u * u.');
    
% compute an offset for the bowl, where bowl_centre = move from bowl_pos
% towards focus by radius
if nargin == 3
    offsetPos = pos1.' + offset * beam_vec.';
else
    offsetPos = 0;
end