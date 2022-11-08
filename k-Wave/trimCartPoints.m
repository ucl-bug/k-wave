function points = trimCartPoints(kgrid, points)
%TRIMCARTPOINTS Remove Cartesian points that are not within a kgrid.
%
% DESCRIPTION:
%     trimCartPoints filters a dim x num_points array of Cartesian points
%     so that only those within the bounds of a given kgrid remain.
%
% USAGE:
%     points = trimCartPoints(kgrid, points)
%
% INPUTS:
%     kgrid       - Object of the kWaveGrid class defining the Cartesian
%                   and k-space grid fields. 
%     points      - dim x num_points array of Cartesian coordinates to trim
%                   [m].
%
% OUTPUTS:
%     points      - dim x num_points array of Cartesian coordinates that
%                   lie within the grid defined by kgrid [m].
%
% ABOUT:
%     author      - Elliott Wise and Bradley Treeby
%     date        - 16th March 2017
%     last update - 4th February 2018
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017-2018 Elliott Wise and Bradley Treeby
%
% See also cart2grid, grid2cart

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

% find indices for points within the simulation domain
ind_x = (points(1, :) >= kgrid.x_vec(1)) & (points(1, :) <= kgrid.x_vec(end));
if kgrid.dim > 1
    ind_y = (points(2, :) >= kgrid.y_vec(1)) & (points(2, :) <= kgrid.y_vec(end));
end
if kgrid.dim > 2
    ind_z = (points(3, :) >= kgrid.z_vec(1)) & (points(3, :) <= kgrid.z_vec(end));
end

% combine indices
switch kgrid.dim
    case 1
        ind = ind_x;
    case 2
        ind = ind_x & ind_y;
    case 3
        ind = ind_x & ind_y & ind_z;
end

% output only valid points
points = points(:, ind);