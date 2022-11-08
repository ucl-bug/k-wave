function mat3D = revolve2D(mat2D)
%REVOLVE2D Form 3D matrix from revolution of 2D matrix.
%
% DESCRIPTION:
%     revolve2D revolves the values of a 2D matrix about the x axis (or
%     matrix rows) to form a 3D matrix. A single point is taken as the
%     x-axis origin, thus for an input matrix of size m by n, the output
%     matrix will be m by (2n-1) by (2n-1). Values outside the
%     interpolation range are set to zero.
%
% USAGE:
%     mat3D = revolve2D(mat2D)
%
% INPUTS:
%     mat2D       - 2D input matrix
%
% OUTPUTS:
%     mat3D       - 3D output matrix
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 19th April 2011
%     last update - 7th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011-2017 Bradley Treeby

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

% start timer
tic;

% update command line status
disp('Revolving 2D matrix to form a 3D matrix...');

% get size of matrix
[m, n] = size(mat2D);

% create the reference axis for the 2D image
r_axis_one_sided = 0:(n - 1);
r_axis_two_sided = -(n - 1):(n - 1);

% compute the distance from every pixel in the z-y cross-section of the 3D
% matrix to the rotation axis
[z, y] = meshgrid(r_axis_two_sided, r_axis_two_sided);
r = sqrt(y.^2 + z.^2);

% create empty image matrix
mat3D = zeros(m, 2 * n - 1, 2 * n - 1);

% loop through each cross section and create 3D matrix
for x_index = 1:m
    mat3D(x_index, :, :) = interp1(r_axis_one_sided, mat2D(x_index, :), r, [], 0);
end

% update command line status
disp(['  completed in ' scaleTime(toc) 's']);