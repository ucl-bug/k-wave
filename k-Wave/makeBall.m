function ball = makeBall(Nx, Ny, Nz, cx, cy, cz, radius, plot_ball, binary)
%MAKEBALL Create a binary map of a filled ball within a 3D grid.
%
% DESCRIPTION:
%     makeBall creates a binary map of a filled ball within a
%     three-dimensional grid (the ball position is denoted by 1's in the
%     matrix with 0's elsewhere). A single grid point is taken as the ball
%     centre thus the total diameter of the ball will always be an odd
%     number of grid points.
%
% USAGE:
%     ball = makeBall(Nx, Ny, Nz, cx, cy, cz, radius)
%     ball = makeBall(Nx, Ny, Nz, cx, cy, cz, radius, plot_ball)
%     ball = makeBall(Nx, Ny, Nz, cx, cy, cz, radius, [], binary)
%
% INPUTS:
%     Nx, Ny, Nz      - size of the 3D grid [grid points]
%     cx, cy, cz      - centre of the ball [grid points]
%     radius          - ball radius [grid points]
%
% OPTIONAL INPUTS:
%     plot_ball       - Boolean controlling whether the ball is
%                       plotted using voxelPlot (default = false)
%     binary          - Boolean controlling whether the ball map is
%                       returned as a double precision matrix (false) or
%                       a logical matrix (true) (default = false)
%
% OUTPUTS:
%     ball            - 3D binary map of a filled ball
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 1st July 2009
%     last update     - 31st July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2019 Bradley Treeby
%
% See also makeCircle, makeDisc, makeSphere 

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

% define literals
MAGNITUDE = 1;

% check for plot_ball input
if nargin < 8 || isempty(plot_ball)
    plot_ball = false;
end

% check for binary input
if nargin < 9 || isempty(binary)
    binary = false;
end

% force integer values
Nx = round(Nx);
Ny = round(Ny);
Nz = round(Nz);
cx = round(cx);
cy = round(cy);
cz = round(cz);

% check for zero values
if cx == 0
    cx = floor(Nx / 2) + 1;
end
if cy == 0
    cy = floor(Ny / 2) + 1;
end
if cz == 0
    cz = floor(Nz / 2) + 1;
end

% create empty matrix
if binary
    ball = false(Nx, Ny, Nz);
else
    ball = zeros(Nx, Ny, Nz);
end

% define pixel map
r = makePixelMap(Nx, Ny, Nz, 'Shift', [0, 0, 0]);

% create ball
ball(r <= radius) = MAGNITUDE;

% shift centre
cx = cx - ceil(Nx / 2);
cy = cy - ceil(Ny / 2);
cz = cz - ceil(Nz / 2);
ball = circshift(ball, [cx, cy, cz]);

% plot results
if plot_ball
    voxelPlot(double(ball));
end