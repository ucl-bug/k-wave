function disc = makeDisc(Nx, Ny, cx, cy, radius, plot_disc)
%MAKEDISC Create a binary map of a filled disc within a 2D grid.
%
% DESCRIPTION:
%     makeDisc creates a binary map of a filled disc within a
%     two-dimensional grid (the disc position is denoted by 1's in the
%     matrix with 0's elsewhere). A single grid point is taken as the disc
%     centre thus the total diameter of the disc will always be an odd
%     number of grid points. As the returned disc has a constant radius, if
%     used within a k-Wave grid where dx ~= dy, the disc will appear oval
%     shaped. If part of the disc overlaps the grid edge, the rest of the
%     disc will wrap to the grid edge on the opposite side.
%
% USAGE:
%     disc = makeDisc(Nx, Ny, cx, cy, radius)
%     disc = makeDisc(Nx, Ny, cx, cy, radius, plot_disc)
%
% INPUTS:
%     Nx, Ny          - size of the 2D grid [grid points]
%     cx, cy          - centre of the disc [grid points]
%     radius          - disc radius [grid points]
%
% OPTIONAL INPUTS:
%     plot_disc       - Boolean controlling whether the disc is plotted
%                       using imagesc (default = false)
%
% OUTPUTS:
%     disc            - 2D binary map of a filled disc
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 9th June 2009
%     last update     - 8th February 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2019 Bradley Treeby
%
% See also makeCircle, makeBall 

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

% check for plot_disc input
if nargin < 6
    plot_disc = false;
end

% force integer values for grid size
Nx = round(Nx);
Ny = round(Ny);
cx = round(cx);
cy = round(cy);

% check for zero values
if cx == 0
    cx = floor(Nx/2) + 1;
end
if cy == 0
    cy = floor(Ny/2) + 1;
end

% check the inputs
if cx < 1 || cx > Nx || cy < 1 || cy > Ny
    error('Disc center must be within grid.');
end

% create empty matrix
disc = zeros(Nx, Ny);

% define pixel map
r = makePixelMap(Nx, Ny, 'Shift', [0, 0]);

% create disc
disc(r <= radius) = MAGNITUDE;

% shift centre
cx = round(cx) - ceil(Nx/2);
cy = round(cy) - ceil(Ny/2);
disc = circshift(disc, [cx, cy]);

% create the figure
if plot_disc
    figure;
    imagesc(disc, [-1, 1]);
    colormap(getColorMap);
    axis image;
    xlabel('y-position [grid points]');
    ylabel('x-position [grid points]');
end