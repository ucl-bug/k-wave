function beamPlot(varargin)
%BEAMPLOT Plot volumetric data using intersecting planes.
%
% DESCRIPTION:
%     beamPlot plots 3D volumetric data as intersecting planes using slice.
%     The data is assumed to be indexed as (x, y, z). The central x-y and
%     x-z planes are extracted and plotted as intersecting planes. The
%     first y-z plane can also be plotted by setting the optional input
%     plot_yz to true. This orientation is useful for plotting the beam
%     patterns produced by ultrasound transducers facing in the
%     x-direction.
%
%     beamPlot can alternatively be called with two 2D slices which are
%     plotted as intersecting planes indexed as (x, y) and (x, z). The
%     slices must have the same size in the x-direction.
%
%     Examples:
%         beamPlot(makeBall(30, 30, 30, 15, 15, 15, 12));
%         beamPlot(makeDisc(40, 30, 20, 15, 10), makeDisc(40, 20, 20, 10, 5))
%
% USAGE:
%     beamPlot(mat)
%     beamPlot(mat, plot_yz)
%     beamPlot(xy_slice, xz_slice)
%
% INPUTS:
%     mat         - 3D matrix to plot
%     plot_yz     - Boolean controlling whether the first y-z plane is
%                   displayed
%     xy_slice    - slice to plot in the x-y plane
%     xz_slice    - slice to plot in the x-z plane
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 8th May 2012
%     last update - 22nd February 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2017 Bradley Treeby
%
% See also slice

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
VIEW_AZIMUTH    = 160;
VIEW_ELEVATION  = 20;
COLOR_MAP       = jet(256);

if nargin == 1 || (nargin == 2 && islogical(varargin{2}))

    % assign pseudonyms for matrix input
    mat = varargin{1};
    
    % check for plot_yz option
    if nargin == 2
        plot_yz = varargin{2};
    else
        plot_yz = false;
    end
    
    % extract size
    sz = size(mat);
    
    % extract x-y and x-z slices from the input matrix, and then plot
    % slices with the edge line switched off
    if plot_yz
        hz = slice(mat, round(sz(2)/2), 1, round(sz(3)/2));
    else
        hz = slice(mat, round(sz(2)/2), [], round(sz(3)/2));
    end
    set(hz, 'EdgeColor','none');
    
elseif nargin == 2
    
    % assign pseudonyms for the inputs
    xy_slice = varargin{1};
    xz_slice = varargin{2};
    
    % check the dimensions of the slices
    if size(xy_slice, 1) ~= size(xz_slice, 1)
        error('The inputs for xy_slice and xz_slice must be the same size in the x-direction.');
    end
    
    % create 3D data from horizontal slice
    xvol(:, : ,1) = xy_slice; 
    xvol(:, :, 2) = xy_slice; 

    % create 3D data from vertical slice
    yvol(:, 1, :) = xz_slice; 
    yvol(:, 2, :) = xz_slice; 

    % get sizes
    sz1 = size(xvol);
    sz2 = size(yvol);

    % plot slices with the edge line switched off
    hz1 = slice((1:sz1(2)) - ceil(sz1(2)/2), (1:sz1(1)) - ceil(sz1(1)/2), (1:sz1(3)) - ceil(sz1(3)/2), xvol, [], [], 1);
    set(hz1, 'EdgeColor','none');
    hold on;
    hz2 = slice((1:sz2(2)) - ceil(sz2(2)/2), (1:sz2(1)) - ceil(sz2(1)/2), (1:sz2(3)) - ceil(sz2(3)/2), yvol, 1, [], []);
    set(hz2, 'EdgeColor','none');
    
end

% figure settings
colormap(COLOR_MAP);
axis image;
axis tight;
view(VIEW_AZIMUTH, VIEW_ELEVATION);
axis off;