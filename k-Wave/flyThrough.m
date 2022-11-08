function flyThrough(mat, dim, fps, loop)
%FLYTHROUGH Display a three-dimensional matrix slice by slice.
%
% DESCRIPTION:
%     flyThrough produces an animated visualisation of the data within a
%     three-dimensional matrix mat by displaying the matrix slice by slice
%     using imagesc and the colormap returned by getColorMap. A constant
%     scaling parameter is automatically computed using the maximum value
%     within mat.
%
% USAGE:
%     flythrough(mat)
%     flythrough(mat, dim)
%     flythrough(mat, dim, speed)
%     flythrough(mat, dim, speed, loop)
%
% INPUTS:
%     mat         - the three-dimensional matrix to visualise
%
% OPTIONAL INPUTS:
%     dim         - matrix dimension through which the slices are taken
%                   (default = 1)
%     fps         - maximum number of frames per second to display
%                   (default = 5)
%     loop        - number of times to loop the animation (default = 1)
%
%     Note: The maximum achievable fps is determined by the speed at which
%     the image slices can be extracted and displayed.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 3rd July 2009
%     last update - 4th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby
%
% See also getColorMap, imagesc

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

% check the matrix input
if numDim(mat) ~= 3
    error('Input matrix must be three-dimensional.');
end

% extract optional inputs
if nargin < 4
    loop = 1;
end
if nargin < 3
    fps = 5;
end
if nargin < 2
    dim = 1;
end

% update command line status
disp('Running fly through visualisation...');

% extract the size
[a, b, c] = size(mat);

% determine the number of slices to display
switch dim
    case 1
        num_slices = a;
    case 2
        num_slices = b;
    case 3
        num_slices = c;
end

% initialise the figure
figure;

% choose a suitable plot scale
plot_scale = max(mat(:));

% load colormap
cmap = getColorMap();

% loop the animation
for loop_index = 1:loop
   
    % display each matrix slice by slice
    for slice_index = 1:num_slices
        tic
        switch dim
            case 1
                slice = squeeze(mat(slice_index, :, :));
            case 2
                slice = squeeze(mat(:, slice_index, :));
            case 3
                slice = squeeze(mat(:, :, slice_index));
        end
        imagesc(slice, [-plot_scale, plot_scale]);
        axis image;
        colormap(cmap);
        drawnow;
        pause_time = 1 / fps - toc;
        pause_time(pause_time < 0) = 0;
        pause(pause_time);
    end   
end

% update command line status
disp('  visualisation complete');