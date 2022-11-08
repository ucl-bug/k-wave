function cm = getColorMap(num_colors)
%GETCOLORMAP Return default k-Wave color map.
%
% DESCRIPTION:
%     getColorMap returns the default color map used for display and
%     visualisation across the k-Wave Toolbox. Zero values are displayed as
%     white, positive values are displayed as yellow through red to black,
%     and negative values are displayed as light to dark blue-greys. If no
%     value for num_colors is provided, cm will have 256 colors.    
%
% USAGE:
%     cm = getColorMap()
%     cm = getColorMap(num_colors)
%
% OPTIONAL INPUTS:
%     num_colors  - number of colors in the color map (default = 256)
%
% OUTPUTS:
%     cm          - three column color map matrix which can be applied
%                   using colormap
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 3rd July 2009
%     last update - 4th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby
%
% See also colormap

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

% set literals
if nargin == 0
    neg_pad = 48;
    num_colors = 256;   
else
    neg_pad = round(48 * num_colors / 256);
end

% define colour spectrums
neg = bone(num_colors / 2 + neg_pad);
neg = neg(1 + neg_pad:end, :);
pos = flipud(hot(num_colors / 2));

% create custom colour map
cm = [neg; pos];