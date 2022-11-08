function overlayPlot(varargin)
%IMAGEOVERLAY Overlay two images.
%
% DESCRIPTION:
%     overlayPlot overlays two 2D images. The background is displayed using
%     a grayscale map. For a non-zero dynamic range, the foreground image
%     is log compressed (discarding negative values), thresholded to a
%     particular dynamic range, and overlaid using an alpha value of 0.5.
%     If the dynamic range is set to zero, the foreground image is
%     overlaid without additional processing.
%
%     Example:
%         x = rand(128);
%         y = peaks(128);
%         overlayPlot(x, y);
%
%     Note, in earlier versions of MATLAB, the ytick labels on the colorbar
%     do not match the plot scale.
%
% USAGE:
%     overlayPlot(bg, fg)
%     overlayPlot(bg, fg, ...)
%     overlayPlot(x, y, bg, fg)
%     overlayPlot(x, y, bg, fg, ...)
%
% INPUTS:
%     x, y         - vectors describing the position of the pixels in the
%                    image equivalent to image(x, y, c)
%     bg           - background image
%     fg           - foreground image
%
% OPTIONAL INPUTS
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%
%     'ColorBar'   - Boolean controlling whether a colorbar is displayed
%                    (default = false).
%     'ColorBarTitle'
%                  - String defining the title used for the colorbar
%                    (default = 'dB' with log compression, otherwise '').
%     'ColorMap'   - String defining the colormap used for the overlay
%                    (default = 'jet').
%     'LogComp'    - Boolean controlling whether the forergound image is
%                    log compressed before display (default = true).
%     'LogCompRef' - Reference value used in the log compression, where
%                    fg_compressed = 20 * log10(fg ./ fg_ref)
%                    (default = max(fg(:))).
%     'NumColors'  - Number of colors used in the colormaps 
%                    (default = 256).
%     'PlotScale'  - Plot scale used to display the foreground image
%                    (default = [-30, 0] with log compression, otherwise
%                    [min(fg(:)), max(fg(:))]).
%     'Transparency' 
%                  - Transparency used for the foreground image, between 0
%                    and 1 (default = 0.5).
%                     
% ABOUT:
%     author       - Bradley Treeby
%     date         - 17th October 2012
%     last update  - 2nd July 2021
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2021 Bradley Treeby

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

% set the literals
num_colors          = 256;
transparency        = 0.5;
log_compression     = true;
color_map_str       = 'jet';
plot_colorbar       = false;
color_bar_title     = '';

% extract the number of inputs and optional inputs
if nargin == 2 || ischar(varargin{3})
    req_inputs = 2;
else
    req_inputs = 4;
end

% extract the required inputs
if req_inputs == 2
    bg = varargin{1};
    fg = varargin{2};
else
    x_vec = varargin{1};
    y_vec = varargin{2};
    bg = varargin{3};
    fg = varargin{4};
end

% set default reference value for log compression
fg_ref = max(fg(:));

% set flags to track changes
plot_scale_default  = true;
colorbar_title_default  = true;

% extract the optional inputs
opt_inputs = nargin - req_inputs;
if rem(opt_inputs, 2)
    error('Optional input parameters must be given as param, value pairs.');
elseif opt_inputs > 0
    for input_index = req_inputs + 1:2:nargin
        switch varargin{input_index}
            case 'ColorBar'
                plot_colorbar = varargin{input_index + 1};
            case 'ColorBarTitle'
                color_bar_title = varargin{input_index + 1};
                colorbar_title_default = false;
            case 'ColorMap'             
                color_map_str = varargin{input_index + 1};
            case 'LogComp'
                log_compression = varargin{input_index + 1};
            case 'LogCompRef'
                fg_ref = varargin{input_index + 1};
            case 'NumColors'
                num_colors = varargin{input_index + 1};
            case 'PlotScale'
                plot_scale = varargin{input_index + 1};
                plot_scale_default = false;
            case 'Transparency'
                transparency = varargin{input_index + 1};
            otherwise
                error(['Unknown optional input ' varargin{input_index} '.']);
        end
    end
end

% set default foreground plot scale if not modified by the user
if plot_scale_default
    if log_compression
        plot_scale = [-30, 0];
    else
        plot_scale = [min(fg(:)), max(fg(:))];
    end
end

% set color bar title if not modified by the user
if colorbar_title_default
    if log_compression
        color_bar_title = 'dB';
    else
        color_bar_title = '';
    end
end

% evaluate foreground colormap
eval(['color_map = ' color_map_str '(' num2str(num_colors) ');']);

% scale the background image from 0 to num_colors
bg = bg - min(bg(:));
bg = round(num_colors * bg / max(bg(:)));

% convert the background image to true color
bg = ind2rgb(bg, gray(num_colors));

% plot the background image
if req_inputs == 4
    image(x_vec, y_vec, bg);
else
    image(bg);
end

% discard negative data, and apply log compression
if log_compression
    fg(fg <= 0) = 0;
    fg = 20 * log10(fg ./ fg_ref);
end

% scale the background image from 1 to num_colors based on plot_scale
fg = fg - plot_scale(1);
fg(fg < 0) = 0;
fg = ceil(fg * (num_colors - 1) / (plot_scale(2) - plot_scale(1)));
fg = fg + 1;

% compute the alpha channel
alpha = transparency * ones(size(fg));
alpha(fg == 1) = 0;

% convert the background image to true color
fg = ind2rgb(fg, color_map);

% plot the foreground image and set the alpha channel
hold on;
if req_inputs == 4
    fg_im = image(x_vec, y_vec, fg);
else
    fg_im = image(fg);
end
set(fg_im, 'AlphaData', alpha);

% color bar
if plot_colorbar
    colormap(color_map);
    cb = colorbar;
    caxis(plot_scale);
    title(cb, color_bar_title);
end