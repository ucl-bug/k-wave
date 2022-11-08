function stackedPlot(varargin)
%STACKEDPLOT Stacked linear plot.
%
% DESCRIPTION:
%     stackedPlot produces a series of stacked linear plots of the rows in
%     data (a 2D matrix) against the vector x. The vector y defines the
%     y-axis label for each linear plot. The plot scaling is defined using
%     the global maximum and minimum of the data. The plot can be annotated
%     in the normal fashion after display.
%
%     Examples:
%         stackedPlot(rand(5, 100));
%         stackedPlot(0.1:0.1:10, rand(5, 100));
%         stackedPlot(0.1:0.1:10, {'a', 'b', 'c', 'd', 'e'}, rand(5, 100));
%
% USAGE:
%         stackedPlot(data)
%         stackedPlot(x, data)
%         stackedPlot(x, y, data)
%
%         stackedPlot(data, ...)
%         stackedPlot(x, data, ...)
%         stackedPlot(x, y, data, ...)
%
% INPUTS:
%     x       - vector defining the x-axis values
%     y       - vector defining the y-axis labels used for each plot
%     data    - 2D matrix to plot
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings. 
%
%     'Spacing'   - Spacing between the individual plots as a fraction of
%                   the maximum plot range. (default = 0.1) 
%     'Symmetric' - Boolean controlling whether the max and minimum of
%                   the individual plots are forced to be symmetric. This
%                   forces the y-ticks to be at 0. (default = false) 
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 14th November 2011
%     last update - 7th June 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011-2017 Bradley Treeby
%
% See also plot

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

% set default parameters
symmetric = false;
spacing = 0.1;

% check which input combination is given
switch nargin
    case 1
        
        % stackedPlot(data)
        inputs = 1;
        
    case 2
        
        % stackedPlot(x, data)
        inputs = 2;
        
    case 3
        if numel(varargin{3}) == 1
            
            % stackedPlot(data, 'string', value)
            inputs = 1;
            
        else
            
            % stackedPlot(x, y, data)
            inputs = 3;
            
        end
    case 4
        
        % stackedPlot(x, data, 'string', value)
        inputs = 2;
        
    otherwise
        if numel(varargin{3}) == 1
            
            % stackedPlot(data, 'string', value, ...)
            inputs = 1;
            
        elseif numel(varargin{4}) == 1
            
            % stackedPlot(x, data, 'string', value, ...)
            inputs = 2;
            
        else
            
            % stackedPlot(x, y, data, 'string', value, ...)
            inputs = 3;
            
        end
end
    
% extract user inputs
switch inputs
    case 1
        data = varargin{1};
        [Ny, Nx] = size(data);
        x = 1:Nx;
        y = 1:Ny;        
    case 2
        x = varargin{1};
        data = varargin{2};
        [Ny, Nx] = size(data);
        y = 1:Ny;        
    case 3
        x = varargin{1};
        y = varargin{2};
        data = varargin{3};
        [Ny, Nx] = size(data);        
end
        
% replace defaults with user defined values if provided
for input_index = (inputs + 1):2:length(varargin)
    switch varargin{input_index}
        case 'Symmetric'
            symmetric = varargin{input_index + 1};
        case 'Spacing'
            spacing = varargin{input_index + 1};
        otherwise
            error('Unknown optional input.');
    end
end

% get the maximum range of the input data
if symmetric
    mx = max(abs(data(:)));
    mn = -mx;
else
    mx = max(data(:));
    mn = min(data(:));
end
rng = (1 + spacing) * (mx - mn);

% create a series of linear plots of each row in data offset by the
% maximum data range
plot(x, data + repmat(rng * (Ny:-1:1).', 1, Nx), 'k-');

% calculate the y-limits based on the data range and spacing
axis tight;
ylims = get(gca, 'YLim');
ylims(1) = ylims(1) - spacing * rng;
ylims(2) = ylims(2) + spacing * rng;

% define the location for the y-ticks
if symmetric
    yticks = (1:Ny) * rng;
else
    yticks = (1:Ny) * rng + (mx + mn) / 2;
end

% set the y-axis labels to y
try
    % use new name for YTickLabel
    set(gca, 'YTick', yticks, 'YTickLabel', fliplr(y), 'YLim', ylims);

catch
    % use old name for YTickLabel
    set(gca, 'YTick', yticks, 'YTickLabels', fliplr(y), 'YLim', ylims); 
end