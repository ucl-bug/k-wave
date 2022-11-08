function [fwhm_val, fwhm_pos] = fwhm(f, x, plot_fwhm)
%FWHM Compute the full width at half maximum.
%
% DESCRIPTION:
%     fwhm calculates the Full Width at Half Maximum (FWHM) of a positive
%     1D input function f(x) with spacing given by x.
%
% USAGE:
%     fwhm_val = fwhm(f, x)
%     fwhm_val = fwhm(f, x, plot_fwhm)
%
% INPUTS:
%     f           - f(x) data
%     x           - x data or dx
%
% OPTIONAL INPUTS:
%     plot_fwhm   - Boolean controlling whether a plot of the function
%                   and FWHM is produced.
%
% OUTPUTS:
%     fwhm_val    - FWHM of f(x)
%     fwhm_pos    - position of leading and trailing edges
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 5th May 2009
%     last update - 19th October 2021
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2021 Bradley Treeby

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

% check for optional inputs
if nargin == 2
    plot_fwhm = false;
end

% check if dx is given in place of an x array
if length(x) == 1
    x = 0:x:(length(f) - 1) * x;
end

% check the input is 1D
if numDim(f) ~= 1
    error('Input function must be 1 dimensional.');
end

% find the maximum value of f(x)
[f_max, i_max] = max(f);

% setup the indexing variables for finding the leading edge
index = i_max;

% loop until the index at half maximum is found
while f(index) > (0.5 * f_max)
    index = index - 1;
    if index < 1
        error('Left half maximum not found.');
    end
end

% linearly interpolate between the previous values to find the value of x
% at the leading edge at half maximum
m = (f(index+1) - f(index)) / (x(index+1) - x(index));
c = f(index) - x(index) * m;
i_leading = (0.5 * f_max  - c) / m;

% setup the indexing variables for finding the trailing edge
index = i_max;

% loop until the index at half maximum is found
while f(index) > (0.5 * f_max)
    index = index + 1;
    if index > length(f)
        error('Right half maximum not found.');
    end
end

% linearly interpolate between the previous values to find the value of x
% at the trailing edge at half maximum
m = (f(index - 1) - f(index)) / (x(index - 1) - x(index));
c = f(index) - x(index) * m;
i_trailing = (0.5 * f_max  - c) / m;
    
% compute the FWHM
fwhm_val = abs(i_trailing - i_leading);

% return position values
if nargout == 2
    fwhm_pos = [i_leading, i_trailing];
end

% plot the function and the FWHM if required
if plot_fwhm   
    figure;
    plot(x, f, 'b-');
    xlabel('x');
    ylabel('f(x)');    
    hold on;
    plot([i_leading, i_trailing], [0.5 * f_max, 0.5 * f_max], 'r*-');
    title(['Full Width At Half Maximum (FWHM): ' num2str(fwhm_val) ]);
end