function gauss_distr = gaussian(x, magnitude, mean, variance)
%GAUSSIAN Create a Gaussian distribution.
% 
% DESCRIPTION:
%     gaussian returns a Gaussian distribution f(x) with the specified
%     magnitude, mean, and variance. If these values are not specified, the
%     magnitude is normalised and values of variance = 1 and mean = 0 are
%     used. For example running 
%
%         x = -3:0.05:3;
%         plot(x, gaussian(x));
%
%     will plot a normalised Gaussian distribution.
%
%     Note, the full width at half maximum of the resulting distribution
%     can be calculated by FWHM = 2 * sqrt(2 * log(2) * variance).
%
% USAGE:
%     gauss_distr = gaussian(x)
%     gauss_distr = gaussian(x, magnitude)
%     gauss_distr = gaussian(x, magnitude, mean)
%     gauss_distr = gaussian(x, magnitude, mean, variance)
%
% INPUTS:
%     x           - x-axis variable
%
% OPTIONAL INPUTS:
%     magnitude   - bell height (default = normalised)
%     mean        - mean or expected value (default = 0)
%     variance    - variance ~ bell width (default = 1)
%
% OUTPUTS:
%     gauss_distr - Gaussian distribution
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 4th December 2009
%     last update - 4th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby

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

% if the variance is not given, define as unity
if nargin < 4
    variance = 1;
end

% if the mean is not given, center the distribution
if nargin < 3
    mean = 0;
end

% if the magnitude is not given, normalise
if nargin < 2
    magnitude = (2 * pi * variance)^-0.5;
end

% compute the gaussian
gauss_distr = magnitude * exp(-(x - mean).^2 ./ (2 * variance));