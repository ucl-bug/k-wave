function alpha = db2neper(alpha, y)
%DB2NEPER Convert decibels to nepers.
%
% DESCRIPTION:
%     db2neper converts an attenuation coefficient in units of
%     dB/(MHz^y cm) to units of Nepers/((rad/s)^y m).
%
% USAGE:
%     alpha = db2neper(alpha)
%     alpha = db2neper(alpha, y)
%
% INPUTS:
%     alpha       - attenuation in dB/(MHz^y cm)
%
% OPTIONAL INPUTS:
%     y           - power law exponent (default = 1)
%
% OUTPUTS:
%     alpha       - attenuation in Nepers/((rad/s)^y m)
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 27th March 2009
%     last update - 4th June 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby
%
% See also neper2db

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

% set default y value if not given by user
if nargin == 1
    y = 1;
end

% calculate conversion
alpha = 100 * alpha .* (1e-6 / (2*pi)).^y ./ (20 * log10(exp(1)));