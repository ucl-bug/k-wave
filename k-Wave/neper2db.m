function alpha = neper2db(alpha, y)
%NEPER2DB Convert nepers to decibels.
%
% DESCRIPTION:
%     neper2db converts an attenuation coefficient in units of
%     Nepers/((rad/s)^y m)to units of dB/(MHz^y cm).
%
% USAGE:
%     alpha = neper2db(alpha)
%     alpha = neper2db(alpha, y)
%
% INPUTS:
%     alpha       - attenuation in Nepers/((rad/s)^y m)
%
% OPTIONAL INPUTS:
%     y           - power law exponent (default = 1)
%
% OUTPUTS:
%     alpha       - attenuation in dB/(MHz^y cm)
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 3rd December 2009
%     last update - 7th June 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2017 Bradley Treeby
%
% See also db2neper

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

if nargin == 1
    y = 1;
end

alpha = 20 * log10(exp(1)) * alpha * (2 * pi * 1e6)^y / 100;