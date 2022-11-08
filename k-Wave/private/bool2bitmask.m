function bitmask = bool2bitmask(boolean_array)
%BOOL2BITMASK Convert Boolean values to a 64-bit bitmask.
%
% DESCRIPTION:
%     bool2bitmask takes an input array of Boolean values, and returns a
%     64-bit bitmask.
%
% USAGE:
%     bitmask = bool2bitmask(boolean_array)
%
% INPUTS:
%     boolean_array - input array
%
% OUTPUTS:
%     bitmask       - output bitmask
%
% ABOUT:
%     author        - Bradley Treeby and Jiri Jaros
%     date          - 6th August 2012
%     last update   - 8th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2017 Bradley Treeby and Jiri Jaros

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

% create 64-bit unsigned integer output with enough bits to hold all the
% Boolean values converted to a bitmask
bitmask = zeros(ceil(length(boolean_array) / 64), 'uint64');

% loop through each element in the input array, convert to a 64-bit
% unsigned integer (which will give either ..0001 or ..0000), do a
% bitwise shift to move it to the appropriate element, then combine with
% output array using a bitwise or
for index = 1:length(boolean_array)
    bitmask = bitor(bitmask, bitshift(uint64(boolean_array(index)), index - 1));
end