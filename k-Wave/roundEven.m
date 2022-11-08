function x = roundEven(x)
%ROUNDEVEN Round towards the nearest even number.
%
% DESCRIPTION:
%     roundEven rounds each element of x towards the nearest even integer.
%
% USAGE:
%     x = roundEven(x)
%      
% INPUTS:
%     x               - input value/s
%
% OUTPUTS:
%     x               - value/s of x rounded to nearest even integer
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 6th August 2017
%     last update     - 29th April 2018
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017-2018 Bradley Treeby
%
% See also round, roundOdd

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

% round towards the nearest even number
x = 2 .* round(x ./ 2);