function test_pass = roundEven_roundOdd(~, ~)
% DESCRIPTION:
%       Unit test to test the functions roundEven and roundOdd.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 29th April 2018
%       last update - 29th April 2018
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018 Bradley Treeby

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

% set pass variable
test_pass = true;

% generate matrix of input values
x = rand(50, 50);

% round odd
x_odd = roundOdd(x);

% check odd
if any(~rem(x_odd(:), 2))
    test_pass = false;
end

% round even
x_even = roundEven(x);

% check odd
if any(rem(x_even(:), 2))
    test_pass = false;
end