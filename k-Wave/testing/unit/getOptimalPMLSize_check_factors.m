function test_pass = getOptimalPMLSize_check_factors(~, ~)
% DESCRIPTION:
%       Unit test to check getOptimalPMLsize returns the size with the
%       smallest prime factors. 
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

% set search range
pml_range = [10, 20];

% set maximum grid size to test
max_grid_size = 300;

% create search vector
pml_vector = pml_range(1):pml_range(2);

% -----------------
% 1D
% -----------------

% generate random grid size
grid_sz = round(max_grid_size * rand(1));

% get optimal pml size
pml_size = getOptimalPMLSize(grid_sz, pml_range);

% check the smallest prime factor in the range
smallest_factor_x = max_grid_size;
for ind = 1:length(pml_vector)
    smallest_factor_x = min(smallest_factor_x, max(factor(grid_sz(1) + 2 * pml_vector(ind))));
end

% check it matches
if smallest_factor_x ~= max(factor(grid_sz + 2 * pml_size))
    test_pass = false;
end

% -----------------
% 2D
% -----------------

% generate random grid size
grid_sz = round(max_grid_size * rand(2, 1));

% get optimal pml size
pml_size = getOptimalPMLSize(grid_sz, pml_range);

% check the smallest prime factor in the range
smallest_factor_x = max_grid_size;
smallest_factor_y = max_grid_size;
for ind = 1:length(pml_vector)
    smallest_factor_x = min(smallest_factor_x, max(factor(grid_sz(1) + 2 * pml_vector(ind))));
    smallest_factor_y = min(smallest_factor_y, max(factor(grid_sz(2) + 2 * pml_vector(ind))));
end

% check it matches
if (smallest_factor_x ~= max(factor(grid_sz(1) + 2 * pml_size(1)))) || ...
   (smallest_factor_y ~= max(factor(grid_sz(2) + 2 * pml_size(2))))
    test_pass = false;
end

% --------------------
% Axisymmetric - WSWA
% --------------------

% generate random grid size
grid_sz = round(max_grid_size * rand(2, 1));

% get optimal pml size
pml_size = getOptimalPMLSize(grid_sz, pml_range, 'WSWA');

% check the smallest prime factor in the range
smallest_factor_x = max_grid_size;
smallest_factor_y = max_grid_size;
for ind = 1:length(pml_vector)
    smallest_factor_x = min(smallest_factor_x, max(factor(grid_sz(1) + 2 * pml_vector(ind))));
    smallest_factor_y = min(smallest_factor_y, max(factor( 4 * (grid_sz(2) + pml_vector(ind)) )));
end

% check it matches
if (smallest_factor_x ~= max(factor(grid_sz(1) + 2 * pml_size(1)))) || ...
   (smallest_factor_y ~= max(factor( 4 * (grid_sz(2) + pml_size(2)) )))
    test_pass = false;
end

% --------------------
% Axisymmetric - WSWS
% --------------------

% generate random grid size
grid_sz = round(max_grid_size * rand(2, 1));

% get optimal pml size
pml_size = getOptimalPMLSize(grid_sz, pml_range, 'WSWS');

% check the smallest prime factor in the range
smallest_factor_x = max_grid_size;
smallest_factor_y = max_grid_size;
for ind = 1:length(pml_vector)
    smallest_factor_x = min(smallest_factor_x, max(factor(grid_sz(1) + 2 * pml_vector(ind))));
    smallest_factor_y = min(smallest_factor_y, max(factor( 2 * (grid_sz(2) + pml_vector(ind)) - 2 )));
end

% check it matches
if (smallest_factor_x ~= max(factor(grid_sz(1) + 2 * pml_size(1)))) || ...
   (smallest_factor_y ~= max(factor( 2 * (grid_sz(2) + pml_size(2)) - 2 )))
    test_pass = false;
end

% -----------------
% 3D
% -----------------

% generate random grid size
grid_sz = round(max_grid_size * rand(3, 1));

% get optimal pml size
pml_size = getOptimalPMLSize(grid_sz, pml_range);

% check the smallest prime factor in the range
smallest_factor_x = max_grid_size;
smallest_factor_y = max_grid_size;
smallest_factor_z = max_grid_size;
for ind = 1:length(pml_vector)
    smallest_factor_x = min(smallest_factor_x, max(factor(grid_sz(1) + 2 * pml_vector(ind))));
    smallest_factor_y = min(smallest_factor_y, max(factor(grid_sz(2) + 2 * pml_vector(ind))));
    smallest_factor_z = min(smallest_factor_z, max(factor(grid_sz(3) + 2 * pml_vector(ind))));
end

% check it matches
if (smallest_factor_x ~= max(factor(grid_sz(1) + 2 * pml_size(1)))) || ...
   (smallest_factor_y ~= max(factor(grid_sz(2) + 2 * pml_size(2)))) || ...
   (smallest_factor_z ~= max(factor(grid_sz(3) + 2 * pml_size(3))))
    test_pass = false;
end