function test_pass = smooth_test_DC(~, ~)
% DESCRIPTION:
%       Unit test to check DC component is maintained when using smooth.
%
% ABOUT:
%       author      - Bradley Treeby
%       date        - 29th July 2019
%       last update - 29th July 2019
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019 Bradley Treeby

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

% comparison threshold
comparison_threshold = 1e-15;

% set size
odd_size = 15;
even_size = 16;

% generate a random DC value
dc_val = rand(1);

% ------------------------------
% 1D (odd, 1)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(odd_size, 1);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 1D (1, odd)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(1, odd_size);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 1D (even, 1)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(even_size, 1);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 1D (1, even)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(1, even_size);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 2D (odd, odd)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(odd_size, odd_size);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 2D (odd, even)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(odd_size, even_size);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 2D (even, odd)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(even_size, odd_size);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 2D (even, even)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(even_size, even_size);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (odd, odd, odd)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(odd_size, odd_size, odd_size);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (odd, odd, even)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(odd_size, odd_size, even_size);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (odd, even, odd)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(odd_size, even_size, odd_size);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (even, odd, odd)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(even_size, odd_size, odd_size);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (odd, even, even)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(odd_size, even_size, even_size);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (even, odd, even)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(even_size, odd_size, even_size);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (even, even, odd)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(even_size, even_size, odd_size);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (even, even, even)
% ------------------------------

% generate input and smooth
func = dc_val .* ones(even_size, even_size, even_size);
func_sm = smooth(func);

% check values
err = max(abs(func(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end