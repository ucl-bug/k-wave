function test_pass = smooth_test_delta(~, ~)
% DESCRIPTION:
%       Unit test to check using a delta function input returns the FFT of
%       the window when using smooth.
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

% window type
win_type = 'Blackman';
restore_max = false;

% ------------------------------
% 1D (odd, 1)
% ------------------------------

% generate input and smooth
func = zeros(odd_size, 1);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(odd_size, win_type, 'Symmetric', true);

% take fft, normalise and shift
func_sm = fftshift(fft(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 1D (1, odd)
% ------------------------------

% generate input and smooth
func = zeros(1, odd_size);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(odd_size, win_type, 'Symmetric', true).';

% take fft, normalise and shift
func_sm = fftshift(fft(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 1D (even, 1)
% ------------------------------

% generate input and smooth
func = zeros(even_size, 1);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(even_size, win_type, 'Symmetric', false);

% take fft, normalise and shift
func_sm = fftshift(fft(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 1D (1, even)
% ------------------------------

% generate input and smooth
func = zeros(1, even_size);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(even_size, win_type, 'Symmetric', false);

% take fft, normalise and shift
func_sm = fftshift(fft(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 2D (odd, odd)
% ------------------------------

% generate input and smooth
sz = [odd_size, odd_size];
sym = [true, true];
func = zeros(sz);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(sz, win_type, 'Rotation', true, 'Symmetric', sym);

% take fft, normalise and shift
func_sm = fftshift(fft2(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 2D (odd, even)
% ------------------------------

% generate input and smooth
sz = [odd_size, even_size];
sym = [true, false];
func = zeros(sz);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(sz, win_type, 'Rotation', true, 'Symmetric', sym);

% take fft, normalise and shift
func_sm = fftshift(fft2(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 2D (even, odd)
% ------------------------------

% generate input and smooth
sz = [even_size, odd_size];
sym = [false, true];
func = zeros(sz);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(sz, win_type, 'Rotation', true, 'Symmetric', sym);

% take fft, normalise and shift
func_sm = fftshift(fft2(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 2D (even, even)
% ------------------------------

% generate input and smooth
sz = [even_size, even_size];
sym = [false, false];
func = zeros(sz);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(sz, win_type, 'Rotation', true, 'Symmetric', sym);

% take fft, normalise and shift
func_sm = fftshift(fft2(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (odd, odd, odd)
% ------------------------------

% generate input and smooth
sz = [odd_size, odd_size, odd_size];
sym = [true, true, true];
func = zeros(sz);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(sz, win_type, 'Rotation', true, 'Symmetric', sym);

% take fft, normalise and shift
func_sm = fftshift(fftn(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (odd, odd, even)
% ------------------------------

% generate input and smooth
sz = [odd_size, odd_size, even_size];
sym = [true, true, false];
func = zeros(sz);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(sz, win_type, 'Rotation', true, 'Symmetric', sym);

% take fft, normalise and shift
func_sm = fftshift(fftn(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (odd, even, odd)
% ------------------------------

% generate input and smooth
sz = [odd_size, even_size, odd_size];
sym = [true, false, true];
func = zeros(sz);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(sz, win_type, 'Rotation', true, 'Symmetric', sym);

% take fft, normalise and shift
func_sm = fftshift(fftn(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (even, odd, odd)
% ------------------------------

% generate input and smooth
sz = [even_size, odd_size, odd_size];
sym = [false, true, true];
func = zeros(sz);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(sz, win_type, 'Rotation', true, 'Symmetric', sym);

% take fft, normalise and shift
func_sm = fftshift(fftn(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (odd, even, even)
% ------------------------------

% generate input and smooth
sz = [odd_size, even_size, even_size];
sym = [true, false, false];
func = zeros(sz);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(sz, win_type, 'Rotation', true, 'Symmetric', sym);

% take fft, normalise and shift
func_sm = fftshift(fftn(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (even, odd, even)
% ------------------------------

% generate input and smooth
sz = [even_size, odd_size, even_size];
sym = [false, true, false];
func = zeros(sz);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(sz, win_type, 'Rotation', true, 'Symmetric', sym);

% take fft, normalise and shift
func_sm = fftshift(fftn(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (even, even, odd)
% ------------------------------

% generate input and smooth
sz = [even_size, even_size, odd_size];
sym = [false, false, true];
func = zeros(sz);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(sz, win_type, 'Rotation', true, 'Symmetric', sym);

% take fft, normalise and shift
func_sm = fftshift(fftn(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end

% ------------------------------
% 3D (even, even, even)
% ------------------------------

% generate input and smooth
sz = [even_size, even_size, even_size];
sym = [false, false, false];
func = zeros(sz);
func(1) = 1;
func_sm = smooth(func, restore_max, win_type);

% get the corresponding window
win = getWin(sz, win_type, 'Rotation', true, 'Symmetric', sym);

% take fft, normalise and shift
func_sm = fftshift(fftn(func_sm));

% check values
err = max(abs(win(:) - func_sm(:)));
if err > comparison_threshold
    test_pass = false;
end