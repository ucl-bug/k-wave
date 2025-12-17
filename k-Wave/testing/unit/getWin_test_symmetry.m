function test_pass = getWin_test_symmetry(~, ~)
% DESCRIPTION:
%     Unit test to make sure the returned windows have the correct
%     symmetry.
%
% ABOUT:
%     author           - Bradley Treeby
%     date             - 25th July 2019
%     last update      - 25th July 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2019- Bradley Treeby

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
zero_val = 1e-15;

% set grid sizes
N_odd = 7;
N_even = 8;

% =========================================================================
% 1D
% =========================================================================

% create a symmetric, odd length window
win = getWin(N_odd, 'Blackman', 'Symmetric', true);

% check it contains 1 in the middle
if abs(win((N_odd + 1)/2) - 1) > zero_val
    test_pass = false;
end

% check values are symmetric
first_half = win(1:(N_odd + 1)/2 - 1);
second_half = win((N_odd + 1)/2 + 1:end);
if max(abs(first_half(:) - flip(second_half(:)))) > zero_val
    test_pass = false;
end

% ----

% create a non-symmetric, odd length window
win = getWin(N_odd, 'Blackman', 'Symmetric', false);

% check values are symmetric
first_half = win(2:(N_odd + 1)/2);
second_half = win((N_odd + 1)/2 + 1:end);
if max(abs(first_half(:) - flip(second_half(:)))) > zero_val
    test_pass = false;
end

% ----

% create a symmetric, even length window
win = getWin(N_even, 'Blackman', 'Symmetric', true);


% check values are symmetric
first_half = win(1:N_even/2);
second_half = win(N_even/2 + 1:end);
if max(abs(first_half(:) - flip(second_half(:)))) > zero_val
    test_pass = false;
end

% ----

% create a non-symmetric, even length window
win = getWin(N_even, 'Blackman', 'Symmetric', false);

% check it contains 1 in the middle
if abs(win(N_even/2 + 1) - 1) > zero_val
    test_pass = false;
end

% check values are symmetric
first_half = win(2:N_even/2);
second_half = win(N_even/2 + 2:end);
if max(abs(first_half(:) - flip(second_half(:)))) > zero_val
    test_pass = false;
end

% =========================================================================
% 2D
% =========================================================================

% create a [sym, sym], [odd, odd] window
win = getWin([N_odd, N_odd], 'Blackman', 'Symmetric', [true, true]);
win1D = getWin(N_odd, 'Blackman', 'Symmetric', true);

% check it contains 1 in the middle
if abs(win((N_odd + 1)/2, (N_odd + 1)/2) - 1) > zero_val
    test_pass = false;
end

% check the central row is the same as the 1D reference
if max(abs(win((N_odd + 1)/2, :) - win1D.')) > zero_val
    test_pass = false;
end

% check the central row is the same as the 1D reference
if max(abs(win(:, (N_odd + 1)/2) - win1D)) > zero_val
    test_pass = false;
end

% ----

% create a [asym, asym], [even, even] window
win = getWin([N_even, N_even], 'Blackman', 'Symmetric', [false, false]);
win1D = getWin(N_even, 'Blackman', 'Symmetric', false);

% check it contains 1 in the middle
if abs(win(N_even/2 + 1, N_even/2 + 1) - 1) > zero_val
    test_pass = false;
end

% check the central row is the same as the 1D reference
if max(abs(win(N_even/2 + 1, :) - win1D.')) > zero_val
    test_pass = false;
end

% check the central row is the same as the 1D reference
if max(abs(win(:, N_even/2 + 1) - win1D)) > zero_val
    test_pass = false;
end

% ----

% create a [sym, asym], [odd, even] window
win = getWin([N_odd, N_even], 'Blackman', 'Symmetric', [true, false]);
win1D = getWin(N_even, 'Blackman', 'Symmetric', false);

% check it contains 1 in the middle
if abs(win((N_odd + 1)/2, N_even/2 + 1) - 1) > zero_val
    test_pass = false;
end

% check the central row is the same as the 1D reference
if max(abs(win((N_odd + 1)/2, :) - win1D.')) > zero_val
    test_pass = false;
end

% ----

% create a [asym, sym], [even, odd] window
win = getWin([N_even, N_odd], 'Blackman', 'Symmetric', [false, true]);
win1D = getWin(N_even, 'Blackman', 'Symmetric', false);

% check it contains 1 in the middle
if abs(win(N_even/2 + 1, (N_odd + 1)/2) - 1) > zero_val
    test_pass = false;
end

% check the central column is the same as the 1D reference
if max(abs(win(:, (N_odd + 1)/2) - win1D)) > zero_val
    test_pass = false;
end

% =========================================================================
% 3D
% =========================================================================

% create a [sym, sym, sym], [odd, odd, odd] window
win = getWin([N_odd, N_odd, N_odd], 'Blackman', 'Symmetric', [true, true, true]);
win1D = getWin(N_odd, 'Blackman', 'Symmetric', true);

% check it contains 1 in the middle
if abs(win((N_odd + 1)/2, (N_odd + 1)/2, (N_odd + 1)/2) - 1) > zero_val
    test_pass = false;
end

% check centre of dim 3 is the same as the 1D reference
if max(abs(win((N_odd + 1)/2, (N_odd + 1)/2, :) - reshape(win1D, 1, 1, []))) > zero_val
    test_pass = false;
end

% check centre of dim 2 is the same as the 1D reference
if max(abs(win((N_odd + 1)/2, :, (N_odd + 1)/2) - reshape(win1D, 1, [], 1))) > zero_val
    test_pass = false;
end

% check centre of dim 1 is the same as the 1D reference
if max(abs(win(:, (N_odd + 1)/2, (N_odd + 1)/2) - reshape(win1D, [], 1, 1))) > zero_val
    test_pass = false;
end

% ----

% create a [asym, asym, asym], [even, even, even] window
win = getWin([N_even, N_even, N_even], 'Blackman', 'Symmetric', [false, false, false]);
win1D = getWin(N_even, 'Blackman', 'Symmetric', false);

% check it contains 1 in the middle
if abs(win(N_even/2 + 1, N_even/2 + 1, N_even/2 + 1) - 1) > zero_val
    test_pass = false;
end

% check centre of dim 3 is the same as the 1D reference
if max(abs(win(N_even/2 + 1, N_even/2 + 1, :) - reshape(win1D, 1, 1, []))) > zero_val
    test_pass = false;
end

% check centre of dim 2 is the same as the 1D reference
if max(abs(win(N_even/2 + 1, :, N_even/2 + 1) - reshape(win1D, 1, [], 1))) > zero_val
    test_pass = false;
end

% check centre of dim 1 is the same as the 1D reference
if max(abs(win(:, N_even/2 + 1, N_even/2 + 1) - reshape(win1D, [], 1, 1))) > zero_val
    test_pass = false;
end

% ----

% create a [sym, asym, asym], [odd, even, even] window
win = getWin([N_odd, N_even, N_even], 'Blackman', 'Symmetric', [true, false, false]);
win1D = getWin(N_even, 'Blackman', 'Symmetric', false);

% check it contains 1 in the middle
if abs(win((N_odd + 1)/2, N_even/2 + 1, N_even/2 + 1) - 1) > zero_val
    test_pass = false;
end

% check centre of longest dim is the same as the 1D reference
if max(abs(win((N_odd + 1)/2, N_even/2 + 1, :) - reshape(win1D, 1, 1, []))) > zero_val
    test_pass = false;
end

% ----

% create a [sym, asym, sym], [odd, even, odd] window
win = getWin([N_odd, N_even, N_odd], 'Blackman', 'Symmetric', [true, false, true]);
win1D = getWin(N_even, 'Blackman', 'Symmetric', false);

% check it contains 1 in the middle
if abs(win((N_odd + 1)/2, N_even/2 + 1, (N_odd + 1)/2) - 1) > zero_val
    test_pass = false;
end

% check centre of longest dim is the same as the 1D reference
if max(abs(win((N_odd + 1)/2, :, (N_odd + 1)/2) - reshape(win1D, 1, [], 1))) > zero_val
    test_pass = false;
end

% ----

% create a [sym, sym, asym], [odd, odd, even] window
win = getWin([N_odd, N_odd, N_even], 'Blackman', 'Symmetric', [true, true, false]);
win1D = getWin(N_even, 'Blackman', 'Symmetric', false);

% check it contains 1 in the middle
if abs(win((N_odd + 1)/2, (N_odd + 1)/2, N_even/2 + 1) - 1) > zero_val
    test_pass = false;
end

% check centre of longest dim is the same as the 1D reference
if max(abs(win((N_odd + 1)/2, (N_odd + 1)/2, :) - reshape(win1D, 1, 1, []))) > zero_val
    test_pass = false;
end

% ----

% create a [asym, sym, sym], [even, odd, odd] window
win = getWin([N_even, N_odd, N_odd], 'Blackman', 'Symmetric', [false, true, true]);
win1D = getWin(N_even, 'Blackman', 'Symmetric', false);

% check it contains 1 in the middle
if abs(win(N_even/2 + 1, (N_odd + 1)/2, (N_odd + 1)/2) - 1) > zero_val
    test_pass = false;
end

% check centre of longest dim is the same as the 1D reference
if max(abs(win(:, (N_odd + 1)/2, (N_odd + 1)/2) - reshape(win1D, [], 1, 1))) > zero_val
    test_pass = false;
end

% ----

% create a [asym, asym, sym], [even, even, odd] window
win = getWin([N_even, N_even, N_odd], 'Blackman', 'Symmetric', [false, false, true]);
win1D = getWin(N_even, 'Blackman', 'Symmetric', false);

% check it contains 1 in the middle
if abs(win(N_even/2 + 1, N_even/2 + 1, (N_odd + 1)/2) - 1) > zero_val
    test_pass = false;
end

% check centre of longest dim is the same as the 1D reference
if max(abs(win(:, N_even/2 + 1, (N_odd + 1)/2) - reshape(win1D, [], 1, 1))) > zero_val
    test_pass = false;
end

% ----

% create a [asym, sym, asym], [even, odd, even] window
win = getWin([N_even, N_odd, N_even], 'Blackman', 'Symmetric', [false, true, false]);
win1D = getWin(N_even, 'Blackman', 'Symmetric', false);

% check it contains 1 in the middle
if abs(win(N_even/2 + 1, (N_odd + 1)/2, N_even/2 + 1) - 1) > zero_val
    test_pass = false;
end

% check centre of longest dim is the same as the 1D reference
if max(abs(win(:, (N_odd + 1)/2, N_even/2 + 1) - reshape(win1D, [], 1, 1))) > zero_val
    test_pass = false;
end