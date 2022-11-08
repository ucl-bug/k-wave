function test_pass = fourierShift_test_cos_shift(plot_comparisons, ~)
%FOURIERSHIFT_TEST_COS_SHIFT Test fourierShift using cosine wave.
%
% DESCRIPTION:
%     fourierShift_test_sin_shift checks the output of fourierShift using a
%     periodic cosine wave.
%
% INPUTS:
%     plot_comparisons - Boolean controlling whether any comparisons (e.g.,
%                        error norms) are plotted 
%
% OUTPUTS:
%     test_pass        - Boolean denoting whether the test passed
%
% ABOUT:
%     author           - Bradley Treeby
%     date             - 19th February 2017
%     last update      - 19th February 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2017 Bradley Treeby

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

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
end

% set pass variable
test_pass = true;

% set comparison threshold
COMPARISON_THRESH   = 1e-13;

% specify periodic cosine
dt      = pi/5;
t       = 0:dt:(2*pi - dt);
y       = cos(t);

% shift, and calculate reference
y_shift = fourierShift(y, 1/2);
y_ref   = cos(t + dt/2);

% compare values
if max(abs(y_ref(:) - y_shift(:))) > COMPARISON_THRESH
    test_pass = false;
end

% plot
if plot_comparisons
    figure;
    plot(t, y, 'k-s', t, y_shift, 'b-s', t, y_ref, 'rx');
    legend('cos', 'shifted cos', 'reference');
end