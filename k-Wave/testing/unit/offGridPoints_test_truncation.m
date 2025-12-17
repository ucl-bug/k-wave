function test_pass = offGridPoints_test_truncation(plot_comparisons, ~)
%OFFGRIDPOINTS_TEST_TRUNCATION Check truncated sinc approximation.
%
% DESCRIPTION:
%     offGridPoints_test_truncation creates an off-grid point source using 
%     a sinc approximation cutoff using a number of different thresholds.
%
% INPUTS:
%     plot_comparisons - Boolean controlling whether any comparisons (e.g.,
%                        error norms) are plotted 
%
% OUTPUTS:
%     test_pass        - Boolean denoting whether the test passed
%
% ABOUT:
%     author           - Elliott Wise and Bradley Treeby
%     date             - 31st May 2018
%     last update      - 8th June 2021
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018- Elliott Wise and Bradley Treeby

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

% specify the computational grid (large enough to encompass smalled
% tolerance value)
Nx = 512;
dx = 1;
Ny = Nx;
dy = dx;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% create list of tolerance values
tolerance_array = logspace(-1, -3, 16).';

% unit scale to avoid needing to multiply tol by this
point = [dx/2; dy/2]; % maximally off-grid
scale = 1;

% generate un-truncated point-source grid weights
source_ref = offGridPoints(kgrid, point, scale, 'BLIType', 'exact');

% create an empty array for the source grid weight differences
source_diff = NaN(size(tolerance_array));

% loop through tolerances
for tol_index = 1:length(tolerance_array)
    
    % get current tolerance
    tol = tolerance_array(tol_index);
    
    % generate truncated point-source grid weights
    source_trunc = offGridPoints(kgrid, point, scale, 'BLITolerance', tol);
    
    % compute the difference between the truncated and reference source
    % grid weights
    source_diff(tol_index) = max(abs(source_trunc(:) - source_ref(:)));
    
end

% check if any differences in the source weights exceed the tolerance
if any(source_diff > tolerance_array)
    test_pass = false;
end

% plot
if plot_comparisons
    figure;
    plot(tolerance_array, source_diff, '.');
    hold on;
    plot(tolerance_array, tolerance_array, 'k--');
    set(gca, 'XScale', 'log', 'YScale', 'log', 'XDir', 'reverse');
    xlim([min(tolerance_array), max(tolerance_array)]);
    xlabel('Tolerance');
    ylabel('Error');
    legend('Difference in source', 'Tolerance');
end