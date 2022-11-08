function test_pass = kWaveArray_test_single_precision(~, plot_comparisons)
%KWAVEARRAY_TEST_SINGLE_PRECISION Compare results for single and double precision.
%
% DESCRIPTION:
%     Compares mask and grid weights for large bowl in single and double
%     precision.
%
% INPUTS:
%     plot_comparisons - Boolean controlling whether any comparisons (e.g.,
%                        error norms) are plotted.
%
% OUTPUTS:
%     test_pass        - Boolean denoting whether the test passed.
%
% ABOUT:
%     author           - Bradley Treeby
%     date             - 14th July 2021
%     last update      - 14th July 2021
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2021 Bradley Treeby

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

% set comparison threshold
comparison_thresh = 1e-5;

% set pass variable
test_pass = true;

% test literals
source_f0          = 500e3;    % source frequency [Hz]
bowl_roc           = 64e-3;    % bowl radius of curvature [m]
bowl_diameter      = 64e-3;    % bowl aperture diameter [m]
bowl_pos           = [0, 0, 0];
focus_pos          = [1, 0, 0];

bli_tolerance      = 0.05;
upsampling_rate    = 8;
cfl                = 0.1;
Nt                 = 500;

% define grid
sc = 2;
Nx = 96 * sc;
Ny = 256 * sc;
Nz = 256 * sc;
dx = 0.5e-3 / sc;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% define input signals
dt = cfl * dx / 1500;
t_array = (0:Nt - 1) * dt;
input_signal = createCWSignals(t_array, source_f0, 1, 0);

% create empty kWaveArray in single precision
karray_single = kWaveArray(...
    'BLITolerance', bli_tolerance, ...
    'UpsamplingRate', upsampling_rate, ...
    'SinglePrecision', true);

% create empty kWaveArray in double precision
karray_double = kWaveArray(...
    'BLITolerance', bli_tolerance, ...
    'UpsamplingRate', upsampling_rate, ...
    'SinglePrecision', false);

% add bowl shaped element
karray_single.addBowlElement(bowl_pos, bowl_roc, bowl_diameter, focus_pos);
karray_double.addBowlElement(bowl_pos, bowl_roc, bowl_diameter, focus_pos);

% get masks
tic; fprintf('Computing masks... ');
karray_single_mask      = karray_single.getArrayBinaryMask(kgrid);
karray_double_mask      = karray_double.getArrayBinaryMask(kgrid);
toc;

% compare masks
mask_err = sum(abs(karray_single_mask(:) - karray_double_mask(:)));
if mask_err > 0
    test_pass = false;
end

% store plot masks and clear memory
if plot_comparisons
    plot_mask_single = karray_single_mask(:, :, ceil(end/2));
    plot_mask_double = karray_double_mask(:, :, ceil(end/2));
end
clear karray_single_mask karray_double_mask

% get weights
tic; fprintf('Computing weights... ');
karray_single_weights   = karray_single.getArrayGridWeights(kgrid);
karray_double_weights   = karray_double.getArrayGridWeights(kgrid);
toc
    
% compare weights
weights_err = max(abs(karray_single_weights(:) - karray_double_weights(:)));
if weights_err > comparison_thresh
    test_pass = false;
end

% store plot masks and clear memory
if plot_comparisons
    plot_weights_single = karray_single_weights(:, :, ceil(end/2));
    plot_weights_double = karray_double_weights(:, :, ceil(end/2));
end
clear karray_single_weights karray_double_weights

% get source input
karray_single_p         = karray_single.getDistributedSourceSignal(kgrid, input_signal);
karray_double_p         = karray_double.getDistributedSourceSignal(kgrid, input_signal);

% compare signals
signals_err = max(abs(karray_single_p(:) - karray_double_p(:)));
if signals_err > comparison_thresh
    test_pass = false;
end

if plot_comparisons
    
    figure;
    
    subplot(2, 3, 1);
    imagesc(plot_mask_single);
    axis image;
    title('Mask - Single');
    colorbar;
    
    subplot(2, 3, 2);
    imagesc(plot_mask_double);
    axis image;
    title('Mask - Double');
    colorbar;
    
    subplot(2, 3, 3);
    imagesc(abs(plot_mask_single - plot_mask_double));
    axis image;
    title('Mask - Difference');
    colorbar;
    
    subplot(2, 3, 4);
    imagesc(plot_weights_single);
    axis image;
    title('Weights - Single');
    colorbar;
    
    subplot(2, 3, 5);
    imagesc(plot_weights_double);
    axis image;
    title('Weights - Double');
    colorbar;
    
    subplot(2, 3, 6);
    imagesc(abs(plot_weights_single - plot_weights_double));
    axis image;
    title('Weights - Difference');
    colorbar;
    
end