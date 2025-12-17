function test_pass = acousticFieldPropagator_compare_kspaceFirstOrder2D(plot_comparisons, ~)
% DESCRIPTION:
%     Unit test to compare acousticFieldPropagator with a k-Wave simulation
%     using kspaceFirstOrder2D.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 8th January 2019
%     last update - 8th January 2019
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

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
end

% set comparison threshold
comparison_threshold = 5e-4;

% set pass variable
test_pass = true;

% =========================================================================
% RUN TEST
% =========================================================================

% define grid size
PML_size = 20;
Nx = 256 - 2*PML_size;
Ny = 128 - 2*PML_size;
dx = 0.1e-3;

% define apeture
aperture_width = 10;
aperture = zeros(Nx, Ny);
aperture(1, Ny/2 - aperture_width + 1: Ny/2 + aperture_width, 1) = 1;

% define source frequency and sound speed
f = 1e6;
c0 = 1500;

% run AFP simulation
[amplitude_afp, phase_afp] = acousticFieldPropagator(aperture, 0, dx, f, c0);
    
% run k-Wave simulation
[amplitude_kw, phase_kw] = acousticFieldPropagator(aperture, 0, dx, f, c0, 'UseFirstOrder', true);
    
% calculate phase different, accounting for constant offset
phase_diff = phase_kw - phase_afp;
phase_diff(phase_diff < 0) = phase_diff(phase_diff < 0) + 2*pi;
phase_diff = phase_diff - min(phase_diff(:));

% compute error
err_amp = max(abs(amplitude_afp(:) - amplitude_kw(:))) / max(abs(amplitude_afp(:)));
err_phase = max(phase_diff(:)) / (2*pi);

% check for test pass
if (err_amp > comparison_threshold) || (err_phase > comparison_threshold)
    test_pass = false;
end

% =========================================================================
% PLOTTING
% =========================================================================

% plot comparison
if plot_comparisons

    figure;

    subplot(2, 3, 1);
    imagesc(amplitude_afp);
    axis image;
    title('AFP - Amplitude');
    colorbar

    subplot(2, 3, 2);
    imagesc(amplitude_kw);
    axis image
    title('k-Wave - Amplitude');
    colorbar

    subplot(2, 3, 3);
    imagesc(abs(amplitude_kw - amplitude_afp) ./ max(abs(amplitude_afp(:))));
    axis image
    title('Difference - Amplitude')
    colorbar

    subplot(2, 3, 4);
    imagesc(phase_afp);
    axis image;
    title('AFP -  Phase');
    colorbar

    subplot(2, 3, 5);
    imagesc(phase_kw);
    axis image;
    title('k-Wave -  Phase');
    colorbar

    subplot(2, 3, 6);
    imagesc(phase_diff);
    axis image
    title('Difference - Phase');
    colorbar

    colormap(getColorMap);
    
end