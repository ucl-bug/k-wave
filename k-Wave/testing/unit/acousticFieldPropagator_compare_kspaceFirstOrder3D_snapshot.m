function test_pass = acousticFieldPropagator_compare_kspaceFirstOrder3D_snapshot(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to compare acousticFieldPropagator with a k-Wave simulation
%     using kspaceFirstOrder3D when calculating a snapshot of the field
%     (rather than amplitude and phase).
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 8th January 2019
%     last update - 8th January 2019
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

% check for plot inputs, and set to true if nargin is zero (to allow the
% test to be run independent of runUnitTests)
if nargin == 0
    plot_comparisons = true;
    plot_simulations = true;
end

% set comparison threshold
comparison_threshold = 0.1;

% set pass variable
test_pass = true;

% =========================================================================
% AFP SIMULATION
% =========================================================================

% define grid size
Nx = 64;
Ny = 64;
Nz = 64;
dx = 0.1e-3;
dy = 0.1e-3;
dz = 0.1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% define apeture
aperture_width = 10;
aperture = zeros(Nx, Ny, Nz);
aperture(Nx/2 - aperture_width + 1: Nx/2 + aperture_width, Ny/2 - aperture_width + 1: Ny/2 + aperture_width, Nz/2) = 1;

% define source frequency and sound speed
f_Hz = 1.5e6;
c0 = 1500;
w0 = 2 * pi * f_Hz;

% choose time
period = 1 / f_Hz;
t_end = 1.5 * period;    

% compute field using AFP
pressure = acousticFieldPropagator(aperture, 0, dx, f_Hz, c0, ...
    'Time', t_end, ...
    'UseRamp', false, ...
    'GridExpansionFactor', 1);

% taking the imaginary part (gives sin source)
beam_pattern_afp = imag(pressure);

% =========================================================================
% K-WAVE SIMULATION
% =========================================================================

% create time array
dt = period / 20;
Nt = round(t_end / dt);
dt = t_end / Nt;
kgrid.t_array = 0:dt:t_end;

% assign medium properties
medium.sound_speed = c0;

% create source
source.p_mask = aperture;
source.p = sin(w0*(kgrid.t_array - kgrid.dt/2));
source.p_mode = 'additive';

% create sensor mask
sensor.mask = zeros(Nx, Ny, Nz);
sensor.record = {'p_final'};

% run simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
    'PMLAlpha', 0, ...
    'DisplayMask', 'off', ...
    'PlotSim', plot_simulations);

% assign output
beam_pattern_kwave = sensor_data.p_final;

% =========================================================================
% ERROR
% =========================================================================

% extract just the middle plane
beam_pattern_afp = squeeze(beam_pattern_afp(:, :, Nz/2));
beam_pattern_kwave = squeeze(beam_pattern_kwave(:, :, Nz/2));

% compute error
err = max(abs(beam_pattern_afp(:) - beam_pattern_kwave(:))) / max(abs(beam_pattern_afp(:)));

% check for test pass
if (err > comparison_threshold)
    test_pass = false;
end

% =========================================================================
% PLOTTING
% =========================================================================

% plot comparison
if plot_comparisons

    figure;
    
    subplot(2, 3, 1)
    mx_kwave = max(abs(beam_pattern_kwave(:)));
    imagesc(beam_pattern_kwave, [-1, 1]*mx_kwave);
    colormap(getColorMap);
    colorbar;
    axis image;
    title('k-Wave');

    subplot(2, 3, 2)
    mx_fast = max(abs(beam_pattern_afp(:)));
    imagesc(beam_pattern_afp, [-1, 1]*mx_fast);
    colormap(getColorMap);
    colorbar;
    axis image;
    title('AFP');

    subplot(2, 3, 3)
    imagesc(abs(beam_pattern_kwave - beam_pattern_afp));
    colormap(getColorMap);
    colorbar;
    axis image;
    title('Difference');

    subplot(2, 3, 4)
    plot(beam_pattern_kwave(end/2, :), 'k-');
    hold on;
    plot(beam_pattern_afp(end/2, :), 'r:');
    legend('k-Wave', 'AFP');
    title('Profile');

    subplot(2, 3, 5)
    plot(abs(beam_pattern_kwave(end/2, :) - beam_pattern_afp(end/2, :)), 'k-');
    title('Difference in profile');
    
    scaleFig(1.5, 1);
    
end