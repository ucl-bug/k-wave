% Filtering A Delta Function Input Signal Example Part 3
%
% This example illustrates how to temporally filter an input signal.
%
% author: Bradley Treeby
% date: 19th January 2010
% last update: 4th June 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2017 Bradley Treeby

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

clearvars;

% modify this parameter to run the different examples
example_number = 1;
% 1: default causal filter
% 2: zero phase filter
% 3: causal filter with modified input settings

% =========================================================================
% FILTERING
% =========================================================================

% create the computational grid
Nx = 256;                   % [grid points]
dx = 10e-3 / Nx;            % [m]
kgrid = kWaveGrid(Nx, dx);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]

% create time array
dt = 7e-9;
Nt = 1024;
kgrid.setTime(Nt, dt);

% define a delta function input pulse
temporal_offset = 100;      % [time steps]
source_magnitude = 2;       % [Pa]
source_func = zeros(size(kgrid.t_array));
source_func(temporal_offset) = source_magnitude;

% filter the input signal
switch example_number
    case 1
        source_func_filtered = filterTimeSeries(kgrid, medium, source_func);
    case 2
        source_func_filtered = filterTimeSeries(kgrid, medium, source_func, 'ZeroPhase', true);
    case 3
        source_func_filtered = filterTimeSeries(kgrid, medium, source_func, 'PPW', 4, 'TransitionWidth', 0.05);
end

% compute the amplitude spectra of the original and filtered input signals
[f, source_func_as] = spect(source_func, 1/dt);
[~, source_func_filtered_as] = spect(source_func_filtered, 1/dt);

% extract the maximum frequency supported by the grid (two points per
% wavelength)
f_max = kgrid.k_max * min(medium.sound_speed(:)) / (2*pi);

% =========================================================================
% SIMULATION
% =========================================================================

% define a single element source
source_offset = 50;
source.p_mask = zeros(Nx, 1);
source.p_mask(1 + source_offset, 1) = 1;

% assigned filtered input signal
source.p = source_func_filtered;

% define a single element sensor
sensor.mask = zeros(Nx, 1);
sensor.mask(end - source_offset, 1) = 1;

% increase the size of the grid using the PML, and turn the PML off
input_args = {'PMLInside', false, 'PMLSize', 128, 'PMLAlpha', 0, 'PlotPML', false};

% run the simulation
sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});

% compute the amplitude spectra of the recorded time series
[~, output_as] = spect(sensor_data, 1/dt);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot original and filtered input signals
figure;
[t_sc, scale, prefix] = scaleSI(max(kgrid.t_array(:)));
subplot(2, 1, 1);
plot(kgrid.t_array * scale, source_func, 'k-', ...
     kgrid.t_array * scale, source_func_filtered, 'b-');
xlabel(['Time [' prefix 's]']);
ylabel('Pressure [au]');
legend('Original input', 'Filtered input');
set(gca, 'XLim', [0, 3], 'YLim', [-0.1, 2.1]);

% plot the amplitude spectra
[f_sc, scale, prefix] = scaleSI(max(f));
subplot(2, 1, 2);
plot(f * scale, source_func_as, 'k-', ...
     f * scale, source_func_filtered_as, 'b-');
xlabel(['Frequency [' prefix 'Hz]']);
ylabel('Amplitude [au]');

% plot the maximum frequency supported by the grid
ylim = get(gca, 'YLim');
hold on;
subplot(2, 1, 2), line([f_max * scale, f_max * scale], [0, ylim(2)], 'LineStyle', '--', 'Color', 'k');
set(gca, 'XLim', [0, 50]);

% plot the input and recorded time series
figure;
[~, scale, prefix] = scaleSI(max(kgrid.t_array(:)));
plot(kgrid.t_array * scale, source.p, 'k-', ...
     kgrid.t_array * scale, sensor_data, 'b-');
xlabel(['Time [' prefix 's]']);
ylabel('Pressure [au]');
legend('Input pulse', 'Recorded pulse');

% plot the amplitude spectra
[~, scale, prefix] = scaleSI(max(f));
figure;
plot(f * scale, source_func_filtered_as, 'k-', ...
     f * scale, output_as, 'b-');
xlabel(['Frequency [' prefix 'Hz]']);
ylabel('Amplitude [au]');

% plot the maximum frequency supported by the grid
ylim = get(gca, 'YLim');
hold on;
line([f_max * scale, f_max * scale], [0, ylim(2)], 'LineStyle','--', 'Color', 'k');
legend('Amplitude spectrum of input pulse', 'Amplitude spectrum of recorded pulse', ...
    'Maximum frequency supported by grid', 'Location', 'SouthEast');