% Defining An Ultrasound Transducer Example
%
% In principle, simulations using diagnostic ultrasound transducers can be
% run following the examples given under Time Varying Source Problems.
% However, assigning the grid points that belong to each transducer
% element, and then assigning the correctly delayed input signals to each
% point of each element can be a laborious task. For this purpose, a
% special input object created using the kWaveTransducer class can be
% substituted for either the source or sensor inputs (or both). This
% example illustrates how this object is created and can be used to
% simulate the field produced by an ultrasound transducer.
%
% Note, transducer inputs can only be used in 3D simulations and thus these
% examples are inherently memory and CPU intensive. Whilst the grid sizes
% and source frequencies used in the examples have been scaled down for the
% purpose of demonstrating the capabilities of the toolbox (the inputs do
% not necessarily represent realistic ultrasound settings), they still
% require a comparatively large amount of computational resources. To
% reduce this load, it is advised to run the simulations in single
% precision by setting the optional input 'DataCast' to 'single'.
% Similarly, if you have access to a recent model GPU and the MATLAB
% Parallel Computing Toolbox R2012a or later, the simulation times can be
% significantly reduced by setting 'DataCast' to 'gpuArray-single'.
% Alternatively, the simulations can be run using the optimised C++ code.
% See the k-Wave Manual for more information. 
%
% The creation of a kWaveTransducer object will only work in versions of
% MATLAB recent enough to support user defined classes. 
%
% author: Bradley Treeby
% date: 20th July 2011
% last update: 11th June 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2011-2017 Bradley Treeby

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

% simulation settings
DATA_CAST = 'single';

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 20;            % [grid points]
PML_Y_SIZE = 10;            % [grid points]
PML_Z_SIZE = 10;            % [grid points]

% set total number of grid points not including the PML
Nx = 128 - 2*PML_X_SIZE;    % [grid points]
Ny = 128 - 2*PML_Y_SIZE;    % [grid points]
Nz = 64 - 2*PML_Z_SIZE;     % [grid points]

% set desired grid size in the x-direction not including the PML
x = 40e-3;                  % [m]

% calculate the spacing between the grid points
dx = x/Nx;                  % [m]
dy = dx;                    % [m]
dz = dx;                    % [m]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
medium.sound_speed = 1540;      % [m/s]
medium.density = 1000;          % [kg/m^3]
medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
medium.BonA = 6;

% create the time array
t_end = 40e-6;                  % [s]
kgrid.makeTime(medium.sound_speed, [], t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 0.5e6;        % [Hz]
tone_burst_cycles = 5;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength ./ (medium.sound_speed * medium.density)) .* input_signal;

% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% physical properties of the transducer
transducer.number_elements = 72;    % total number of transducer elements
transducer.element_width = 1;       % width of each element [grid points]
transducer.element_length = 12;     % length of each element [grid points]
transducer.element_spacing = 0;     % spacing (kerf width) between the elements [grid points]
transducer.radius = inf;            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = 1540;                  % sound speed [m/s]
transducer.focus_distance = 20e-3;              % focus distance [m]
transducer.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
transducer.steering_angle = 0;                  % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Rectangular';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
transducer.active_elements = zeros(transducer.number_elements, 1);
transducer.active_elements(21:52) = 1;

% append input signal used to drive the transducer
transducer.input_signal = input_signal;

% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties;

% =========================================================================
% DEFINE SENSOR MASK
% =========================================================================

% create a binary sensor mask with four detection positions
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask([Nx/4, Nx/2, 3*Nx/4], Ny/2, Nz/2) = 1;

% =========================================================================
% RUN THE SIMULATION
% =========================================================================

% set the input settings
input_args = {'DisplayMask', transducer.all_elements_mask | sensor.mask, ...
    'PMLInside', false, 'PlotPML', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
    'DataCast', DATA_CAST, 'PlotScale', [-1/2, 1/2] * source_strength};

% run the simulation
[sensor_data] = kspaceFirstOrder3D(kgrid, medium, transducer, sensor, input_args{:});

% calculate the amplitude spectrum of the input signal and the signal
% recorded each of the sensor positions 
[f_input, as_input] = spect([input_signal, zeros(1, 2 * length(input_signal))], 1/kgrid.dt);
[~, as_1] = spect(sensor_data(1, :), 1/kgrid.dt);
[~, as_2] = spect(sensor_data(2, :), 1/kgrid.dt);
[f, as_3] = spect(sensor_data(3, :), 1/kgrid.dt);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the input signal and its frequency spectrum
figure;
subplot(2, 1, 1);
plot((0:kgrid.dt:(length(input_signal) - 1) * kgrid.dt) * 1e6, input_signal, 'k-');
xlabel('Time [\mus]');
ylabel('Particle Velocity [m/s]');

subplot(2, 1, 2);
plot(f_input .* 1e-6, as_input./max(as_input(:)), 'k-');
hold on;
line([tone_burst_freq, tone_burst_freq] .* 1e-6, [0 1], 'Color', 'k', 'LineStyle', '--');
xlabel('Frequency [MHz]');
ylabel('Amplitude Spectrum [au]');
f_max = medium.sound_speed / (2 * dx);
set(gca, 'XLim', [0, f_max .* 1e-6]);

% plot the recorded time series
figure;
stackedPlot(kgrid.t_array * 1e6, {'Sensor Position 1', 'Sensor Position 2', 'Sensor Position 3'}, sensor_data);
xlabel('Time [\mus]');

% plot the corresponding amplitude spectrums
figure;
plot(f .* 1e-6, as_1 ./ max(as_1(:)), 'k-', ...
     f .* 1e-6, as_2 ./ max(as_1(:)), 'b-', ...
     f .* 1e-6, as_3 ./ max(as_1(:)), 'r-');
legend('Sensor Position 1', 'Sensor Position 2', 'Sensor Position 3');
xlabel('Frequency [MHz]');
ylabel('Normalised Amplitude Spectrum [au]');
f_max = medium.sound_speed / (2 * dx);
set(gca, 'XLim', [0, f_max .* 1e-6]);