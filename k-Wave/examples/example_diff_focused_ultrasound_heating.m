% Heating By A Focused Ultrasound Transducer
%
% This example demonstrates how to combine acoustic and thermal simulations
% in k-Wave to calculate the heating by a focused ultrasound transducer. It
% builds on the Simulating Transducer Field Patterns and Using A Binary
% Sensor Mask examples.
%
% author: Bradley Treeby
% date: 27th April 2017
% last update: 28th April 2017
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

clearvars;

% =========================================================================
% ACOUSTIC SIMULATION
% =========================================================================

% define the PML size
pml_size = 20;              % [grid points]

% define the grid parameters
Nx = 256 - 2 * pml_size;    % [grid points]
Ny = 256 - 2 * pml_size;    % [grid points]
dx = 0.25e-3;               % [m]
dy = 0.25e-3;               % [m]

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1510;  % [m/s]
medium.density     = 1020;  % [kg/m^3]
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;

% define the source parameters
diameter = 45e-3;           % [m]
radius   = 35e-3;           % [m]
freq     = 1e6;             % [Hz]
amp      = 0.5e6;           % [Pa]

% define a focused ultrasound transducer
source.p_mask = makeArc([Nx, Ny], [1, Ny/2], round(radius / dx), round(diameter / dx) + 1, [Nx/2, Ny/2]);

% calculate the time step using an integer number of points per period
ppw = medium.sound_speed / (freq * dx); % points per wavelength
cfl = 0.3;                              % cfl number
ppp = ceil(ppw / cfl);                  % points per period
T   = 1 / freq;                         % period [s]
dt  = T / ppp;                          % time step [s]

% calculate the number of time steps to reach steady state
t_end = sqrt( kgrid.x_size.^2 + kgrid.y_size.^2 ) / medium.sound_speed; 
Nt = round(t_end / dt);

% create the time array
kgrid.setTime(Nt, dt);

% define the input signal
source.p = createCWSignals(kgrid.t_array, freq, amp, 0);

% set the sensor mask to cover the entire grid
sensor.mask = ones(Nx, Ny);
sensor.record = {'p', 'p_max_all'};

% record the last 3 cycles in steady state
num_periods = 3;
T_points = round(num_periods * T / kgrid.dt);
sensor.record_start_index = Nt - T_points + 1;

% set the input arguements
input_args = {'PMLInside', false, 'PlotPML', false, 'DisplayMask', ...
    'off', 'PlotScale', [-1, 1] * amp};

% run the acoustic simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% CALCULATE HEATING
% =========================================================================

% convert the absorption coefficient to nepers/m
alpha_np = db2neper(medium.alpha_coeff, medium.alpha_power) * ...
    (2 * pi * freq).^medium.alpha_power;

% extract the pressure amplitude at each position
p = extractAmpPhase(sensor_data.p, 1/kgrid.dt, freq);

% reshape the data, and calculate the volume rate of heat deposition
p = reshape(p, Nx, Ny);
Q = alpha_np .* p.^2 ./ (medium.density .* medium.sound_speed);

% =========================================================================
% THERMAL SIMULATION
% =========================================================================

% clear the input structures
clear medium source sensor;

% set the background temperature and heating term
source.Q = Q;
source.T0 = 37;

% define medium properties related to diffusion
medium.density              = 1020;     % [kg/m^3]
medium.thermal_conductivity = 0.5;      % [W/(m.K)]
medium.specific_heat        = 3600;     % [J/(kg.K)]

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium, source, []);

% set source on time and off time
on_time  = 10;  % [s]
off_time = 20;  % [s]

% set time step size
dt = 0.1;

% take time steps
kdiff.takeTimeStep(round(on_time / dt), dt);

% store the current temperature field
T1 = kdiff.T;

% turn off heat source and take time steps
kdiff.Q = 0;
kdiff.takeTimeStep(round(off_time / dt), dt);

% store the current temperature field
T2 = kdiff.T;

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the thermal dose and lesion map
figure;

% plot the acoustic pressure
subplot(2, 3, 1);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p * 1e-6);
h = colorbar;
xlabel(h, '[MPa]');
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Acoustic Pressure Amplitude');

% plot the volume rate of heat deposition
subplot(2, 3, 2);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, Q * 1e-7);
h = colorbar;
xlabel(h, '[kW/cm^2]');
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Volume Rate Of Heat Deposition');

% plot the temperature after heating
subplot(2, 3, 3);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, T1);
h = colorbar;
xlabel(h, '[degC]');
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Temperature After Heating');

% plot the temperature after cooling
subplot(2, 3, 4);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, T2);
h = colorbar;
xlabel(h, '[degC]');
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Temperature After Cooling');

% plot thermal dose
subplot(2, 3, 5);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, kdiff.cem43, [0, 1000]);
h = colorbar;
xlabel(h, '[CEM43]');
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Thermal Dose');

% plot lesion map
subplot(2, 3, 6);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, kdiff.lesion_map, [0, 1]);
colorbar;
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Ablated Tissue');

% set colormap and enlarge figure window
colormap(jet(256));
scaleFig(1.5, 1);