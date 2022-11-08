% Modelling Nonlinear Wave Propagation Example
%
% This example describes the characteristics of the nonlinearity
% encapsulated by the first-order k-Wave simulation functions. 
%
% author: Bradley Treeby
% date: 8th December 2011
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

% =========================================================================
% DEFINE SIMULATION PROPERTIES
% =========================================================================

% define the properties used in the simulation
p0 = 10e6;                     % source pressure [Pa]
c0 = 1500;                     % sound speed [m/s]
rho0 = 1000;                   % density [kg/m^3]
alpha_0 = 0.25;                % absorption coefficient [dB/(MHz^2 cm)]
sigma = 2;                     % shock parameter
source_freq = 1e6;             % frequency [Hz]
points_per_wavelength = 100;   % number of grid points per wavelength at f0
wavelength_separation = 15;    % separation between the source and detector
pml_size = 80;                 % PML size
pml_alpha = 1.5;               % PML absorption coefficient [Np/grid point]
CFL = 0.25;                    % CFL number

% =========================================================================
% RUN SIMULATION
% =========================================================================

% compute grid spacing 
dx = c0 / (points_per_wavelength * source_freq);            % [m]

% compute grid size 
Nx = wavelength_separation * points_per_wavelength + 20;    % [grid points]

% create the computational grid
kgrid = kWaveGrid(Nx, dx);

% assign the properties of the propagation medium
medium.sound_speed = c0;
medium.density = rho0;
medium.alpha_power = 2;
medium.alpha_coeff = alpha_0;

% extract the maximum frequency supported by the grid
f_max = min(medium.sound_speed) / (2 * dx);     % [Hz]

% define a single source element
source_pos = 10;
source.p_mask = zeros(Nx, 1);
source.p_mask(source_pos) = 1;

% define a single sensor position an integer number of wavelengths away
sensor.mask = zeros(Nx, 1);
x_px = wavelength_separation * points_per_wavelength;
sensor.mask(source_pos + x_px) = 1;
x = x_px * dx;

% compute the nonlinearity coefficient required to give the correct shock
% parameter
mach_num = p0 / (rho0 * c0.^2);
k = 2 * pi * source_freq / c0;
BonA = 2 * (sigma / (mach_num * k * x) - 1);
medium.BonA = BonA;

% set the simulation options
input_args = {'PlotFreq', 20, 'PlotScale', [-1, 1] * p0 * 1.05, ...
    'PMLInside', false, 'PMLSize', pml_size, 'PMLAlpha', pml_alpha};

% compute points per temporal period
points_per_period = round(points_per_wavelength / CFL);

% compute corresponding time spacing
dt = 1 / (points_per_period * source_freq);    

% create the time array using an integer number of points per period
t_end = 25e-6;
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% set the start time to only record the last three periods
sensor.record_start_index = kgrid.Nt - 3 * points_per_period + 1;

% create the source term
source.p = p0 * sin(2 * pi * source_freq * kgrid.t_array);

% run the simulation
sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});

% create time axis for mendousse solution
t_axis = (0:(length(sensor_data) - 1)) .* dt; 

% compute mendousse solution for comparison
p_mendousse = mendousse(x * ones(size(t_axis)), t_axis, source_freq, p0, c0, rho0, BonA, alpha_0);

% get the amplitude spectra
[f, as_mendousse] = spect(p_mendousse, 1/dt);
f = f * 1e-6;
as_kspace = spect(sensor_data, 1/dt);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the time series
figure; 
subplot(2, 1, 1);
plot(t_axis * 1e6, sensor_data * 1e-6, 'r-');
hold on;
plot(t_axis * 1e6, p_mendousse * 1e-6, 'k--');
xlabel('Time [\mus]');
ylabel('Pressure [MPa]');
legend('k-Wave', 'Mendousse');

% get the number of harmonics
num_harm = floor(f_max / source_freq);

% extract and plot the data at the harmonics
subplot(2, 1, 2);
for harm = 1:num_harm

    % find index of frequency
    [f_val, f_index] = findClosest(f, harm);

    % plot simulation
    plot(f_val, as_kspace(f_index) * 1e-6, 'rx');
    hold on;
    
    % plot reference
    stem(f_val, as_mendousse(f_index) * 1e-6, 'Color', 'k', 'Marker', 'o');

end

% annotate the plot
xlabel('Frequency [MHz]');
ylabel('Amplitude [MPa]');
legend('k-Wave', 'Mendousse');
set(gca, 'XLim', [0.5, num_harm + 0.5]);

% plot the continuous frequency spectra and error as a percetange of the
% maximum amplitude
figure;
subplot(2, 1, 1);
plot(f(2:end), as_kspace(2:end) * 1e-6, 'r-');
hold on;
plot(f(2:end), as_mendousse(2:end) * 1e-6, 'k--');
xlabel('Frequency [MHz]');
ylabel('Amplitude [MPa]');
legend('k-Wave', 'Mendousse');
set(gca, 'XLim', [0, 2 * num_harm]);

subplot(2, 1, 2);
plot(f(2:end), 100 * abs(as_mendousse(2:end) - as_kspace(2:end)) ./ max(as_mendousse(2:end)), 'k-');
xlabel('Frequency [MHz]');
ylabel('Error [% of max]');
set(gca, 'XLim', [0, 2 * num_harm]);