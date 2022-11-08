% The Doppler Effect
%
% This example demonstrates the doppler effect in which a stationary sensor
% point records a shift in frequency as a moving source travels past. It
% builds on the Monopole Point Source In A Homogeneous Propagation Medium
% Example.   
%
% author: Bradley Treeby
% date: 23rd December 2010
% last update: 4th May 2017
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

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 64;            % number of grid points in the x (row) direction
Ny = Nx*2;          % number of grid points in the y (column) direction
dy = 20e-3/Ny;    	% grid point spacing in the y direction [m]
dx = dy;            % grid point spacing in the x direction [m]
pml_size = 20;      % [grid points]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500;      % [m/s]
medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5; 

% set the velocity of the moving source
source_vel = 150;               % [m/s]

% set the relative x-position between the source and sensor
source_sensor_x_distance = 5;   % [grid points]

% manually create the time array
Nt = 4500;
dt = 20e-9;                     % [s]
kgrid.setTime(Nt, dt);

% define a single time varying sinusoidal source
source_freq = 0.75e6;           % [MHz]
source_mag = 3;                 % [Pa]
source_pressure = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

% filter the source to remove high frequencies not supported by the grid
source_pressure = filterTimeSeries(kgrid, medium, source_pressure);

% define a line of source points
source_x_pos = 5;               % [grid points]
source.p_mask = zeros(Nx, Ny);
source.p_mask(end - pml_size - source_x_pos, 1 + pml_size:end - pml_size) = 1;

% preallocate an empty pressure source matrix
num_source_positions = sum(source.p_mask(:));
source.p = zeros(num_source_positions, length(kgrid.t_array));

% move the source along the source mask by interpolating the pressure
% series between the source elements
sensor_index = 1;
t_index = 1;
while t_index < length(kgrid.t_array) && sensor_index < num_source_positions - 1
    
    % check if the source has moved to the next pair of grid points
    if kgrid.t_array(t_index) > (sensor_index*dy/source_vel)
        sensor_index = sensor_index + 1;
    end    
    
    % calculate the position of source in between the two current grid
    % points
    exact_pos = source_vel * kgrid.t_array(t_index);
    discrete_pos = sensor_index * dy;
    pos_ratio = (discrete_pos - exact_pos) ./ dy;
    
    % update the pressure at the two current grid points using linear
    % interpolation
    source.p(sensor_index, t_index) = pos_ratio * source_pressure(t_index);
    source.p(sensor_index + 1, t_index) = (1 - pos_ratio) * source_pressure(t_index);
    
    % update the time index
    t_index = t_index + 1;
    
end

% define a single sensor point
sensor.mask = zeros(Nx, Ny);
sensor.mask(end - pml_size - source_x_pos - source_sensor_x_distance, Ny/2) = 1;

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PlotPML', false);

% compute the number of time steps before any sound is recorded at the
% sensor position [time steps]
source_lag = round(sqrt( (Ny/2 - pml_size)^2 + source_sensor_x_distance^2 ) * ...
    dx / (medium.sound_speed * kgrid.dt)) + 50;

% calculate the observed source frequency during the approach and depart
time_steps = 1500;
Fs = 1/kgrid.dt;
[~, approach_as] = spect(sensor_data(source_lag:source_lag + time_steps), Fs, 'Window', 'Hanning');
[f, retreat_as] = spect(sensor_data(end - time_steps:end), Fs, 'Window', 'Hanning');

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the moving pressure field
figure;
subplot(2, 1, 1);
imagesc(source.p);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
title('Input Pressure Signal');

subplot(2, 1, 2);
plot(sum(source.p));
ylabel('Pressure [au]');
xlabel('Time Step');
title('Sum Of Input Pressure Signal Across All Sensor Positions');

% plot the simulated sensor data
figure;
[t_sc, scale, prefix] = scaleSI(max(kgrid.t_array(:)));
plot(kgrid.t_array*scale, sensor_data, 'r-');
xlabel(['Time [' prefix 's]']);
ylabel('Signal Amplitude');
axis tight;
title('Sensor Pressure Signal');

% plot the frequency content of the signal as the source approaches the
% sensor and as it retreats
figure;
[f_sc, scale, prefix] = scaleSI(max(f(:)));
plot(f*scale, approach_as, 'b-', f*scale, retreat_as, 'r-');
xlabel(['Frequency [' prefix 'Hz]']);
ylabel('Amplitude Spectrum')
axis([0, 4, 0, 0.2])
title('Spectra During Approach and Retreat')
legend('approach', 'retreat')