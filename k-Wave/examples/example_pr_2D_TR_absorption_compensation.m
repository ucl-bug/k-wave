% Attenuation Compensation Using Time Reversal Example
%
% This example demonstrates how the acoustic attenuation present in the
% photoacoustic forward problem can be compensated for using time reversal
% image reconstruction. It builds on the 2D Time Reversal Reconstruction
% For A Circular Sensor Example. 
%
% For a more detailed discussion of this example and the underlying
% techniques, see B. E. Treeby, E. Z. Zhang, and B. T. Cox, "Photoacoustic
% tomography in absorbing acoustic media using time reversal," Inverse
% Problems, vol. 26, no. 11, p. 115003, 2010.   
%
% author: Bradley Treeby
% date: 6th September 2010
% last update: 25th July 2019
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2010-2019 Bradley Treeby

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
% FORWARD SIMULATION
% =========================================================================

% define the size of the simulation grid and the PML
simulation_size = 512;                      % [grid points]
PML_size_forw = 20;                         % [grid points]
x = 52e-3;                                  % [m]
y = x;                                      % [m]

% reduce the number of grid points in Nx and Ny by the size of the PML so
% that using 'PMLInside' set to false will still give the correct
% simulation size 
Nx = simulation_size - 2 * PML_size_forw;   % [grid points]    
Ny = Nx;                                    % [grid points]
dx = x / Nx;                                % [m]
dy = dx;                                    % [m]

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium_non_absorbing.sound_speed = 1510;	% [m/s]

% create a duplicate of the propagation medium structure and append the
% absorption properties
medium = medium_non_absorbing;
medium.alpha_power = 1.5;      
medium.alpha_coeff = 3;                     % [dB/(MHz^y cm)]

% store maximum supported frequency
f_max = kgrid.k_max * medium.sound_speed / (2 * pi);

% load the shepp logan phantom (note if the simulation or PML sizes are
% changed, the loaded data will need to be resized)
load EXAMPLE_shepp_logan;

% smooth the phantom and assign it to the initial pressure
shepp_logan = smooth(shepp_logan, true);
source.p0 = shepp_logan;

% define a circular Cartesian sensor mask
sensor_radius = 25e-3;                      % [m]
sensor_points = 200;
cart_sensor_mask = makeCartCircle(sensor_radius, sensor_points);
sensor.mask = cart_sensor_mask;

% create the time array used for the simulation, with t_max defined using
% Huygens' principle to avoid artifact trapping in the reconstruction
t_max = 2 * sensor_radius / medium.sound_speed;
kgrid.makeTime(medium.sound_speed, [], t_max);

% set the input options, switching off the smoothing (the input has already
% been smoothed), setting the PML to be outside the defined grid, casting
% to 'single' to speed up the example, and switching off visualisation
input_args = {'Smooth', false, 'PMLInside', false, ...
    'PMLSize', PML_size_forw, 'DataCast', 'single', 'PlotSim', false};

% run the forward simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% add noise to the recorded sensor data
signal_to_noise_ratio = 40;                 % [dB]
sensor_data = addNoise(sensor_data, signal_to_noise_ratio, 'peak');

% =========================================================================
% IMAGE RECONSTRUCTION WITHOUT ABSORPTION COMPENSATION
% =========================================================================

% create a second computation grid for the reconstruction to avoid the
% inverse crime
PML_size_recon = 25;                        % [grid points]
Nx = simulation_size - 2 * PML_size_recon;  % [grid points]
Ny = Nx;                                    % [grid points]
dx = x / Nx;                                % [m]
dy = dx;                                    % [m]
kgrid_recon = kWaveGrid(Nx, dx, Ny, dy);

% use the same time array for the reconstruction
kgrid_recon.setTime(kgrid.Nt, kgrid.dt);

% remove the initial pressure field from the source structure
source = rmfield(source, 'p0');

% create a continuous binary sensor mask with the same radius as the
% Cartesian sensor mask used in the forward simulation
radius_gp = round(sensor_radius / kgrid_recon.dx);
binary_sensor_mask = makeCircle(kgrid_recon.Nx, kgrid_recon.Ny, floor(kgrid_recon.Nx/2) + 1, floor(kgrid_recon.Ny/2) + 1, radius_gp);

% assign the sensor mask to the sensor structure
sensor.mask = binary_sensor_mask;

% interpolate the simulated sensor data onto the continuous binary sensor
% mask to remove any gaps and assign to the time reversal field
sensor.time_reversal_boundary_data = interpCartData(kgrid_recon, sensor_data, cart_sensor_mask, binary_sensor_mask, 'linear');

% re-assign the input options (the PML_size has changed)
input_args = {'Smooth', false, 'PMLInside', false, ...
    'PMLSize', PML_size_recon, 'DataCast', 'single', 'PlotSim', false};

% run the time-reversal reconstruction using the non absorbing medium
p0_recon = kspaceFirstOrder2D(kgrid_recon, medium_non_absorbing, source, sensor, input_args{:});

% =========================================================================
% CHOOSING THE CUTOFF FREQUENCY
% =========================================================================

% get the average frequency spectrum of the simulated sensor data
num_signals = length(sensor_data(:, 1));
[as_f, as] = spect(sensor_data(1, :), 1/kgrid.dt);
for index = 2:num_signals
    [as_f, sp] = spect(sensor_data(index, :), 1/kgrid.dt);
    as = as + sp;
end
as = as/num_signals;

% compute the relative power spectrum
ps = 10 * log10(as.^2);
offset = max(ps(:));
ps = ps - offset;

% get the frequency spectrum of a single measurement
as_sing = spect(sensor_data(1, :), 1/kgrid.dt);
ps_sing = 10 * log10(as_sing.^2) - offset;

% scale the frequency variable
[f_sc, scale, prefix] = scaleSI(as_f(end));

% define the cutoff frequency for the filter
f_cutoff = 3e6;                     % [Hz]

% =========================================================================
% IMAGE RECONSTRUCTION WITH ABSORPTION COMPENSATION
% =========================================================================

% create the filter to regularise the absorption parameters
medium.alpha_filter = getAlphaFilter(kgrid_recon, medium, f_cutoff);

% reverse the sign of the absorption proportionality coefficient
medium.alpha_sign = [-1, 1];        % [absorption, dispersion];

% run the time-reversal reconstruction
p0_recon_compensated = kspaceFirstOrder2D(kgrid_recon, medium, source, sensor, input_args{:});

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the amplitude spectrum
figure;
plot(as_f * scale, ps_sing, 'r-');
ylabel('Power Spectrum [dB]');
xlabel(['Frequency [' prefix 'Hz]']);
hold on;
plot(as_f * scale, ps, 'k-');

% plot the cutoff frequency and maximum supported frequency
ylim = [-70, 0];
plot([f_cutoff * scale, f_cutoff * scale], [ylim(1), ylim(2)], 'k--');   
plot([f_max * scale, f_max * scale], [ylim(1), ylim(2)], 'k--');  
set(gca, 'XLim', [0, 15], 'Ylim', ylim);

% plot the filter, reducing the number of plotting points by a factor of 10
figure;
mesh(medium.alpha_filter(1:10:end, 1:10:end), 'EdgeColor', 'black');
axis tight; 
view([-29, 46]);

% set the pressure outside the sensor mask to be zero
binary_sensor_map = makeDisc(kgrid_recon.Nx, kgrid_recon.Ny, kgrid_recon.Nx/2, kgrid_recon.Ny/2, radius_gp - 50);
p0_recon(binary_sensor_map ~= 1) = 0;
p0_recon_compensated(binary_sensor_map ~= 1) = 0;

% set the plot limits
x_lim = [-20, 20];
y_lim = [-20, 20];

% plot the initial pressure distribution
figure;
subplot(2, 2, 1);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec  * 1e3, shepp_logan, [0, 1]);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
set(gca, 'XLim', x_lim, 'YLim', y_lim);

% plot the reconstruction without attenuation compensation
subplot(2, 2, 2);
imagesc(kgrid_recon.y_vec * 1e3, kgrid_recon.x_vec * 1e3, p0_recon, [0, 1]);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
set(gca, 'XLim', x_lim, 'YLim', y_lim);

% plot the reconstruction with attenuation compensation
subplot(2, 2, 3);
imagesc(kgrid_recon.y_vec * 1e3, kgrid_recon.x_vec * 1e3, p0_recon_compensated, [0, 1]);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
set(gca, 'XLim', x_lim, 'YLim', y_lim);
colormap(flipud(gray));
scaleFig(1, 1.25);

% plot a reconstruction profile
figure;
plot(kgrid.y_vec * 1e3, shepp_logan(end/2, :), 'k');
hold on;
plot(kgrid_recon.y_vec * 1e3, p0_recon(end/2, :), 'r-');
plot(kgrid_recon.y_vec * 1e3, p0_recon_compensated(end/2, :), 'b-');
set(gca, 'XLim', x_lim, 'YLim', [-0.1, 1.1]);
ylabel('Pressure Magnitude [au]');
xlabel('y-position [mm]');
legend('Original', 'No Compensation', 'With Attenuation Compensation', 'Location', 'Best');